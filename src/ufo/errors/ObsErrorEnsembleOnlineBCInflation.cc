/*
 * (C) Copyright 2026 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/errors/ObsErrorEnsembleOnlineBCInflation.h"

#include <utility>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

namespace ufo {

// -----------------------------------------------------------------------------

static ObsErrorMaker<ObsErrorEnsembleOnlineBCInflation>
    makerEnsembleOnlineBCInflation_("ensemble online bc inflation");

// -----------------------------------------------------------------------------

ObsErrorEnsembleOnlineBCInflation::ObsErrorEnsembleOnlineBCInflation(
    const Parameters_ & params, ioda::ObsSpace & obsgeom, const eckit::mpi::Comm & timeComm)
  : ObsErrorBase(timeComm),
    base_(ObsErrorFactory::create(params.baseError.value().errorParams(), obsgeom)),
    method_(params.inflationMethod),
    obsSpace_(obsgeom),
    nvars_(obsgeom.assimvariables().size())
{
  if (params.bias.value().BiasCorrectionByRecord.value()) {
    throw eckit::NotImplemented(
        "ObsErrorEnsembleOnlineBCInflation: 'bc by record' is not yet supported; "
        "the bias coefficient prior variance indexing assumes a single record.", Here());
  }
  if (method_ != "diagonal" && method_ != "full_matrix") {
    throw eckit::BadParameter(
        "ObsErrorEnsembleOnlineBCInflation: 'inflation method' must be 'diagonal' or "
        "'full_matrix', got '" + method_ + "'", Here());
  }

  // Reuses ObsBiasCovariance's existing "read prior file + apply cycle-to-cycle inflation"
  // logic to obtain B_b, instead of re-implementing it here.
  biasCov_ = std::make_unique<ObsBiasCovariance>(obsgeom, params.bias.value().toConfiguration());
  predictorNames_ = biasCov_->predictorNames();

  oops::Log::info() << "ObsErrorEnsembleOnlineBCInflation: inflation method = " << method_
                    << ", " << predictorNames_.size()
                    << " variational-BC predictor(s): ";
  for (const auto & name : predictorNames_) oops::Log::info() << name << " ";
  oops::Log::info() << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorEnsembleOnlineBCInflation::ensurePredictorsLoaded() const {
  if (predictorsLoaded_) return;
  // Deferred until first use: the "<predictor>Predictor" groups are only written by
  // ObsBiasOperator as part of computeHofX(), which runs strictly after this object is
  // constructed (R_ is built before computeHofX4D() in LocalEnsembleSolver::computeHofX),
  // but strictly before localize() is ever called (measurementUpdate() runs after
  // computeHofX() completes). So it is safe -- and necessary -- to read them here rather
  // than in the constructor.
  predictorVectors_.reserve(predictorNames_.size());
  for (const auto & name : predictorNames_) {
    predictorVectors_.push_back(
        std::make_unique<ioda::ObsVector>(obsSpace_, name + "Predictor"));
  }
  predictorsLoaded_ = true;
}

// -----------------------------------------------------------------------------

void ObsErrorEnsembleOnlineBCInflation::ensureGlobalInvVarEffComputed() const {
  if (globalInvVarEffComputed_ || method_ != "diagonal") return;
  ensurePredictorsLoaded();

  // GLOBAL (pre-localization) inverse variance of the wrapped R -- deliberately NOT the
  // localized base_->local_invVarR(), so that R_eff = R + P B_b P^T is formed before any
  // localization weight is applied (see the ordering note in the class-level comment in the
  // header: localizing R first and adding a fixed inflation afterwards is NOT equivalent).
  const std::unique_ptr<ioda::ObsVector> baseInvVar = base_->getInverseVariance();
  const size_t npred = predictorNames_.size();
  const Eigen::VectorXd & priorVar = biasCov_->priorVariances();  // B_b diagonal, nrecs_ == 1
  const size_t ntotal = baseInvVar->size();

  globalInvVarEff_.resize(ntotal);
  for (size_t jj = 0; jj < ntotal; ++jj) {
    const int jvar = static_cast<int>(jj % nvars_);
    double inflation = 0.0;
    for (size_t p = 0; p < npred; ++p) {
      const double pval = (*predictorVectors_[p])[jj];
      const double bb = priorVar[jvar * static_cast<int>(npred) + static_cast<int>(p)];
      inflation += pval * pval * bb;
    }
    // R_eff_jj = R_jj + inflation ;  (R_eff^-1)_jj = 1 / R_eff_jj
    const double rjj = 1.0 / (*baseInvVar)[jj];
    globalInvVarEff_[jj] = 1.0 / (rjj + inflation);
  }
  globalInvVarEffComputed_ = true;
}

// -----------------------------------------------------------------------------

void ObsErrorEnsembleOnlineBCInflation::localize(ioda::ObsVector & locvector) const {
  base_->localize(locvector);
  ensurePredictorsLoaded();
  ensureGlobalInvVarEffComputed();

  const double missing = util::missingValue<double>();
  const size_t npred = predictorNames_.size();

  size_t nlocal = 0;
  for (size_t jj = 0; jj < locvector.size(); ++jj) {
    if (locvector[jj] != missing) ++nlocal;
  }

  localPredictors_.resize(npred, nlocal);
  localChannelIdx_.resize(nlocal);
  if (method_ == "diagonal") localInvVarEff_.resize(nlocal);

  size_t col = 0;
  for (size_t jj = 0; jj < locvector.size(); ++jj) {
    if (locvector[jj] != missing) {
      // "<predictor>Predictor" ObsVectors share the same (nlocs x nvars), variable-fastest
      // serialization as locvector itself (see ObsBiasOperator::computeObsBias), so jj % nvars_
      // recovers the channel index for this entry.
      localChannelIdx_[col] = static_cast<int>(jj % nvars_);
      for (size_t p = 0; p < npred; ++p) {
        localPredictors_(p, col) = static_cast<float>((*predictorVectors_[p])[jj]);
      }
      if (method_ == "diagonal") {
        // Apply the SAME localization weight to the already-inflated GLOBAL R_eff^-1 that
        // base_->localize() would apply to the un-inflated R (compare to
        // ObsErrorDiagonal::localize(): local_inverseVariance_ = locvector[jj] / stddev_[jj]^2).
        localInvVarEff_[col] = locvector[jj] * globalInvVarEff_[jj];
      }
      ++col;
    }
  }

  // Defensive consistency check: our masking (locvector != missing) must select exactly the
  // same local subset, in the same order, as base_->localize() did -- otherwise P_local and
  // R_local would not correspond to the same observations.
  if (static_cast<int>(nlocal) != base_->localDim()) {
    throw eckit::BadValue(
        "ObsErrorEnsembleOnlineBCInflation::localize: local dimension mismatch between "
        "predictor masking (" + std::to_string(nlocal) + ") and base_->localDim() (" +
        std::to_string(base_->localDim()) + ").", Here());
  }
}

// -----------------------------------------------------------------------------

Eigen::MatrixXf ObsErrorEnsembleOnlineBCInflation::localInverseMultiply(
    const Eigen::MatrixXf & zz) const {
  if (method_ == "full_matrix") {
    throw eckit::NotImplemented(
        "ObsErrorEnsembleOnlineBCInflation: 'full_matrix' (Woodbury) inflation method is "
        "not yet implemented; use 'diagonal' for now.", Here());
  }

  Eigen::MatrixXf zzRinv(zz.rows(), zz.cols());
  for (int ii = 0; ii < zz.rows(); ++ii) {
    zzRinv(ii, Eigen::placeholders::all) =
        zz(ii, Eigen::placeholders::all).cwiseProduct(localInvVarEff_.cast<float>().transpose());
  }
  return zzRinv;
}

// -----------------------------------------------------------------------------

Eigen::VectorXd ObsErrorEnsembleOnlineBCInflation::local_invVarR() const {
  if (method_ == "full_matrix") {
    throw eckit::NotImplemented(
        "ObsErrorEnsembleOnlineBCInflation: local_invVarR() for 'full_matrix' inflation "
        "method is not yet implemented.", Here());
  }
  return localInvVarEff_;
}

// -----------------------------------------------------------------------------

int ObsErrorEnsembleOnlineBCInflation::localDim() const {
  return base_->localDim();
}

// -----------------------------------------------------------------------------
// Everything below operates on the full (non-localized) R and is delegated unchanged to
// base_: the online-BC correction only applies to the local analysis-space operator used
// by LETKF/GETKF (localize()/localInverseMultiply()/local_invVarR() above).

void ObsErrorEnsembleOnlineBCInflation::update(const ioda::ObsVector & stddev) {
  base_->update(stddev);
}

void ObsErrorEnsembleOnlineBCInflation::multiply(ioda::ObsVector & dy) const {
  base_->multiply(dy);
}

void ObsErrorEnsembleOnlineBCInflation::inverseMultiply(ioda::ObsVector & dy) const {
  base_->inverseMultiply(dy);
}

void ObsErrorEnsembleOnlineBCInflation::randomize(ioda::ObsVector & dy) const {
  base_->randomize(dy);
}

void ObsErrorEnsembleOnlineBCInflation::save(const std::string & name) const {
  base_->save(name);
}

double ObsErrorEnsembleOnlineBCInflation::getRMSE() const {
  return base_->getRMSE();
}

std::unique_ptr<ioda::ObsVector> ObsErrorEnsembleOnlineBCInflation::getObsErrors() const {
  return base_->getObsErrors();
}

std::unique_ptr<ioda::ObsVector> ObsErrorEnsembleOnlineBCInflation::getInverseVariance() const {
  return base_->getInverseVariance();
}

// -----------------------------------------------------------------------------

void ObsErrorEnsembleOnlineBCInflation::print(std::ostream & os) const {
  os << "ObsErrorEnsembleOnlineBCInflation wrapping: " << *base_
     << " with inflation method: " << method_
     << ", " << predictorNames_.size() << " variational-BC predictor(s)" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
