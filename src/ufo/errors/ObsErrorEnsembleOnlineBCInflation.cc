/*
 * (C) Copyright 2026 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/errors/ObsErrorEnsembleOnlineBCInflation.h"

#include <algorithm>
#include <cmath>
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

  if (method_ == "full_matrix") {
    // The Woodbury correction needs B_{b,c}^{-1} for every (channel, predictor); a
    // non-positive prior variance would silently produce a NaN/Inf inverse deep inside
    // localInverseMultiply(), once per grid point. Check it once, up front, instead.
    const Eigen::VectorXd & priorVar = biasCov_->priorVariances();
    const int npred = static_cast<int>(predictorNames_.size());
    for (int jvar = 0; jvar < static_cast<int>(nvars_); ++jvar) {
      for (int jpred = 0; jpred < npred; ++jpred) {
        const double bb = priorVar[jvar * npred + jpred];
        if (!(bb > 0.0)) {
          throw eckit::BadValue(
              "ObsErrorEnsembleOnlineBCInflation: 'full_matrix' inflation method requires a "
              "strictly positive bias coefficient prior variance for every (channel, "
              "predictor); got " + std::to_string(bb) + " for predictor '" +
              predictorNames_[jpred] + "', channel index " + std::to_string(jvar) +
              ". Check 'bias.covariance.prior'.", Here());
        }
      }
    }
  }

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
  if (globalInvVarEffComputed_) return;
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
    // R_eff_jj = R_jj + inflation ;  (R_eff^-1)_jj = 1 / R_eff_jj  --- global, pre-localization.
    const double rjj = 1.0 / (*baseInvVar)[jj];
    globalInvVarEff_[jj] = 1.0 / (rjj + inflation);
  }
  globalInvVarEffComputed_ = true;
}

// -----------------------------------------------------------------------------

void ObsErrorEnsembleOnlineBCInflation::ensureGlobalInvVarRawComputed() const {
  if (globalInvVarRawComputed_ || method_ != "full_matrix") return;
  const std::unique_ptr<ioda::ObsVector> baseInvVar = base_->getInverseVariance();
  globalInvVarRaw_.resize(baseInvVar->size());
  for (size_t jj = 0; jj < baseInvVar->size(); ++jj) {
    globalInvVarRaw_[jj] = (*baseInvVar)[jj];
  }
  globalInvVarRawComputed_ = true;
}

// -----------------------------------------------------------------------------

void ObsErrorEnsembleOnlineBCInflation::localize(ioda::ObsVector & locvector) const {
  base_->localize(locvector);
  ensurePredictorsLoaded();
  ensureGlobalInvVarEffComputed();
  ensureGlobalInvVarRawComputed();

  const double missing = util::missingValue<double>();
  const size_t npred = predictorNames_.size();

  size_t nlocal = 0;
  for (size_t jj = 0; jj < locvector.size(); ++jj) {
    if (locvector[jj] != missing) ++nlocal;
  }

  localPredictors_.resize(npred, nlocal);
  localChannelIdx_.resize(nlocal);
  if (method_ == "diagonal") localInvVarEff_.resize(nlocal);
  if (method_ == "full_matrix") {
    localInvVarRaw_.resize(nlocal);
    localSqrtWeight_.resize(nlocal);
  }

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
      if (method_ == "full_matrix") {
        // Raw (un-inflated) R_jj^-1 and sqrt(w_j) are cached here, un-combined, because the
        // Woodbury correction in localInverseMultiply() needs them separately (see the
        // derivation in the comment above that function): the localization congruence
        // transform sqrt(w) R_eff^-1 sqrt(w) must be applied OUTSIDE the Woodbury correction,
        // not folded into P or R beforehand.
        localInvVarRaw_[col] = globalInvVarRaw_[jj];
        localSqrtWeight_[col] = std::sqrt(locvector[jj]);
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
// Woodbury correction for "full_matrix", derived per channel c (B_b has no cross-channel
// correlation, so channels are independent). Using
//   R_c        : n_c x n_c diagonal, RAW (un-localized) inverse variance of the wrapped R
//   P_c        : npred x n_c, RAW (un-weighted) predictor values for channel c's local obs
//   B_{b,c}^-1 : npred x npred diagonal, inverse prior variance of channel c's coefficients
//   D_w        : n_c x n_c diagonal, localization weights w_j for channel c's local obs
// the GLOBAL (pre-localization) Woodbury identity is
//   R_{c,eff}^-1 = R_c^-1 - R_c^-1 P_c^T (B_{b,c}^-1 + P_c R_c^-1 P_c^T)^-1 P_c R_c^-1
// and, per the ordering principle in the class-level comment, the LOCALIZED R_eff is obtained
// by the SAME congruence transform base_->localize() applies to a diagonal R, generalized to a
// matrix: R_{c,eff,local}^-1 = D_w^{1/2} R_{c,eff}^-1 D_w^{1/2}. Applying this to zz_c and
// simplifying (D_w^{1/2} must be applied OUTSIDE the Woodbury correction, not folded into P_c
// beforehand -- doing the latter does not reduce to the validated "diagonal" formula) gives,
// writing M_c = B_{b,c}^-1 + P_c R_c^-1 P_c^T and PR_c = P_c * diag(R_c^-1):
//   Z    = zz_c * diag(sqrt(w))
//   ZR   = Z * diag(R_c^-1)
//   W    = (ZR * P_c^T) * M_c^-1
//   zzRinv_c = (ZR - W * PR_c) * diag(sqrt(w))
// This reduces exactly to the validated "diagonal" formula w_j / (R_jj + inflation_j) in the
// n_c == npred == 1 case, which is how the derivation was checked.
Eigen::MatrixXf ObsErrorEnsembleOnlineBCInflation::localInverseMultiply(
    const Eigen::MatrixXf & zz) const {
  if (method_ != "full_matrix") {
    Eigen::MatrixXf zzRinv(zz.rows(), zz.cols());
    for (int ii = 0; ii < zz.rows(); ++ii) {
      zzRinv(ii, Eigen::placeholders::all) =
          zz(ii, Eigen::placeholders::all).cwiseProduct(localInvVarEff_.cast<float>().transpose());
    }
    return zzRinv;
  }

  const size_t npred = predictorNames_.size();
  const size_t nlocal = localChannelIdx_.size();
  const Eigen::VectorXd & priorVar = biasCov_->priorVariances();
  const Eigen::MatrixXd zzD = zz.cast<double>();

  // Group local columns by channel; channel indices are 0..nvars_-1.
  std::vector<std::vector<size_t>> channelCols(nvars_);
  for (size_t col = 0; col < nlocal; ++col) {
    channelCols[localChannelIdx_[col]].push_back(col);
  }

  Eigen::MatrixXd zzRinv = Eigen::MatrixXd::Zero(zz.rows(), zz.cols());
  for (size_t channel = 0; channel < nvars_; ++channel) {
    const std::vector<size_t> & cols = channelCols[channel];
    const size_t nc = cols.size();
    if (nc == 0) continue;

    Eigen::MatrixXd Pc(npred, nc);
    Eigen::VectorXd RcInv(nc);
    Eigen::VectorXd sqrtWc(nc);
    for (size_t k = 0; k < nc; ++k) {
      const size_t col = cols[k];
      for (size_t p = 0; p < npred; ++p) Pc(p, k) = localPredictors_(p, col);
      RcInv[k] = localInvVarRaw_[col];
      sqrtWc[k] = localSqrtWeight_[col];
    }

    Eigen::VectorXd BcInvDiag(npred);
    for (size_t p = 0; p < npred; ++p) {
      BcInvDiag[p] = 1.0 / priorVar[static_cast<int>(channel * npred + p)];
    }

    const Eigen::MatrixXd PRc = Pc * RcInv.asDiagonal();                     // npred x nc
    Eigen::MatrixXd Mc = BcInvDiag.asDiagonal();
    Mc.noalias() += Pc * PRc.transpose();                                   // npred x npred
    const Eigen::LDLT<Eigen::MatrixXd> McSolver(Mc);

    Eigen::MatrixXd Z(zz.rows(), nc);
    for (size_t k = 0; k < nc; ++k) Z.col(k) = zzD.col(cols[k]) * sqrtWc[k];

    const Eigen::MatrixXd ZR = Z * RcInv.asDiagonal();                       // rows x nc
    const Eigen::MatrixXd V = ZR * Pc.transpose();                          // rows x npred
    const Eigen::MatrixXd W = McSolver.solve(V.transpose()).transpose();    // rows x npred
    const Eigen::MatrixXd correction = W * PRc;                            // rows x nc
    const Eigen::MatrixXd resultC = ZR - correction;                       // rows x nc

    for (size_t k = 0; k < nc; ++k) zzRinv.col(cols[k]) = resultC.col(k) * sqrtWc[k];
  }
  return zzRinv.cast<float>();
}

// -----------------------------------------------------------------------------

Eigen::VectorXd ObsErrorEnsembleOnlineBCInflation::local_invVarR() const {
  if (method_ == "full_matrix") {
    throw eckit::NotImplemented(
        "ObsErrorEnsembleOnlineBCInflation: local_invVarR() is not implemented for the "
        "'full_matrix' inflation method, because R_eff_local^-1 is not diagonal in general "
        "(observations of the same channel are correlated through the Woodbury correction). "
        "local_invVarR() is only ever called by the GSI Fortran LETKF/GETKF weight-computation "
        "path ('fortran ETKF: true'); use the default C++ path (which calls "
        "localInverseMultiply() instead) with 'full_matrix'.", Here());
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

  // Also write out sqrt(R_eff) and sqrt(R_eff - R) (i.e. the stddev equivalent of the
  // online-BC inflation alone) as extra groups in the obs file, so the effect of the
  // inflation can be inspected/compared directly with ncdump/h5dump/python.
  // For "diagonal" this is exactly the R_eff actually used by the solver. For "full_matrix"
  // it is DIAGONAL-ONLY: it reflects R_jj + inflation_j (the diagonal of R + P B_b P^T) but
  // omits the off-diagonal Woodbury correction that the solver actually applies, so it is a
  // diagnostic view only, not the true R_eff for that method.
  ensureGlobalInvVarEffComputed();

  const std::unique_ptr<ioda::ObsVector> baseInvVar = base_->getInverseVariance();
  std::unique_ptr<ioda::ObsVector> effStddev = base_->getObsErrors();
  std::unique_ptr<ioda::ObsVector> inflationStddev = base_->getObsErrors();
  for (size_t jj = 0; jj < effStddev->size(); ++jj) {
    const double rbase = 1.0 / (*baseInvVar)[jj];
    const double reff = 1.0 / globalInvVarEff_[jj];
    (*effStddev)[jj] = std::sqrt(reff);
    (*inflationStddev)[jj] = std::sqrt(std::max(reff - rbase, 0.0));
  }
  effStddev->save(name + "EnsembleOnlineBCInflated");
  inflationStddev->save(name + "EnsembleOnlineBCInflationAmount");
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
