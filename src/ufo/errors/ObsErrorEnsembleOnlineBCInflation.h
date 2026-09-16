/*
 * (C) Copyright 2026 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_ERRORS_OBSERRORENSEMBLEONLINEBCINFLATION_H_
#define UFO_ERRORS_OBSERRORENSEMBLEONLINEBCINFLATION_H_

#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "ioda/ObsVector.h"

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/parameters/Parameters.h"

#include "ufo/errors/ObsErrorBase.h"
#include "ufo/errors/ObsErrorParametersBase.h"
#include "ufo/ObsBiasCovariance.h"
#include "ufo/ObsBiasParameters.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {

/// \brief Parameters for ObsErrorEnsembleOnlineBCInflation.
class ObsErrorEnsembleOnlineBCInflationParameters : public ObsErrorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsErrorEnsembleOnlineBCInflationParameters, ObsErrorParametersBase)
 public:
  /// "diagonal": R_eff = R + diag(P * B_b * P^T)  (only the diagonal terms; no matrix inversion).
  /// "full_matrix": Woodbury-based non-diagonal correction (not yet implemented).
  oops::RequiredParameter<std::string> inflationMethod{"inflation method", this};

  /// The R covariance model this inflation is applied on top of (e.g. "diagonal").
  oops::Parameter<ObsErrorParametersWrapper> baseError{"base error", {}, this};

  /// The full "obs bias" configuration for this obs space (typically referenced via a YAML
  /// anchor from the same observer's "obs bias" block, so predictors/covariance settings are
  /// written once and not duplicated). Must include a "covariance.prior" so that the bias
  /// coefficient prior error variances (B_b) can be read.
  oops::RequiredParameter<ObsBiasParameters> bias{"bias", this};
};

// -----------------------------------------------------------------------------
/// \brief Observation error covariance with an additional inflation term that accounts for
/// the prior uncertainty of online (ensemble) variational bias-correction coefficients.
///
/// \details Implements
///   R_eff = R + P * B_b * P^T
/// where P is the matrix of bias-correction predictor values for the (variational-BC)
/// predictors actually used, and B_b is the (diagonal) prior error variance of the
/// corresponding bias coefficients. Because B_b is diagonal and channel-specific (no
/// cross-channel prior correlation), R_eff is block-diagonal by channel.
///
/// Two methods are supported, selected via "inflation method":
///  - "diagonal": only the diagonal of P B_b P^T is added to R (cheap, no inversion needed).
///  - "full_matrix": the full non-diagonal correction, computed per-channel via the Woodbury
///    identity. Not yet implemented.
///
/// \note On ordering: R_eff = R + P B_b P^T is formed at the GLOBAL (pre-localization) level,
/// and the localization weight is then applied to R_eff exactly the way base_ would apply it
/// to R alone (i.e. localize(R_eff), not localize(R) + inflation). These are NOT equivalent
/// in general: localizing first and adding a fixed inflation afterwards would leave the
/// inflation term unaffected by distance from the analysis point, which is not the intended
/// behaviour -- the bias-correction-uncertainty contribution to the observation error should
/// be down-weighted by localization exactly like the rest of R.
///
/// This class wraps an existing ObsErrorBase ("base error") and only overrides localize()/
/// localInverseMultiply() (and local_invVarR() for the diagonal method); all other operations
/// (multiply, inverseMultiply, randomize, save, update, ...) are delegated unchanged to the
/// base R, since the online-BC correction only applies to the local analysis-space operator
/// used by LETKF/GETKF.
class ObsErrorEnsembleOnlineBCInflation : public ObsErrorBase {
 public:
  /// The type of parameters for this class.
  typedef ObsErrorEnsembleOnlineBCInflationParameters Parameters_;

  static const std::string classname() {return "ufo::ObsErrorEnsembleOnlineBCInflation";}

  ObsErrorEnsembleOnlineBCInflation(const Parameters_ &, ioda::ObsSpace &,
                                    const eckit::mpi::Comm &);

  // --- Delegated unchanged to the wrapped base_ error covariance. ---
  void update(const ioda::ObsVector &) override;
  void multiply(ioda::ObsVector &) const override;
  void inverseMultiply(ioda::ObsVector &) const override;
  void randomize(ioda::ObsVector &) const override;
  void save(const std::string &) const override;
  double getRMSE() const override;
  std::unique_ptr<ioda::ObsVector> getObsErrors() const override;
  std::unique_ptr<ioda::ObsVector> getInverseVariance() const override;
  int localDim() const override;

  // --- Where the R_eff correction actually happens. ---
  void localize(ioda::ObsVector &) const override;
  Eigen::MatrixXf localInverseMultiply(const Eigen::MatrixXf & zz) const override;
  Eigen::VectorXd local_invVarR() const override;

 private:
  void print(std::ostream &) const override;

  /// Ensures predictorVectors_ has been read from the ObsSpace (lazily, since the
  /// "<predictor>Predictor" groups are only written by ObsBiasOperator during computeHofX,
  /// which runs after this object is constructed but before localize() is ever called).
  void ensurePredictorsLoaded() const;

  /// Computes globalInvVarEff_ = 1 / (R_jj + inflation_j) for every observation in this obs
  /// space (NOT just the local subset), where inflation_j = P_j^T B_{b,channel(j)} P_j.
  /// Deliberately done at the GLOBAL (pre-localization) level -- see the class-level comment
  /// on ordering -- so that localize() can apply the localization weight to R_eff exactly the
  /// same way base_ would apply it to the un-inflated R. No-op after the first call, and a
  /// no-op entirely when method_ != "diagonal".
  void ensureGlobalInvVarEffComputed() const;

  std::unique_ptr<ObsErrorBase> base_;
  std::unique_ptr<ObsBiasCovariance> biasCov_;
  std::vector<std::string> predictorNames_;   ///< variational-BC predictor names (from biasCov_)
  std::string method_;                        ///< "diagonal" or "full_matrix"

  ioda::ObsSpace & obsSpace_;
  size_t nvars_;                               ///< number of channels/variables (assimvariables())

  // Lazily-loaded, read-only for the lifetime of the object once loaded.
  mutable bool predictorsLoaded_ = false;
  mutable std::vector<std::unique_ptr<ioda::ObsVector>> predictorVectors_;

  // Lazily-computed once (method_ == "diagonal" only); GLOBAL (pre-localization) effective
  // inverse variance 1/(R_jj + inflation_j), indexed the same way as any ioda::ObsVector
  // (size nlocs x nvars, channel-fastest).
  mutable bool globalInvVarEffComputed_ = false;
  mutable Eigen::VectorXd globalInvVarEff_;

  // Filled in by localize().
  //  - localPredictors_/localChannelIdx_: local subset of P, kept for the future "full_matrix"
  //    (Woodbury) method, which is not yet implemented.
  //  - localInvVarEff_: the localization weight applied to globalInvVarEff_ (diagonal method);
  //    consumed by localInverseMultiply()/local_invVarR().
  mutable Eigen::MatrixXf localPredictors_;    ///< (n_pred x n_local)
  mutable std::vector<int> localChannelIdx_;   ///< size n_local; channel index of each local obs
  mutable Eigen::VectorXd localInvVarEff_;     ///< size n_local; w_j * globalInvVarEff_[j]
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_ERRORS_OBSERRORENSEMBLEONLINEBCINFLATION_H_
