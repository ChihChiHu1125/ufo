/*
 * Use the clear-sky BT to define model & obs cloud, use abs(model cloud - obs cloud) as criteria to filter
 * Chih-Chi Hu, Sep 2026
 */

#ifndef UFO_QCBIASFILTERS_CLOUDDIFFFILTER_H_
#define UFO_QCBIASFILTERS_CLOUDDIFFFILTER_H_

#include <string>
#include <vector>

#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/qcbiasfilters/QcBiasFilterBase.h"

namespace ufo {

class CloudDiffFilterParameters : public QcBiasFilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(CloudDiffFilterParameters, QcBiasFilterParametersBase)
 public:
  /// One threshold per channel, in the same order as this obs space's simulated channels.
  oops::RequiredParameter<std::vector<float>> threshold{"threshold", this};
};

/// \brief Flags an observation's contribution to VarBC coefficient update when the
/// model- and observation-space cloud signals (estimated from the clear-sky minus
/// all-sky brightness temperature) disagree by more than a channel-specific threshold.
class CloudDiffFilter : public QcBiasFilterBase {
 public:
  typedef CloudDiffFilterParameters Parameters_;

  CloudDiffFilter(const Parameters_ &, const oops::ObsVariables &);

  void compute(const ioda::ObsSpace &, const GeoVaLs &, const ObsDiagnostics &,
              const ioda::ObsVector & hofx, ioda::ObsDataVector<int> & qc_bias) const override;

 private:
  std::vector<float> threshold_;
};

}  // namespace ufo

#endif  // UFO_QCBIASFILTERS_CLOUDDIFFFILTER_H_
