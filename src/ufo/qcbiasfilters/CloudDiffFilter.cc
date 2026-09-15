/*
 * Chih-Chi Hu, Sep 2026
 */

#include "ufo/qcbiasfilters/CloudDiffFilter.h"

#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "ufo/ObsDiagnostics.h"

namespace ufo {

static QcBiasFilterMaker<CloudDiffFilter> makerCloudDiffFilter_("cloud_diff_filter");

// -----------------------------------------------------------------------------

CloudDiffFilter::CloudDiffFilter(const Parameters_ & parameters, const oops::ObsVariables & vars)
  : QcBiasFilterBase(parameters, vars)
{
  threshold_ = parameters.threshold.value();
  if (threshold_.size() != selectedChannels_.size()) {
    throw eckit::BadValue("cloud_diff_filter: 'threshold' must have one entry per channel "
                          "selected via 'channels' (or per channel of this obs space, "
                          "if 'channels' is not set)", Here());
  }

  hdiags_ += oops::ObsVariables({"brightness_temperature_assuming_clear_sky"}, selectedChannels_);
}

// -----------------------------------------------------------------------------

void CloudDiffFilter::compute(const ioda::ObsSpace & odb, const GeoVaLs &,
                              const ObsDiagnostics & ydiags, const ioda::ObsVector & hofx,
                              ioda::ObsDataVector<int> & qc_bias) const {
  const std::size_t nlocs = hofx.nlocs();

  std::vector<float> obsBT(nlocs);
  std::vector<float> btClr(nlocs);

  for (std::size_t jsel = 0; jsel < selectedChannels_.size(); ++jsel) {
    const int channel = selectedChannels_[jsel];
    const std::size_t jch = channelIndices_[jsel];   // position in the FULL channel list

    // observed brightness temperature for this channel, straight from the obs database
    odb.get_db("ObsValue", "brightnessTemperature", obsBT, {channel});
    // model clear-sky brightness temperature diagnostic for this channel
    ydiags.get(btClr, "brightness_temperature_assuming_clear_sky_" + std::to_string(channel));

    for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
      const std::size_t idx = jloc * hofx.nvars() + jch;
      const double modelBT = hofx[idx];   // pre-bias-correction model BT, per ObsOperator.cc

      const double modelCloud = std::max(0.0, static_cast<double>(btClr[jloc]) - modelBT);
      const double obsCloud   = std::max(0.0, static_cast<double>(btClr[jloc]) -
                                              static_cast<double>(obsBT[jloc]));

      qc_bias[jch][jloc] =
          (std::abs(modelCloud - obsCloud) > threshold_[jsel]) ? 1 : 0;
    }
  }
}

}  // namespace ufo
