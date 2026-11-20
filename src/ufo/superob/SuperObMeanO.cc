/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/superob/SuperObMeanO.h"

#include <optional>
#include <cassert>
#include <map>
#include <string>

#include "eckit/exception/Exceptions.h"

namespace ufo {

static SuperObMaker<SuperObMeanO> makerSuperObMeanO_("mean obs");

SuperObMeanO::SuperObMeanO(const Parameters_ &params,
                           const ObsFilterData &obsdb,
                           const std::vector<bool> &apply,
                           const Variables &filtervars,
                           const ioda::ObsDataVector<int> &flags,
                           std::vector<std::vector<bool>> &flagged)
    : SuperObBase(params, obsdb, apply, filtervars, flags, flagged),
      params_(params) {}

bool SuperObMeanO::computeSuperOb(const std::vector<std::size_t> &locs,
                                  const std::vector<float> &obs,
                                  const std::vector<float> & /*hofx*/,
                                  const ioda::ObsDataRow<int> &flags,
                                  std::vector<float> &superobs,
                                  std::vector<bool> &flagged) const {
  const float missing = util::missingValue<float>();

  const bool assignToAllValues = params_.assignToAllValues.value();

  // Deduplicate locations if grouping values are provided
  const std::vector<std::size_t> &locsForComputation =
      params_.groupingVariable.value()
          ? getUniqueLocations(locs, *params_.groupingVariable.value(), obs,
                               flags)
          : locs;

  bool superObComputed = false;
  if (locsForComputation.empty()) return superObComputed;

  float sumO = 0.0f;
  size_t count = 0;
  // Set on first valid contributing location; optional avoids implying a
  // default location.
  std::optional<size_t> superobloc;

  for (std::size_t jloc : locsForComputation) {
    const float obsValue = obs[jloc];
    // Only consider locations which have valid observation values
    // and are either passing QC or have been marked as passive (i.e. H(x) is
    // computed but the observation is not assimilated).
    if ((flags[jloc] != QCflags::pass && flags[jloc] != QCflags::passive) ||
        obsValue == missing) {
      continue;
    }
    sumO += obsValue;
    count++;
    if (!superObComputed) {
      superObComputed = true;
      superobloc = jloc;
    }
  }

  if (superObComputed) {
    const float meanO = sumO / static_cast<float>(count);
    assert(superobloc.has_value());
    assignSuperObToLocations(locs, meanO, *superobloc, assignToAllValues,
                             superobs, flagged);
  }
  return superObComputed;
}

}  // namespace ufo
