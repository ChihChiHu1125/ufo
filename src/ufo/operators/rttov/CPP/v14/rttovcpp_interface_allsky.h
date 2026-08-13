/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_RTTOV_CPP_V14_RTTOVCPP_INTERFACE_ALLSKY_H_
#define UFO_OPERATORS_RTTOV_CPP_V14_RTTOVCPP_INTERFACE_ALLSKY_H_

#include <string>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsDiagnostics.h"

#include "rttov/wrapper/RttovProfile.h"
#include "rttov/wrapper/RttovSafe.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;
}

namespace ufo {

/// \brief Maps a yaml "Absorbers"/"MWClouds" entry to its GeoVaLs variable name and the
/// RTTOV profile item id used to read back the Jacobian (rttov::Profile::getItemK).
/// Single source of truth shared by the nonlinear profile setup (rttovcpp_interface_allsky.cc)
/// and the TL/AD Jacobian lookup (ObsRadianceRTTOVCPPTLAD_allsky.cc).
struct RttovCppItemSpec {
  std::string yamlName;      // e.g. "O3", "Water"
  std::string geovalsName;   // e.g. "mole_fraction_of_ozone_in_air",
                              // "cloud_liquid_water_mixing_ratio_wrt_moist_air"
  rttov::itemIdType item;    // e.g. rttov::O3, rttov::MW_CLW
};

/// Gas absorbers supported by the operator. H2O is mandatory (always in Absorbers and always
/// set on the profile in rttovcpp_interface_allsky.cc); O3 is optional. Both are listed here so
/// that the TL/AD (ObsRadianceRTTOVCPPTLAD_allsky.cc) can resolve either one generically by name.
const std::vector<RttovCppItemSpec>& rttovcppAbsorberSpecs();

/// The 5 hydrometeors of RTTOV's built-in MW scattering scheme
const std::vector<RttovCppItemSpec>& rttovcppMWCloudSpecs();

/// Looks up a spec by its yaml name in either table; throws eckit::BadParameter if not found.
const RttovCppItemSpec& rttovcppFindSpec(const std::string& yamlName);

/// The "Water"/"Ice" species of RTTOV's IR/VIS OPAC(liquid)/Baum(ice) full-layer scattering
/// scheme. Unlike rttovcppMWCloudSpecs()/rttovcppAbsorberSpecs(), the item id for "Water"
/// depends on which OPAC liquid-water-cloud regime is selected (yaml "LiquidCloudType"), so
/// this cannot be a fixed table; it is rebuilt on each call (cheap: 2 entries).
/// Throws eckit::BadParameter if liquidCloudType is not one of stco/stma/cucc/cucp/cuma/clwde.
std::vector<RttovCppItemSpec> rttovcppIRCloudSpecs(const std::string& liquidCloudType);

/// Looks up a spec by its yaml name ("Water" or "Ice") in rttovcppIRCloudSpecs(liquidCloudType);
/// throws eckit::BadParameter if not found.
const RttovCppItemSpec rttovcppFindIRCloudSpec(const std::string& yamlName,
                                                const std::string& liquidCloudType);

/// Maps an IRClouds yaml name ("Water"/"Ice") to the GeoVaLs holding its effective radius
/// (e.g. MPAS's own diagnosed particle size), used only when UseModelEffectiveRadius is
/// enabled. Note: this GeoVaLs is a *radius* in microns; RTTOV's setClwDeff/setBaumIceDeff
/// expect an effective *diameter* in microns (rttovcpp_interface.cc doubles it).
/// Throws eckit::BadParameter if yamlName is not "Water" or "Ice".
std::string rttovcppIRCloudEffectiveRadiusGeovalsName(const std::string& yamlName);

void rttovcpp_interface(const GeoVaLs &, const ioda::ObsSpace & odb_,
                        rttov::RttovSafe & aRttov_, const std::string CoefFileName,
                        const std::vector<int> channels_, std::size_t & nlevels,
                        std::vector<bool> & skip_profile,
                        const bool doO3, const std::vector<std::string> & mwCloudSpecies,
                        const std::vector<std::string> & irCloudSpecies,
                        const std::string & liquidCloudType,
                        const bool useModelEffectiveRadius,
                        const std::string & thermalSolver,
                        const std::string & overlapParam,
                        const std::string & hydrotablePath);

void rttovcpp_setYdiag(const GeoVaLs &, rttov::RttovSafe & aRttov_,
                      ObsDiagnostics & d,
                      const std::vector<int> channels_, std::size_t & nlevels,
                      const std::vector<std::string> & mwCloudSpecies,
                      const std::vector<std::string> & irCloudSpecies,
                      const std::string & liquidCloudType);

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_RTTOV_CPP_V14_RTTOVCPP_INTERFACE_ALLSKY_H_
