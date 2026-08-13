/*
 * (C) Copyright 2021 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_RTTOV_CPP_OBSRADIANCERTTOVCPPPARAMETERS_ALLSKY_H_
#define UFO_OPERATORS_RTTOV_CPP_OBSRADIANCERTTOVCPPPARAMETERS_ALLSKY_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/ObsOperatorParametersBase.h"

namespace ufo {

/// \brief Sub-parameters for the linear obs operator (TL/AD active variable selection).
/// Mirrors ufo::CRTMLinearObsOperatorParameters.
class RTTOVCPPLinearObsOperatorParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(RTTOVCPPLinearObsOperatorParameters, ObsOperatorParametersBase)

 public:
  /// Gas absorbers perturbed/adjointed by the TL/AD (e.g. [H2O] or [H2O, O3])
  oops::RequiredParameter<std::vector<std::string>> Absorbers{"Absorbers", this};
  /// MW hydrometeors perturbed/adjointed by the TL/AD (e.g. [Water, Ice, Rain, Snow, Graupel])
  oops::OptionalParameter<std::vector<std::string>> MWClouds{"MWClouds", this};
  /// IR/VIS cloud species perturbed/adjointed by the TL/AD (subset of {Water, Ice})
  oops::OptionalParameter<std::vector<std::string>> IRClouds{"IRClouds", this};
};

/// Configuration options recognized by the ObsRadianceRTTOVCPP operator.
class ObsRadianceRTTOVCPPParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsRadianceRTTOVCPPParameters, ObsOperatorParametersBase)

 public:
  oops::RequiredParameter<std::string> CoefPath
    {"CoefPath",
     "Path to optical depth coefficients file",
     this};

  oops::RequiredParameter<std::string> SensorID
    {"SensorID",
     "Name of optical depth coefficients file",
     this};

  /// Gas absorbers to include in the RTTOV profile. "H2O" is always required; "O3" may be
  /// added to read ozone from GeoVaLs instead of using RTTOV's climatology.
  oops::RequiredParameter<std::vector<std::string>> Absorbers{"Absorbers", this};

  /// Surface wind GeoVaLs convention. Only "uv" (eastward/northward wind components) is
  /// currently supported; present for parity with the CRTM operator's yaml.
  oops::Parameter<std::string> SurfaceWindGeoVars{"SurfaceWindGeoVars", "uv", this};

  /// Hydrometeors to pass to RTTOV's built-in MW scattering scheme, e.g.
  /// [Water, Ice, Rain, Snow, Graupel]. Omit (or leave empty) for clear-sky.
  /// Mutually exclusive with IRClouds (a single RttovSafe instance loads one hydrotable file).
  oops::OptionalParameter<std::vector<std::string>> MWClouds{"MWClouds", this};

  /// Cloud species to pass to RTTOV's IR/VIS OPAC(liquid)/Baum(ice) full-layer scattering
  /// scheme, e.g. [Water, Ice]. RTTOV's OPAC/Baum scheme has no rain/snow/graupel categories,
  /// so only "Water" and "Ice" are supported here. Omit (or leave empty) for clear-sky.
  /// Mutually exclusive with MWClouds.
  oops::OptionalParameter<std::vector<std::string>> IRClouds{"IRClouds", this};

  /// Which liquid-water-cloud scheme to assume for the "Water" species in IRClouds: one of the
  /// 5 discrete OPAC regimes -- stco (stratus continental), stma (stratus maritime), cucc
  /// (cumulus continental clean), cucp (cumulus continental polluted), cuma (cumulus maritime)
  /// -- each of which RTTOV tabulates at a single fixed effective diameter (setClwDeff /
  /// UseModelEffectiveRadius has no effect for these); or clwde, RTTOV's continuous
  /// effective-diameter-indexed liquid water scheme (hydrotable type "clw_deff"), which is the
  /// only "Water" option that actually varies with a supplied/model effective radius. Only used
  /// if "Water" is in IRClouds.
  oops::Parameter<std::string> LiquidCloudType{"LiquidCloudType", "clwde", this};

  /// If true, effective radius for IRClouds "Water"/"Ice" is read from GeoVaLs
  /// (effective_radius_of_cloud_liquid_water_particle / effective_radius_of_cloud_ice_particle
  /// -- e.g. MPAS's own diagnosed particle size, when config_microp_re is enabled in the MPAS
  /// namelist) and passed to RTTOV via setClwDeff/setBaumIceDeff, instead of letting RTTOV's
  /// own Martin(1994)/Ou-Liou(1995)-family parameterization compute it internally. Useful for
  /// an apples-to-apples comparison against CRTM, which always uses the model-diagnosed
  /// particle size (or a fixed fallback) for its "Water"/"Ice" cloud types. NOTE: for "Water",
  /// this is only meaningful when LiquidCloudType is "clwde" -- RTTOV's 5 discrete OPAC regimes
  /// each use a single fixed tabulated diameter regardless of this setting (see hydrotable
  /// "Effective diameters" entries: stco/stma/cucc/cucp/cuma each have exactly 1 tabulated
  /// value), so combining UseModelEffectiveRadius=true with "Water" in IRClouds and a
  /// non-"clwde" LiquidCloudType throws eckit::BadParameter rather than silently doing nothing.
  /// For "Ice" (Baum), this always has an effect since the Baum table is continuously indexed
  /// by diameter. Only used if IRClouds is non-empty; defaults to false (RTTOV's internal
  /// parameterization), matching the operator's original all-sky behaviour.
  oops::Parameter<bool> UseModelEffectiveRadius{"UseModelEffectiveRadius", true, this};

  /// Which thermal (IR) multiple-scattering solver RTTOV uses when MWClouds or IRClouds is
  /// non-empty: "dom" (Discrete Ordinates, most general-purpose), "chou" (Chou-scaling, a
  /// faster approximation), or "delta_edd" (delta-Eddington 2-stream, RTTOV's historical
  /// MW/RTTOV-SCATT method). Only used if MWClouds or IRClouds is non-empty; defaults to
  /// "delta_edd", matching RTTOV's own internal default -- this operator never overrode it
  /// before this option existed, so omitting it preserves the previous all-sky behaviour.
  oops::Parameter<std::string> ThermalSolver{"ThermalSolver", "delta_edd", this};

  /// Cloud overlap scheme RTTOV uses to combine per-layer cloud fractions (HydroFrac) into an
  /// effective column cloud cover: "auto_select", "max_random", "2col_max", "2col_weighted",
  /// or "2col_user". Only used if MWClouds or IRClouds is non-empty; defaults to "auto_select",
  /// matching RTTOV's own internal default -- this operator never overrode it before this
  /// option existed, so omitting it preserves the previous all-sky behaviour.
  oops::Parameter<std::string> OverlapParam{"OverlapParam", "auto_select", this};

  /// Path to the RTTOV hydrotable file. Required if MWClouds or IRClouds is non-empty.
  oops::OptionalParameter<std::string> HydrotablePath{"HydrotablePath", this};

  /// Selects which of the Absorbers/MWClouds/IRClouds are active (perturbed) in the TL/AD. If
  /// omitted, defaults to H2O only (clear-sky TL/AD), matching the operator's previous
  /// behaviour.
  oops::OptionalParameter<RTTOVCPPLinearObsOperatorParameters>
      LinearObsOperator{"linear obs operator", this};
};

}  // namespace ufo
#endif  // UFO_OPERATORS_RTTOV_CPP_OBSRADIANCERTTOVCPPPARAMETERS_ALLSKY_H_
