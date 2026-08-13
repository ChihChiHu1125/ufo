/*
 * (C) Copyright 2017-2021 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <algorithm>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "ioda/ObsVector.h"

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/operators/rttov/CPP/ObsRadianceRTTOVCPP_allsky.h"
#include "ufo/operators/rttov/CPP/v14/rttovcpp_interface_allsky.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsRadianceRTTOVCPP> makerRTTOVCPP_("RTTOVCPP");

// -----------------------------------------------------------------------------

ObsRadianceRTTOVCPP::ObsRadianceRTTOVCPP(const ioda::ObsSpace & odb,
                                         const Parameters_ & parameters)
  : ObsOperatorBase(odb), odb_(odb), varin_()
{
  // Fields to be requested from getvalues and stored in geovals
  // need to be consistent with those defined in ufo_variables_mod.F90
  //-----------------------------------------------------------------------------
  const std::vector<std::string> vv{
    "air_pressure_levels",
    "air_temperature",
    "water_vapor_mixing_ratio_wrt_moist_air",
    "air_pressure_at_surface",
    "air_temperature_at_2m",   // this is actually var_sfc_t2m
    "water_vapor_mixing_ratio_wrt_moist_air_at_2m",
    "eastward_wind_at_10m",
    "northward_wind_at_10m",
    "skin_temperature_at_surface",
    "seaice_fraction",   // this is actually var_sfc_seaicefrac
    "landmask",   // this is actually var_sfc_landmask
    "geopotential_height_at_surface"
  };

  for (size_t jvar = 0; jvar < vv.size(); ++jvar) {
     varin_.push_back(vv[jvar]);  // set private data member varin_
  }

  // Absorbers: H2O is implicit (already in vv above); O3 is optional
  const std::vector<std::string> absorbers = parameters.Absorbers;
  if (std::find(absorbers.begin(), absorbers.end(), "H2O") == absorbers.end()) {
    throw eckit::BadParameter("ObsRadianceRTTOVCPP: Absorbers must include 'H2O'");
  }
  for (const std::string & name : absorbers) {
    if (name == "H2O") continue;
    doO3_ = doO3_ || (name == "O3");
    varin_.push_back(rttovcppFindSpec(name).geovalsName);  // validates the name too
  }

  // Only the uv surface wind convention is currently implemented
  if (parameters.SurfaceWindGeoVars.value() != "uv") {
    throw eckit::BadParameter("ObsRadianceRTTOVCPP: SurfaceWindGeoVars only supports 'uv'");
  }

  // MWClouds: RTTOV's built-in MW scattering scheme (empty => clear-sky)
  if (parameters.MWClouds.value() != boost::none) {
    mwCloudSpecies_ = *parameters.MWClouds.value();
    for (const std::string & name : mwCloudSpecies_) {
      varin_.push_back(rttovcppFindSpec(name).geovalsName);  // validates the name too
    }
  }

  // IRClouds: RTTOV's built-in IR/VIS OPAC(liquid)/Baum(ice) scattering scheme
  // (empty => clear-sky). Mutually exclusive with MWClouds (see rttovcpp_interface).
  liquidCloudType_ = parameters.LiquidCloudType.value();
  useModelEffectiveRadius_ = parameters.UseModelEffectiveRadius.value();
  thermalSolver_ = parameters.ThermalSolver.value();
  overlapParam_ = parameters.OverlapParam.value();
  if (parameters.IRClouds.value() != boost::none) {
    irCloudSpecies_ = *parameters.IRClouds.value();
    for (const std::string & name : irCloudSpecies_) {
      // validates the name and, for "Water", the LiquidCloudType too
      varin_.push_back(rttovcppFindIRCloudSpec(name, liquidCloudType_).geovalsName);
      if (useModelEffectiveRadius_) {
        // For "Water", RTTOV only actually varies with a supplied effective diameter when
        // LiquidCloudType is "clwde" -- the 5 discrete OPAC regimes are each tabulated at a
        // single fixed diameter, so UseModelEffectiveRadius would silently have no effect.
        // Reject this combination outright rather than let it silently do nothing.
        if (name == "Water" && liquidCloudType_ != "clwde") {
          throw eckit::BadParameter("ObsRadianceRTTOVCPP: UseModelEffectiveRadius has no "
                                     "effect on \"Water\" unless LiquidCloudType is \"clwde\" "
                                     "(RTTOV's 5 discrete OPAC regimes each use a single fixed "
                                     "tabulated diameter) -- got LiquidCloudType '" +
                                     liquidCloudType_ + "'");
        }
        varin_.push_back(rttovcppIRCloudEffectiveRadiusGeovalsName(name));
      }
    }
  }

  if (!mwCloudSpecies_.empty() && !irCloudSpecies_.empty()) {
    throw eckit::BadParameter("ObsRadianceRTTOVCPP: MWClouds and IRClouds are mutually "
                               "exclusive");
  }
  if (!mwCloudSpecies_.empty() || !irCloudSpecies_.empty()) {
    varin_.push_back("cloud_area_fraction_in_atmosphere_layer");
    if (parameters.HydrotablePath.value() == boost::none) {
      throw eckit::BadParameter("ObsRadianceRTTOVCPP: HydrotablePath is required when "
                                 "MWClouds or IRClouds is non-empty");
    }
    hydrotablePath_ = *parameters.HydrotablePath.value();
  }

  // get channels from observations
  const oops::ObsVariables & observed = odb.assimvariables();
  channels_ = observed.channels();  // set private data member channels_

  // get optical depth coef file name from yaml
  const std::string CoefPath = parameters.CoefPath;
  const std::string SensorID = parameters.SensorID;
  CoefFileName = CoefPath + "rtcoef_" + SensorID + ".dat";
  oops::Log::info() << CoefFileName << std::endl;

  oops::Log::trace() << "ObsRadianceRTTOVCPP created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadianceRTTOVCPP::~ObsRadianceRTTOVCPP() {
  oops::Log::trace() << "ObsRadianceRTTOVCPP destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOVCPP::simulateObs(const GeoVaLs & geovals, ioda::ObsVector & hofx,
                                   ObsDiagnostics & d, const QCFlags_t & qc_flags) const {
//
  std::vector<bool>  skip_profile;
  ufo::rttovcpp_interface(geovals, odb_, aRttov_, CoefFileName, channels_,
                          nlevels, skip_profile, doO3_, mwCloudSpecies_, irCloudSpecies_,
                          liquidCloudType_, useModelEffectiveRadius_,
                          thermalSolver_, overlapParam_, hydrotablePath_);
  ufo::rttovcpp_setYdiag(geovals, aRttov_, d, channels_, nlevels,
                         mwCloudSpecies_, irCloudSpecies_, liquidCloudType_);

// ------------------------------------------------------------------------
// Obtain calculated brightness temperature for all profiles/channels
// ------------------------------------------------------------------------
  std::size_t nprofiles = geovals.nlocs();
  std::size_t nchannels = aRttov_.getNchannels();

  ASSERT(geovals.nlocs() == hofx.nlocs());
  hofx.zero();  // this may not be necessary

  const double missing = util::missingValue<double>();

  std::vector <double> bt;

  for (size_t p = 0; p < nprofiles; p++) {
      for (size_t c = 0; c < nchannels; c++) hofx[p*nchannels+c] = missing;
      if (skip_profile[p]) continue;
      bt = aRttov_.getBtRefl(p);
      for (size_t c = 0; c < nchannels; c++) hofx[p*nchannels+c] = bt[c];
  }

  oops::Log::trace() << "ObsRadianceRTTOVCPP::simulateObs done." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOVCPP::print(std::ostream & os) const {
  os << "ObsRadianceRTTOVCPP::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
