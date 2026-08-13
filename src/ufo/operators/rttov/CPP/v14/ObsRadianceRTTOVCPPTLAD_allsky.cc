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

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/operators/rttov/CPP/ObsRadianceRTTOVCPPTLAD_allsky.h"
#include "ufo/operators/rttov/CPP/v14/rttovcpp_interface_allsky.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsRadianceRTTOVCPPTLAD> makerRTTOVCPPTL_("RTTOVCPP");
// -----------------------------------------------------------------------------

ObsRadianceRTTOVCPPTLAD::ObsRadianceRTTOVCPPTLAD(const ioda::ObsSpace & odb,
                                           const Parameters_ & parameters)
  : LinearObsOperatorBase(odb), varin_()
{
  // Increment fields to be requested from getvalues and stored in geovals.
  // Temperature is always active; H2O/O3/MWClouds depend on "linear obs operator" below.
  varin_.push_back("air_temperature");

  // Trajectory configuration: must match the nonlinear operator's Absorbers/MWClouds so the
  // Jacobians retrieved in simulateObsTL/AD actually exist on the RTTOV K-matrix.
  const std::vector<std::string> absorbers = parameters.Absorbers;
  if (std::find(absorbers.begin(), absorbers.end(), "H2O") == absorbers.end()) {
    throw eckit::BadParameter("ObsRadianceRTTOVCPPTLAD: Absorbers must include 'H2O'");
  }
  for (const std::string & name : absorbers) {
    if (name == "H2O") continue;
    trajDoO3_ = trajDoO3_ || (name == "O3");
    rttovcppFindSpec(name);  // validates the name
  }

  if (parameters.MWClouds.value() != boost::none) {
    trajMwCloudSpecies_ = *parameters.MWClouds.value();
    for (const std::string & name : trajMwCloudSpecies_) rttovcppFindSpec(name);  // validate
  }

  trajLiquidCloudType_ = parameters.LiquidCloudType.value();
  trajUseModelEffectiveRadius_ = parameters.UseModelEffectiveRadius.value();
  trajThermalSolver_ = parameters.ThermalSolver.value();
  trajOverlapParam_ = parameters.OverlapParam.value();
  if (parameters.IRClouds.value() != boost::none) {
    trajIRCloudSpecies_ = *parameters.IRClouds.value();
    for (const std::string & name : trajIRCloudSpecies_) {
      rttovcppFindIRCloudSpec(name, trajLiquidCloudType_);  // validate
      if (trajUseModelEffectiveRadius_) {
        // See matching check in ObsRadianceRTTOVCPP_allsky.cc: UseModelEffectiveRadius has no
        // effect on "Water" unless LiquidCloudType is "clwde".
        if (name == "Water" && trajLiquidCloudType_ != "clwde") {
          throw eckit::BadParameter("ObsRadianceRTTOVCPPTLAD: UseModelEffectiveRadius has no "
                                     "effect on \"Water\" unless LiquidCloudType is \"clwde\" "
                                     "(RTTOV's 5 discrete OPAC regimes each use a single fixed "
                                     "tabulated diameter) -- got LiquidCloudType '" +
                                     trajLiquidCloudType_ + "'");
        }
        varin_.push_back(rttovcppIRCloudEffectiveRadiusGeovalsName(name));
      }
    }
  }

  if (!trajMwCloudSpecies_.empty() && !trajIRCloudSpecies_.empty()) {
    throw eckit::BadParameter("ObsRadianceRTTOVCPPTLAD: MWClouds and IRClouds are mutually "
                               "exclusive");
  }
  if (!trajMwCloudSpecies_.empty() || !trajIRCloudSpecies_.empty()) {
    if (parameters.HydrotablePath.value() == boost::none) {
      throw eckit::BadParameter("ObsRadianceRTTOVCPPTLAD: HydrotablePath is required when "
                                 "MWClouds or IRClouds is non-empty");
    }
    hydrotablePath_ = *parameters.HydrotablePath.value();
  }

  // "linear obs operator": which of the above Absorbers/MWClouds/IRClouds are active in the
  // TL/AD. Default (block omitted) preserves the operator's previous clear-sky, H2O-only
  // behaviour.
  if (parameters.LinearObsOperator.value() != boost::none) {
    const auto & linParams = *parameters.LinearObsOperator.value();
    linearAbsorbers_ = linParams.Absorbers;
    if (linParams.MWClouds.value() != boost::none) {
      linearMWClouds_ = *linParams.MWClouds.value();
    }
    if (linParams.IRClouds.value() != boost::none) {
      linearIRClouds_ = *linParams.IRClouds.value();
    }
  } else {
    linearAbsorbers_ = {"H2O"};
  }

  // Every linear-active variable must also be part of the trajectory, otherwise its
  // Jacobian was never computed by RTTOV.
  for (const std::string & name : linearAbsorbers_) {
    if (name != "H2O" && std::find(absorbers.begin(), absorbers.end(), name) == absorbers.end()) {
      throw eckit::BadParameter("ObsRadianceRTTOVCPPTLAD: linear obs operator Absorbers '" +
                                 name + "' is not in the (trajectory) Absorbers list");
    }
    varin_.push_back(rttovcppFindSpec(name).geovalsName);
  }
  for (const std::string & name : linearMWClouds_) {
    if (std::find(trajMwCloudSpecies_.begin(), trajMwCloudSpecies_.end(), name) ==
        trajMwCloudSpecies_.end()) {
      throw eckit::BadParameter("ObsRadianceRTTOVCPPTLAD: linear obs operator MWClouds '" +
                                 name + "' is not in the (trajectory) MWClouds list");
    }
    varin_.push_back(rttovcppFindSpec(name).geovalsName);
  }
  for (const std::string & name : linearIRClouds_) {
    if (std::find(trajIRCloudSpecies_.begin(), trajIRCloudSpecies_.end(), name) ==
        trajIRCloudSpecies_.end()) {
      throw eckit::BadParameter("ObsRadianceRTTOVCPPTLAD: linear obs operator IRClouds '" +
                                 name + "' is not in the (trajectory) IRClouds list");
    }
    varin_.push_back(rttovcppFindIRCloudSpec(name, trajLiquidCloudType_).geovalsName);
  }

  // get channels from observations
  const oops::ObsVariables & observed = odb.assimvariables();
  channels_ = observed.channels();  // set private data member channels_

  // get optical depth coef file name from yaml
  const std::string CoefPath = parameters.CoefPath;
  const std::string SensorID = parameters.SensorID;
  CoefFileName = CoefPath + "rtcoef_" + SensorID + ".dat";

  oops::Log::trace() << "ObsRadianceRTTOVCPPTLAD created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadianceRTTOVCPPTLAD::~ObsRadianceRTTOVCPPTLAD() {
  oops::Log::trace() << "ObsRadianceRTTOVCPPTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOVCPPTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics & d,
                                            const QCFlags_t & qc_flags) {
//
  ufo::rttovcpp_interface(geovals, obsspace(), aRttov_, CoefFileName, channels_,
                          nlevels, skip_profile, trajDoO3_, trajMwCloudSpecies_,
                          trajIRCloudSpecies_, trajLiquidCloudType_,
                          trajUseModelEffectiveRadius_,
                          trajThermalSolver_, trajOverlapParam_, hydrotablePath_);
  ufo::rttovcpp_setYdiag(geovals, aRttov_, d, channels_, nlevels,
                         trajMwCloudSpecies_, trajIRCloudSpecies_, trajLiquidCloudType_);

  oops::Log::trace() << "ObsRadianceRTTOVCPPTLAD::setTrajectory done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOVCPPTLAD::simulateObsTL(const GeoVaLs & dx, ioda::ObsVector & dy) const {
//
  std::size_t nprofiles = dy.nlocs();
  std::size_t nchannels = aRttov_.getNchannels();

  std::vector<double>  tmpvar2d(nprofiles, 0.0);  // one single level field
  std::vector<std::vector<double>> dT;  // [nlevels][nprofiles]

  // Retrieve temperature increment in K (always active)
  for (std::size_t i = 0; i < nlevels; ++i) {
      dx.getAtLevel(tmpvar2d, oops::Variable{"air_temperature"}, i);  // get one level T
      dT.push_back(tmpvar2d);  // push one level T into 3D T
  }

  // Retrieve increments for every active Absorber/Cloud, keyed by RTTOV item id
  std::vector<rttov::itemIdType> items;
  std::vector<std::vector<std::vector<double>>> dItems;  // [item][nlevels][nprofiles]
  for (const std::string & name : linearAbsorbers_) {
    const RttovCppItemSpec & spec = rttovcppFindSpec(name);
    std::vector<std::vector<double>> dVar;
    for (std::size_t i = 0; i < nlevels; ++i) {
      dx.getAtLevel(tmpvar2d, oops::Variable{spec.geovalsName}, i);
      dVar.push_back(tmpvar2d);
    }
    items.push_back(spec.item);
    dItems.push_back(dVar);
  }
  for (const std::string & name : linearMWClouds_) {
    const RttovCppItemSpec & spec = rttovcppFindSpec(name);
    std::vector<std::vector<double>> dVar;
    for (std::size_t i = 0; i < nlevels; ++i) {
      dx.getAtLevel(tmpvar2d, oops::Variable{spec.geovalsName}, i);
      dVar.push_back(tmpvar2d);
    }
    items.push_back(spec.item);
    dItems.push_back(dVar);
  }
  for (const std::string & name : linearIRClouds_) {
    const RttovCppItemSpec spec = rttovcppFindIRCloudSpec(name, trajLiquidCloudType_);
    std::vector<std::vector<double>> dVar;
    for (std::size_t i = 0; i < nlevels; ++i) {
      dx.getAtLevel(tmpvar2d, oops::Variable{spec.geovalsName}, i);
      dVar.push_back(tmpvar2d);
    }
    items.push_back(spec.item);
    dItems.push_back(dVar);
  }

//-------------------------------------------
  ASSERT(dx.nlocs() == dy.nlocs());
  ASSERT(nchannels == dy.nvars());
  dy.zero();

  std::vector <double> var_k(nlevels, 0.0);

  for (size_t p = 0; p < nprofiles; p++) {
    if (skip_profile[p]) continue;
    for (size_t c = 0; c < nchannels; c++) {
      var_k = aRttov_.getTK(p, c);              // T Jacobian for a single profile/channel
      for (size_t l = 0; l < nlevels; l++)
          dy[p*nchannels+c] += var_k[l]*dT[l][p];

      for (size_t j = 0; j < items.size(); j++) {
        var_k = aRttov_.getItemK(items[j], p, c);
        for (size_t l = 0; l < nlevels; l++)
            dy[p*nchannels+c] += var_k[l]*dItems[j][l][p];
      }
    }
  }

  oops::Log::trace() << "ObsRadianceRTTOVCPPTLAD::simulateObsTL done" << std::endl;
}

// -----------------------------------------------------------------------------

  void ObsRadianceRTTOVCPPTLAD::simulateObsAD(GeoVaLs & dx, const ioda::ObsVector & dy) const {
  std::size_t nprofiles = dy.nlocs();
  std::size_t nchannels = aRttov_.getNchannels();

  std::vector<double>  tmpvar2d(nprofiles, 0.0);  // one single level field
  std::vector<std::vector<double>> dT;  // [nlevels][nprofiles]

  // Adjoint accumulation must add onto whatever dx already holds (other obs operators may
  // contribute to the same Increment), so seed each accumulator from dx's current content.
  for (std::size_t i = 0; i < nlevels; ++i) {
      dx.getAtLevel(tmpvar2d, oops::Variable{"air_temperature"}, i);
      dT.push_back(tmpvar2d);
  }

  // Accumulators for every active Absorber/Cloud, keyed by RTTOV item id
  std::vector<std::string> geovalsNames;
  std::vector<rttov::itemIdType> items;
  std::vector<std::vector<std::vector<double>>> dItems;  // [item][nlevels][nprofiles]
  for (const std::string & name : linearAbsorbers_) {
    const RttovCppItemSpec & spec = rttovcppFindSpec(name);
    geovalsNames.push_back(spec.geovalsName);
    items.push_back(spec.item);
    std::vector<std::vector<double>> dVar;
    for (std::size_t i = 0; i < nlevels; ++i) {
      dx.getAtLevel(tmpvar2d, oops::Variable{spec.geovalsName}, i);
      dVar.push_back(tmpvar2d);
    }
    dItems.push_back(dVar);
  }
  for (const std::string & name : linearMWClouds_) {
    const RttovCppItemSpec & spec = rttovcppFindSpec(name);
    geovalsNames.push_back(spec.geovalsName);
    items.push_back(spec.item);
    std::vector<std::vector<double>> dVar;
    for (std::size_t i = 0; i < nlevels; ++i) {
      dx.getAtLevel(tmpvar2d, oops::Variable{spec.geovalsName}, i);
      dVar.push_back(tmpvar2d);
    }
    dItems.push_back(dVar);
  }
  for (const std::string & name : linearIRClouds_) {
    const RttovCppItemSpec spec = rttovcppFindIRCloudSpec(name, trajLiquidCloudType_);
    geovalsNames.push_back(spec.geovalsName);
    items.push_back(spec.item);
    std::vector<std::vector<double>> dVar;
    for (std::size_t i = 0; i < nlevels; ++i) {
      dx.getAtLevel(tmpvar2d, oops::Variable{spec.geovalsName}, i);
      dVar.push_back(tmpvar2d);
    }
    dItems.push_back(dVar);
  }

//-------------------------------------------
  ASSERT(dx.nlocs() == dy.nlocs());

  const double missing = util::missingValue<double>();

  std::vector <double> var_k(nlevels, 0.0);

  for (size_t p = 0; p < nprofiles; p++) {
    if (skip_profile[p]) continue;
    for (size_t c = 0; c < nchannels; c++) {
      if (dy[p*nchannels+c] == missing) continue;

      var_k = aRttov_.getTK(p, c);               // T Jacobian, nlevels
      for (size_t l = 0; l < nlevels; ++l) {
          dT[l][p] += dy[p*nchannels+c] * var_k[l];
      }

      for (size_t j = 0; j < items.size(); j++) {
        var_k = aRttov_.getItemK(items[j], p, c);
        for (size_t l = 0; l < nlevels; l++) {
            dItems[j][l][p] += dy[p*nchannels+c] * var_k[l];
        }
      }
    }
  }

  // Put temperature increment in K
  for (std::size_t l = 0; l < nlevels; ++l) {
      for (size_t p = 0; p < nprofiles; p++) {
          tmpvar2d[p] = dT[l][p];
      }
      dx.putAtLevel(tmpvar2d, oops::Variable{"air_temperature"}, l);  // put one level T
  }

  // Put back increments for every active Absorber/Cloud
  for (size_t j = 0; j < items.size(); j++) {
    for (std::size_t l = 0; l < nlevels; ++l) {
        for (size_t p = 0; p < nprofiles; p++) {
            tmpvar2d[p] = dItems[j][l][p];
        }
        dx.putAtLevel(tmpvar2d, oops::Variable{geovalsNames[j]}, l);
    }
  }

  oops::Log::trace() << "ObsRadianceRTTOVCPPTLAD::simulateObsAD done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceRTTOVCPPTLAD::print(std::ostream & os) const {
  os << "ObsRadianceRTTOVCPPTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
