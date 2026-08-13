/*
 * (C) Copyright 2017-2021 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_RTTOV_CPP_OBSRADIANCERTTOVCPP_ALLSKY_H_
#define UFO_OPERATORS_RTTOV_CPP_OBSRADIANCERTTOVCPP_ALLSKY_H_

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/rttov/CPP/ObsRadianceRTTOVCPPParameters_allsky.h"

#include "rttov/wrapper/RttovSafe.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// RadianceRTTOVCPP observation operator class
class ObsRadianceRTTOVCPP : public ObsOperatorBase,
                   private util::ObjectCounter<ObsRadianceRTTOVCPP> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the ObsOperatorFactory.
  typedef ObsRadianceRTTOVCPPParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsRadianceRTTOVCPP";}

  ObsRadianceRTTOVCPP(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsRadianceRTTOVCPP();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

// Other: declare variable function with return type of oops:Variables
  const oops::Variables & requiredVars() const override {return varin_;}

 private:
  void print(std::ostream &) const override;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
  std::string        CoefFileName;
  std::vector<int>   channels_;
  mutable std::size_t        nlevels;  // need this in order to allocate dx

// All-sky / trace gas configuration (see ObsRadianceRTTOVCPPParameters)
  bool                     doO3_ = false;
  std::vector<std::string> mwCloudSpecies_;  // MW hydrometeors, e.g. {"Water","Rain"}
  std::vector<std::string> irCloudSpecies_;  // IR/VIS cloud species, subset of {"Water","Ice"}
  std::string              liquidCloudType_;  // OPAC regime for "Water" in irCloudSpecies_
  bool                     useModelEffectiveRadius_ = false;
  std::string              thermalSolver_ = "delta_edd";
  std::string              overlapParam_ = "auto_select";
  std::string              hydrotablePath_;

// Declare a RttovSafe object for one single sensor
  mutable rttov::RttovSafe  aRttov_ = rttov::RttovSafe();
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_RTTOV_CPP_OBSRADIANCERTTOVCPP_ALLSKY_H_
