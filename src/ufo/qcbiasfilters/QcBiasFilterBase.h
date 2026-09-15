/*
 * A new type of QC for bias coefficients
 * the idea is to block the gradient of the observations whose qcbiasflag is non-zero,
 * while their gradient to model state remains the same 
 * This mimics the "predictorbase"
 * Chih-Chi Hu, Sep 2026
 */

#ifndef UFO_QCBIASFILTERS_QCBIASFILTERBASE_H_
#define UFO_QCBIASFILTERS_QCBIASFILTERBASE_H_

#include <map>
#include <memory>
#include <string>

#include <boost/noncopyable.hpp>

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/AssociativeContainers.h"


namespace ioda {
class ObsSpace;
class ObsVector;
template <typename DATATYPE> class ObsDataVector;
}

namespace ufo {
class GeoVaLs;
class ObsDiagnostics;

/// \brief Base class for the parameters of a single qc-bias filter.
class QcBiasFilterParametersBase : public oops::Parameters {
  OOPS_ABSTRACT_PARAMETERS(QcBiasFilterParametersBase, Parameters)
 public:
  /// Filter name, used to look up the filter implementation in QcBiasFilterFactory.
  oops::RequiredParameter<std::string> name{"name", this};
  /// Optional subset of channels (e.g. "9-14") this filter entry applies to. If omitted,
  /// applies to every channel simulated by this obs space.
  oops::OptionalParameter<std::string> channels{"channels", this};
};

// -----------------------------------------------------------------------------

/// \brief Base class for filters deciding whether an observation should contribute
/// to the update of the VarBC bias coefficients (independent of whether it
/// contributes to Jo / the state increment).
///
/// The implementer should define a subclass of QcBiasFilterParametersBase, typedef
/// it to `Parameters_`, and provide a constructor with the signature
///     SubclassName(const Parameters_ &, const oops::ObsVariables &);
class QcBiasFilterBase : private boost::noncopyable {
 public:
  QcBiasFilterBase(const QcBiasFilterParametersBase &, const oops::ObsVariables &);
  virtual ~QcBiasFilterBase() = default;

  /// Compute this filter's own qc_bias flag (0 = keep, nonzero = exclude from bias update)
  /// for every location/channel in \p qc_bias. \p qc_bias is assumed to already be
  /// correctly sized and does not need to be zeroed by the implementation beforehand
  /// (ObsBias::computeQcBias() takes care of that).
  virtual void compute(const ioda::ObsSpace &,
                       const GeoVaLs &,
                       const ObsDiagnostics &,
                       const ioda::ObsVector & hofx,
                       ioda::ObsDataVector<int> & qc_bias) const = 0;

  /// hdiags names required to compute this filter (e.g. clear-sky brightness temperature).
  const oops::ObsVariables & requiredHdiagnostics() const {return hdiags_;}

  const std::string & name() const {return func_name_;}

 protected:
  oops::ObsVariables vars_;      ///< full variable/channel list of this obs space
  oops::ObsVariables hdiags_;    ///< required ObsDiagnostics -- filled in by subclass ctor

  /// Channel numbers this filter entry actually applies to (parsed from 'channels', or
  /// = vars_.channels() if 'channels' was not set). Per-channel Parameters (e.g. a filter's
  /// own 'threshold') should be sized against this, not against vars_.channels().
  std::vector<int> selectedChannels_;
  /// For each entry of selectedChannels_, its position within vars_.channels(). hofx and
  /// qc_bias in compute() are always indexed by the FULL channel list, so subclasses must use
  /// channelIndices_[jsel] (not jsel) whenever indexing into them.
  std::vector<std::size_t> channelIndices_;

 private:
  std::string func_name_;
};

typedef std::vector<std::shared_ptr<QcBiasFilterBase>> QcBiasFilters;

// -----------------------------------------------------------------------------

/// \brief Factory for qc-bias filters. Mirrors ufo::PredictorFactory.
class QcBiasFilterFactory {
 public:
  static std::unique_ptr<QcBiasFilterBase> create(const QcBiasFilterParametersBase &,
                                                  const oops::ObsVariables &);
  static std::unique_ptr<QcBiasFilterParametersBase> createParameters(const std::string &);
  static bool filterExists(const std::string &);
  static std::vector<std::string> getMakerNames() {         
    return oops::keys(getMakers());
  }
  virtual ~QcBiasFilterFactory() = default;

 protected:
  explicit QcBiasFilterFactory(const std::string &);

 private:
  virtual std::unique_ptr<QcBiasFilterBase> make(const QcBiasFilterParametersBase &,
                                                 const oops::ObsVariables &) = 0;
  virtual std::unique_ptr<QcBiasFilterParametersBase> makeParameters() const = 0;

  static std::map<std::string, QcBiasFilterFactory *> & getMakers() {
    static std::map<std::string, QcBiasFilterFactory *> makers_;
    return makers_;
  }
};

template <class T>
class QcBiasFilterMaker : public QcBiasFilterFactory {
  typedef typename T::Parameters_ Parameters_;

  std::unique_ptr<QcBiasFilterBase> make(const QcBiasFilterParametersBase & parameters,
                                         const oops::ObsVariables & vars) override {
    const auto & stronglyTypedParameters = dynamic_cast<const Parameters_ &>(parameters);
    return std::make_unique<T>(stronglyTypedParameters, vars);
  }
  std::unique_ptr<QcBiasFilterParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit QcBiasFilterMaker(const std::string & name) : QcBiasFilterFactory(name) {}
};

}  // namespace ufo

#endif  // UFO_QCBIASFILTERS_QCBIASFILTERBASE_H_
