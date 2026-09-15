/*
 * Chih-Chi Hu, Sep 2026
 */

#include "ufo/qcbiasfilters/QcBiasFilterBase.h"

#include <algorithm>
#include <map>
#include <set>

#include "eckit/exception/Exceptions.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"

namespace ufo {

QcBiasFilterBase::QcBiasFilterBase(const QcBiasFilterParametersBase & parameters,
                                   const oops::ObsVariables & vars)
  : vars_(vars), hdiags_(), func_name_(parameters.name) {
  const std::vector<int> & allChannels = vars_.channels();

  if (parameters.channels.value() == boost::none || allChannels.empty()) {
    // No restriction requested (or this obs space has no channels at all): apply to everything.
    selectedChannels_ = allChannels;
    channelIndices_.resize(allChannels.size());
    for (std::size_t j = 0; j < channelIndices_.size(); ++j) channelIndices_[j] = j;
  } else {
    const std::set<int> selected = oops::parseIntSet(*parameters.channels.value());
    selectedChannels_.assign(selected.begin(), selected.end());   // ascending order
    channelIndices_.clear();
    for (const int ch : selectedChannels_) {
      const auto it = std::find(allChannels.begin(), allChannels.end(), ch);
      if (it == allChannels.end()) {
        throw eckit::BadValue("qc bias filter '" + func_name_ + "': channel " +
                              std::to_string(ch) + " listed in 'channels' is not simulated "
                              "by this obs space", Here());
      }
      channelIndices_.push_back(static_cast<std::size_t>(it - allChannels.begin()));
    }
  }
}

QcBiasFilterFactory::QcBiasFilterFactory(const std::string & name) {
  if (filterExists(name)) {
    oops::Log::error() << name << " already registered in ufo::QcBiasFilterFactory."
                       << std::endl;
    ABORT("Element already registered in ufo::QcBiasFilterFactory.");
  }
  getMakers()[name] = this;
}

std::unique_ptr<QcBiasFilterBase> QcBiasFilterFactory::create(
    const QcBiasFilterParametersBase & parameters, const oops::ObsVariables & vars) {
  const std::string name = parameters.name;
  if (!filterExists(name)) {
    oops::Log::error() << name << " does not exist in ufo::QcBiasFilterFactory." << std::endl;
    ABORT("Element does not exist in ufo::QcBiasFilterFactory.");
  }
  return getMakers().find(name)->second->make(parameters, vars);
}

std::unique_ptr<QcBiasFilterParametersBase> QcBiasFilterFactory::createParameters(
    const std::string & name) {
  auto it = getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in ufo::QcBiasFilterFactory");
  }
  return it->second->makeParameters();
}

bool QcBiasFilterFactory::filterExists(const std::string & name) {
  return (getMakers().find(name) != getMakers().end());
}

}  // namespace ufo
