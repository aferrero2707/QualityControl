// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \file   TrendCheck.h
/// \author Andrea Ferrero
///

#ifndef QC_MODULE_GLO_MEANVERTEXCHECK_H
#define QC_MODULE_GLO_MEANVERTEXCHECK_H

#include "QualityControl/CheckInterface.h"

#include <optional>
#include <array>
#include <unordered_map>

class TGraph;
class TObject;

using namespace o2::quality_control::core;

namespace o2::quality_control_modules::common
{

/// \brief  Check whether the matching efficiency is within some configurable limits
///
/// \author Andrea Ferrero
class TrendCheck : public o2::quality_control::checker::CheckInterface
{
 public:
  struct ThresholdsParameters
  {
    struct ThresholdValues
    {
      /// \brief min and max threshold values for Bad quality
      std::pair<double, double> mThresholdsBad;
      /// \brief min and max threshold values for Medium quality
      std::optional<std::pair<double, double>> mThresholdsMedium;
    };
    struct ThresholdsElement
    {
      /// \brief min and max threshold values for Bad quality
      std::pair<double, double> mThresholdsBad;
      /// \brief min and max threshold values for Medium quality
      std::optional<std::pair<double, double>> mThresholdsMedium;
      /// \brief nominal interaction rate for which the thresholds are valid
      std::optional<double> mNominalRate;
    };

    std::vector<ThresholdsElement> mThresholds;

    struct ThresholdsElement_
    {
      /// \brief min and max threshold values
      std::pair<double, double> mThresholds;
      /// \brief nominal interaction rate for which the thresholds are valid
      std::optional<double> mNominalRate;
    };

    /// vectors of [min,max] threshold pairs, each with the associated reference interaction rate (optional)
    /// the array index=0 corresponds to the Bad thresholds, while the index=1 corresponds to the Medium thresholds
    std::array<std::vector<ThresholdsElement_>, 2> mThresholds_;

    /// \brief function to retrieve the optimal thresholds (Medium and Bad) for a given interaction rate
    std::array<std::optional<std::pair<double, double>>, 2> getThresholds_(double rate);
    std::array<std::optional<std::pair<double, double>>, 2> getThresholds(double rate);

    void initFromConfiguration_(const CustomParameters& customParameters, const std::string& plotName, const Activity& activity);
    void initFromConfiguration(const CustomParameters& customParameters, const std::string& plotName, const Activity& activity);
  };

  /// Default constructor
  TrendCheck() = default;
  /// Destructor
  ~TrendCheck() override = default;

  void configure() override;
  Quality check(std::map<std::string, std::shared_ptr<MonitorObject>>* moMap) override;
  void beautify(std::shared_ptr<MonitorObject> mo, Quality checkResult = Quality::Null) override;
  std::string getAcceptedType() override;

  void startOfActivity(const Activity& activity) override;
  void endOfActivity(const Activity& activity) override;

  ClassDefOverride(TrendCheck, 1);

 private:
  enum ThresholdsMode
  {
    Fixed,
    Mean,
    StdDeviation
  };
  void initThresholds(std::string plotName);
  std::array<std::optional<std::pair<double, double>>, 2> getThresholds(std::string key, TGraph* graph);
  void getGraphs(TObject* object, std::vector<TGraph*>& graphs);

  Activity mActivity;
  bool mSliceTrend{ false };
  bool mThresholdsAreAbsolute{ true };
  ThresholdsMode mThresholdsMode{ Fixed };
  int mNPointsForAverage{ 0 };
  std::unordered_map<std::string, ThresholdsParameters> mThresholdsParameters;
  //std::unordered_map<std::string, std::vector<std::pair<double, std::array<std::optional<std::pair<double, double>>, 2>>>> mThresholdsTrends;
  std::unordered_map<std::string, std::vector<std::pair<double, std::pair<double, double>>>> mAverageTrend;
  std::unordered_map<std::string, std::vector<std::pair<double, std::pair<double, double>>>> mThresholdsBadTrend;
  std::unordered_map<std::string, std::vector<std::pair<double, std::pair<double, double>>>> mThresholdsMediumTrend;
  std::unordered_map<std::string, Quality> mQualities;
};

} // namespace o2::quality_control_modules::common

#endif // QC_MODULE_GLO_MEANVERTEXCHECK_H
