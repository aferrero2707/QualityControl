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
/// \file   TrendCheck.cxx
/// \author Andrea Ferrero
///

#include "Common/TrendCheck.h"
#include "QualityControl/MonitorObject.h"
#include "QualityControl/Quality.h"
#include "QualityControl/QcInfoLogger.h"
#include <CommonUtils/StringUtils.h>
#include <TMath.h>
#include <TGraph.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TPolyLine.h>

namespace o2::quality_control_modules::common
{
static std::string getCustomParameter(const CustomParameters& customParameters, const std::string& key, const Activity& activity)
{
  auto parOpt = customParameters.atOptional(key, activity);
  if (parOpt.has_value()) {
    return parOpt.value();
  } else {
    // search in standard parameters if key not found in extended ones
    parOpt = customParameters.atOptional(key);
    if (parOpt.has_value()) {
      return parOpt.value();
    }
  }
  return std::string();
}

void TrendCheck::ThresholdsParameters::initFromConfiguration(const CustomParameters& customParameters, const std::string& plotName, const Activity& activity)
{
  std::array<std::string, 2> qualityLabels{ "Bad", "Medium" };

  for (size_t index = 0; index < 2; index++) {
    std::cout << "Initializing " << qualityLabels[index] << " thresholds for \"" << plotName << "\"" << std::endl;

    // get configuration parameter associated to the key
    std::string parKey = std::string("thresholds") + qualityLabels[index] + ":" + plotName;
    std::string parValue;
    // search in extended parameters first
    auto parOpt = customParameters.atOptional(parKey, activity);
    if (!parOpt) {
      // search in standard parameters if key not found in extended ones
      parOpt = customParameters.atOptional(parKey);
    }
    if (parOpt.has_value()) {
      parValue = parOpt.value();
    }

    // extract information for each nominal interaction rate value
    auto interactionRatePoints = o2::utils::Str::tokenize(parValue, ';', false, true);
    std::cout << "  interactionRatePoints.size(): " << interactionRatePoints.size() << std::endl;
    for (const auto& interactionRatePoint : interactionRatePoints) {
      ThresholdsElement_ element;
      // extract interaction rate and threshold pairs
      auto rateAndThresholds = o2::utils::Str::tokenize(interactionRatePoint, ':', false, true);
      std::cout << "    rateAndThresholds.size(): " << rateAndThresholds.size() << std::endl;
      std::cout << "    rateAndThresholds[0]: \"" << rateAndThresholds[0] << "\"" << std::endl;
      std::string thresholdValues;
      if (rateAndThresholds.size() == 2) {
        thresholdValues = rateAndThresholds[1];
        element.mNominalRate = std::stod(rateAndThresholds[0]) * 1000;
      } else {
        thresholdValues = rateAndThresholds[0];
      }

      std::vector<std::string> thresholdMinMax = o2::utils::Str::tokenize(thresholdValues, ',', false, true);
      if (thresholdMinMax.size() == 2) {
        try {
          element.mThresholds.first = std::stod(thresholdMinMax[0]);
          element.mThresholds.second = std::stod(thresholdMinMax[1]);
          //ILOG(Info, Support) << "TrendCheck: range for " << key << " set to [" << min << "," << max << "]" << ENDM;
        } catch (std::exception& exception) {
          ILOG(Error, Support) << "Cannot convert values from string to double for " << qualityLabels[index] << " thresholds of plot \"" << plotName << "\", string is " << interactionRatePoint << ENDM;
          element.mThresholds.first = 0;
          element.mThresholds.second = 0;
        }
      }

      mThresholds_[index].emplace_back(element);
    }
  }

  std::cout << "Thresholds for \"" << plotName << "\"" << std::endl;
  for (size_t index = 0; index < 2; index++) {
    std::cout << "  " << qualityLabels[index] << " thresholds:" << std::endl;
    for (auto & element : mThresholds_[index]) {
      if (element.mNominalRate.has_value()) {
        std::cout << "  nominal rate: " << element.mNominalRate.value() << " Hz" << std::endl;
      }
      std::cout << "    min:" << element.mThresholds.first << std::endl;
      std::cout << "    max:" << element.mThresholds.second << std::endl;
    }
  }
}

void TrendCheck::ThresholdsParameters::initFromConfiguration_(const CustomParameters& customParameters, const std::string& plotName, const Activity& activity)
{
  // get configuration parameter associated to the key
  std::string parKey = std::string("thresholds:") + plotName;
  std::string parValue;
  // search in extended parameters first
  auto parOpt = customParameters.atOptional(parKey, activity);
  if (parOpt.has_value()) {
    parValue = parOpt.value();
  } else {
    // search in standard parameters if key not found in extended ones
    parOpt = customParameters.atOptional(parKey);
    if (parOpt.has_value()) {
      parValue = parOpt.value();
    }
  }

  // extract information for each nominal interaction rate value
  auto interactionRatePoints = o2::utils::Str::tokenize(parValue, '|', false, true);
  std::cout << "  interactionRatePoints.size(): " << interactionRatePoints.size() << std::endl;
  for (const auto& interactionRatePoint : interactionRatePoints) {
    ThresholdsElement element;
    // extract interaction rate and threshold pairs
    auto rateAndThresholds = o2::utils::Str::tokenize(interactionRatePoint, ':', false, true);
    std::cout << "    rateAndThresholds.size(): " << rateAndThresholds.size() << std::endl;
    std::cout << "    rateAndThresholds[0]: \"" << rateAndThresholds[0] << "\"" << std::endl;
    std::string thresholds;
    if (rateAndThresholds.size() == 2) {
      thresholds = rateAndThresholds[1];
      element.mNominalRate = std::stod(rateAndThresholds[0]) * 1000;
    } else {
      thresholds = rateAndThresholds[0];
    }

    auto thresholdBadMedium = o2::utils::Str::tokenize(thresholds, ';', false, true);
    if (thresholdBadMedium.empty()) {
      continue;
    }

    std::vector<std::string> thresholdBad = o2::utils::Str::tokenize(thresholdBadMedium[0], ',', false, true);
    if (thresholdBad.size() == 2) {
      try {
        element.mThresholdsBad.first = std::stod(thresholdBad[0]);
        element.mThresholdsBad.second = std::stod(thresholdBad[1]);
        //ILOG(Info, Support) << "TrendCheck: range for " << key << " set to [" << min << "," << max << "]" << ENDM;
      } catch (std::exception& exception) {
        ILOG(Error, Support) << "Cannot convert values from string to double for Bad thresholds of plot \"" << plotName << "\", string is " << thresholdBadMedium[0] << ENDM;
        element.mThresholdsBad.first = 0;
        element.mThresholdsBad.second = 0;
      }
    }

    if (thresholdBadMedium.size() == 2) {
      auto thresholdMedium = o2::utils::Str::tokenize(thresholdBadMedium[1], ',', false, true);
      if (thresholdMedium.size() == 2) {
        try {
          element.mThresholdsMedium = std::make_pair<double, double>(std::stod(thresholdMedium[0]), std::stod(thresholdMedium[1]));
        } catch (std::exception& exception) {
          ILOG(Error, Support) << "Cannot convert values from string to double for Medium thresholds of plot \"" << plotName << "\", string is " << thresholdBadMedium[1] << ENDM;
        }
      }
    }

    mThresholds.emplace_back(element);
  }

  std::cout << "Thresholds for \"" << plotName << "\"" << std::endl;
  for (auto & element : mThresholds) {
    if (element.mNominalRate.has_value()) {
      std::cout << "  nominal rate: " << element.mNominalRate.value() << " Hz" << std::endl;
    }
    std::cout << "  Bad thresholds:" << std::endl;
    std::cout << "    min:" << element.mThresholdsBad.first << std::endl;
    std::cout << "    max:" << element.mThresholdsBad.second << std::endl;
    if (element.mThresholdsMedium.has_value()) {
      std::cout << "  Medium thresholds:" << std::endl;
      std::cout << "    min:" << element.mThresholdsMedium->first << std::endl;
      std::cout << "    max:" << element.mThresholdsMedium->second << std::endl;
    }
  }
}

static std::pair<double, double> interpolateThresholds(double fraction, const std::pair<double, double>& thresholdsLow, const std::pair<double, double>& thresholdsHigh)
{
  double thresholdMin = thresholdsLow.first * (1.0 - fraction) + thresholdsHigh.first * fraction;
  double thresholdMax = thresholdsLow.second * (1.0 - fraction) + thresholdsHigh.second * fraction;
  return std::make_pair(thresholdMin, thresholdMax);
}

std::array<std::optional<std::pair<double, double>>, 2> TrendCheck::ThresholdsParameters::getThresholds_(double rate)
{
  std::array<std::optional<std::pair<double, double>>, 2> result;

  int indexLow = -1;
  int indexHigh = -1;

  double rateLow = -1;
  double rateHigh = -1;

  std::cout << "    getThresholds(): getting thresholds for " << rate << " Hz" << std::endl;
  std::cout << "    getThresholds(): mThresholds.size() = " << mThresholds.size() << std::endl;

  for (int index = 0; index < mThresholds.size(); index++) {
    const auto& element = mThresholds[index];
    if (!element.mNominalRate) {
      continue;
    }

    double rateCurrent = element.mNominalRate.value();

    if (rateCurrent <= rate) {
      if (indexLow < 0 || rateLow < rateCurrent) {
        indexLow = index;
        rateLow = rateCurrent;
      }
    }

    if (rateCurrent >= rate) {
      if (indexHigh < 0 || rateHigh > rateCurrent) {
        indexHigh = index;
        rateHigh = rateCurrent;
      }
    }
  }

  std::cout << "    getThresholds(): indexLow=" << indexLow << "  indexHigh=" << indexHigh << std::endl;
  std::cout << "    getThresholds(): rateLow=" << rateLow << "  rateHigh=" << rateHigh << std::endl;

  if (indexLow >= 0 && indexHigh >= 0 && indexLow != indexHigh) {
    // interpolation
    double rateLow = mThresholds[indexLow].mNominalRate.value();
    double rateHigh = mThresholds[indexHigh].mNominalRate.value();

    double fraction = (rate - rateLow) / (rateHigh - rateLow);

    result[0] = interpolateThresholds(fraction, mThresholds[indexLow].mThresholdsBad, mThresholds[indexHigh].mThresholdsBad);
    std::cout << "    getThresholds(): interpolated values for " << rate << " Hz: " << result[0].value().first << ", " << result[0].value().second << std::endl;

    if (mThresholds[indexLow].mThresholdsMedium.has_value() && mThresholds[indexHigh].mThresholdsMedium.has_value()) {
      result[1] = interpolateThresholds(fraction, mThresholds[indexLow].mThresholdsMedium.value(), mThresholds[indexHigh].mThresholdsMedium.value());
    }
  } else if (indexLow >= 0) {
    result[0] = mThresholds[indexLow].mThresholdsBad;
    result[1] = mThresholds[indexLow].mThresholdsMedium;
  } else if (indexHigh >= 0) {
    result[0] = mThresholds[indexHigh].mThresholdsBad;
    result[1] = mThresholds[indexHigh].mThresholdsMedium;
  } else if (mThresholds.size() > 0) {
    result[0] = mThresholds[0].mThresholdsBad;
    result[1] = mThresholds[0].mThresholdsMedium;
  }

  return result;
}

std::array<std::optional<std::pair<double, double>>, 2> TrendCheck::ThresholdsParameters::getThresholds(double rate)
{
  std::array<std::optional<std::pair<double, double>>, 2> result;

  for (size_t qualityIndex = 0; qualityIndex < 2; qualityIndex++) {

    int indexLow = -1;
    int indexHigh = -1;

    double rateLow = -1;
    double rateHigh = -1;

    std::cout << "    getThresholds(): getting thresholds for " << rate << " Hz" << std::endl;
    std::cout << "    getThresholds(): mThresholds.size() = " << mThresholds_[qualityIndex].size() << std::endl;

    for (int index = 0; index < mThresholds_[qualityIndex].size(); index++) {
      const auto& element = mThresholds_[qualityIndex][index];
      if (!element.mNominalRate) {
        continue;
      }

      double rateCurrent = element.mNominalRate.value();

      if (rateCurrent <= rate) {
        if (indexLow < 0 || rateLow < rateCurrent) {
          indexLow = index;
          rateLow = rateCurrent;
        }
      }

      if (rateCurrent >= rate) {
        if (indexHigh < 0 || rateHigh > rateCurrent) {
          indexHigh = index;
          rateHigh = rateCurrent;
        }
      }
    }

    std::cout << "    getThresholds(): indexLow=" << indexLow << "  indexHigh=" << indexHigh << std::endl;
    std::cout << "    getThresholds(): rateLow=" << rateLow << "  rateHigh=" << rateHigh << std::endl;

    if (indexLow >= 0 && indexHigh >= 0 && indexLow != indexHigh) {
      // interpolation
      std::cout << "    getThresholds(): interpolating values for " << rate << " Hz" << std::endl;
      double rateLow = mThresholds_[qualityIndex][indexLow].mNominalRate.value();
      double rateHigh = mThresholds_[qualityIndex][indexHigh].mNominalRate.value();
      std::cout << "    getThresholds(): rates min/max " << rateLow << " / " << rateHigh << std::endl;

      double fraction = (rate - rateLow) / (rateHigh - rateLow);

      result[qualityIndex] = interpolateThresholds(fraction, mThresholds_[qualityIndex][indexLow].mThresholds, mThresholds_[qualityIndex][indexHigh].mThresholds);
      std::cout << "    getThresholds(): interpolated values for " << rate << " Hz: " << result[0].value().first << ", " << result[0].value().second << std::endl;
    } else if (indexLow >= 0) {
      result[qualityIndex] = mThresholds_[qualityIndex][indexLow].mThresholds;
    } else if (indexHigh >= 0) {
      result[qualityIndex] = mThresholds_[qualityIndex][indexHigh].mThresholds;
    } else if (mThresholds.size() > 0) {
      result[qualityIndex] = mThresholds_[qualityIndex][0].mThresholds;
    }
  }

  return result;
}

void TrendCheck::configure() {}

void TrendCheck::startOfActivity(const Activity& activity)
{
  mActivity = activity;

  std::string parKey = "thresholdsMode";
  // search in extended parameters first
  auto parOpt = mCustomParameters.atOptional(parKey, mActivity);
  if (!parOpt) {
    // search in standard parameters if key not found in extended ones
    parOpt = mCustomParameters.atOptional(parKey);
  }
  if (parOpt.has_value()) {
    if (parOpt.value() == "Mean") {
      mThresholdsMode = Mean;
    } else if (parOpt.value() == "StdDeviation") {
      mThresholdsMode = StdDeviation;
    } else if (parOpt.value() != "Fixed") {
      ILOG(Warning, Support) << "unrecognized threshold mode \"" << parOpt.value() << "\", using default \"Fixed\" mode" << ENDM;
    }
  }
  switch (mThresholdsMode) {
  case Fixed:
    ILOG(Info, Support) << "thresholds mode set to \"Fixed\"" << ENDM;
    break;
  case Mean:
    ILOG(Info, Support) << "thresholds mode set to \"Mean\"" << ENDM;
    break;
  case StdDeviation:
    ILOG(Info, Support) << "thresholds mode set to \"StdDeviation\"" << ENDM;
    break;
  }

  parKey = "nPointsForAverage";
  // search in extended parameters first
  parOpt = mCustomParameters.atOptional(parKey, mActivity);
  if (!parOpt) {
    // search in standard parameters if key not found in extended ones
    parOpt = mCustomParameters.atOptional(parKey);
  }
  if (parOpt.has_value()) {
    if (parOpt.value() == "all") {
      mNPointsForAverage = 0;
    } else {
      mNPointsForAverage = std::stoi(parOpt.value());
    }
  }

  if (mNPointsForAverage == 0) {
    ILOG(Info, Support) << "using all points for statistics calculation" << ENDM;
  } else {
    ILOG(Info, Support) << "using at most " << mNPointsForAverage << " points for statistics calculation" << ENDM;
  }


}

void TrendCheck::endOfActivity(const Activity& activity)
{
  mActivity = Activity{};
  mThresholdsParameters.clear();
  mAverageTrend.clear();
  mThresholdsBadTrend.clear();
  mThresholdsMediumTrend.clear();
  mQualities.clear();
}

static std::string getBaseName(std::string name)
{
  auto pos = name.rfind("/");
  return ((pos < std::string::npos) ? name.substr(pos + 1) : name);
}

void TrendCheck::initThresholds(std::string plotName)
{
  // Get acceptable range for this histogram
  if (mThresholdsParameters.count(plotName) > 0) {
    return;
  }

  mThresholdsParameters.emplace(std::unordered_map<std::string, ThresholdsParameters>::value_type (plotName, ThresholdsParameters()));
  mThresholdsParameters[plotName].initFromConfiguration(mCustomParameters, plotName, mActivity);
}

static std::optional<std::pair<double, double>> getGraphStatistics(TGraph* graph, int nPointsForAverage)
{
  std::optional<std::pair<double, double>> result;

  int nPoints = graph->GetN();
  std::cout << "      getGraphAverage(): nPoints = " << nPoints << std::endl;
  std::cout << "      getGraphAverage(): nPointsForAverage = " << nPointsForAverage << std::endl;
  if (nPoints < 2) {
    return result;
  }

  int pointIndexMax = nPoints - 2;
  int pointIndexMin = (nPointsForAverage > 0 && pointIndexMax >= nPointsForAverage) ? pointIndexMax - nPointsForAverage + 1 : 0;
  int N = pointIndexMax - pointIndexMin + 1;
  std::cout << "      getGraphAverage(): indexes min = " << pointIndexMin << "  max = " << pointIndexMax << std::endl;

  if (N < 2) {
    return result;
  }

  double mean = 0;
  for (int pointIndex = pointIndexMin; pointIndex <= pointIndexMax; pointIndex++) {
    mean += graph->GetPointY(pointIndex);
  }
  std::cout << "      getGraphAverage(): integral = " << mean << std::endl;

  mean /= N;

  double stdDev = 0;
  for (int pointIndex = pointIndexMin; pointIndex <= pointIndexMax; pointIndex++) {
    double delta = graph->GetPointY(pointIndex) - mean;
    stdDev += delta * delta;
  }
  stdDev /= (N - 1) * N;
  stdDev = std::sqrt(stdDev);

  result = std::make_pair(mean, stdDev);

  return result;
}

std::array<std::optional<std::pair<double, double>>, 2> TrendCheck::getThresholds(std::string plotName, TGraph* graph)
{
  double rateScale = 1;
  for (int i = 0; i < graph->GetN() - 1; i++) {
    rateScale *= 0.9;
  }
  double rate = 1000000 * rateScale;
  auto result = mThresholdsParameters[plotName].getThresholds(rate);

  if (mThresholdsMode != Fixed) {
    // the thresholds retrieved from the configuration are relative to the trend graphStatistics, we need to convert them into absolute values
    auto graphStatistics = getGraphStatistics(graph, mNPointsForAverage);
    if (!graphStatistics) {
      result[0].reset();
      result[1].reset();
      return result;
    }

    if (mThresholdsMode == Mean) {
      double mean = graphStatistics.value().first;
      // the thresholds retrieved from the configuration are relative to the mean value of the last N points, we need to convert them into absolute values
      if (result[0].has_value()) {
        result[0].value().first = mean + result[0].value().first * std::fabs(mean);
        result[0].value().second = mean + result[0].value().second * std::fabs(mean);
      }

      if (result[1].has_value()) {
        result[1].value().first = mean + result[1].value().first * std::fabs(mean);
        result[1].value().second = mean + result[1].value().second * std::fabs(mean);
      }
    } else if (mThresholdsMode == StdDeviation) {
      // the thresholds retrieved from the configuration are expressed as number of sigmas from the mean value of the last N points, we need to convert them into absolute values
      double mean = graphStatistics.value().first;
      double stdDevOfMean = graphStatistics.value().second;
      double lastPointValue = graph->GetPointY(graph->GetN() - 1);
      double lastPointError = graph->GetErrorY(graph->GetN() - 1);
      if (lastPointError < 0) {
        lastPointError = 0;
      }
      const double totalError = sqrt(stdDevOfMean * stdDevOfMean + lastPointError * lastPointError);

      result[0].value().first = mean + result[0].value().first * totalError;
      result[0].value().second = mean + result[0].value().second * totalError;

      if (result[1].has_value()) {
        result[1].value().first = mean + result[1].value().first * totalError;
        result[1].value().second = mean + result[1].value().second * totalError;
      }
    }
  }

  return result;
}

void TrendCheck::getGraphs(TObject* object, std::vector<TGraph*>& graphs)
{
  TCanvas* canvas = dynamic_cast<TCanvas*>(object);
  if (canvas) {
    if (mSliceTrend) {
      TList* padList = (TList*)canvas->GetListOfPrimitives();
      padList->SetOwner(kTRUE);
      const int numberPads = padList->GetEntries();
      for (int iPad = 0; iPad < numberPads; iPad++) {
        auto pad = static_cast<TPad*>(padList->At(iPad));
        graphs.push_back(static_cast<TGraph*>(pad->GetPrimitive("Graph")));
      }
    } else {
      // If we have a standard trending with a TGraphErrors, then there will be a TGraph and a TGraphErrors with the same name ("Graph")
      // with the current configuration the TGraphErrors will always be added after the TGraph
      // loop from the back and find the last element with the name "Graph"
      TGraph* g;
      int jList = canvas->GetListOfPrimitives()->LastIndex();
      for (; jList > 0; jList--) {
        if (!strcmp(canvas->GetListOfPrimitives()->At(jList)->GetName(), "Graph")) {
          g = (TGraph*)canvas->GetListOfPrimitives()->At(jList);
          break;
        }
      }
      if (!g) {
        // if the upper loop somehow fails, log a warning and get the object the old fashioned way
        ILOG(Warning, Support) << "No TGraph found in the List of primitives." << ENDM;
        g = (TGraph*)canvas->GetListOfPrimitives()->FindObject("Graph");
      }
      graphs.push_back(g);
    }
  } else {
    TGraph* graph = dynamic_cast<TGraph*>(object);
    if (graph) {
      graphs.push_back(graph);
    }
  }
}

Quality TrendCheck::check(std::map<std::string, std::shared_ptr<MonitorObject>>* moMap)
{
  for (auto& [moKey, mo] : *moMap) {

    std::vector<TGraph*> graphs;
    getGraphs(mo->getObject(), graphs);
    std::cout << "Checking MO " << mo->getName() << " with " << graphs.size() << " graphs" << std::endl;
    if (graphs.empty()) {
      continue;
    }

    auto moName = mo->getName();
    auto key = getBaseName(moName);

    for (size_t graphIndex = 0; graphIndex < graphs.size(); graphIndex++) {

      TGraph* graph = graphs[graphIndex];
      if (!graph) {
        continue;
      }

      auto graphName = moName + "_" + std::to_string(graphIndex);

      // check that the graph is not empty
      int nPoints = graph->GetN();
      std::cout << "Checking graph " << graphName << " (" << graph->GetN() << " points)" << std::endl;
      if (nPoints < 1) {
        continue;
      }

      // get the value for the last point
      double value = graph->GetPointY(nPoints - 1);

      // get acceptable range for the current plot
      initThresholds(key);
      auto thresholds = getThresholds(key, graph);
      // check that at least the thresholds for Bad quality are available
      if (!thresholds[0]) {
        std::cout << "  Thresholds cannot be retrieved" << std::endl;
        continue;
      }

      std::cout << "  Thresholds Bad: " << thresholds[0].value().first << ", " << thresholds[0].value().second << std::endl;
      std::cout << "  Value: " << value << std::endl;

      // Quality is good by default, unless the last point is outside the acceptable range
      mQualities[graphName] = Quality::Good;

      mThresholdsBadTrend[graphName].emplace_back(graph->GetPointX(nPoints - 1), thresholds[0].value());
      if (thresholds[1].has_value()) {
        mThresholdsMediumTrend[graphName].emplace_back(graph->GetPointX(nPoints - 1), thresholds[1].value());
      }

      if (value < thresholds[0]->first || value > thresholds[0]->second) {
        mQualities[graphName] = Quality::Bad;
      } else if (thresholds[1].has_value()) {
        if (value < thresholds[1]->first || value > thresholds[1]->second) {
          mQualities[graphName] = Quality::Medium;
        }
      }
    }
  }

  // compute overall quality
  Quality result = mQualities.empty() ? Quality::Null : Quality::Good;
  for (auto& [key, q] : mQualities) {
    if (q.isWorseThan(result)) {
      result = q;
    }
  }

  return result;
}

std::string TrendCheck::getAcceptedType() { return "TObject"; }

static void drawThresholds(TGraph* graph, const std::vector<std::pair<double, std::pair<double, double>>>& thresholds, int lineColor, int lineStyle)
{
  if (thresholds.empty()) {
    return;
  }

  double rangeMin = TMath::MinElement(graph->GetN(), graph->GetY());
  double rangeMax = TMath::MaxElement(graph->GetN(), graph->GetY());

  double* xValues = new double[thresholds.size()];
  double* yValuesMin = new double[thresholds.size()];
  double* yValuesMax = new double[thresholds.size()];

  for (size_t index = 0; index < thresholds.size(); index++) {
    const auto& thresholdsPoint = thresholds[index];
    xValues[index] = thresholdsPoint.first;
    yValuesMin[index] = thresholdsPoint.second.first;
    yValuesMax[index] = thresholdsPoint.second.second;

    std::cout << "    threshold trend point [" << index << "]: " << xValues[index] << " -> " << yValuesMin[index] << "," <<yValuesMax[index] << std::endl;

    rangeMin = TMath::Min(rangeMin, yValuesMin[index]);
    rangeMax = TMath::Max(rangeMax, yValuesMax[index]);
  }

  TPolyLine* lineMin = new TPolyLine(thresholds.size(), xValues, yValuesMin);
  lineMin->SetLineColor(lineColor);
  lineMin->SetLineStyle(lineStyle);
  graph->GetListOfFunctions()->Add(lineMin);

  TPolyLine* lineMax = new TPolyLine(thresholds.size(), xValues, yValuesMax);
  lineMax->SetLineColor(lineColor);
  lineMax->SetLineStyle(lineStyle);
  graph->GetListOfFunctions()->Add(lineMax);

  auto delta = rangeMax - rangeMin;
  rangeMin -= 0.1 * delta;
  rangeMax += 0.1 * delta;

  graph->SetMinimum(rangeMin);
  graph->SetMaximum(rangeMax);
}

void TrendCheck::beautify(std::shared_ptr<MonitorObject> mo, Quality checkResult)
{
  std::vector<TGraph*> graphs;
  getGraphs(mo->getObject(), graphs);
  std::cout << "Beautifying MO " << mo->getName() << " with " << graphs.size() << " graphs" << std::endl;
  if (graphs.empty()) {
    return;
  }

  auto moName = mo->getName();

  for (size_t graphIndex = 0; graphIndex < graphs.size(); graphIndex++) {

    TGraph* graph = graphs[graphIndex];
    if (!graph || graph->GetN() < 1) {
      return;
    }

    auto graphName = moName + "_" + std::to_string(graphIndex);

    Quality quality = mQualities[graphName];

    std::cout << "Beautifying " << graphName << std::endl;

    // draw the graph in red if the quality is Bad
    if (quality == Quality::Bad) {
      graph->SetLineColor(kRed);
      graph->SetMarkerColor(kRed);
    }

    if (mThresholdsMediumTrend.count(graphName) > 0) {
      const auto& thresholds = mThresholdsMediumTrend[graphName];
      drawThresholds(graph, thresholds, kOrange, kDotted);
    }

    if (mThresholdsBadTrend.count(graphName) > 0) {
      const auto& thresholds = mThresholdsBadTrend[graphName];

      std::cout << "  Bad thresholds trend size: " << thresholds.size() << std::endl;

      drawThresholds(graph, thresholds, kRed, kDashed);
    }
  }
}

} // namespace o2::quality_control_modules::glo
