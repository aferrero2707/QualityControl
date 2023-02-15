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
/// \file    QualityPostProcessing.cxx
/// \author  Andrea Ferrero andrea.ferrero@cern.ch
/// \brief   Post-processing of the MCH quality flags
/// \since   21/06/2022
///

#include "MCH/QualityPostProcessing.h"
#include "QualityControl/QcInfoLogger.h"
#include "QualityControl/MonitorObject.h"
#include "QualityControl/DatabaseInterface.h"
#include <TDatime.h>

using namespace o2::quality_control;
using namespace o2::quality_control::core;
using namespace o2::quality_control::postprocessing;
using namespace o2::quality_control_modules::muonchambers;

void QualityPostProcessing::configure(std::string name, const boost::property_tree::ptree& config)
{
  mConfig = PostProcessingConfigMCH(name, config);

  mCcdbObjects2.emplace(digitsSourceName(), CcdbObjectHelper());
  mCcdbObjects2.emplace(preclustersSourceName(), CcdbObjectHelper());
  mCcdbObjects2.emplace(combinedSourceName(), CcdbObjectHelper());

  for (auto source : mConfig.dataSources) {
    std::string sourceType, sourceName;
    splitDataSourceName(source.name, sourceType, sourceName);
    if (sourceType.empty()) {
      continue;
    }

    mCcdbObjects.emplace_back(source.path, sourceName);

    auto obj = mCcdbObjects2.find(sourceType);
    if (obj != mCcdbObjects2.end()) {
      obj->second.mPath = source.path;
      obj->second.mName = sourceName;
    }
  }
}

//_________________________________________________________________________________________

static void setQualityLabels(TH1F* h)
{
  TAxis* ax = h->GetXaxis();
  ax->SetBinLabel(1, "Null");
  ax->ChangeLabel(1, 45);
  ax->SetBinLabel(2, "Bad");
  ax->ChangeLabel(2, 45);
  ax->SetBinLabel(3, "Good");
  ax->ChangeLabel(3, 45);
}

void QualityPostProcessing::initialize(Trigger t, framework::ServiceRegistryRef services)
{
  for (auto& obj : mCcdbObjects) {
    auto iter1 = mHistogramsQuality.emplace(std::make_pair(obj.mName, std::make_unique<TH1F>(obj.mName.c_str(), obj.mName.c_str(), 3, 0, 3)));
    setQualityLabels(iter1.first->second.get());
    publishHisto(iter1.first->second.get(), false, "hist", "");

    auto iter2 = mTrendsQuality.emplace(std::make_pair(obj.mName, std::make_unique<QualityTrendGraph>((std::string("Trends/") + obj.mName).c_str(), obj.mName.c_str())));
    getObjectsManager()->startPublishing(iter2.first->second.get());
    getObjectsManager()->setDisplayHint(iter2.first->second.get(), "gridy");
  }
}

//_________________________________________________________________________________________

void QualityPostProcessing::update(Trigger t, framework::ServiceRegistryRef services)
{
  auto& qcdb = services.get<repository::DatabaseInterface>();

  for (auto& obj : mCcdbObjects) {
    if (obj.update(&qcdb, t.timestamp, t.activity)) {
      auto qo = obj.get<QualityObject>();
      if (qo) {
        auto time = obj.getTimeStamp() / 1000; // ROOT expects seconds since epoch
        int Q = 0;
        if (qo->getQuality() == Quality::Bad) {
          Q = 1;
        }
        if (qo->getQuality() == Quality::Good) {
          Q = 2;
        }

        auto iter1 = mHistogramsQuality.find(obj.mName);
        if (iter1 != mHistogramsQuality.end()) {
          iter1->second->Fill(Q);
        }

        auto iter2 = mTrendsQuality.find(obj.mName);
        if (iter2 != mTrendsQuality.end()) {
          iter2->second->update(time, qo->getQuality());
        }
      }
    }
  }
}

//_________________________________________________________________________________________

void QualityPostProcessing::finalize(Trigger t, framework::ServiceRegistryRef)
{
}
