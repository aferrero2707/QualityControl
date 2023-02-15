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
/// \file   HeartBeatPacketsPlotter.cxx
/// \author Andrea Ferrero
///

#include "MCH/HeartBeatPacketsPlotter.h"
#include "MCH/Helpers.h"
#include "MCHMappingInterface/Segmentation.h"
#include <iostream>
#include <fmt/format.h>

using namespace o2::mch::raw;

namespace o2
{
namespace quality_control_modules
{
namespace muonchambers
{

static void setYAxisLabels(TH2F* hErrors)
{
  TAxis* ay = hErrors->GetYaxis();
  for (int i = 1; i <= 10; i++) {
    auto label = fmt::format("CH{}", i);
    ay->SetBinLabel(i, label.c_str());
  }
}

HeartBeatPacketsPlotter::HeartBeatPacketsPlotter(std::string path, int hbExpectedBc) : mHBExpectedBc(hbExpectedBc)
{
  bool fullPlots = true;
  mElec2DetMapper = createElec2DetMapper<ElectronicMapperGenerated>();
  mFeeLink2SolarMapper = createFeeLink2SolarMapper<ElectronicMapperGenerated>();

  //--------------------------------------------
  // Decoding errors per chamber, DE and FEEID
  //--------------------------------------------
  const uint32_t nElecXbins = FecId::max();

  mSyncStatusFEC = std::make_unique<TH2F>(TString::Format("%sSyncStatusFEC", path.c_str()), "Heart-beat status vs. FEC ID", nElecXbins, 0, nElecXbins, 3, 0, 3);
  mSyncStatusFEC->GetYaxis()->SetBinLabel(1, "OK");
  mSyncStatusFEC->GetYaxis()->SetBinLabel(2, "Out-of-sync");
  mSyncStatusFEC->GetYaxis()->SetBinLabel(3, "Missing");
  addHisto(mSyncStatusFEC.get(), false, "col", "");

  mHistogramSynchErrorsPerDE = std::make_unique<TH1F>(TString::Format("%sSynchErrorsPerDE", path.c_str()), "Out-of-sync boards fraction per DE", getNumDE(), 0, getNumDE());
  addHisto(mHistogramSynchErrorsPerDE.get(), false, "hist", "");

  mHistogramSynchErrorsPerChamber = std::make_unique<TH1F>(TString::Format("%sSynchErrorsPerChamber", path.c_str()), "Out-of-sync boards fraction per chamber", 10, 0, 10);
  addHisto(mHistogramSynchErrorsPerChamber.get(), false, "hist", "");

  //--------------------------------------------------
  // Rates histograms in global detector coordinates
  //--------------------------------------------------

  mHistogramHBRateGlobal[0] = std::make_shared<GlobalHistogram>(fmt::format("{}Rate_ST12", path), "ST12 Rate", 0, 5);
  mHistogramHBRateGlobal[0]->init();
  addHisto(mHistogramHBRateGlobal[0]->getHist(), false, "colz", "colz");

  mHistogramHBRateGlobal[1] = std::make_shared<GlobalHistogram>(fmt::format("{}Rate_ST345", path), "ST345 Rate", 1, 10);
  mHistogramHBRateGlobal[1]->init();
  addHisto(mHistogramHBRateGlobal[1]->getHist(), false, "colz", "colz");

  //--------------------------------------------------
  // Rates histograms in detector coordinates
  //--------------------------------------------------

  for (auto de : o2::mch::constants::deIdsForAllMCH) {
    auto h = std::make_shared<DetectorHistogram>(TString::Format("%s%sRate_XY_B_%03d", path.c_str(), getHistoPath(de).c_str(), de),
                                                 TString::Format("Hit Rate (DE%03d B)", de), de, int(0));
    mHistogramHBRateDE[0].insert(make_pair(de, h));
    if (fullPlots) {
      addHisto(h->getHist(), false, "colz", "colz");
    }

    h = std::make_shared<DetectorHistogram>(TString::Format("%s%sRate_XY_NB_%03d", path.c_str(), getHistoPath(de).c_str(), de),
                                            TString::Format("Hit Rate (DE%03d NB)", de), de, int(1));
    mHistogramHBRateDE[1].insert(make_pair(de, h));
    if (fullPlots) {
      addHisto(h->getHist(), false, "colz", "colz");
    }
  }
}

//_____________________________________________________________________________

int HeartBeatPacketsPlotter::getDeId(uint16_t feeId, uint8_t linkId, uint8_t eLinkId)
{
  uint16_t solarId = -1;
  int deId = -1;
  int dsIddet = -1;
  int padId = -1;

  o2::mch::raw::FeeLinkId feeLinkId{ feeId, linkId };

  if (auto opt = mFeeLink2SolarMapper(feeLinkId); opt.has_value()) {
    solarId = opt.value();
  }
  if (solarId < 0 || solarId > 1023) {
    return -1;
  }

  o2::mch::raw::DsElecId dsElecId{ solarId, static_cast<uint8_t>(eLinkId / 5), static_cast<uint8_t>(eLinkId % 5) };

  if (auto opt = mElec2DetMapper(dsElecId); opt.has_value()) {
    o2::mch::raw::DsDetId dsDetId = opt.value();
    dsIddet = dsDetId.dsId();
    deId = dsDetId.deId();
  }

  if (deId < 0 || dsIddet < 0) {
    return -1;
  }

  return deId;
}

//_________________________________________________________________________________________

void HeartBeatPacketsPlotter::update(TH2F* h)
{
  if (!h) {
    return;
  }
  int mBcMin{ mHBExpectedBc - 2 };
  int mBcMax{ mHBExpectedBc + 2 };

  std::vector<double> deNum(static_cast<size_t>(getNumDE()));
  std::vector<double> deDen(static_cast<size_t>(getNumDE()));
  std::array<double, 10> chNum;
  std::array<double, 10> chDen;

  std::fill(deNum.begin(), deNum.end(), 0);
  std::fill(deDen.begin(), deDen.end(), 0);
  std::fill(chNum.begin(), chNum.end(), 0);
  std::fill(chDen.begin(), chDen.end(), 0);

  int nbinsx = h->GetXaxis()->GetNbins();
  int nbinsy = h->GetYaxis()->GetNbins();
  for (int i = 1; i <= nbinsx; i++) {
    FecId fecId(i - 1);
    uint16_t feeId = fecId.getFeeId();
    uint8_t linkId = fecId.getLinkId();
    uint8_t eLinkId = fecId.getDsAddr();

    uint16_t solarId = -1;
    int deId = -1;
    int dsIddet = -1;

    o2::mch::raw::FeeLinkId feeLinkId{ feeId, linkId };

    if (auto opt = mFeeLink2SolarMapper(feeLinkId); opt.has_value()) {
      solarId = opt.value();
    }
    if (solarId < 0 || solarId > 1023) {
      continue;
    }

    o2::mch::raw::DsElecId dsElecId{ solarId, static_cast<uint8_t>(eLinkId / 5), static_cast<uint8_t>(eLinkId % 5) };

    if (auto opt = mElec2DetMapper(dsElecId); opt.has_value()) {
      o2::mch::raw::DsDetId dsDetId = opt.value();
      dsIddet = dsDetId.dsId();
      deId = dsDetId.deId();
    }

    if (deId < 0 || dsIddet < 0) {
      continue;
    }

    int deIndex = getDEindex(deId);
    if (deIndex < 0) {
      continue;
    }

    int chamber = deId / 100 - 1;
    if (chamber < 0 || chamber >= 10) {
      continue;
    }

    deDen[deIndex] += 1;
    chDen[chamber] += 1;

    int ybinmin = h->GetYaxis()->FindBin(mBcMin);
    int ybinmax = h->GetYaxis()->FindBin(mBcMax);
    // number of HB packets in the good bc range, normalized to the number of processed time-frames
    // we expect 2 HB packets per TF and per DS (one per SAMPA chip)
    auto ngood = h->Integral(i, i, ybinmin, ybinmax);

    // total number of HB packets received, including underflow/overflow
    auto total = h->Integral(i, i, 1, nbinsy); // integral over all the bc values
    total += h->GetBinContent(i, 0);           // add underflow
    total += h->GetBinContent(i, nbinsy + 1);  // add overflow

    bool isOutOfSync = false;
    auto nbad = total - ngood;
    if (nbad > 0) {
      isOutOfSync = true;
    }

    bool isMissing = false;
    if (ngood < 1.5) {
      isMissing = true;
    }

    if (isOutOfSync || isMissing) {
      deNum[deIndex] += 1;
      chNum[chamber] += 1;
    }

    if (!isOutOfSync && !isMissing) {
      mSyncStatusFEC->Fill(i - 1, 0);
    }
    if (isOutOfSync) {
      mSyncStatusFEC->Fill(i - 1, 1);
    }
    if (isMissing) {
      mSyncStatusFEC->Fill(i - 1, 2);
    }

    for (int channel = 0; channel < 64; channel++) {
      int padId = -1;

      const o2::mch::mapping::Segmentation& segment = o2::mch::mapping::segmentation(deId);
      padId = segment.findPadByFEE(dsIddet, int(channel));

      if (padId < 0) {
        continue;
      }

      double rate = total;

      double padX = segment.padPositionX(padId);
      double padY = segment.padPositionY(padId);
      float padSizeX = segment.padSizeX(padId);
      float padSizeY = segment.padSizeY(padId);
      int cathode = segment.isBendingPad(padId) ? 0 : 1;

      // Fill 2D rate histograms
      auto hRate = mHistogramHBRateDE[cathode].find(deId);
      if ((hRate != mHistogramHBRateDE[cathode].end()) && (hRate->second != NULL)) {
        hRate->second->Set(padX, padY, padSizeX, padSizeY, rate);
      }
    }
  }

  // update the average number of out-of-sync boards
  for (size_t i = 0; i < deDen.size(); i++) {
    if (deDen[i] > 0) {
      mHistogramSynchErrorsPerDE->SetBinContent(i + 1, deNum[i] / deDen[i]);
    } else {
      mHistogramSynchErrorsPerDE->SetBinContent(i + 1, 0);
    }
  }
  for (size_t ch = 0; ch < chDen.size(); ch++) {
    if (chDen[ch] > 0) {
      mHistogramSynchErrorsPerChamber->SetBinContent(ch + 1, chNum[ch] / chDen[ch]);
    } else {
      mHistogramSynchErrorsPerChamber->SetBinContent(ch + 1, 0);
    }
  }

  mHistogramHBRateGlobal[0]->set(mHistogramHBRateDE[0], mHistogramHBRateDE[1]);
  mHistogramHBRateGlobal[1]->set(mHistogramHBRateDE[0], mHistogramHBRateDE[1]);
}

} // namespace muonchambers
} // namespace quality_control_modules
} // namespace o2
