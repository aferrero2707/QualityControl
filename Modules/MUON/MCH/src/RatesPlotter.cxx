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
/// \file   RatesPlotter.cxx
/// \author Andrea Ferrero
///

#include "MCH/RatesPlotter.h"
#include "MCH/Helpers.h"
#include "MCHMappingInterface/Segmentation.h"
#include <fmt/format.h>

using namespace o2::mch::raw;

namespace o2
{
namespace quality_control_modules
{
namespace muonchambers
{

RatesPlotter::RatesPlotter(std::string path, TH2F* hRef, float rateMin, float rateMax, bool fullPlots)
{
  // mappers used for filling the histograms in detector coordinates
  mElec2DetMapper = createElec2DetMapper<ElectronicMapperGenerated>();
  mFeeLink2SolarMapper = createFeeLink2SolarMapper<ElectronicMapperGenerated>();

  // reductor for the rates plot in electronics coordinates
  mElecMapReductor = std::make_unique<TH2ElecMapReductor>(rateMin, rateMax);

  //----------------------------------
  // Reference mean rates histogram
  //----------------------------------

  if (hRef) {
    TH2ElecMapReductor elecMapReductorRef(rateMin, rateMax);
    elecMapReductorRef.update(hRef);

    mHistogramMeanRatePerDERef =
      std::make_unique<TH1F>(TString::Format("%sMeanRateRef", path.c_str()), "Mean Rate vs DE, reference",
                             getNumDE(), 0, getNumDE());
    mHistogramMeanRatePerDERef->SetLineColor(kRed);
    mHistogramMeanRatePerDERef->SetLineStyle(kDashed);
    mHistogramMeanRatePerDERef->SetLineWidth(2);

    for (size_t de = 0; de < mHistogramMeanRatePerDERef->GetXaxis()->GetNbins(); de++) {
      mHistogramMeanRatePerDERef->SetBinContent(de + 1, elecMapReductorRef.getDeValue(de));
      mHistogramMeanRatePerDERef->SetBinError(de + 1, 0);
    }
  }

  //----------------------------------
  // Mean rates histograms
  //----------------------------------

  mHistogramMeanRatePerDE = std::make_unique<TH1F>(TString::Format("%sMeanRateHist", path.c_str()), "Mean Rate vs DE",
                                                   getNumDE(), 0, getNumDE());

  mHistogramMeanRateRefRatio = std::make_unique<TH1F>(TString::Format("%sMeanRateRefRatio", path.c_str()), "Mean Rate vs DE, ratio wrt reference",
                                                      getNumDE(), 0, getNumDE());
  addHisto(mHistogramMeanRateRefRatio.get(), false, "histo", "histo");

  mCanvasMeanRatePerDE = std::make_unique<TCanvas>(TString::Format("%sMeanRate", path.c_str()), "Mean Rate vs DE", 800, 600);
  mCanvasMeanRatePerDE->SetLogy(kTRUE);
  addCanvas(mCanvasMeanRatePerDE.get(), mHistogramMeanRatePerDE.get(), false, "histo", "histo");

  mHistogramBadChannelsFractionPerDE = std::make_unique<TH1F>(TString::Format("%sBadChannelsFraction", path.c_str()),
                                                              "Bad channels fraction", getNumDE(), 0, getNumDE());
  addHisto(mHistogramBadChannelsFractionPerDE.get(), false, "histo", "histo");

  //--------------------------------------------------
  // Rates histograms in global detector coordinates
  //--------------------------------------------------

  mHistogramRateGlobal[0] = std::make_shared<GlobalHistogram>(fmt::format("{}Rate_ST12", path), "ST12 Rate", 0, 5);
  mHistogramRateGlobal[0]->init();
  addHisto(mHistogramRateGlobal[0]->getHist(), false, "colz", "colz");

  mHistogramRateGlobal[1] = std::make_shared<GlobalHistogram>(fmt::format("{}Rate_ST345", path), "ST345 Rate", 1, 10);
  mHistogramRateGlobal[1]->init();
  addHisto(mHistogramRateGlobal[1]->getHist(), false, "colz", "colz");

  //--------------------------------------------------
  // Rates histograms in detector coordinates
  //--------------------------------------------------

  for (auto de : o2::mch::constants::deIdsForAllMCH) {
    auto h = std::make_shared<DetectorHistogram>(TString::Format("%s%sRate_XY_B_%03d", path.c_str(), getHistoPath(de).c_str(), de),
                                                 TString::Format("Hit Rate (DE%03d B)", de), de, int(0));
    mHistogramRateDE[0].insert(make_pair(de, h));
    if (fullPlots) {
      addHisto(h->getHist(), false, "colz", "colz");
    }

    h = std::make_shared<DetectorHistogram>(TString::Format("%s%sRate_XY_NB_%03d", path.c_str(), getHistoPath(de).c_str(), de),
                                            TString::Format("Hit Rate (DE%03d NB)", de), de, int(1));
    mHistogramRateDE[1].insert(make_pair(de, h));
    if (fullPlots) {
      addHisto(h->getHist(), false, "colz", "colz");
    }
  }
}

//_________________________________________________________________________________________

void RatesPlotter::fillAverageHistos(TH2F* hRates)
{
  // extract the integrated average occupancies
  mElecMapReductor->update(hRates);
  for (size_t de = 0; de < mHistogramMeanRatePerDE->GetXaxis()->GetNbins(); de++) {
    mHistogramMeanRatePerDE->SetBinContent(de + 1, mElecMapReductor->getDeValue(de));
    mHistogramMeanRatePerDE->SetBinError(de + 1, 0.1);
  }

  mCanvasMeanRatePerDE->Clear();
  mCanvasMeanRatePerDE->cd();
  mHistogramMeanRatePerDE->Draw();

  if (mHistogramMeanRatePerDERef) {
    mHistogramMeanRatePerDERef->Draw("histsame");

    mHistogramMeanRateRefRatio->Reset();
    mHistogramMeanRateRefRatio->Add(mHistogramMeanRatePerDE.get());
    mHistogramMeanRateRefRatio->Divide(mHistogramMeanRatePerDERef.get());

    // special handling of bins with zero rate in reference
    int nbinsx = mHistogramMeanRatePerDERef->GetXaxis()->GetNbins();
    for (int b = 1; b <= nbinsx; b++) {
      if (mHistogramMeanRatePerDERef->GetBinContent(b) == 0) {
        mHistogramMeanRateRefRatio->SetBinContent(b, 1);
        mHistogramMeanRateRefRatio->SetBinError(b, 0);
      }
    }
  }

  for (size_t de = 0; de < mHistogramBadChannelsFractionPerDE->GetXaxis()->GetNbins(); de++) {
    float nPads = mElecMapReductor->getNumPads(de);
    float nPadsBad = mElecMapReductor->getNumPadsBad(de) + mElecMapReductor->getNumPadsNoStat(de);
    if (nPads > 0) {
      mHistogramBadChannelsFractionPerDE->SetBinContent(de + 1, nPadsBad / nPads);
      mHistogramBadChannelsFractionPerDE->SetBinError(de + 1, 0.1);
    } else {
      mHistogramBadChannelsFractionPerDE->SetBinContent(de + 1, 0);
      mHistogramBadChannelsFractionPerDE->SetBinError(de + 1, 1);
    }
  }
}

//_________________________________________________________________________________________

void RatesPlotter::fillGlobalHistos(TH2F* h)
{
  if (!h) {
    return;
  }

  // loop over bins in electronics coordinates, and map the channels to the corresponding cathode pads
  int nbinsx = h->GetXaxis()->GetNbins();
  int nbinsy = h->GetYaxis()->GetNbins();
  for (int i = 1; i <= nbinsx; i++) {
    // address of the DS board in FEC representation
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

    for (int j = 1; j <= nbinsy; j++) {
      int channel = j - 1;
      int padId = -1;

      const o2::mch::mapping::Segmentation& segment = o2::mch::mapping::segmentation(deId);
      padId = segment.findPadByFEE(dsIddet, int(channel));

      if (padId < 0) {
        continue;
      }

      double rate = h->GetBinContent(i, j);

      double padX = segment.padPositionX(padId);
      double padY = segment.padPositionY(padId);
      float padSizeX = segment.padSizeX(padId);
      float padSizeY = segment.padSizeY(padId);
      int cathode = segment.isBendingPad(padId) ? 0 : 1;

      // Fill 2D rate histograms
      auto hRate = mHistogramRateDE[cathode].find(deId);
      if ((hRate != mHistogramRateDE[cathode].end()) && (hRate->second != NULL)) {
        hRate->second->Set(padX, padY, padSizeX, padSizeY, rate);
      }
    }
  }

  mHistogramRateGlobal[0]->set(mHistogramRateDE[0], mHistogramRateDE[1]);
  mHistogramRateGlobal[1]->set(mHistogramRateDE[0], mHistogramRateDE[1]);
}

//_________________________________________________________________________________________

void RatesPlotter::update(TH2F* hRates)
{
  fillAverageHistos(hRates);
  fillGlobalHistos(hRates);
}

} // namespace muonchambers
} // namespace quality_control_modules
} // namespace o2
