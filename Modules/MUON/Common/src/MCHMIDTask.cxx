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
/// \file   MCHMIDQcTask.cxx
/// \author Bogdan Vulpescu
/// \author Xavier Lopez
/// \author Diego Stocco
/// \author Guillaume Taillepied
///

#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "QualityControl/QcInfoLogger.h"
#include "MUONCommon/MCHMIDTask.h"
#include <Framework/InputRecord.h>

#include "Framework/DataRefUtils.h"
#include "DataFormatsMID/ColumnData.h"
#include "DataFormatsMID/ROFRecord.h"
#include "DataFormatsMCH/Digit.h"
#include "DataFormatsMCH/ROFRecord.h"
#include "MCHBase/PreCluster.h"
#ifdef MCH_HAS_MAPPING_FACTORY
#include "MCHMappingFactory/CreateSegmentation.h"
#endif
#include "MCHMappingInterface/Segmentation.h"
#include "MCHMappingInterface/CathodeSegmentation.h"
#include "MCHMappingSegContour/CathodeSegmentationContours.h"

namespace o2::quality_control_modules::muon
{

MCHMIDQcTask::~MCHMIDQcTask()
{
}

void MCHMIDQcTask::initialize(o2::framework::InitContext& /*ctx*/)
{
  ILOG(Info) << "initialize MCHMIDQcTask" << ENDM; // QcInfoLogger is used. FairMQ logs will go to there as well.

  // Histograms to be published
  mTimeCorrelation = new TH1F("TimeCorrelation", "Time correlation", 2000, -1000, 1000);
  getObjectsManager()->startPublishing(mTimeCorrelation);

  mColumnSize = new TH2F("ColumnSize", "Column size", 100, 0, 100, 100, 0, 100);
  mColumnSize->SetOption("colz");
  getObjectsManager()->startPublishing(mColumnSize);

  mRofSize = new TH2F("RofSize", "ROF size", 100, 0, 100, 100, 0, 100);
  mRofSize->SetOption("colz");
  getObjectsManager()->startPublishing(mRofSize);

  mRofSizeInTF_MID = new TH1F("mRofSizeInTF_MID", "ROF size in TF - MID", 3600 * 128, 0, 3600 * 128);
  mRofSizeInTF_MCH = new TH1F("mRofSizeInTF_MCH", "ROF size in TF - MCH", 3600 * 128, 0, 3600 * 128);
  mTcRofSizeInTF_MCH = new TH1F("mTcRofSizeInTF_MCH", "TC ROF size in TF - MCH", 3600 * 128, 0, 3600 * 128);
  mRofSizeInTF_MCHms = new TH1F("mRofSizeInTF_MCHms", "ROF size in TF - MCH (ms)", (3564 * 128) * 25 / 10000, 0, (3564 * 128) * 25 / 1000000);
  mTcRofNStationsInTF_MCH = new TH1F("mTcRofNStationsInTF_MCH", "TC ROF # of stations in TF - MCH", 3600 * 128, 0, 3600 * 128);
  mTcRofNStations_MCH = new TH1F("mTcRofNStations_MCH", "TC ROF # of stations - MCH", 12, 0, 12);
  getObjectsManager()->startPublishing(mRofSizeInTF_MID);
  getObjectsManager()->startPublishing(mRofSizeInTF_MCH);
  getObjectsManager()->startPublishing(mTcRofSizeInTF_MCH);
  getObjectsManager()->startPublishing(mRofSizeInTF_MCHms);
  getObjectsManager()->startPublishing(mTcRofNStationsInTF_MCH);
  getObjectsManager()->startPublishing(mTcRofNStations_MCH);

  mRofSizeInOrbit_MID = new TH1F("mRofSizeInOrbit_MID", "ROF size in Orbit - MID", 3564, 0, 3564);
  mRofSizeInOrbit_MCH = new TH2F("mRofSizeInOrbit_MCH", "ROF size in Orbit - MCH", 3564/4, 0, 3564, 10, 1, 11);
  mRofSizeInOrbit_MCH->SetDrawOption("col");
  getObjectsManager()->startPublishing(mRofSizeInOrbit_MID);
  getObjectsManager()->startPublishing(mRofSizeInOrbit_MCH);

  mDigitsInOrbit_MCH = new TH2F("mDigitsInOrbit_MCH", "Digits in Orbit - MCH", 3564/4, 0, 3564, 1200, 0, 1200);
  mDigitsInOrbit_MCH->SetDrawOption("col");
  getObjectsManager()->startPublishing(mDigitsInOrbit_MCH);
}

void MCHMIDQcTask::startOfActivity(Activity& /*activity*/)
{
  ILOG(Info) << "startOfActivity" << ENDM;
}

void MCHMIDQcTask::startOfCycle()
{
  ILOG(Info) << "startOfCycle" << ENDM;
}

static int countColumnDataHits(const o2::mid::ColumnData& digit, int id)
{
  int nHits = 0;
  int mask = 1;
  for (int j = 0; j < 16; j++) {
    if ((digit.patterns[id] & mask) != 0) {
      nHits += 1;
    }
    mask <<= 1;
  }
  return nHits;
}

static int getBendingHits(const o2::mid::ColumnData& digit)
{
  int nHits = 0;
  for (int i = 0; i < 4; i++) {
    nHits += countColumnDataHits(digit, i);
  }
  return nHits;
}

static int getNonBendingHits(const o2::mid::ColumnData& digit)
{
  return countColumnDataHits(digit, 4);
}

static std::pair<uint32_t, uint32_t> getROFSize(const o2::mid::ROFRecord& rof, gsl::span<const o2::mid::ColumnData> digits)
{
  uint32_t nHitsB{ 0 };
  uint32_t nHitsNB{ 0 };

  auto lastEntry = rof.getEndIndex();
  for (auto i = rof.firstEntry; i < lastEntry; i++) {
    const auto& digit = digits[i];
    nHitsB += getBendingHits(digit);
    nHitsNB += getNonBendingHits(digit);
  }

  return std::make_pair(nHitsB, nHitsNB);
}

void MCHMIDQcTask::monitorData(o2::framework::ProcessingContext& ctx)
{
  ILOG(Info) << "startOfDataMonitoring" << ENDM;

  const auto* dh = o2::framework::DataRefUtils::getHeader<o2::header::DataHeader*>(ctx.inputs().getFirstValid(true));
  int firstOrbit = dh->firstTForbit;

  auto mchrofs = ctx.inputs().get<gsl::span<o2::mch::ROFRecord>>("mchrofs");
  auto tcrofs = ctx.inputs().get<gsl::span<o2::mch::ROFRecord>>("tcrofs");
  auto mchdigits = ctx.inputs().get<gsl::span<o2::mch::Digit>>("mchdigits");
  auto midDigits = ctx.inputs().get<gsl::span<o2::mid::ColumnData>>("middigits");
  auto midrofs = ctx.inputs().get<gsl::span<o2::mid::ROFRecord>>("midrofs");
  auto mchpreClusters = ctx.inputs().get<gsl::span<o2::mch::PreCluster>>("preclusters");
  auto mchdigits2 = ctx.inputs().get<gsl::span<o2::mch::Digit>>("preclusterdigits");

  std::cout << fmt::format("MCH digits {}  rofs {}", mchdigits.size(), mchrofs.size()) << std::endl;

  for (const auto& digit : midDigits) {
    // total number of bending plane hits
    int nHitsB = getBendingHits(digit);
    // total number of non-bending plane hits
    int nHitsNB = getNonBendingHits(digit);
    mColumnSize->Fill(nHitsB, nHitsNB);
  }

  //int firstOrbit = 717215;
  //int firstOrbit = 717343;

  for (auto& midrof: midrofs) {
    auto rofsize = getROFSize(midrof, midDigits);
    int64_t orbit = midrof.interactionRecord.orbit;
    int64_t dOrbit = orbit - firstOrbit;

    mRofSizeInTF_MID->SetBinContent(dOrbit * 3564 + midrof.interactionRecord.bc, rofsize.first + rofsize.second);

    mRofSizeInOrbit_MID->Fill(midrof.interactionRecord.bc, rofsize.first + rofsize.second);
  }

  for (auto& mchrof: mchrofs) {
    //auto rofsize = mchrof.getNEntries();
    int rofsize[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int rofsizeTot = 0;
    for (int i = mchrof.getFirstIdx(); i <= mchrof.getLastIdx(); i++) {
      auto& digit = mchdigits[i];
      int chamber = digit.getDetID() / 100;
      //std::cout << digit.getDetID() << "  chamber  " << chamber <<  std::endl;
      //if (digit.getDetID() >= 100 && digit.getDetID() < 700) { rofsize += 1; }
      rofsize[chamber - 1] += 1;
      rofsizeTot += 1;
    }

    int64_t orbit = mchrof.getBCData().orbit;
    int64_t dOrbit = orbit - firstOrbit;

    mRofSizeInTF_MCH->SetBinContent(dOrbit * 3564 + mchrof.getBCData().bc, rofsizeTot);

    for (int chamber = 1; chamber <= 10; chamber++) {
      mRofSizeInOrbit_MCH->Fill(mchrof.getBCData().bc, chamber + 0.1, rofsize[chamber - 1]);
    }

    double ms = dOrbit * 3564 + mchrof.getBCData().bc;
    ms *= 25;
    ms /= 1000000;
    //if(rofsize > 0 && dOrbit < 128) std::cout<<fmt::format("dOrbit {}  BC {}  ms {}", dOrbit, mchrof.getBCData().bc, ms) << std::endl;
    mRofSizeInTF_MCHms->Fill(ms, rofsizeTot);
  }

  for (auto& rof: tcrofs) {
    //auto rofsize = mchrof.getNEntries();
    int rofsize[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int rofsizeTot = 0;
    for (int i = rof.getFirstIdx(); i <= rof.getLastIdx(); i++) {
      auto& digit = mchdigits[i];
      int chamber = digit.getDetID() / 100;
      //std::cout << digit.getDetID() << "  chamber  " << chamber <<  std::endl;
      //if (digit.getDetID() >= 100 && digit.getDetID() < 700) { rofsize += 1; }
      rofsize[chamber - 1] += 1;
      rofsizeTot += 1;
    }

    int64_t orbit = rof.getBCData().orbit;
    int64_t dOrbit = orbit - firstOrbit;

    mTcRofSizeInTF_MCH->SetBinContent(dOrbit * 3564 + rof.getBCData().bc, rofsizeTot);

    int nChambers = 0;
    for (int ch = 0; ch < 10; ch += 1) {
      if (rofsize[ch] > 0) {
        nChambers += 1;
      }
    }

    int nStations = 0;
    for (int ch = 0; ch < 10; ch += 2) {
      if (rofsize[ch] > 0 || rofsize[ch + 1] > 0) {
        nStations += 1;
      }
    }
    mTcRofNStationsInTF_MCH->SetBinContent(dOrbit * 3564 + rof.getBCData().bc, nStations);

    mTcRofNStations_MCH->Fill(nStations);
}

  for (auto& preCluster : mchpreClusters) {

    // get the digits of this precluster
    auto preClusterDigits = mchdigits2.subspan(preCluster.firstDigit, preCluster.nDigits);

    // whether a cathode has digits or not
    bool cathode[2] = { false, false };
    // total charge on each cathode
    float chargeSum[2] = { 0, 0 };
    // largest signal in each cathode
    float chargeMax[2] = { 0, 0 };
    // number of digits in each cathode
    int multiplicity[2] = { 0, 0 };

    int detid = preClusterDigits[0].getDetID();
    const o2::mch::mapping::Segmentation& segment = o2::mch::mapping::segmentation(detid);

    // loop over digits and collect information on charge and multiplicity
    for (const o2::mch::Digit& digit : preClusterDigits) {
      int padid = digit.getPadID();

      // cathode index
      int cid = segment.isBendingPad(padid) ? 0 : 1;
      cathode[cid] = true;
      chargeSum[cid] += digit.getADC();
      multiplicity[cid] += 1;

      if (digit.getADC() > chargeMax[cid]) {
        chargeMax[cid] = digit.getADC();
      }
    }

    if (multiplicity[0] < 2 || multiplicity[1] < 2) { continue; }
    //if (multiplicity[0] < 3) { continue; }
    if ((chargeSum[0]+chargeSum[1]) < 100) { continue; }

    int chamber = detid / 100;
    for (const o2::mch::Digit& digit : preClusterDigits) {
      mDigitsInOrbit_MCH->Fill(digit.getTime() % 3564, detid + 0.1);
    }
  }


  for (auto& midrof: midrofs) {

    //if (midrof.interactionRecord.bc > 10) { continue; }
    auto rofsize = getROFSize(midrof, midDigits);
    mRofSize->Fill(rofsize.first, rofsize.second);

    if (rofsize.first < 5) { continue; }
    if (rofsize.second < 1) { continue; }
    //if(rofsize.second > 5)
    //  std::cout << fmt::format("MID: {},{} {}/{}", rofsize.first, rofsize.second, midrof.interactionRecord.orbit, midrof.interactionRecord.bc) << std::endl;

    for (auto& mchrof: mchrofs) {
      if (mchrof.getNEntries() < 50) { continue; }

      int64_t orbit1 = mchrof.getBCData().orbit;
      int64_t orbit2 = midrof.interactionRecord.orbit;
      int64_t dOrbit = orbit2 - orbit1;
      if (dOrbit < -1 || dOrbit > 1) { continue; }

      auto bcdiff = mchrof.getBCData().differenceInBC(midrof.interactionRecord);
      mTimeCorrelation->Fill(bcdiff);
      //if (rofsize.second > 5) {
      //  std::cout << fmt::format("entries: MCH {}   MID {}/{}   diff {}",
      //      mchrof.getNEntries(), rofsize.first, rofsize.second, bcdiff) << std::endl;
      //}
    }
  }
}

void MCHMIDQcTask::endOfCycle()
{
  ILOG(Info) << "endOfCycle" << ENDM;
}

void MCHMIDQcTask::endOfActivity(Activity& /*activity*/)
{
  ILOG(Info) << "endOfActivity" << ENDM;
}

void MCHMIDQcTask::reset()
{
  // clean all the monitor objects here

  ILOG(Info) << "Resetting the histogram" << ENDM;
}

} // namespace o2::quality_control_modules::mid
