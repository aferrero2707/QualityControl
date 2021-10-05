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
#include "DataFormatsMCH/ROFRecord.h"

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
  getObjectsManager()->startPublishing(mRofSizeInTF_MID);
  getObjectsManager()->startPublishing(mRofSizeInTF_MCH);
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

  auto mchrofs = ctx.inputs().get<gsl::span<o2::mch::ROFRecord>>("mchrofs");
  auto midDigits = ctx.inputs().get<gsl::span<o2::mid::ColumnData>>("middigits");
  auto midrofs = ctx.inputs().get<gsl::span<o2::mid::ROFRecord>>("midrofs");

  for (const auto& digit : midDigits) {
    // total number of bending plane hits
    int nHitsB = getBendingHits(digit);
    // total number of non-bending plane hits
    int nHitsNB = getNonBendingHits(digit);
    mColumnSize->Fill(nHitsB, nHitsNB);
  }

  int firstOrbit = 717215;

  for (auto& midrof: midrofs) {
    auto rofsize = getROFSize(midrof, midDigits);
    int64_t orbit = midrof.interactionRecord.orbit;
    int64_t dOrbit = orbit - firstOrbit;

    mRofSizeInTF_MID->SetBinContent(dOrbit * 3564 + midrof.interactionRecord.bc, rofsize.first + rofsize.second);
  }

  for (auto& mchrof: mchrofs) {
    auto rofsize = mchrof.getNEntries();
    int64_t orbit = mchrof.getBCData().orbit;
    int64_t dOrbit = orbit - firstOrbit;

    mRofSizeInTF_MCH->SetBinContent(dOrbit * 3564 + mchrof.getBCData().bc, rofsize);
  }

  for (auto& midrof: midrofs) {

    //if (midrof.interactionRecord.bc > 10) { continue; }
    auto rofsize = getROFSize(midrof, midDigits);
    mRofSize->Fill(rofsize.first, rofsize.second);

    if (rofsize.first < 5) { continue; }
    if (rofsize.second < 1) { continue; }
    if(rofsize.second > 5)
      std::cout << fmt::format("MID: {},{} {}/{}", rofsize.first, rofsize.second, midrof.interactionRecord.orbit, midrof.interactionRecord.bc) << std::endl;

    for (auto& mchrof: mchrofs) {
      if (mchrof.getNEntries() < 50) { continue; }

      int64_t orbit1 = mchrof.getBCData().orbit;
      int64_t orbit2 = midrof.interactionRecord.orbit;
      int64_t dOrbit = orbit2 - orbit1;
      if (dOrbit < -1 || dOrbit > 1) { continue; }

      auto bcdiff = mchrof.getBCData().differenceInBC(midrof.interactionRecord);
      mTimeCorrelation->Fill(bcdiff);
      if (rofsize.second > 5) {
        std::cout << fmt::format("entries: MCH {}   MID {}/{}   diff {}",
            mchrof.getNEntries(), rofsize.first, rofsize.second, bcdiff) << std::endl;
      }
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
