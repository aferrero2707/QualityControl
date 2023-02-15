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
/// \file   HeartBeatPacketsPlotter.h
/// \author Andrea Ferrero
///

#ifndef QC_MODULE_HBPACKETSPLOTTER_H
#define QC_MODULE_HBPACKETSPLOTTER_H

#include "MUONCommon/HistPlotter.h"
#include "MCH/GlobalHistogram.h"
#include "MCHRawElecMap/Mapper.h"
#include <TH1F.h>
#include <TH2F.h>

using namespace o2::quality_control_modules::muon;

namespace o2
{
namespace quality_control_modules
{
namespace muonchambers
{

class HeartBeatPacketsPlotter : public HistPlotter
{
 public:
  HeartBeatPacketsPlotter(std::string path, int hbExpectedBc = 456190);

  void update(TH2F* h2);

 private:
  int getDeId(uint16_t feeId, uint8_t linkId, uint8_t eLinkId);

  void addHisto(TH1* h, bool statBox, const char* drawOptions, const char* displayHints)
  {
    h->SetOption(drawOptions);
    if (!statBox) {
      h->SetStats(0);
    }
    histograms().emplace_back(HistInfo{ h, drawOptions, displayHints });
  }

  o2::mch::raw::Elec2DetMapper mElec2DetMapper;
  o2::mch::raw::FeeLink2SolarMapper mFeeLink2SolarMapper;

  /// \brief expected bunch-crossing value in heart-beat packets
  int mHBExpectedBc;

  std::unique_ptr<TH1F> mHistogramSynchErrorsPerDE;      ///< fraction of out-of-sync DS boards per detection element
  std::unique_ptr<TH1F> mHistogramSynchErrorsPerChamber; ///< fraction of out-of-sync DS boards per chamber
  std::unique_ptr<TH2F> mSyncStatusFEC;                  ///< time synchronization status of each DS board (OK, out-of-synch, missing good HB)

  std::map<int, std::shared_ptr<DetectorHistogram>> mHistogramHBRateDE[2]; // 2D hit rate map for each DE
  std::shared_ptr<GlobalHistogram> mHistogramHBRateGlobal[2];              // Rate histogram (global XY view)
};

} // namespace muonchambers
} // namespace quality_control_modules
} // namespace o2

#endif // QC_MODULE_HBPACKETSPLOTTER_H
