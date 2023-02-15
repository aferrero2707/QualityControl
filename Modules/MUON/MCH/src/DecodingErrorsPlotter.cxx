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
/// \file   DecodingErrorsPlotter.cxx
/// \author Andrea Ferrero
///

#include "MCH/DecodingErrorsPlotter.h"
#include "MCH/Helpers.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHRawDecoder/ErrorCodes.h"
#include "MCHBase/DecoderError.h"
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
  TAxis* ax = hErrors->GetYaxis();
  for (int i = 0; i < getErrorCodesSize(); i++) {
    ax->SetBinLabel(i + 1, errorCodeAsString(1 << i).c_str());
    ax->ChangeLabel(i + 1, 45);
  }
}

static void setXAxisLabels(TH2F* hErrors)
{
  TAxis* ay = hErrors->GetXaxis();
  for (int i = 1; i <= 10; i++) {
    auto label = fmt::format("CH{}", i);
    ay->SetBinLabel(i, label.c_str());
  }
}

DecodingErrorsPlotter::DecodingErrorsPlotter(std::string path) : mPath(path)
{
  mElec2DetMapper = createElec2DetMapper<ElectronicMapperGenerated>();
  mDet2ElecMapper = createDet2ElecMapper<ElectronicMapperGenerated>();
  mFeeLink2SolarMapper = createFeeLink2SolarMapper<ElectronicMapperGenerated>();
  mSolar2FeeLinkMapper = createSolar2FeeLinkMapper<ElectronicMapperGenerated>();

  //--------------------------------------------
  // Decoding errors per chamber, DE and FEEID
  //--------------------------------------------

  // Number of decoding errors, grouped by FEE ID and normalized to the number of processed TF
  mHistogramErrorsPerFeeId = std::make_unique<TH2F>(TString::Format("%sDecodingErrorsPerFeeId", path.c_str()),
                                                    "FEE ID vs. Error Type", 64, 0, 64,
                                                    getErrorCodesSize(), 0, getErrorCodesSize());
  setYAxisLabels(mHistogramErrorsPerFeeId.get());
  addHisto(mHistogramErrorsPerFeeId.get(), false, "colz", "gridy logz");

  // Number of decoding errors, grouped by DE ID and normalized to the number of processed TF
  mHistogramErrorsPerDE = std::make_unique<TH2F>(TString::Format("%sDecodingErrorsPerDE", path.c_str()),
                                                 "Error Type vs. DE ID",
                                                 getNumDE(), 0, getNumDE(),
                                                 getErrorCodesSize(), 0, getErrorCodesSize());
  setYAxisLabels(mHistogramErrorsPerDE.get());
  // setYAxisLabels(mHistogramErrorsPerDE.get());
  addHisto(mHistogramErrorsPerDE.get(), false, "colz", "gridy logz");

  // Number of decoding errors, grouped by chamber ID and normalized to the number of processed TF
  mHistogramErrorsPerChamber = std::make_unique<TH2F>(TString::Format("%sDecodingErrorsPerChamber", path.c_str()),
                                                      "Chamber Number vs. Error Type", 10, 1, 11,
                                                      getErrorCodesSize(), 0, getErrorCodesSize());
  setXAxisLabels(mHistogramErrorsPerChamber.get());
  setYAxisLabels(mHistogramErrorsPerChamber.get());
  addHisto(mHistogramErrorsPerChamber.get(), false, "colz", "gridx gridy logz");
}

//_________________________________________________________________________________________

void DecodingErrorsPlotter::update(TH2F* h)
{
  if (!h) {
    return;
  }

  auto incrementBin = [&](TH2F* h, int bx, int by, float val) {
    auto entries = h->GetBinContent(bx, by);
    h->SetBinContent(bx, by, entries + val);
  };

  mHistogramErrorsPerFeeId->Reset("ICES");
  mHistogramErrorsPerDE->Reset("ICES");
  mHistogramErrorsPerChamber->Reset("ICES");

  // loop over bins in electronics coordinates, and map the channels to the corresponding cathode pads
  int nbinsx = h->GetXaxis()->GetNbins();
  int nbinsy = h->GetYaxis()->GetNbins();
  for (int i = 1; i <= nbinsx; i++) {
    FecId fecId(i - 1);
    // address of the DS board in FEC representation
    uint16_t feeId = fecId.getFeeId();
    uint8_t linkId = fecId.getLinkId();
    uint8_t eLinkId = fecId.getDsAddr();

    uint16_t solarId = -1;
    int de = -1;

    o2::mch::raw::FeeLinkId feeLinkId{ feeId, linkId };

    if (auto opt = mFeeLink2SolarMapper(feeLinkId); opt.has_value()) {
      solarId = opt.value();
    }

    if (solarId >= 0 && solarId <= 1023) {
      o2::mch::raw::DsElecId dsElecId{ solarId, static_cast<uint8_t>(eLinkId / 5), static_cast<uint8_t>(eLinkId % 5) };

      if (auto opt = mElec2DetMapper(dsElecId); opt.has_value()) {
        o2::mch::raw::DsDetId dsDetId = opt.value();
        de = dsDetId.deId();
      }
    }

    for (int j = 1; j <= nbinsy; j++) {
      auto count = h->GetBinContent(i, j);
      incrementBin(mHistogramErrorsPerFeeId.get(), feeId + 1, j, count);
      // mHistogramErrorsPerFeeId->SetBinContent(j, feeId + 1, mHistogramErrorsPerFeeId->GetBinContent(j, feeId + 1) + count);
      if (count > 0) {
        std::cout << fmt::format("FEE={} L={} DS={}  CODE={}   COUNT={}  ENTRIES={}", feeId, linkId, eLinkId, j, count, mHistogramErrorsPerFeeId->GetBinContent(j, feeId + 1)) << std::endl;
      }

      if (de < 0) {
        continue;
      }

      int deId = getDEindex(de);
      if (deId < 0) {
        continue;
      }
      incrementBin(mHistogramErrorsPerDE.get(), deId + 1, j, count);
      // mHistogramErrorsPerDE->SetBinContent(j, deId, mHistogramErrorsPerDE->GetBinContent(j, deId) + count);

      int chamber = de / 100;
      incrementBin(mHistogramErrorsPerChamber.get(), chamber, j, count);
      // mHistogramErrorsPerChamber->SetBinContent(j, chamber, mHistogramErrorsPerChamber->GetBinContent(j, chamber) + count);
    }
  }
}

} // namespace muonchambers
} // namespace quality_control_modules
} // namespace o2
