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
/// \file   OrbitsPlotter.cxx
/// \author Andrea Ferrero
///

#include "MCH/OrbitsPlotter.h"
#include "MCH/Helpers.h"
#include <fmt/format.h>

using namespace o2::mch::raw;

namespace o2
{
namespace quality_control_modules
{
namespace muonchambers
{

OrbitsPlotter::OrbitsPlotter(std::string path)
{
  mElec2DetMapper = createElec2DetMapper<ElectronicMapperGenerated>();
  mFeeLink2SolarMapper = createFeeLink2SolarMapper<ElectronicMapperGenerated>();

  //----------------------------------
  // Orbits histogram
  //----------------------------------
  mHistogramOrbits = std::make_unique<TH2F>(TString::Format("%sDigitOrbitInTFDE", path.c_str()), "Digit orbits vs DE", getNumDE(), 0, getNumDE(), 768, -384, 384);
  addHisto(mHistogramOrbits.get(), false, "colz", "colz");
}

//_________________________________________________________________________________________

void OrbitsPlotter::update(TH2F* h)
{
  if (!h) {
    return;
  }

  mHistogramOrbits->Reset();

  // loop over bins in electronics coordinates, and map the channels to the corresponding cathode pads
  int nbinsx = h->GetXaxis()->GetNbins();
  int nbinsy = h->GetYaxis()->GetNbins();
  for (int i = 1; i <= nbinsx; i++) {
    int index = i - 1;
    // address of the DS board in FEC representation
    uint16_t feeId = index / (12 * 40);
    uint8_t linkId = (index / 40) % 12;
    uint8_t eLinkId = (index % 40);

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

    if (deId < 0) {
      continue;
    }

    int deIndex = getDEindex(deId);

    for (int j = 1; j <= nbinsy; j++) {
      float entries = h->GetBinContent(i, j);
      int orbit = h->GetYaxis()->GetBinCenter(j);
      for (int i = 0; i < entries; i++) {
        mHistogramOrbits->Fill(deIndex, orbit);
      }
    }
  }
}

} // namespace muonchambers
} // namespace quality_control_modules
} // namespace o2
