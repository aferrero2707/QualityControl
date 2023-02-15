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
/// \file   TH2ElecMapReductor.cxx
/// \author Piotr Konopka, Sebastien Perrin
///

#include <TH1.h>
#include "MCH/TH2ElecMapReductor.h"
#include "MCH/Helpers.h"
#include "MUONCommon/MergeableTH2Ratio.h"
#include "MCHMappingInterface/Segmentation.h"
#include <iostream>
#include <string>
#include <regex>
#include <gsl/gsl>
#include <limits>

using namespace o2::quality_control_modules::muon;

namespace o2::quality_control_modules::muonchambers
{

TH2ElecMapReductor::TH2ElecMapReductor(float min, float max)
  : quality_control::postprocessing::Reductor(),
    mMin(min),
    mMax(max)
{
  mElec2DetMapper = o2::mch::raw::createElec2DetMapper<o2::mch::raw::ElectronicMapperGenerated>();
  mDet2ElecMapper = o2::mch::raw::createDet2ElecMapper<o2::mch::raw::ElectronicMapperGenerated>();
  mFeeLink2SolarMapper = o2::mch::raw::createFeeLink2SolarMapper<o2::mch::raw::ElectronicMapperGenerated>();
  mSolar2FeeLinkMapper = o2::mch::raw::createSolar2FeeLinkMapper<o2::mch::raw::ElectronicMapperGenerated>();
}

void* TH2ElecMapReductor::getBranchAddress()
{
  return nullptr;
}

const char* TH2ElecMapReductor::getBranchLeafList()
{
  return "";
}

float TH2ElecMapReductor::getDeValue(int deid, int cathode)
{
  if (deid < 0 || deid >= sDeNum) {
    return 0;
  }
  if (cathode < 0 || cathode > 2) {
    return 0;
  }
  return deValues[cathode][deid];
}

float TH2ElecMapReductor::getChamberValue(int chid)
{
  if (chid < 0 || chid >= 10) {
    return 0;
  }
  return chValues[chid];
}

int TH2ElecMapReductor::getNumPads(int deid, int cathode)
{
  if (deid < 0 || deid >= sDeNum) {
    return 0;
  }
  if (cathode < 0 || cathode > 1) {
    return 0;
  }
  return deNumPads[cathode][deid];
}

int TH2ElecMapReductor::getNumPadsBad(int deid, int cathode)
{
  if (deid < 0 || deid >= sDeNum) {
    return 0;
  }
  if (cathode < 0 || cathode > 1) {
    return 0;
  }
  return deNumPadsBad[cathode][deid];
}

int TH2ElecMapReductor::getNumPadsNoStat(int deid, int cathode)
{
  if (deid < 0 || deid >= sDeNum) {
    return 0;
  }
  if (cathode < 0 || cathode > 1) {
    return 0;
  }
  return deNumPadsNoStat[cathode][deid];
}

int TH2ElecMapReductor::checkPadMapping(uint16_t feeId, uint8_t linkId, uint8_t eLinkId, o2::mch::raw::DualSampaChannelId channel, int& cid)
{
  uint16_t solarId = -1;
  int deId = -1;
  int dsIddet = -1;
  int padId = -1;
  cid = -1;

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

  const o2::mch::mapping::Segmentation& segment = o2::mch::mapping::segmentation(deId);
  padId = segment.findPadByFEE(dsIddet, int(channel));

  if (padId < 0) {
    return -1;
  }
  cid = segment.isBendingPad(padId) ? 0 : 1;
  return deId;
}

void TH2ElecMapReductor::update(TObject* obj)
{
  if (sDeNum != getNumDE()) {
    std::cout << "wrong sDeNum" << std::endl;
    return;
  }

  auto h = dynamic_cast<TH2*>(obj);
  if (!h) {
    std::cout << "cannot cast to TH2F" << std::endl;
    return;
  }

  auto* hr = dynamic_cast<MergeableTH2Ratio*>(obj);
  if (!hr) {
    std::cout << "cannot cast to MergeableTH2Ratio" << std::endl;
    return;
  }

  // cumulative numerators and denominators for the computation of
  // the average value over detection elements and chambers
  std::vector<double> deValueNum(getNumDE());
  std::vector<double> deValueDen(getNumDE());
  std::fill(deValueNum.begin(), deValueNum.end(), 0);
  std::fill(deValueDen.begin(), deValueDen.end(), 0);

  std::vector<double> deValueNumB(getNumDE());
  std::vector<double> deValueDenB(getNumDE());
  std::fill(deValueNumB.begin(), deValueNumB.end(), 0);
  std::fill(deValueDenB.begin(), deValueDenB.end(), 0);

  std::vector<double> deValueNumNB(getNumDE());
  std::vector<double> deValueDenNB(getNumDE());
  std::fill(deValueNumNB.begin(), deValueNumNB.end(), 0);
  std::fill(deValueDenNB.begin(), deValueDenNB.end(), 0);

  std::vector<double> chValueNum(10);
  std::vector<double> chValueDen(10);
  std::fill(chValueNum.begin(), chValueNum.end(), 0);
  std::fill(chValueDen.begin(), chValueDen.end(), 0);

  for (int i = 0; i < 2; i++) {
    for (int de = 0; de < sDeNum; de++) {
      deNumPads[i][de] = 0;
      deNumPadsBad[i][de] = 0;
      deNumPadsNoStat[i][de] = 0;
    }
  }

  double nOrbits = 0;

  entries = h->GetEntries();

  int nbinsx = h->GetXaxis()->GetNbins();
  int nbinsy = h->GetYaxis()->GetNbins();
  int ngood = 0;
  int nPads = 0;
  for (int i = 1; i <= nbinsx; i++) {
    int index = i - 1;
    int ds_addr = (index % 40);
    int link_id = (index / 40) % 12;
    int fee_id = index / (12 * 40);

    for (int j = 1; j <= nbinsy; j++) {
      int chan_addr = j - 1;

      int cid;
      int de = checkPadMapping(fee_id, link_id, ds_addr, chan_addr, cid);
      if (de < 0 || cid < 0) {
        continue;
      }

      int deIndex = getDEindex(de);
      if (deIndex < 0) {
        continue;
      }

      deNumPads[cid][deIndex] += 1;

      Float_t stat = hr->getDen()->GetBinContent(i, j);
      if (stat == 0) {
        deNumPadsNoStat[cid][deIndex] += 1;
        continue;
      }

      Float_t value = h->GetBinContent(i, j);
      if (value <= mMin || value >= mMax) {
        deNumPadsBad[cid][deIndex] += 1;
      }

      nOrbits += stat;
      nPads += 1;

      deValueNum[deIndex] += value;
      deValueDen[deIndex] += 1;
      if (cid == 0) {
        deValueNumB[deIndex] += value;
        deValueDenB[deIndex] += 1;
      } else {
        deValueNumNB[deIndex] += value;
        deValueDenNB[deIndex] += 1;
      }

      int chamber = de / 100 - 1;
      if (chamber >= 0 && chamber < 10) {
        chValueNum[chamber] += value;
        chValueDen[chamber] += 1;
      }
    }
  }

  // update the average values
  for (size_t de = 0; de < deValueDenB.size(); de++) {
    // integrated values
    if (deValueDenB[de] > 0) {
      deValues[0][de] = deValueNumB[de] / deValueDenB[de];
    } else {
      deValues[0][de] = 0;
    }
  }
  for (size_t de = 0; de < deValueDenNB.size(); de++) {
    // integrated values
    if (deValueDenNB[de] > 0) {
      deValues[1][de] = deValueNumNB[de] / deValueDenNB[de];
    } else {
      deValues[1][de] = 0;
    }
  }
  for (size_t de = 0; de < deValueDen.size(); de++) {
    // integrated values
    if (deValueDen[de] > 0) {
      deValues[2][de] = deValueNum[de] / deValueDen[de];
    } else {
      deValues[2][de] = 0;
    }
  }
  for (size_t ch = 0; ch < chValueDen.size(); ch++) {
    // integrated values
    if (chValueDen[ch] > 0) {
      chValues[ch] = chValueNum[ch] / chValueDen[ch];
    } else {
      chValues[ch] = 0;
    }
  }

  if (nPads > 0) {
    meanOrbits = nOrbits / nPads;
  } else {
    meanOrbits = 0;
  }
}

} // namespace o2::quality_control_modules::muonchambers
