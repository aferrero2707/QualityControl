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
/// \file   TH2ElecHBPacketsReductor.cxx
/// \author Piotr Konopka, Sebastien Perrin
///

#include <TH1.h>
#include "MCH/TH2ElecHBPacketsReductor.h"
#include "MCH/Helpers.h"
#include "MUONCommon/MergeableTH2Ratio.h"
#include "MCHMappingInterface/Segmentation.h"
#include <iostream>
#include <string>
#include <regex>
#include <gsl/gsl>

using namespace o2::quality_control_modules::muon;

namespace o2::quality_control_modules::muonchambers
{

TH2ElecHBPacketsReductor::TH2ElecHBPacketsReductor() : quality_control::postprocessing::Reductor()
{
  mElec2DetMapper = o2::mch::raw::createElec2DetMapper<o2::mch::raw::ElectronicMapperGenerated>();
  mDet2ElecMapper = o2::mch::raw::createDet2ElecMapper<o2::mch::raw::ElectronicMapperGenerated>();
  mFeeLink2SolarMapper = o2::mch::raw::createFeeLink2SolarMapper<o2::mch::raw::ElectronicMapperGenerated>();
  mSolar2FeeLinkMapper = o2::mch::raw::createSolar2FeeLinkMapper<o2::mch::raw::ElectronicMapperGenerated>();
}

void* TH2ElecHBPacketsReductor::getBranchAddress()
{
  return &mStats;
}

const char* TH2ElecHBPacketsReductor::getBranchLeafList()
{
  return "DE100/D:DE101:DE102:DE103:DE200:DE201:DE202:DE203:DE300:DE301:DE302:DE303:DE400:DE401:DE402:DE403:DE500:DE501:DE502:DE503:DE504:DE505:DE506:DE507:DE508:DE509:DE510:DE511:DE512:DE513:DE514:DE515:DE516:DE517:DE600:DE601:DE602:DE603:DE604:DE605:DE606:DE607:DE608:DE609:DE610:DE611:DE612:DE613:DE614:DE615:DE616:DE617:DE700:DE701:DE702:DE703:DE704:DE705:DE706:DE707:DE708:DE709:DE710:DE711:DE712:DE713:DE714:DE715:DE716:DE717:DE718:DE719:DE720:DE721:DE722:DE723:DE724:DE725:DE800:DE801:DE802:DE803:DE804:DE805:DE806:DE807:DE808:DE809:DE810:DE811:DE812:DE813:DE814:DE815:DE816:DE817:DE818:DE819:DE820:DE821:DE822:DE823:DE824:DE825:DE900:DE901:DE902:DE903:DE904:DE905:DE906:DE907:DE908:DE909:DE910:DE911:DE912:DE913:DE914:DE915:DE916:DE917:DE918:DE919:DE920:DE921:DE922:DE923:DE924:DE925:DE1000:DE1001:DE1002:DE1003:DE1004:DE1005:DE1006:DE1007:DE1008:DE1009:DE1010:DE1011:DE1012:DE1013:DE1014:DE1015:DE1016:DE1017:DE1018:DE1019:DE1020:DE1021:DE1022:DE1023:DE1024:DE1025:CH1:CH2:CH3:CH4:CH5:CH6:CH7:CH8:CH9:CH10:entries";
}

Double_t TH2ElecHBPacketsReductor::getDeValue(int deid)
{
  if (deid < 0 || deid >= sDeNum) {
    return 0;
  }
  return mStats.deValues.values[deid];
}

Double_t TH2ElecHBPacketsReductor::getChamberValue(int chid)
{
  if (chid < 0 || chid >= 10) {
    return 0;
  }
  return mStats.chValues.values[chid];
}

int TH2ElecHBPacketsReductor::checkMapping(uint16_t feeId, uint8_t linkId, uint8_t eLinkId)
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

void TH2ElecHBPacketsReductor::update(TObject* obj)
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
  // the average number of out-of-sync FEC boards
  std::vector<double> deNum(getNumDE());
  std::vector<double> deDen(getNumDE());
  std::fill(deNum.begin(), deNum.end(), 0);
  std::fill(deDen.begin(), deDen.end(), 0);

  std::vector<double> chNum(10);
  std::vector<double> chDen(10);
  std::fill(chNum.begin(), chNum.end(), 0);
  std::fill(chDen.begin(), chDen.end(), 0);

  mStats.entries = h->GetEntries();

  int nbinsx = h->GetXaxis()->GetNbins();
  int nbinsy = h->GetYaxis()->GetNbins();
  double nbad = 0;
  int nboards = 0;
  for (int i = 1; i <= nbinsx; i++) {
    int index = i - 1;
    int ds_addr = (index % 40);
    int link_id = (index / 40) % 12;
    int fee_id = index / (12 * 40);

    int de = checkMapping(fee_id, link_id, ds_addr);
    if (de < 0) {
      continue;
    }

    int deIndex = getDEindex(de);
    if (deIndex < 0) {
      continue;
    }

    int chamber = de / 100 - 1;
    if (chamber < 0 || chamber >= 10) {
      continue;
    }

    deDen[deIndex] += 1;
    chDen[chamber] += 1;

    bool isbad = false;

    int ybinmin = h->GetYaxis()->FindBin(mBcMin);
    int ybinmax = h->GetYaxis()->FindBin(mBcMax);
    auto good = hr->getNum()->Integral(i, i, ybinmin, ybinmax);
    auto total = hr->getNum()->Integral(i, i, 1, nbinsy);
    // add underflow
    total += hr->getNum()->GetBinContent(i, 0);
    // add overflow
    total += hr->getNum()->GetBinContent(i, nbinsy + 1);
    auto bad = total - good;
    if (bad > 0) {
      deNum[deIndex] += 1;
      chNum[chamber] += 1;
    }
  }

  // update the average occupancy values
  for (size_t de = 0; de < deDen.size(); de++) {
    // integrated occupancies
    if (deDen[de] > 0) {
      mStats.deValues.values[de] = deNum[de] / deDen[de];
    } else {
      mStats.deValues.values[de] = 0;
    }
  }
  for (size_t ch = 0; ch < chDen.size(); ch++) {
    // integrated occupancies
    if (chDen[ch] > 0) {
      mStats.chValues.values[ch] = chNum[ch] / chDen[ch];
    } else {
      mStats.chValues.values[ch] = 0;
    }
  }
}

} // namespace o2::quality_control_modules::muonchambers
