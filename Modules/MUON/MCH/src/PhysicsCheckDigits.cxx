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
/// \file   PhysicsCheck.cxx
/// \author Andrea Ferrero, Sebastien Perrin
///

#include "MCHMappingInterface/Segmentation.h"
//#include "MCHMappingSegContour/CathodeSegmentationContours.h"
#include "MCH/PhysicsCheckDigits.h"
#include "MCH/GlobalHistogram.h"

// ROOT
#include <fairlogger/Logger.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TLine.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

namespace o2::quality_control_modules::muonchambers
{

PhysicsCheckDigits::PhysicsCheckDigits()
{
  mPrintLevel = 0;
  minOccupancy = 0.05;
  maxOccupancy = 100.00;

  mElec2DetMapper = o2::mch::raw::createElec2DetMapper<o2::mch::raw::ElectronicMapperGenerated>();
  mDet2ElecMapper = o2::mch::raw::createDet2ElecMapper<o2::mch::raw::ElectronicMapperGenerated>();
  mFeeLink2SolarMapper = o2::mch::raw::createFeeLink2SolarMapper<o2::mch::raw::ElectronicMapperGenerated>();
  mSolar2FeeLinkMapper = o2::mch::raw::createSolar2FeeLinkMapper<o2::mch::raw::ElectronicMapperGenerated>();
}

PhysicsCheckDigits::~PhysicsCheckDigits() {}

void PhysicsCheckDigits::configure(std::string)
{
}

bool PhysicsCheckDigits::checkPadMapping(uint16_t feeId, uint8_t linkId, uint8_t eLinkId, o2::mch::raw::DualSampaChannelId channel)
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
    return false;
  }

  o2::mch::raw::DsElecId dsElecId{ solarId, eLinkId / 5, eLinkId % 5 };

  if (auto opt = mElec2DetMapper(dsElecId); opt.has_value()) {
    o2::mch::raw::DsDetId dsDetId = opt.value();
    dsIddet = dsDetId.dsId();
    deId = dsDetId.deId();
  }

  if (deId < 0 || dsIddet < 0) {
    return false;
  }

  const o2::mch::mapping::Segmentation& segment = o2::mch::mapping::segmentation(deId);
  padId = segment.findPadByFEE(dsIddet, int(channel));

  if (padId < 0) {
    return false;
  }
  return true;
}

Quality PhysicsCheckDigits::check(std::map<std::string, std::shared_ptr<MonitorObject>>* moMap)
{
  Quality result = Quality::Null;

  for (auto& [moName, mo] : *moMap) {

    (void)moName;
    if (mo->getName().find("Occupancy_Elec") != std::string::npos) {
      auto* h = dynamic_cast<TH2F*>(mo->getObject());
      if (!h)
        return result;

      if (h->GetEntries() == 0) {
        result = Quality::Medium;
      } else {
        int nbinsx = h->GetXaxis()->GetNbins();
        int nbinsy = h->GetYaxis()->GetNbins();
        int nbad = 0;
        int npads = 0;
        for (int i = 1; i <= nbinsx; i++) {
          int index = i - 1;
          int ds_addr = (index % 40);
          int link_id = (index / 40) % 12;
          int fee_id = index / (12 * 40);

          for (int j = 1; j <= nbinsy; j++) {
            int chan_addr = j - 1;

            if (!checkPadMapping(fee_id, link_id, ds_addr, chan_addr)) {
              continue;
            }
            npads += 1;

            Float_t occupancy = h->GetBinContent(i, j);
            if (occupancy >= minOccupancy && occupancy <= maxOccupancy) {
              continue;
            }

            nbad += 1;

            if (mPrintLevel >= 1) {
              std::cout << "Channel with unusual occupancy read from OccupancyElec histogrm: fee_id = " << fee_id << ", link_id = " << link_id << ", ds_addr = " << ds_addr << " , chan_addr = " << chan_addr << " with an occupancy of " << occupancy << std::endl;
            }
          }
        }
        if (nbad < 0.1 * npads)
          result = Quality::Good;
        else
          result = Quality::Bad;
      }
    }
  }
  return result;
}

std::string PhysicsCheckDigits::getAcceptedType() { return "TH1"; }

void PhysicsCheckDigits::beautify(std::shared_ptr<MonitorObject> mo, Quality checkResult)
{
  std::cout << "PhysicsCheckDigits::beautify(): MO = " << mo->getName() << std::endl;
  if ((mo->getName().find("Occupancy_Elec") != std::string::npos)) {
    auto* h = dynamic_cast<TH2F*>(mo->getObject());
    h->SetDrawOption("colz");
    h->SetMinimum(0);
    h->SetMaximum(10);
    TPaveText* msg = new TPaveText(0.1, 0.9, 0.9, 0.95, "NDC");
    h->GetListOfFunctions()->Add(msg);
    msg->SetName(Form("%s_msg", mo->GetName()));

    if (checkResult == Quality::Good) {
      msg->Clear();
      msg->AddText("All occupancies within limits: OK!!!");
      msg->SetFillColor(kGreen);

      h->SetFillColor(kGreen);
    } else if (checkResult == Quality::Bad) {
      LOG(info) << "Quality::Bad, setting to red";
      //
      msg->Clear();
      msg->AddText("Call MCH on-call.");
      msg->SetFillColor(kRed);

      h->SetFillColor(kRed);
    } else if (checkResult == Quality::Medium) {
      LOG(info) << "Quality::medium, setting to orange";

      msg->Clear();
      msg->AddText("No entries. If MCH in the run, check MCH TWiki");
      msg->SetFillColor(kYellow);
      h->SetFillColor(kOrange);
    }
    h->SetLineColor(kBlack);
  }

  if ((mo->getName().find("MeanOccupancy") != std::string::npos)) {
    auto* h = dynamic_cast<TH1F*>(mo->getObject());
    h->SetDrawOption("hist");
    h->SetMinimum(0);
    h->SetMaximum(10);
    TPaveText* msg = new TPaveText(0.1, 0.9, 0.9, 0.95, "NDC");
    h->GetListOfFunctions()->Add(msg);
    msg->SetName(Form("%s_msg", mo->GetName()));

    if (checkResult == Quality::Good) {
      msg->Clear();
      msg->AddText("All occupancies within limits: OK!!!");
      msg->SetFillColor(kGreen);

      h->SetFillColor(kGreen);
    } else if (checkResult == Quality::Bad) {
      LOG(info) << "Quality::Bad, setting to red";
      //
      msg->Clear();
      msg->AddText("Call MCH on-call.");
      msg->SetFillColor(kRed);

      h->SetFillColor(kRed);
    } else if (checkResult == Quality::Medium) {
      LOG(info) << "Quality::medium, setting to orange";

      msg->Clear();
      msg->AddText("No entries. If MCH in the run, check MCH TWiki");
      msg->SetFillColor(kYellow);
      h->SetFillColor(kOrange);
    }
    h->SetLineColor(kBlack);

    if (mo->getName().find("MeanOccupancy") != std::string::npos) {
      // Draw lines delimiting the chambers
      for (int de = 200; de < 1100; de += 100) {
        int idx = getDEindex(de) - 1;
        TLine* l = new TLine(idx, 0, idx, h->GetMaximum());
        l->SetLineColor(kRed);
        h->GetListOfFunctions()->Add(l);
      }
    }
  }
}

} // namespace o2::quality_control_modules::muonchambers
