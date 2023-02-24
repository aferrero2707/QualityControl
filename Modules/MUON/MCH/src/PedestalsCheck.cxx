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
/// \file   PedestalsCheck.cxx
/// \author Andrea Ferrero
///

#include "MCHMappingInterface/Segmentation.h"
#include "MCHMappingSegContour/CathodeSegmentationContours.h"
#include "MCH/PedestalsCheck.h"

#include <fairlogger/Logger.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TList.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TColor.h>

using namespace std;

namespace o2::quality_control_modules::muonchambers
{

void PedestalsCheck::configure()
{
  if (auto param = mCustomParameters.find("MaxBadDE"); param != mCustomParameters.end()) {
    mMaxBadDE = std::stoi(param->second);
  }
  if (auto param = mCustomParameters.find("MaxBadFractionPerDE"); param != mCustomParameters.end()) {
    mMaxBadFractionPerDE = std::stof(param->second);
  }
  if (auto param = mCustomParameters.find("MaxEmptyFractionPerDE"); param != mCustomParameters.end()) {
    mMaxEmptyFractionPerDE = std::stof(param->second);
  }
  if (auto param = mCustomParameters.find("PedestalsPlotScaleMin"); param != mCustomParameters.end()) {
    mPedestalsPlotScaleMin = std::stof(param->second);
  }
  if (auto param = mCustomParameters.find("PedestalsPlotScaleMax"); param != mCustomParameters.end()) {
    mPedestalsPlotScaleMax = std::stof(param->second);
  }
  if (auto param = mCustomParameters.find("NoisePlotScaleMin"); param != mCustomParameters.end()) {
    mNoisePlotScaleMin = std::stof(param->second);
  }
  if (auto param = mCustomParameters.find("NoisePlotScaleMax"); param != mCustomParameters.end()) {
    mNoisePlotScaleMax = std::stof(param->second);
  }
}

Quality PedestalsCheck::check(std::map<std::string, std::shared_ptr<MonitorObject>>* moMap)
{
  //Quality resultEmpty = Quality::Null;
  Quality resultTable = Quality::Null;
  //Quality resultBad = Quality::Null;
  Quality result = Quality::Null;

  mQualityBadChannels = Quality::Null;
  mQualityEmptyChannels = Quality::Null;

  mErrorMessages.clear();

  for (auto& [moName, mo] : *moMap) {
    if (mo->getName().find("BadChannels_Elec") != std::string::npos) {
      auto* h = dynamic_cast<TH2F*>(mo->getObject());
      if (!h) {
        return result;
      }

      std::cout << "BadChannels_Elec->GetEntries(): " << h->GetEntries() << std::endl;

      if (h->GetEntries() == 0) {
        resultTable = Quality::Bad;
        mErrorMessages.emplace_back("Missing Bad Channels Table");
      } else {
        resultTable = Quality::Good;
      }
    }

    if (mo->getName().find("BadChannelsPerDE") != std::string::npos) {
      auto* h = dynamic_cast<TH1F*>(mo->getObject());
      if (!h)
        return result;

      if (h->GetEntries() == 0) {
        mQualityBadChannels = Quality::Medium;
      } else {
        int nbinsx = h->GetXaxis()->GetNbins();
        int nbad = 0;
        for (int i = 1; i <= nbinsx; i++) {
          if (h->GetBinContent(i) > mMaxBadFractionPerDE) {
            nbad += 1;
          }
        }
        mQualityBadChannels = Quality::Good;
        if (nbad > mMaxBadDE) {
          mQualityBadChannels = Quality::Bad;
          mErrorMessages.emplace_back("Too many bad channels");
        }
      }
    }

    if (mo->getName().find("EmptyChannelsPerDE") != std::string::npos) {
      auto* h = dynamic_cast<TH1F*>(mo->getObject());
      if (!h)
        return result;

      if (h->GetEntries() == 0) {
        mQualityEmptyChannels = Quality::Medium;
      } else {
        int nbinsx = h->GetXaxis()->GetNbins();
        int nbad = 0;
        for (int i = 1; i <= nbinsx; i++) {
          if (h->GetBinContent(i) > mMaxEmptyFractionPerDE) {
            nbad += 1;
          }
        }
        mQualityEmptyChannels = Quality::Good;
        if (nbad > mMaxBadDE) {
          mQualityEmptyChannels = Quality::Bad;
          mErrorMessages.emplace_back("Too many empty channels");
        }
      }
    }
  }

  result = Quality::Good;
  if (resultTable == Quality::Bad) {
    result = Quality::Bad;
  }
  if (mQualityBadChannels == Quality::Bad) {
    result = Quality::Bad;
  }
  if (mQualityEmptyChannels == Quality::Bad) {
    result = Quality::Bad;
  }

  if (result == Quality::Good) {
    mErrorMessages.insert(mErrorMessages.begin(), "Quality: GOOD\n");
  }
  if (result == Quality::Medium) {
    mErrorMessages.insert(mErrorMessages.begin(), "Quality: MEDIUM\n");
  }
  if (result == Quality::Bad) {
    mErrorMessages.insert(mErrorMessages.begin(), "Quality: BAD\n");
  }
  if (result == Quality::Null) {
    mErrorMessages.insert(mErrorMessages.begin(), "Quality: NULL\n");
  }

  return result;
}

std::string PedestalsCheck::getAcceptedType() { return "TH1"; }

static void updateTitle(TH1* hist, std::string suffix)
{
  if (!hist) {
    return;
  }
  TString title = hist->GetTitle();
  title.Append(" ");
  title.Append(suffix.c_str());
  hist->SetTitle(title);
}

static std::string getCurrentTime()
{
  time_t t;
  time(&t);

  struct tm* tmp;
  tmp = localtime(&t);

  char timestr[500];
  strftime(timestr, sizeof(timestr), "(%x - %X)", tmp);

  std::string result = timestr;
  return result;
}

void PedestalsCheck::beautify(std::shared_ptr<MonitorObject> mo, Quality checkResult)
{
  auto currentTime = getCurrentTime();
  updateTitle(dynamic_cast<TH1*>(mo->getObject()), currentTime);

  if (mo->getName().find("CheckerMessages") != std::string::npos) {
    auto* canvas = dynamic_cast<TCanvas*>(mo->getObject());
    if (!canvas) {
      return;
    }
    canvas->cd();

    TPaveText* msg = new TPaveText(0.2, 0.3, 0.8, 0.7, "NDC");
    for (auto s : mErrorMessages) {
      msg->AddText(s.c_str());
    }
    if (checkResult == Quality::Good) {
      msg->SetTextColor(kGreen + 2);
    }
    if (checkResult == Quality::Medium) {
      msg->SetTextColor(kOrange);
    }
    if (checkResult == Quality::Bad) {
      msg->SetTextColor(kRed);
    }
    if (checkResult == Quality::Null) {
      msg->SetTextColor(kBlack);
    }
    msg->SetBorderSize(0);
    msg->SetFillColor(kWhite);
    msg->Draw();
  }

  if (mo->getName().find("EmptyChannelsPerDE") != std::string::npos) {
    auto* h = dynamic_cast<TH1F*>(mo->getObject());

    h->SetMinimum(0);
    h->SetMaximum(1.1);

    TLine* delimiter = new TLine(h->GetXaxis()->GetXmin(), mMaxEmptyFractionPerDE, h->GetXaxis()->GetXmax(), mMaxEmptyFractionPerDE);
    delimiter->SetLineColor(kBlack);
    delimiter->SetLineStyle(kDashed);
    h->GetListOfFunctions()->Add(delimiter);

    if (mQualityEmptyChannels == Quality::Good) {
      h->SetFillColor(kGreen);
    } else if (mQualityEmptyChannels == Quality::Bad) {
      h->SetFillColor(kRed);
    } else if (mQualityEmptyChannels == Quality::Medium) {
      h->SetFillColor(kOrange);
    }
    h->SetLineColor(kBlack);
  }

  if (mo->getName().find("BadChannelsPerDE") != std::string::npos) {
    auto* h = dynamic_cast<TH1F*>(mo->getObject());

    h->SetMinimum(0);
    h->SetMaximum(1.1);

    TLine* delimiter = new TLine(h->GetXaxis()->GetXmin(), mMaxBadFractionPerDE, h->GetXaxis()->GetXmax(), mMaxBadFractionPerDE);
    delimiter->SetLineColor(kBlack);
    delimiter->SetLineStyle(kDashed);
    h->GetListOfFunctions()->Add(delimiter);

    if (mQualityBadChannels == Quality::Good) {
      h->SetFillColor(kGreen);
    } else if (mQualityBadChannels == Quality::Bad) {
      h->SetFillColor(kRed);
    } else if (mQualityBadChannels == Quality::Medium) {
      h->SetFillColor(kOrange);
    }
    h->SetLineColor(kBlack);
  }

  if (mo->getName().find("Pedestals_Elec") != std::string::npos) {
    auto* h = dynamic_cast<TH2F*>(mo->getObject());

    h->SetMinimum(mPedestalsPlotScaleMin);
    h->SetMaximum(mPedestalsPlotScaleMax);
  }

  if (mo->getName().find("Noise_Elec") != std::string::npos) {
    auto* h = dynamic_cast<TH2F*>(mo->getObject());

    h->SetMinimum(mNoisePlotScaleMin);
    h->SetMaximum(mNoisePlotScaleMax);
  }

  if ((mo->getName().find("Pedestals_ST12") != std::string::npos) ||
      (mo->getName().find("Pedestals_ST345") != std::string::npos)) {
    auto* h = dynamic_cast<TH2F*>(mo->getObject());
    h->SetMinimum(mPedestalsPlotScaleMin);
    h->SetMaximum(mPedestalsPlotScaleMax);
    h->GetXaxis()->SetTickLength(0.0);
    h->GetXaxis()->SetLabelSize(0.0);
    h->GetYaxis()->SetTickLength(0.0);
    h->GetYaxis()->SetLabelSize(0.0);
  }

  if ((mo->getName().find("Noise_ST12") != std::string::npos) ||
      (mo->getName().find("Noise_ST345") != std::string::npos)) {
    auto* h = dynamic_cast<TH2F*>(mo->getObject());
    h->SetMinimum(mNoisePlotScaleMin);
    h->SetMaximum(mNoisePlotScaleMax);
    h->GetXaxis()->SetTickLength(0.0);
    h->GetXaxis()->SetLabelSize(0.0);
    h->GetYaxis()->SetTickLength(0.0);
    h->GetYaxis()->SetLabelSize(0.0);
  }

  if ((mo->getName().find("BadChannels_ST12") != std::string::npos) ||
      (mo->getName().find("BadChannels_ST345") != std::string::npos)) {
    auto* h = dynamic_cast<TH2F*>(mo->getObject());
    double min = 1;
    int nbinsx = h->GetXaxis()->GetNbins();
    int nbinsy = h->GetYaxis()->GetNbins();
    for (int i = 1; i <= nbinsx; i++) {
      for (int j = 1; j <= nbinsy; j++) {
        auto value = h->GetBinContent(i, j);
        if (value > 0 && value < min) {
          min = value;
        }
      }
    }
    std::cout << "[PedestalsCheck] min " << min << std::endl;
    h->SetMinimum(0.99 * min);
    //h->SetMinimum(0);
    //h->SetMaximum(2);
    h->GetXaxis()->SetTickLength(0.0);
    h->GetXaxis()->SetLabelSize(0.0);
    h->GetYaxis()->SetTickLength(0.0);
    h->GetYaxis()->SetLabelSize(0.0);
  }

  if (mo->getName().find("Pedestals_XY") != std::string::npos) {
    auto* h = dynamic_cast<TH2F*>(mo->getObject());
    h->SetMinimum(mPedestalsPlotScaleMin);
    h->SetMaximum(mPedestalsPlotScaleMax);
    h->GetXaxis()->SetTickLength(0.0);
    h->GetXaxis()->SetLabelSize(0.0);
    h->GetYaxis()->SetTickLength(0.0);
    h->GetYaxis()->SetLabelSize(0.0);
  }

  if (mo->getName().find("Noise_XY") != std::string::npos) {
    auto* h = dynamic_cast<TH2F*>(mo->getObject());
    h->SetMinimum(mNoisePlotScaleMin);
    h->SetMaximum(mNoisePlotScaleMax);
    h->GetXaxis()->SetTickLength(0.0);
    h->GetXaxis()->SetLabelSize(0.0);
    h->GetYaxis()->SetTickLength(0.0);
    h->GetYaxis()->SetLabelSize(0.0);
  }

  if (mo->getName().find("BadChannels_XY") != std::string::npos) {
    auto* h = dynamic_cast<TH2F*>(mo->getObject());
    h->SetMinimum(0);
    h->SetMaximum(1);
    h->GetXaxis()->SetTickLength(0.0);
    h->GetXaxis()->SetLabelSize(0.0);
    h->GetYaxis()->SetTickLength(0.0);
    h->GetYaxis()->SetLabelSize(0.0);
  }
}

} // namespace o2::quality_control_modules::muonchambers
