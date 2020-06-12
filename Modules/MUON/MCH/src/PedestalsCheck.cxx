// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
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

// ROOT
#include <fairlogger/Logger.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TMath.h>
#include <TPaveText.h>
#include <iostream>
#include <fstream>
#include <string>


using namespace std;

namespace o2::quality_control_modules::muonchambers
{

PedestalsCheck::PedestalsCheck() : minMCHpedestal(50.f), maxMCHpedestal(100.f)
{
    ofstream myfile;
    noisyfilename = "/tmp/NoisyChannelsPedestalsCheck.txt";
    myfile.open(noisyfilename.c_str(), ios::out | ios::trunc);
    myfile << "fee_id    link_id   ds_addr  chan_addr\n";
    myfile.close();
}

PedestalsCheck::~PedestalsCheck() {}

void PedestalsCheck::configure(std::string)
{
  // if (AliRecoParam::ConvertIndex(specie) == AliRecoParam::kCosmic) {
  //   minTOFrawTime = 150.; //ns
  //   maxTOFrawTime = 250.; //ns
  // }
}

Quality PedestalsCheck::check(std::map<std::string, std::shared_ptr<MonitorObject>>* moMap)
{
  std::cout<<"================================="<<std::endl;
  std::cout<<"PedestalsCheck::check() called"<<std::endl;
  std::cout<<"================================="<<std::endl;
  Quality result = Quality::Null;
  ofstream myfile;

  for (auto& [moName, mo] : *moMap) {

    (void)moName;
    if (mo->getName().find("QcMuonChambers_PedestalsA") != std::string::npos) {
      auto* h = dynamic_cast<TH2F*>(mo->getObject());
      if (!h)
        return result;

      if (h->GetEntries() == 0) {
        result = Quality::Medium;
      } else {
        int nbinsx = 6; //h->GetXaxis()->GetNbins();
        int nbinsy = h->GetYaxis()->GetNbins();
        int nbad = 0;
        for (int i = 1; i <= nbinsx; i++) {
          for (int j = 1; j <= nbinsy; j++) {
            Float_t ped = h->GetBinContent(i, j);
            if (ped < minMCHpedestal || ped > maxMCHpedestal)
              nbad += 1;
          }
        }
        if (nbad < 1)
          result = Quality::Good;
        else
          result = Quality::Bad;
      }
    }
      
        if (mo->getName().find("QcMuonChambers_NoiseA") != std::string::npos) {
            auto* h = dynamic_cast<TH2F*>(mo->getObject());
            if (!h)
                continue;
            if (h->GetEntries() == 0) {
                std::cout << "Nothing to check, Noise histogram is empty" << std::endl;
            } else {
              int nbinsx = h->GetXaxis()->GetNbins();
              int nbinsy = h->GetYaxis()->GetNbins();
              int nnoisy = 0;
              for (int i = 1; i <= nbinsx; i++) {
                for (int j = 1; j <= nbinsy; j++) {
                  Float_t noise = h->GetBinContent(i, j);
                    if (noise > 1.2){
                    nnoisy += 1;
                    int ds_addr =  (i%40)-1;
                    int link_id = ( (i-1-ds_addr) / 40 ) % 12;
                    int fee_id = (i-1-ds_addr-40*link_id) / (12*40);
                    int chan_addr = j-1;
                    
                    std::cout << "Noisy channel read from Noise histogram: fee_id = "<< fee_id << ", link_id = "<< link_id << ", ds_addr = "<< ds_addr << " , chan_addr = " << chan_addr <<" with a noise of " << noise << std::endl;
                    
                    myfile.open(noisyfilename.c_str(), ios_base::out | ios_base::app);
                    myfile << fee_id << "   " << link_id << "    " << ds_addr << "   " << chan_addr << "\n";
                    myfile.close();
                    }
                }
              }
                if (nnoisy == 0){
                    std::cout << "No noisy channel detected" << std::endl;
                }
              else
                std::cout << nnoisy << " noisy channel(s) detected" << std::endl;
            }
        }
      
  }
  return result;
}

std::string PedestalsCheck::getAcceptedType() { return "TH1"; }

void PedestalsCheck::beautify(std::shared_ptr<MonitorObject> mo, Quality checkResult)
{
  //std::cout<<"===================================="<<std::endl;
  //std::cout<<"PedestalsCheck::beautify() called"<<std::endl;
  //std::cout<<"===================================="<<std::endl;
  if (mo->getName().find("QcMuonChambers_PedestalsA") != std::string::npos) {
    auto* h = dynamic_cast<TH2F*>(mo->getObject());
    h->SetDrawOption("colz");
    TPaveText* msg = new TPaveText(0.1, 0.9, 0.9, 0.95, "NDC");
    h->GetListOfFunctions()->Add(msg);
    msg->SetName(Form("%s_msg", mo->GetName()));

    if (checkResult == Quality::Good) {
      msg->Clear();
      msg->AddText("All pedestals within limits: OK!!!");
      msg->SetFillColor(kGreen);
      //
      h->SetFillColor(kGreen);
    } else if (checkResult == Quality::Bad) {
      LOG(INFO) << "Quality::Bad, setting to red";
      //
      msg->Clear();
      msg->AddText("Call MCH on-call.");
      msg->SetFillColor(kRed);
      //
      h->SetFillColor(kRed);
    } else if (checkResult == Quality::Medium) {
      LOG(INFO) << "Quality::medium, setting to orange";
      //
      msg->Clear();
      msg->AddText("No entries. If MCH in the run, check MCH TWiki");
      msg->SetFillColor(kYellow);
      // text->Clear();
      // text->AddText(Form("Raw time peak/total integral = %5.2f%%", peakIntegral * 100. / totIntegral));
      // text->AddText(Form("Mean = %5.2f ns", timeMean));
      // text->AddText(Form("Allowed range: %3.0f-%3.0f ns", minTOFrawTime, maxTOFrawTime));
      // text->AddText("If multiple peaks, check filling scheme");
      // text->AddText("See TOF TWiki.");
      // text->SetFillColor(kYellow);
      //
      h->SetFillColor(kOrange);
    }
    h->SetLineColor(kBlack);
  }

  //____________________________________________________________________________
  // Noise histograms
  if (mo->getName().find("QcMuonChambers_NoiseA") != std::string::npos) {
    auto* h = dynamic_cast<TH2F*>(mo->getObject());
    if (!h)
      return;
    h->SetDrawOption("colz");
    h->SetMaximum(1.5);

    if (mo->getName().find("QcMuonChambers_Noise_XYb") != std::string::npos) {
      int deid;
      sscanf(mo->getName().c_str(), "QcMuonChambers_Noise_XYb_%d", &deid);

      try {
        o2::mch::mapping::Segmentation segment(deid);
        const o2::mch::mapping::CathodeSegmentation& csegment = segment.bending();
        //std::vector<std::vector<o2::mch::contour::Polygon<double>>>poly = o2::mch::mapping::getPadPolygons(csegment);
        o2::mch::contour::Contour<double> envelop = o2::mch::mapping::getEnvelop(csegment);
        std::vector<o2::mch::contour::Vertex<double>> vertices = envelop.getVertices();
        for (unsigned int vi = 0; vi < vertices.size(); vi++) {
          const o2::mch::contour::Vertex<double> v1 = vertices[vi];
          const o2::mch::contour::Vertex<double> v2 = (vi < (vertices.size() - 1)) ? vertices[vi + 1] : vertices[0];
          TLine* line = new TLine(v1.x, v1.y, v2.x, v2.y);
          h->GetListOfFunctions()->Add(line);
        }
      } catch (std::exception& e) {
        return;
      }
    }

    if (mo->getName().find("QcMuonChambers_Noise_XYnb") != std::string::npos) {
      int deid;
      sscanf(mo->getName().c_str(), "QcMuonChambers_Noise_XYnb_%d", &deid);

      try {
        o2::mch::mapping::Segmentation segment(deid);
        const o2::mch::mapping::CathodeSegmentation& csegment = segment.nonBending();
        //std::vector<std::vector<o2::mch::contour::Polygon<double>>>poly = o2::mch::mapping::getPadPolygons(csegment);
        o2::mch::contour::Contour<double> envelop = o2::mch::mapping::getEnvelop(csegment);
        std::vector<o2::mch::contour::Vertex<double>> vertices = envelop.getVertices();
        for (unsigned int vi = 0; vi < vertices.size(); vi++) {
          const o2::mch::contour::Vertex<double> v1 = vertices[vi];
          const o2::mch::contour::Vertex<double> v2 = (vi < (vertices.size() - 1)) ? vertices[vi + 1] : vertices[0];
          TLine* line = new TLine(v1.x, v1.y, v2.x, v2.y);
          h->GetListOfFunctions()->Add(line);
          std::cout << "v1=" << v1.x << "," << v1.y << "  v2=" << v2.x << "," << v2.y << std::endl;
        }
      } catch (std::exception& e) {
        return;
      }
    }
  }

  //____________________________________________________________________________
  // Pedestals histograms
  if (mo->getName().find("QcMuonChambers_Pedestals_XYb") != std::string::npos) {
    auto* h = dynamic_cast<TH2F*>(mo->getObject());
    if (!h)
      return;

    if (mo->getName().find("QcMuonChambers_Pedestals_XYb") != std::string::npos) {
      int deid;
      sscanf(mo->getName().c_str(), "QcMuonChambers_Pedestals_XYb_%d", &deid);

      try {
        o2::mch::mapping::Segmentation segment(deid);
        const o2::mch::mapping::CathodeSegmentation& csegment = segment.bending();
        //std::vector<std::vector<o2::mch::contour::Polygon<double>>>poly = o2::mch::mapping::getPadPolygons(csegment);
        o2::mch::contour::Contour<double> envelop = o2::mch::mapping::getEnvelop(csegment);
        std::vector<o2::mch::contour::Vertex<double>> vertices = envelop.getVertices();
        for (unsigned int vi = 0; vi < vertices.size(); vi++) {
          const o2::mch::contour::Vertex<double> v1 = vertices[vi];
          const o2::mch::contour::Vertex<double> v2 = (vi < (vertices.size() - 1)) ? vertices[vi + 1] : vertices[0];
          TLine* line = new TLine(v1.x, v1.y, v2.x, v2.y);
          h->GetListOfFunctions()->Add(line);
        }
      } catch (std::exception& e) {
        return;
      }
    }

    if (mo->getName().find("QcMuonChambers_Pedestals_XYnb") != std::string::npos) {
      int deid;
      sscanf(mo->getName().c_str(), "QcMuonChambers_Pedestals_XYnb_%d", &deid);

      try {
        o2::mch::mapping::Segmentation segment(deid);
        const o2::mch::mapping::CathodeSegmentation& csegment = segment.nonBending();
        //std::vector<std::vector<o2::mch::contour::Polygon<double>>>poly = o2::mch::mapping::getPadPolygons(csegment);
        o2::mch::contour::Contour<double> envelop = o2::mch::mapping::getEnvelop(csegment);
        std::vector<o2::mch::contour::Vertex<double>> vertices = envelop.getVertices();
        for (unsigned int vi = 0; vi < vertices.size(); vi++) {
          const o2::mch::contour::Vertex<double> v1 = vertices[vi];
          const o2::mch::contour::Vertex<double> v2 = (vi < (vertices.size() - 1)) ? vertices[vi + 1] : vertices[0];
          TLine* line = new TLine(v1.x, v1.y, v2.x, v2.y);
          h->GetListOfFunctions()->Add(line);
        }
      } catch (std::exception& e) {
        return;
      }
    }
  }
}

} // namespace o2::quality_control_modules::muonchambers
