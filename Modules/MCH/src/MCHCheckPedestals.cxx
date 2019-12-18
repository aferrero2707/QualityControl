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
/// \file   MCHCheckPedestals.cxx
/// \author Andrea Ferrero
///

#include "MCHMappingInterface/Segmentation.h"
#include "MCHMappingSegContour/CathodeSegmentationContours.h"
#include "MCH/MCHCheckPedestals.h"

// ROOT
#include <fairlogger/Logger.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TMath.h>
#include <TPaveText.h>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

namespace o2
{
namespace quality_control_modules
{
namespace muonchambers
{

  MCHCheckPedestals::MCHCheckPedestals() : minMCHpedestal(50.f), maxMCHpedestal(100.f), missing(0)
  {
  }

  MCHCheckPedestals::~MCHCheckPedestals() {}

  void MCHCheckPedestals::configure(std::string)
  {
    // if (AliRecoParam::ConvertIndex(specie) == AliRecoParam::kCosmic) {
    //   minTOFrawTime = 150.; //ns
    //   maxTOFrawTime = 250.; //ns
    // }
  }

  Quality MCHCheckPedestals::check(const MonitorObject* mo)
  {
    //std::cout<<"================================="<<std::endl;
    //std::cout<<"MCHCheckPedestals::check() called"<<std::endl;
    //std::cout<<"================================="<<std::endl;
    Quality result = Quality::Null;

    // const Double_t binWidthTOFrawTime = 2.44;

    // if ((histname.EndsWith("RawsTime")) || (histname.Contains("RawsTime") && suffixTrgCl)) {
    if (mo->getName().find("QcMuonChambers_Pedestals") != std::string::npos) {
        
      int deid;
      sscanf(mo->getName().c_str(), "QcMuonChambers_Pedestals_DE%d", &deid);
        
      auto* h = dynamic_cast<TH2F*>(mo->getObject());
      if( !h ) return result;

      if (h->GetEntries() == 0) {
        result = Quality::Medium;
        // flag = AliQAv1::kWARNING;
      } else {
          
        int nbinsx = h->GetXaxis()->GetNbins();
        int nbinsy = h->GetYaxis()->GetNbins();
        int nbad = 0;
        int nDSmiss = 0;
        std::vector<int> DSIdTestedb, DSIdTestednb;
          
        o2::mch::mapping::Segmentation segment(deid);
        const o2::mch::mapping::CathodeSegmentation& csegmentb = segment.bending();
        const o2::mch::mapping::CathodeSegmentation& csegmentnb = segment.nonBending();
          
          for(int i=0; i< csegmentb.nofDualSampas(); i++){
              DSIdTestedb[i] = csegmentb.dualSampaId(i);
          }
          
          for(int i=0; i< csegmentnb.nofDualSampas(); i++){
              DSIdTestednb[i] = csegmentnb.dualSampaId(i);
          }
          
        for(int i = 1; i <= nbinsx; i++) {
          float sum = 0;
            
          for(int j = 1; j <= nbinsy; j++) {
            Float_t ped = h->GetBinContent(i, j);
            sum += ped;
            if((ped < minMCHpedestal || ped > maxMCHpedestal) && ((std::find(DSIdTestedb.begin(), DSIdTestedb.end(), i) != DSIdTestedb.end()) ||  (std::find(DSIdTestednb.begin(), DSIdTestednb.end(), i) != DSIdTestednb.end()) )) nbad += 1;
          }
            
            if((sum == 0) && ((std::find(DSIdTestedb.begin(), DSIdTestedb.end(), i) != DSIdTestedb.end()) ||  (std::find(DSIdTestednb.begin(), DSIdTestednb.end(), i) != DSIdTestednb.end()) ) ){
                nDSmiss += 1;
                missing.push_back(i);
                int J = 1 + ((i-(i%5))/5);
                int DSnum = i%5;
                LOG(INFO) << " Missing J" << J << " DS" << DSnum;
            }
            
        }
        if( nbad < 1 ) result = Quality::Good;
        else if( nbad != 0 && nDSmiss == 0) result = Quality::Bad;
        else if( nbad != 0 && nDSmiss != 0) result = Quality::BadAndMissing;
      }
    }
    return result;
  }

  std::string MCHCheckPedestals::getAcceptedType() { return "TH1"; }

  void MCHCheckPedestals::beautify(MonitorObject* mo, Quality checkResult)
  {
    //std::cout<<"===================================="<<std::endl;
    //std::cout<<"MCHCheckPedestals::beautify() called"<<std::endl;
    //std::cout<<"===================================="<<std::endl;
    if (mo->getName().find("QcMuonChambers_Pedestals") != std::string::npos) {
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
      } else if (checkResult == Quality::BadAndMissing) {
        LOG(INFO) << "Quality::BadAndMissing, setting to black";
        //
        msg->Clear();
        msg->AddText("There are missing DSs here !!!");
        msg->SetFillColor(kBlack);
        //
        h->SetFillColor(kBlack);
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
    if (mo->getName().find("QcMuonChambers_Noise") != std::string::npos) {
      auto* h = dynamic_cast<TH2F*>(mo->getObject());
      if( !h ) return;
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
          for(unsigned int vi = 0; vi < vertices.size(); vi++) {
            const o2::mch::contour::Vertex<double> v1 = vertices[vi];
            const o2::mch::contour::Vertex<double> v2 = (vi < (vertices.size()-1)) ? vertices[vi+1] : vertices[0];
            TLine* line = new TLine(v1.x, v1.y, v2.x, v2.y);
            h->GetListOfFunctions()->Add(line);
          }
        
            for(unsigned long num = 0; num < missing.size(); num++){
                int dualSampaId = missing[num];
          o2::mch::contour::Contour<double> dscontour = o2::mch::mapping::getDualSampaContour(csegment, dualSampaId);
          std::vector<o2::mch::contour::Vertex<double>> vertices = dscontour.getVertices();
          for(unsigned int vi = 0; vi < vertices.size(); vi++) {
            const o2::mch::contour::Vertex<double> v1 = vertices[vi];
            const o2::mch::contour::Vertex<double> v2 = (vi < (vertices.size()-1)) ? vertices[vi+1] : vertices[0];
            TLine* line = new TLine(v1.x, v1.y, v2.x, v2.y);
              line->SetLineColor(kRed);
            h->GetListOfFunctions()->Add(line);
          }
            }
              
        } catch(std::exception& e) {
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
          for(unsigned int vi = 0; vi < vertices.size(); vi++) {
            const o2::mch::contour::Vertex<double> v1 = vertices[vi];
            const o2::mch::contour::Vertex<double> v2 = (vi < (vertices.size()-1)) ? vertices[vi+1] : vertices[0];
            TLine* line = new TLine(v1.x, v1.y, v2.x, v2.y);
            h->GetListOfFunctions()->Add(line);
            std::cout<<"v1="<<v1.x<<","<<v1.y<<"  v2="<<v2.x<<","<<v2.y<<std::endl;
          }
            
            for(unsigned long num = 0; num < missing.size(); num++){
                int dualSampaId = missing[num];
                o2::mch::contour::Contour<double> dscontour = o2::mch::mapping::getDualSampaContour(csegment, dualSampaId);
                std::vector<o2::mch::contour::Vertex<double>> vertices = dscontour.getVertices();
                for(unsigned int vi = 0; vi < vertices.size(); vi++) {
                  const o2::mch::contour::Vertex<double> v1 = vertices[vi];
                  const o2::mch::contour::Vertex<double> v2 = (vi < (vertices.size()-1)) ? vertices[vi+1] : vertices[0];
                  TLine* line = new TLine(v1.x, v1.y, v2.x, v2.y);
                    line->SetLineColor(kRed);
                  h->GetListOfFunctions()->Add(line);
                }
            }
        } catch(std::exception& e) {
          return;
        }
      }
    }


    //____________________________________________________________________________
    // Pedestals histograms
    if (mo->getName().find("QcMuonChambers_Pedestals_XYb") != std::string::npos) {
      auto* h = dynamic_cast<TH2F*>(mo->getObject());
      if( !h ) return;

      if (mo->getName().find("QcMuonChambers_Pedestals_XYb") != std::string::npos) {
        int deid;
        sscanf(mo->getName().c_str(), "QcMuonChambers_Pedestals_XYb_%d", &deid);

        try {
          o2::mch::mapping::Segmentation segment(deid);
          const o2::mch::mapping::CathodeSegmentation& csegment = segment.bending();
          //std::vector<std::vector<o2::mch::contour::Polygon<double>>>poly = o2::mch::mapping::getPadPolygons(csegment);
          o2::mch::contour::Contour<double> envelop = o2::mch::mapping::getEnvelop(csegment);
          std::vector<o2::mch::contour::Vertex<double>> vertices = envelop.getVertices();
          for(unsigned int vi = 0; vi < vertices.size(); vi++) {
            const o2::mch::contour::Vertex<double> v1 = vertices[vi];
            const o2::mch::contour::Vertex<double> v2 = (vi < (vertices.size()-1)) ? vertices[vi+1] : vertices[0];
            TLine* line = new TLine(v1.x, v1.y, v2.x, v2.y);
            h->GetListOfFunctions()->Add(line);
          }
        } catch(std::exception& e) {
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
          for(unsigned int vi = 0; vi < vertices.size(); vi++) {
            const o2::mch::contour::Vertex<double> v1 = vertices[vi];
            const o2::mch::contour::Vertex<double> v2 = (vi < (vertices.size()-1)) ? vertices[vi+1] : vertices[0];
            TLine* line = new TLine(v1.x, v1.y, v2.x, v2.y);
            h->GetListOfFunctions()->Add(line);
          }
        } catch(std::exception& e) {
          return;
        }
      }
    }
  }

} // namespace muonchambers
} // namespace quality_control_modules
} // namespace o2
