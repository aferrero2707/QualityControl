// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file TestMerging_cxx
/// \brief Tests of QC Merger
///
/// \author Sebastien Perrin


#include "Mergers/MergerAlgorithm.h"
#include "MCH/CustomMergeableTH2Quotient.h"

#include <TObjArray.h>
#include <TObjString.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THn.h>
#include <TTree.h>
#include <THnSparse.h>
#include <TF1.h>
#include <TGraph.h>
#include <TProfile.h>
#include <TCanvas.h>

//using namespace o2::framework;
using namespace o2::mergers;
  
int main(){
    
    {
    TH2F* obj1n(NULL);
    TH2F* obj2n(NULL);
    TH2F* obj1d(NULL);
    TH2F* obj2d(NULL);
    
    obj1n = new TH2F("obj1n",
                      "obj1n",
                      10,0,10,5,0,5);
    obj1n->SetXTitle("X obj1n");
    obj1n->SetYTitle("Y obj1n");
        obj1d = new TH2F("obj1d",
                          "obj1d",
                          10,0,10,5,0,5);
        obj1d->SetXTitle("X obj1d");
        obj1d->SetYTitle("Y obj1d");
    
    obj2n = new TH2F("obj2n",
                      "obj2n",
                      10,0,10,5,0,5);
    obj2n->SetXTitle("X obj2n");
    obj2n->SetYTitle("Y obj2n");
        
        obj2d = new TH2F("obj2d",
                          "obj2d",
                          10,0,10,5,0,5);
        obj2d->SetXTitle("X obj2d");
        obj2d->SetYTitle("Y obj2d");
    
    obj1n->Fill(1.5,1.5,5);
        obj1d->Fill(1.5,1.5,10);
    obj2n->Fill(4.5,4.5,2);
        obj2d->Fill(4.5,4.5,2);
        obj2n->Fill(1.5,1.5,1);
        obj2d->Fill(1.5,1.5,8);
        
        std::cout << "Will create objects : "<<std::endl;
  
    auto* target = new o2::quality_control_modules::muonchambers::CustomMergeableTH2Quotient(obj1n, obj1d);
    auto* other1 = new o2::quality_control_modules::muonchambers::CustomMergeableTH2Quotient(obj2n, obj2d);
        
        std::cout << "Before merge, ratio obj1 = 5/10 : " << target->GetBinContent(2,2) <<std::endl;
        
        std::cout << "Before merge, ratio obj2 5,5 = 2/2 : " << other1->GetBinContent(5,5) <<std::endl;
        std::cout << "Before merge, ratio obj2 2,2 = 1/8 : " << other1->GetBinContent(2,2) <<std::endl;
        
    target->merge(other1);
    
    std::cout << "result ratio 5,5 = 2/2 = " << target->GetBinContent(5,5) <<std::endl;
        
        std::cout << "result ratio 2,2 = 6/18 = " << target->GetBinContent(2,2) <<std::endl;

    delete other1;
    delete target;
    }
    
    
  
    return 0;
}
