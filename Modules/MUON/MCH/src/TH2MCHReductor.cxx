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
/// \file   TH2MCHReductor.cxx
/// \author Piotr Konopka, Sebastien Perrin
///

#include <TH2.h>
#include "MCH/TH2MCHReductor.h"
#include <iostream>
#include <string>
#include <regex>
#include "MCHMappingInterface/Segmentation.h"
#include "MCHMappingInterface/CathodeSegmentation.h"
#ifdef MCH_HAS_MAPPING_FACTORY
#include "MCHMappingFactory/CreateSegmentation.h"
#endif
#include "MCHMappingSegContour/CathodeSegmentationContours.h"

namespace o2::quality_control_modules::muonchambers
{

void* TH2MCHReductor::getBranchAddress()
{
  return &mStats;
}

const char* TH2MCHReductor::getBranchLeafList()
{
  return "sumw/D:sumw2:sumwx:sumwx2:sumwy:sumwy2:sumwxy:entries:mean_nonull:mean_all";
}

void TH2MCHReductor::update(TObject* obj)
{
  auto histo = dynamic_cast<TH2*>(obj);
  if (histo) {
      
      //Get DE number from HistoName
      
      std::string input = histo->GetName();
      std::string output = std::regex_replace(
          input,
          std::regex("[^0-9]*([0-9]+).*"),
          std::string("$1")
          );
//      std::cout << "Name input histo " << input << std::endl;
//      std::cout << "Int output " << output << std::endl;
      
      int de = stoi(output);
      
      const o2::mch::mapping::Segmentation& segment = o2::mch::mapping::segmentation(de);
      const o2::mch::mapping::CathodeSegmentation& csegmentB = segment.bending();
      o2::mch::contour::BBox<double> bboxB = o2::mch::mapping::getBBox(csegmentB);
        
        double xmin = bboxB.xmin();
        double xmax = bboxB.xmax();
        double ymin = bboxB.ymin();
        double ymax = bboxB.ymax();
      
      double sum_bins = 0;
      int n_bins_all_DE = 0;
      int n_bins_nonull_DE = 0;
    histo->GetStats(mStats.sums.array);
    mStats.entries = histo->GetEntries();
      
      int nbinsx = histo->GetXaxis()->GetNbins();
      int nbinsy = histo->GetYaxis()->GetNbins();
      for(int i=0; i < nbinsx; i++){
          for(int j=0; j < nbinsy; j++){
              if((histo->GetXaxis()->GetBinCenter(i) < xmax && histo->GetXaxis()->GetBinCenter(i) > xmin) && (histo->GetYaxis()->GetBinCenter(j) < ymax && histo->GetYaxis()->GetBinCenter(j) > ymin)){
                  
                  if(histo->GetBinContent(i,j) != 0){
                      sum_bins += histo->GetBinContent(i,j);
                      n_bins_nonull_DE++;
                  }
                  
                  n_bins_all_DE++;
              }
          }
      }
      mStats.mean_nonull = sum_bins/n_bins_nonull_DE;
      
      mStats.mean_all = sum_bins/n_bins_all_DE;
      
      std::cout << "sum_bins " << sum_bins <<std::endl;
      std::cout << "n_bins_nonull_DE " << n_bins_nonull_DE <<std::endl;
      std::cout << "n_bins_all_DE " << n_bins_all_DE <<std::endl;
  }
}

} // namespace o2::quality_control_modules::muonchambers

