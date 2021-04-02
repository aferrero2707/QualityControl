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
/// \file   TH1MCHReductor.cxx
/// \author Piotr Konopka, Sebastien Perrin
///

#include <TH1.h>
#include "MCH/TH1MCHReductor.h"
#include <iostream>
#include <string>
#include <regex>

namespace o2::quality_control_modules::muonchambers
{

void* TH1MCHReductor::getBranchAddress()
{
  return &mStats;
}

const char* TH1MCHReductor::getBranchLeafList()
{
  return "occ819/D:occ5I:occ5O:occ6I:occ6O:occ7I:occ7O:occ8I:occ8O:occ9I:occ9O:occ10I:occ10O:entries";
}

void TH1MCHReductor::update(TObject* obj)
{
 auto histo = dynamic_cast<TH1*>(obj);
  if (histo) {
      double mean = 0;
    mStats.entries = histo->GetEntries();
    mStats.indiv_occs.indiv[0] = histo->GetBinContent(819+1);
      
      //5I
      for(int i : detCH5I){
          mean += histo->GetBinContent(i+1);
      }
      mean /= 9;
      mStats.halfch_occs.halfch[0] = mean;
      mean = 0;
      
      //5O
      for(int i : detCH5O){
          mean += histo->GetBinContent(i+1);
      }
      mean /= 9;
      mStats.halfch_occs.halfch[1] = mean;
      mean = 0;
      
      //6I
      for(int i : detCH6I){
          mean += histo->GetBinContent(i+1);
      }
      mean /= 9;
      mStats.halfch_occs.halfch[2] = mean;
      mean = 0;
      
      //6O
      for(int i : detCH6O){
          mean += histo->GetBinContent(i+1);
      }
      mean /= 9;
      mStats.halfch_occs.halfch[3] = mean;
      mean = 0;
      
      //7I
      for(int i : detCH7I){
          mean += histo->GetBinContent(i+1);
      }
      mean /= 13;
      mStats.halfch_occs.halfch[4] = mean;
      mean = 0;
      
      //7O
      for(int i : detCH7O){
          mean += histo->GetBinContent(i+1);
      }
      mean /= 13;
      mStats.halfch_occs.halfch[5] = mean;
      mean = 0;
      
      //8I
      for(int i : detCH8I){
          mean += histo->GetBinContent(i+1);
      }
      mean /= 13;
      mStats.halfch_occs.halfch[6] = mean;
      mean = 0;
      
      //8O
      for(int i : detCH8O){
          mean += histo->GetBinContent(i+1);
      }
      mean /= 13;
      mStats.halfch_occs.halfch[7] = mean;
      mean = 0;
      
      //9I
      for(int i : detCH9I){
          mean += histo->GetBinContent(i+1);
      }
      mean /= 13;
      mStats.halfch_occs.halfch[8] = mean;
      mean = 0;
      
      //9O
      for(int i : detCH9O){
          mean += histo->GetBinContent(i+1);
      }
      mean /= 13;
      mStats.halfch_occs.halfch[9] = mean;
      mean = 0;
      
      //10I
      for(int i : detCH10I){
          mean += histo->GetBinContent(i+1);
      }
      mean /= 13;
      mStats.halfch_occs.halfch[10] = mean;
      mean = 0;
      
      //10O
      for(int i : detCH10O){
          mean += histo->GetBinContent(i+1);
      }
      mean /= 13;
      mStats.halfch_occs.halfch[11] = mean;
      mean = 0;
      
//      for(int i=807; i<=819; i++){
//          mean += histo->GetBinContent(i+1);
//      }
//      mean /= 13;
//      mStats.halfch_occs.halfch[7] = mean;
    
      std::cout << "MeanOccupancy DE819 obtained from reductor " << mStats.indiv_occs.indiv[0] << std::endl;
      std::cout << "MeanOccupancy Ch5I obtained from reductor " << mStats.halfch_occs.halfch[0] << std::endl;
      std::cout << "MeanOccupancy Ch8O obtained from reductor " << mStats.halfch_occs.halfch[7] << std::endl;
      
  }
}

} // namespace o2::quality_control_modules::muonchambers


