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
/// \file   GlobalHistogramReductor.cxx
/// \author Piotr Konopka, Sebastien Perrin
///

#include <TH2.h>
#include "MCH/GlobalHistogramReductor.h"
#include <iostream>

namespace o2::quality_control_modules::muonchambers
{

void* GlobalHistogramReductor::getBranchAddress()
{
  return &mStats;
}

const char* GlobalHistogramReductor::getBranchLeafList()
{
  return "occ819/D:entries";
}

void GlobalHistogramReductor::update(TObject* obj)
{
  auto histo = dynamic_cast<TH2*>(obj);
  if (histo) {

    mStats.entries = histo->GetEntries();
    mStats.occs.array[0] = histo->GetBinContent(853,406);
    
      std::cout << "Mean_nonull occ DE819 " << mStats.occs.array[0] << std::endl;
      
  }
}

} // namespace o2::quality_control_modules::muonchambers


