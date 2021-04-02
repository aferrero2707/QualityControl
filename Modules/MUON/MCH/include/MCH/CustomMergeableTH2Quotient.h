// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CustomMergeableTH2Quotient.h
/// \brief An example of a custom TH2Quotient inheriting MergeInterface
///
/// \author Piotr Konopka, piotr.jan.konopka@cern.ch

#ifndef O2_CUSTOMMERGEABLETH2QUOTIENT_H
#define O2_CUSTOMMERGEABLETH2QUOTIENT_H

#include <sstream>
#include <iostream>
#include <TObject.h>
#include <TH2.h>
#include <TList.h>
#include "Mergers/MergeInterface.h"

using namespace std;
namespace o2::quality_control_modules::muonchambers
{

class CustomMergeableTH2Quotient : public TH2F, public o2::mergers::MergeInterface
{
 public:
  CustomMergeableTH2Quotient() = default;
  CustomMergeableTH2Quotient(TH2F* histonum, TH2F* histoden)
    : TH2F(*histonum), o2::mergers::MergeInterface(), mhistoNum(histonum), mhistoDen(histoden)
    {
        cout << "Object created" <<endl;
        this->SetTitle("Quotient");
        this->Divide(mhistoNum,mhistoDen);
    }

  ~CustomMergeableTH2Quotient() override = default;

   // using MergeInterface::merge;
    void merge(MergeInterface* const other) override
  {
      mhistoNum->Add(dynamic_cast<const CustomMergeableTH2Quotient* const>(other)->getNum());

      mhistoDen->Add(dynamic_cast<const CustomMergeableTH2Quotient* const>(other)->getDen());
      
      this->Divide(mhistoNum,mhistoDen);
  }
          
  TH2F* getNum() const
  {
      return mhistoNum;
  }
    
  TH2F* getDen() const
  {
      return mhistoDen;
  }


    
  private:

  TH2F* mhistoNum;
  TH2F* mhistoDen;

  ClassDefOverride(CustomMergeableTH2Quotient, 1);
};

} // namespace o2::mergers

#endif //O2_CUSTOMMERGEABLETH2QUOTIENT_H

