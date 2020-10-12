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
  return "occ500/D:occ501:occ502:occ503:occ504:occ505:occ506:occ507:occ508:occ509:occ510:occ511:occ512:occ513:occ514:occ515:occ516:occ517:occ600:occ601:occ602:occ603:occ604:occ605:occ606:occ607:occ608:occ609:occ610:occ611:occ612:occ613:occ614:occ615:occ616:occ617:occ700:occ701:occ702:occ703:occ704:occ705:occ706:occ707:occ708:occ709:occ710:occ711:occ712:occ713:occ714:occ715:occ716:occ717:occ718:occ719:occ720:occ721:occ722:occ723:occ724:occ725:occ800:occ801:occ802:occ803:occ804:occ805:occ806:occ807:occ808:occ809:occ810:occ811:occ812:occ813:occ814:occ815:occ816:occ817:occ818:occ819:occ820:occ821:occ822:occ823:occ824:occ825:occ900:occ901:occ902:occ903:occ904:occ905:occ906:occ907:occ908:occ909:occ910:occ911:occ912:occ913:occ914:occ915:occ916:occ917:occ918:occ919:occ920:occ921:occ922:occ923:occ924:occ925:occ1000:occ1001:occ1002:occ1003:occ1004:occ1005:occ1006:occ1007:occ1008:occ1009:occ1010:occ1011:occ1012:occ1013:occ1014:occ1015:occ1016:occ1017:occ1018:occ1019:occ1020:occ1021:occ1022:occ1023:occ1024:occ1025:occ5I:occ5O:occ6I:occ6O:occ7I:occ7O:occ8I:occ8O:occ9I:occ9O:occ10I:occ10O:entries";
}

void TH1MCHReductor::update(TObject* obj)
{
  auto histo = dynamic_cast<TH1*>(obj);
  if (histo) {
    double mean = 0;
    mStats.entries = histo->GetEntries();
    for (int i = 0; i <= 17; i++) {
      mStats.indiv_occs.indiv[i] = histo->GetBinContent(500 + i + 1);
    }
    for (int i = 18; i <= 35; i++) {
      mStats.indiv_occs.indiv[i] = histo->GetBinContent(600 + (i - 18) + 1);
    }
    for (int i = 36; i <= 61; i++) {
      mStats.indiv_occs.indiv[i] = histo->GetBinContent(700 + (i - 36) + 1);
    }
    for (int i = 62; i <= 87; i++) {
      mStats.indiv_occs.indiv[i] = histo->GetBinContent(800 + (i - 62) + 1);
    }
    for (int i = 88; i <= 113; i++) {
      mStats.indiv_occs.indiv[i] = histo->GetBinContent(900 + (i - 88) + 1);
    }
    for (int i = 114; i <= 139; i++) {
      mStats.indiv_occs.indiv[i] = histo->GetBinContent(1000 + (i - 114) + 1);
    }

    //5I
    for (int i : detCH5I) {
      mean += histo->GetBinContent(i + 1);
    }
    mean /= 9;
    mStats.halfch_occs.halfch[0] = mean;
    mean = 0;

    //5O
    for (int i : detCH5O) {
      mean += histo->GetBinContent(i + 1);
    }
    mean /= 9;
    mStats.halfch_occs.halfch[1] = mean;
    mean = 0;

    //6I
    for (int i : detCH6I) {
      mean += histo->GetBinContent(i + 1);
    }
    mean /= 9;
    mStats.halfch_occs.halfch[2] = mean;
    mean = 0;

    //6O
    for (int i : detCH6O) {
      mean += histo->GetBinContent(i + 1);
    }
    mean /= 9;
    mStats.halfch_occs.halfch[3] = mean;
    mean = 0;

    //7I
    for (int i : detCH7I) {
      mean += histo->GetBinContent(i + 1);
    }
    mean /= 13;
    mStats.halfch_occs.halfch[4] = mean;
    mean = 0;

    //7O
    for (int i : detCH7O) {
      mean += histo->GetBinContent(i + 1);
    }
    mean /= 13;
    mStats.halfch_occs.halfch[5] = mean;
    mean = 0;

    //8I
    for (int i : detCH8I) {
      mean += histo->GetBinContent(i + 1);
    }
    mean /= 13;
    mStats.halfch_occs.halfch[6] = mean;
    mean = 0;

    //8O
    for (int i : detCH8O) {
      mean += histo->GetBinContent(i + 1);
    }
    mean /= 13;
    mStats.halfch_occs.halfch[7] = mean;
    mean = 0;

    //9I
    for (int i : detCH9I) {
      mean += histo->GetBinContent(i + 1);
    }
    mean /= 13;
    mStats.halfch_occs.halfch[8] = mean;
    mean = 0;

    //9O
    for (int i : detCH9O) {
      mean += histo->GetBinContent(i + 1);
    }
    mean /= 13;
    mStats.halfch_occs.halfch[9] = mean;
    mean = 0;

    //10I
    for (int i : detCH10I) {
      mean += histo->GetBinContent(i + 1);
    }
    mean /= 13;
    mStats.halfch_occs.halfch[10] = mean;
    mean = 0;

    //10O
    for (int i : detCH10O) {
      mean += histo->GetBinContent(i + 1);
    }
    mean /= 13;
    mStats.halfch_occs.halfch[11] = mean;
    mean = 0;

    //      for(int i=807; i<=819; i++){
    //          mean += histo->GetBinContent(i+1);
    //      }
    //      mean /= 13;
    //      mStats.halfch_occs.halfch[7] = mean;

    std::cout << "Value DE500 obtained from reductor " << mStats.indiv_occs.indiv[0] << std::endl;
    std::cout << "Value Ch5I obtained from reductor " << mStats.halfch_occs.halfch[0] << std::endl;
    std::cout << "Value Ch8O obtained from reductor " << mStats.halfch_occs.halfch[7] << std::endl;
  }
}

} // namespace o2::quality_control_modules::muonchambers
