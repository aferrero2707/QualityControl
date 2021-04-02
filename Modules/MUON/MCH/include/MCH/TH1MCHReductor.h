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
/// \file   TH1MCHReductor.h
/// \author Piotr Konopka, Sebastien Perrin
///
#ifndef QUALITYCONTROL_TH1MCHREDUCTOR_H
#define QUALITYCONTROL_TH1MCHREDUCTOR_H

#include "QualityControl/Reductor.h"

namespace o2::quality_control_modules::muonchambers
{

/// \brief A Reductor which obtains the most popular characteristics of TH1.
///
/// A Reductor which obtains the most popular characteristics of TH1.

class TH1MCHReductor : public quality_control::postprocessing::Reductor
{
 public:
  TH1MCHReductor() = default;
  ~TH1MCHReductor() = default;

  void* getBranchAddress() override;
  const char* getBranchLeafList() override;
  void update(TObject* obj) override;

 private:
 struct {
   union {
     struct {
         Double_t occ819;
     } named;
     Double_t indiv[1];
   } indiv_occs;
     union {
       struct {
           Double_t occ5I;
           Double_t occ5O;
           Double_t occ6I;
           Double_t occ6O;
           Double_t occ7I;
           Double_t occ7O;
           Double_t occ8I;
           Double_t occ8O;
           Double_t occ9I;
           Double_t occ9O;
           Double_t occ10I;
           Double_t occ10O;
       } named;
         Double_t halfch[12];
 } halfch_occs;
   Double_t entries;
 } mStats;
    int detCH5I[9] = {500,501,502,503,504,514,515,516,517};
    int detCH5O[9] = {505,506,507,508,509,510,511,512,513};
    int detCH6I[9] = {600,601,602,603,604,614,615,616,617};
    int detCH6O[9] = {605,606,607,608,609,610,611,612,613};
    int detCH7I[13] = {700,701,702,703,704,705,706,720,721,722,723,724,725};
    int detCH7O[13] = {707,708,709,710,711,712,713,714,715,716,717,718,719};
    int detCH8I[13] = {800,801,802,803,804,805,806,820,821,822,823,824,825};
    int detCH8O[13] = {807,808,809,810,811,812,813,814,815,816,817,818,819};
    int detCH9I[13] = {900,901,902,903,904,905,906,920,921,922,923,924,925};
    int detCH9O[13] = {907,908,909,910,911,912,913,914,915,916,917,918,919};
    int detCH10I[13] = {1000,1001,1002,1003,1004,1005,1006,1020,1021,1022,1023,1024,1025};
    int detCH10O[13] = {1007,1008,1009,1010,1011,1012,1013,1014,1015,1016,1017,1018,1019};
};

} // namespace o2::quality_control_modules::muonchambers

#endif //QUALITYCONTROL_TH1mchREDUCTOR_H


