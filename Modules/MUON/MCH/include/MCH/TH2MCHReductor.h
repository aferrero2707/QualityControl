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
/// \file   TH2MCHReductor.h
/// \author Piotr Konopka, Sebastien Perrin
///
#ifndef QUALITYCONTROL_TH2MCHREDUCTOR_H
#define QUALITYCONTROL_TH2MCHREDUCTOR_H

#include "QualityControl/Reductor.h"

namespace o2::quality_control_modules::muonchambers
{

/// \brief A Reductor which obtains the most popular characteristics of TH2.
///
/// A Reductor which obtains the most popular characteristics of TH2.
/// It produces a branch in the format: "sumw/D:sumw2:sumwx:sumwx2:sumwy:sumwy2:sumwxy:entries"
class TH2MCHReductor : public quality_control::postprocessing::Reductor
{
 public:
  TH2MCHReductor() = default;
  ~TH2MCHReductor() = default;

  void* getBranchAddress() override;
  const char* getBranchLeafList() override;
  void update(TObject* obj) override;

 private:
  struct {
    union {
      struct {
        Double_t sumw;
        Double_t sumw2;
        Double_t sumwx;
        Double_t sumwx2;
        Double_t sumwy;
        Double_t sumwy2;
        Double_t sumwxy;
      } named;
      Double_t array[7];
    } sums;
    Double_t entries; // is sumw == entries always? maybe not for values which land into the edge bins?
    Double_t mean_nonull;
    Double_t mean_all;
  } mStats;
};

} // namespace o2::quality_control_modules::muonchambers

#endif //QUALITYCONTROL_TH2mchREDUCTOR_H

