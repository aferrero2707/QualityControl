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
/// \file   GlobalHistogramReductor.h
/// \author Piotr Konopka, Sebastien Perrin
///
#ifndef QUALITYCONTROL_GLOBALHISTOGRAMREDUCTOR_H
#define QUALITYCONTROL_GLOBALHISTOGRAMREDUCTOR_H

#include "QualityControl/Reductor.h"

namespace o2::quality_control_modules::muonchambers
{

/// \brief A Reductor which obtains the most popular characteristics of TH2.
///
/// A Reductor which obtains the most popular characteristics of TH2.
/// It produces a branch in the format: "sumw/D:sumw2:sumwx:sumwx2:sumwy:sumwy2:sumwxy:entries"
class GlobalHistogramReductor : public quality_control::postprocessing::Reductor
{
 public:
  GlobalHistogramReductor() = default;
  ~GlobalHistogramReductor() = default;

  void* getBranchAddress() override;
  const char* getBranchLeafList() override;
  void update(TObject* obj) override;

 private:
  struct {
    union {
      struct {
          Double_t occ819;
      } named;
      Double_t array[1];
    } occs;
    Double_t entries; // is sumw == entries always? maybe not for values which land into the edge bins?
  } mStats;
};

} // namespace o2::quality_control_modules::muonchambers

#endif //QUALITYCONTROL_GLOBALHISTOGRAMREDUCTOR_H


