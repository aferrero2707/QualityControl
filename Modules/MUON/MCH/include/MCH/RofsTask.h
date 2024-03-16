// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \file   RofsTask.h
/// \author Andrea Ferrero
/// \author Sebastien Perrin
///

#ifndef QC_MODULE_MCH_ROFSTASK_H
#define QC_MODULE_MCH_ROFSTASK_H

#include <TH2.h>
#include <gsl/span>

#include "DataFormatsMCH/Digit.h"
#include "DataFormatsMCH/ROFRecord.h"
#include "MCHDigitFiltering/DigitFilter.h"
#include "QualityControl/QcInfoLogger.h"
#include "QualityControl/TaskInterface.h"

namespace o2::quality_control_modules::muonchambers
{

/// \brief Quality Control Task for the analysis of raw data decoding errors
/// \author Andrea Ferrero
/// \author Sebastien Perrin
class RofsTask /*final*/ : public o2::quality_control::core::TaskInterface
{
 public:
  /// \brief Constructor
  RofsTask() = default;
  /// Destructor
  ~RofsTask() override = default;

  void initialize(o2::framework::InitContext& ctx) override;
  void startOfActivity(const o2::quality_control::core::Activity& activity) override;
  void startOfCycle() override;
  void monitorData(o2::framework::ProcessingContext& ctx) override;
  void endOfCycle() override;
  void endOfActivity(const o2::quality_control::core::Activity& activity) override;
  void reset() override;

 private:
  template <typename T>
  void publishObject(std::shared_ptr<T> histo, std::string drawOption, bool statBox);

  void plotROF(const o2::mch::ROFRecord& rof, gsl::span<const o2::mch::Digit> digits);

  o2::mch::DigitFilter mIsSignalDigit; ///< functor to select signal-like digits

  std::shared_ptr<TH1F> mHistRofSize;             ///< number of digits per ROF
  std::shared_ptr<TH1F> mHistRofSize_Signal;      ///< number of signal-like digits per ROF
  std::shared_ptr<TH1F> mHistRofNStations;        ///< number of stations per ROF
  std::shared_ptr<TH1F> mHistRofNStations_Signal; ///< number of stations per ROF from signal-like digits
  std::shared_ptr<TH1F> mHistRofTime;             ///< average ROF time in orbit
  std::shared_ptr<TH1F> mHistRofTime_Signal;      ///< average ROF time in orbit from signal-like digits
  std::shared_ptr<TH1F> mHistRofWidth;            ///< ROF width in BC

  std::vector<TH1*> mAllHistograms;
};

} // namespace o2::quality_control_modules::muonchambers

#endif // QC_MODULE_MCH_ROFSTASK_H
