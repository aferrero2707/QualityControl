///
/// \file   PedestalsTask.h
/// \author Andrea Ferrero
///

#ifndef QC_MODULE_MUONCHAMBERS_PEDESTALSTASK_H
#define QC_MODULE_MUONCHAMBERS_PEDESTALSTASK_H

#include "QualityControl/TaskInterface.h"
#include "MCH/Mapping.h"
#include "MCH/Decoding.h"
#include "MCHBase/Digit.h"
#include "MCH/GlobalHistogram.h"
#include "Framework/DataRef.h"

class TH1F;
class TH2F;

#define MCH_MAX_DE 1100

using namespace o2::quality_control::core;

namespace o2::quality_control_modules::muonchambers
{

/// \brief Quality Control Task for the analysis of MCH pedestal data
/// \author Andrea Ferrero
/// \author Sebastien Perrin
class PedestalsTask final : public TaskInterface
{
 public:
  /// \brief Constructor
  PedestalsTask();
  /// Destructor
  ~PedestalsTask() override;

  // Definition of the methods for the template method pattern
  void initialize(o2::framework::InitContext& ctx) override;
  void startOfActivity(Activity& activity) override;
  void startOfCycle() override;
  void monitorDataTF(o2::framework::ProcessingContext& ctx);
  void monitorDataReadout(const o2::framework::DataRef& input);
  void monitorDataDigits(const o2::framework::DataRef& input);
  void monitorData(o2::framework::ProcessingContext& ctx) override;
  void endOfCycle() override;
  void endOfActivity(Activity& activity) override;
  void reset() override;

 private:
  Decoder mDecoder;
  uint64_t nhits[MCH_NCRU][24][40][64];
  double pedestal[MCH_NCRU][24][40][64];
  double noise[MCH_NCRU][24][40][64];

  int64_t bxc[MCH_NCRU][24][40][64];
  double bxcMean[MCH_NCRU][24];

  //Matrices [de][padid], stated an upper value for de# and padid#

  uint64_t nhitsDigits[1100][1500];
  double pedestalDigits[1100][1500];
  double noiseDigits[1100][1500];

  MapCRU mMapCRU[MCH_NCRU];
  TH2F* mHistogramPedestals;
  TH2F* mHistogramNoise;
  TH2F* mHistogramBunchCrossing;

  std::vector<int> DEs;
  //MapFEC mMapFEC;
  std::map<int, TH2F*> mHistogramPedestalsDE;
  std::map<int, TH2F*> mHistogramNoiseDE;
  std::map<int, TH1F*> mHistogramDeltaDE[2];
  std::map<int, TH2F*> mHistogramPedestalsXY[2];
  std::map<int, TH2F*> mHistogramNoiseXY[2];
  std::map<int, TH2F*> mHistogramDsIDXY[2];

  std::map<int, TH1F*> mHistogramNoiseDistributionDE[5][2];

  // mFlipDE[DE][0] -> horizontal flipping
  // mFlipDE[DE][1] -> vertical flipping
  bool mFlipDE[MCH_MAX_DE][2];

  GlobalHistogram* mHistogramPedestalsMCH;
  GlobalHistogram* mHistogramNoiseMCH;

  int mPrintLevel;

  void fill_noise_distributions();
  void save_histograms();
};

} // namespace o2::quality_control_modules::muonchambers

#endif // QC_MODULE_MUONCHAMBERS_PEDESTALSTASK_H
