///
/// \file   PedestalsTask.h
/// \author Andrea Ferrero
///

#ifndef QC_MODULE_MUONCHAMBERS_PEDESTALSTASK_H
#define QC_MODULE_MUONCHAMBERS_PEDESTALSTASK_H

#include "QualityControl/TaskInterface.h"
#include "MCHRawElecMap/Mapper.h"
#include "MCH/GlobalHistogram.h"
#include "Framework/DataRef.h"
#include "MCHCalibration/MCHChannelCalibrator.h"

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
  struct PadCoordinates
  {
    double padX;
    double padSizeX;
    double padY;
    double padSizeY;
    int padId;
    int cathode;
    int density;
  };
  o2::mch::raw::Solar2FeeLinkMapper mSolar2FeeLinkMapper;
  o2::mch::raw::Elec2DetMapper mElec2DetMapper;

  /// helper class that performs the actual computation of the pedestals from the input digits
  o2::mch::calibration::PedestalProcessor mPedestalProcessor;

  //int64_t bxc[MCH_NCRU][24][40][64];
  double bxcMean;

  TH2F* mHistogramPedestals;
  TH2F* mHistogramNoise;
  TH2F* mHistogramBunchCrossing;

  std::vector<int> DEs;
  //MapFEC mMapFEC;
  std::map<int, TH2F*> mHistogramPedestalsDE;
  std::map<int, TH2F*> mHistogramPedMinDE;
  std::map<int, TH2F*> mHistogramPedMaxDE;
  std::map<int, TH2F*> mHistogramNoiseDE;
  std::map<int, TH1F*> mHistogramDeltaDE[5][2];
  std::map<int, TH2F*> mHistogramPedestalsXY[2];
  std::map<int, TH2F*> mHistogramNoiseXY[2];
  std::map<int, TH2F*> mHistogramPedMinXY[2];
  std::map<int, TH2F*> mHistogramPedMaxXY[2];
  std::map<int, TH2F*> mHistogramDeltaXY[2];
  std::map<int, TH2F*> mHistogramDsIDXY[2];

  std::map<int, TH1F*> mHistogramNoiseDistributionDE[5][2];

  GlobalHistogram* mHistogramPedestalsMCH;
  GlobalHistogram* mHistogramNoiseMCH;

  // mFlipDE[DE][0] -> horizontal flipping
  // mFlipDE[DE][1] -> vertical flipping
  bool mFlipDE[MCH_MAX_DE][2];

  int mPrintLevel;

  bool getDeMapping(uint16_t solarID, uint8_t dsID, int& deId, int& dsIddet);
  bool getPadMapping(uint16_t solarID, uint8_t dsID, uint8_t channel, int& deId, int& dsIddet, PadCoordinates& coords);
  bool getFeeMapping(uint16_t solarID, int& feeId, int& linkId);
  int getFeeBin(int feeId, int linkId, int dsId)
  {
    return (feeId * 12 * 40 + (linkId % 12) * 40 + dsId + 1);
  }
  int getPadDensity(float padSizeX, float padSizeY)
  {
    float szmax = padSizeX;
    if (szmax < padSizeY) {
      szmax = padSizeY;
    }

    int szid = 0;
    if (fabs(szmax - 2.5) < 0.001) {
      szid = 1;
    } else if (fabs(szmax - 5.0) < 0.001) {
      szid = 2;
    } else if (fabs(szmax - 10.0) < 0.001) {
      szid = 3;
    }

    return szid;
  }

  void monitorDataDigits(o2::framework::ProcessingContext& ctx);
  void monitorDataPedestals(o2::framework::ProcessingContext& ctx);

  void PlotBunchCrossingAndDeltas(uint16_t solarID, uint8_t dsID, uint8_t channel, uint32_t bxc, int16_t pedMin, int16_t pedMax);
  void PlotPedestal(uint16_t solarID, uint8_t dsID, uint8_t channel, double mean, double rms);
  void PlotPedestalDE(uint16_t solarID, uint8_t dsID, uint8_t channel, double mean, double rms);
  void fill_noise_distributions();
  void save_histograms();
};

} // namespace o2::quality_control_modules::muonchambers

#endif // QC_MODULE_MUONCHAMBERS_PEDESTALSTASK_H
