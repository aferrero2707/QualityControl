///
/// \file   PhysicsTaskPreclusters.h
/// \author Barthelemy von Haller
/// \author Piotr Konopka
/// \author Andrea Ferrero
///

#ifndef QC_MODULE_MUONCHAMBERS_PHYSICSTASKPRECLUSTERS_H
#define QC_MODULE_MUONCHAMBERS_PHYSICSTASKPRECLUSTERS_H

#include <vector>

#include "QualityControl/TaskInterface.h"
#include "MCH/Mapping.h"
#include "MCH/Decoding.h"
#include "MCH/GlobalHistogram.h"
#include "MCHBase/Digit.h"
#include "MCHBase/PreCluster.h"

class TH1F;
class TH2F;


#define MCH_FFEID_MAX (31*2 + 1)

using namespace o2::quality_control::core;

namespace o2
{
namespace quality_control_modules
{
namespace muonchambers
{

/// \brief Quality Control Task for the analysis of MCH physics data
/// \author Andrea Ferrero
/// \author Sebastien Perrin
class PhysicsTaskPreclusters /*final*/ : public TaskInterface // todo add back the "final" when doxygen is fixed
{
 public:
  /// \brief Constructor
  PhysicsTaskPreclusters();
  /// Destructor
  ~PhysicsTaskPreclusters() override;

  // Definition of the methods for the template method pattern
  void initialize(o2::framework::InitContext& ctx) override;
  void startOfActivity(Activity& activity) override;
  void startOfCycle() override;
  void monitorDataPreclusters(o2::framework::ProcessingContext& ctx);
  void monitorData(o2::framework::ProcessingContext& ctx) override;
  void endOfCycle() override;
  void endOfActivity(Activity& activity) override;
  void reset() override;

  bool plotPrecluster(const o2::mch::PreCluster& preCluster, gsl::span<const o2::mch::Digit> digits);
  void printPreclusters(gsl::span<const o2::mch::PreCluster> preClusters, gsl::span<const o2::mch::Digit> digits);

 private:
  int count;
  Decoder mDecoder;
    
    double MeanPseudoeffDE[1100];
    double MeanPseudoeffDECycle[1100];
    double LastPreclBNBDE[1100];
    double NewPreclBNBDE[1100];
    double LastPreclNumDE[1100];
    double NewPreclNumDE[1100];

  std::vector<std::unique_ptr<mch::Digit>> digits;
    
    // TH1 de la pseudoeff moyenne par DE (intégré ou sur le cycle écoulé)
    TH1F* mMeanPseudoeffPerDE;
    TH1F* mMeanPseudoeffPerDECycle;

    
//  std::map<int, TH2F*> mHistogramMeanNhitsPerDE;
//  std::map<int, TH2F*> mHistogramMeanNorbitsPerDE;

  std::map<int, TH1F*> mHistogramClchgDE;
  std::map<int, TH1F*> mHistogramClchgDEOnCycle;
  std::map<int, TH1F*> mHistogramClsizeDE;

  std::map<int, TH2F*> mHistogramPreclustersXY[4];
  std::map<int, TH2F*> mHistogramPseudoeffXY[3];

  GlobalHistogram* mHistogramPseudoeff[3];
    
};

} // namespace muonchambers
} // namespace quality_control_modules
} // namespace o2

#endif // QC_MODULE_MUONCHAMBERS_PHYSICSDATAPROCESSOR_H

