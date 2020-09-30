///
/// \file   PhysicsTaskPreclusters.cxx
/// \author Barthelemy von Haller
/// \author Piotr Konopka
/// \author Andrea Ferrero
/// \author Sebastien Perrin
///

#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TFile.h>
#include <algorithm>

#include "DPLUtils/DPLRawParser.h"
#include "Headers/RAWDataHeader.h"
#include "MCH/PhysicsTaskPreclusters.h"
#ifdef MCH_HAS_MAPPING_FACTORY
#include "MCHMappingFactory/CreateSegmentation.h"
#endif
#include "MCHMappingInterface/Segmentation.h"
#include "MCHMappingInterface/CathodeSegmentation.h"
#include "MCHMappingSegContour/CathodeSegmentationContours.h"
#include "MCHRawDecoder/PageDecoder.h"
#include "QualityControl/QcInfoLogger.h"

using namespace std;

#define QC_MCH_SAVE_TEMP_ROOTFILE 1

static FILE* flog = NULL;

struct CRUheader {
  uint8_t header_version;
  uint8_t header_size;
  uint16_t block_length;
  uint16_t fee_id;
  uint8_t priority_bit;
  uint8_t reserved_1;
  uint16_t next_packet_offset;
  uint16_t memory_size;
  uint8_t link_id;
  uint8_t packet_counter;
  uint16_t source_id;
  uint32_t hb_orbit;
};

enum decode_state_t {
  DECODE_STATE_UNKNOWN,
  DECODE_STATE_SYNC_FOUND,
  DECODE_STATE_HEADER_FOUND,
  DECODE_STATE_CSIZE_FOUND,
  DECODE_STATE_CTIME_FOUND,
  DECODE_STATE_SAMPLE_FOUND
};

namespace o2
{
namespace quality_control_modules
{
namespace muonchambers
{
PhysicsTaskPreclusters::PhysicsTaskPreclusters() : TaskInterface(), count(1) {}

PhysicsTaskPreclusters::~PhysicsTaskPreclusters() { fclose(flog); }

void PhysicsTaskPreclusters::initialize(o2::framework::InitContext& /*ctx*/)
{
  QcInfoLogger::GetInstance() << "initialize PhysicsTaskPreclusters" << AliceO2::InfoLogger::InfoLogger::endm;
  fprintf(stdout, "initialize PhysicsTaskPreclusters\n");

  mDecoder.initialize();

  flog = stdout;

  for(int de = 0; de < 1100; de++) {
    MeanPseudoeffDE[de] = MeanPseudoeffDECycle[de] = LastPreclBNBDE[de] = NewPreclBNBDE[de] = LastPreclNumDE[de] = NewPreclNumDE[de] = 0;
  }
    
    // 1D histograms for mean pseudoeff per DE (integrated or per elapsed cycle)
    mMeanPseudoeffPerDE = new TH1F("QcMuonChambers_MeanPseudoeff", "Mean Pseudoeff of each DE", 1100, -0.5, 1099.5);
    getObjectsManager()->startPublishing(mMeanPseudoeffPerDE);
    mMeanPseudoeffPerDECycle = new TH1F("QcMuonChambers_MeanPseudoeff_OnCycle", "Mean Pseudoeff of each DE during the cycle", 1100, -0.5, 1099.5);
    getObjectsManager()->startPublishing(mMeanPseudoeffPerDECycle);

  for (int de = 0; de < 1030; de++){
    const o2::mch::mapping::Segmentation* segment = &(o2::mch::mapping::segmentation(de));
    if (segment == nullptr) continue;

    TH1F* h = new TH1F(TString::Format("QcMuonChambers_Cluster_Charge_DE%03d", de),
        TString::Format("QcMuonChambers - cluster charge (DE%03d)", de), 1000, 0, 50000);
    mHistogramClchgDE.insert(make_pair(de, h));
    h = new TH1F(TString::Format("QcMuonChambers_Cluster_Charge_OnCycle_DE%03d", de),
        TString::Format("QcMuonChambers - cluster charge on cycle (DE%03d)", de), 1000, 0, 50000);
    mHistogramClchgDEOnCycle.insert(make_pair(de, h));

    float Xsize = 40 * 5;
    float Xsize2 = Xsize / 2;
    float Ysize = 50;
    float Ysize2 = Ysize / 2;
    float scale = 0.5;

    // Histograms using the XY Mapping

    {
      TH2F* hXY = new TH2F(TString::Format("QcMuonChambers_Preclusters_Number_XY_%03d", de),
          TString::Format("QcMuonChambers - Preclusters Number XY (DE%03d B)", de), Xsize / scale, -Xsize2, Xsize2, Ysize / scale, -Ysize2, Ysize2);
      mHistogramPreclustersXY[0].insert(make_pair(de, hXY));
      hXY = new TH2F(TString::Format("QcMuonChambers_Preclusters_B_XY_%03d", de),
          TString::Format("QcMuonChambers - Preclusters XY (DE%03d B)", de), Xsize / scale, -Xsize2, Xsize2, Ysize / scale, -Ysize2, Ysize2);
      mHistogramPreclustersXY[1].insert(make_pair(de, hXY));
      hXY = new TH2F(TString::Format("QcMuonChambers_Preclusters_NB_XY_%03d", de),
          TString::Format("QcMuonChambers - Preclusters XY (DE%03d NB)", de), Xsize / scale, -Xsize2, Xsize2, Ysize / scale, -Ysize2, Ysize2);
      mHistogramPreclustersXY[2].insert(make_pair(de, hXY));
      hXY = new TH2F(TString::Format("QcMuonChambers_Preclusters_BNB_XY_%03d", de),
          TString::Format("QcMuonChambers - Preclusters XY (DE%03d B+NB)", de), Xsize / scale, -Xsize2, Xsize2, Ysize / scale, -Ysize2, Ysize2);
      mHistogramPreclustersXY[3].insert(make_pair(de, hXY));
      hXY = new TH2F(TString::Format("QcMuonChambers_Pseudoeff_B_XY_%03d", de),
          TString::Format("QcMuonChambers - Pseudo-efficiency XY (DE%03d B)", de), Xsize / scale, -Xsize2, Xsize2, Ysize / scale, -Ysize2, Ysize2);
      mHistogramPseudoeffXY[0].insert(make_pair(de, hXY));
      hXY = new TH2F(TString::Format("QcMuonChambers_Pseudoeff_NB_XY_%03d", de),
          TString::Format("QcMuonChambers - Pseudo-efficiency XY (DE%03d NB)", de), Xsize / scale, -Xsize2, Xsize2, Ysize / scale, -Ysize2, Ysize2);
      mHistogramPseudoeffXY[1].insert(make_pair(de, hXY));
      hXY = new TH2F(TString::Format("QcMuonChambers_Pseudoeff_BNB_XY_%03d", de),
          TString::Format("QcMuonChambers - Pseudo-efficiency XY (DE%03d B+NB)", de), Xsize / scale, -Xsize2, Xsize2, Ysize / scale, -Ysize2, Ysize2);
      mHistogramPseudoeffXY[2].insert(make_pair(de, hXY));
    }
  }

  mHistogramPseudoeff[0] = new GlobalHistogram("QcMuonChambers_Pseudoeff_den", "Pseudo-efficiency");
  mHistogramPseudoeff[0]->init();
  mHistogramPseudoeff[1] = new GlobalHistogram("QcMuonChambers_Pseudoeff", "Pseudo-efficiency");
  mHistogramPseudoeff[1]->init();
  mHistogramPseudoeff[2] = new GlobalHistogram("QcMuonChambers_Pseudoeff_BNB", "Pseudo-efficiency - B+NB");
  mHistogramPseudoeff[2]->init();
}

void PhysicsTaskPreclusters::startOfActivity(Activity& /*activity*/)
{
  QcInfoLogger::GetInstance() << "startOfActivity" << AliceO2::InfoLogger::InfoLogger::endm;
}

void PhysicsTaskPreclusters::startOfCycle()
{
  QcInfoLogger::GetInstance() << "startOfCycle" << AliceO2::InfoLogger::InfoLogger::endm;
}

void PhysicsTaskPreclusters::monitorDataPreclusters(o2::framework::ProcessingContext& ctx)
{
    fprintf(flog, "\n================\nmonitorDataPreClusters\n================\n");
  // get the input preclusters and associated digits
  auto preClusters = ctx.inputs().get<gsl::span<o2::mch::PreCluster>>("preclusters");
  auto digits = ctx.inputs().get<gsl::span<o2::mch::Digit>>("preclusterdigits");

  bool print = false;
  for(auto& p : preClusters) {
    if (!plotPrecluster(p, digits)) {
      print = true;
    }
  }

  if (print) {
    printPreclusters(preClusters, digits);
  }
}


void PhysicsTaskPreclusters::monitorData(o2::framework::ProcessingContext& ctx)
{
  bool preclustersFound = false;
  bool preclusterDigitsFound = false;
  for (auto&& input : ctx.inputs()) {
    if (input.spec->binding == "preclusters") {
      preclustersFound = true;
    }
    if (input.spec->binding == "preclusterdigits") {
      preclusterDigitsFound = true;
    }
  }
  if(preclustersFound && preclusterDigitsFound) {
    monitorDataPreclusters(ctx);
  }
}


//_________________________________________________________________________________________________
static void CoG(gsl::span<const o2::mch::Digit> precluster, double& Xcog, double& Ycog, bool isWide[2])
{
    fprintf(flog, "\n================\nCoG\n================\n");
    
  double xmin = 1E9;
  double ymin = 1E9;
  double xmax = -1E9;
  double ymax = -1E9;
  double charge[] = { 0.0, 0.0 };
  int multiplicity[] = { 0, 0 };
    double padXPos[] = {0,0};
    double padYPos[] = {0,0};
    isWide[0] = isWide[1] = false;

  double x[] = { 0.0, 0.0 };
  double y[] = { 0.0, 0.0 };

  double xsize[] = { 0.0, 0.0 };
  double ysize[] = { 0.0, 0.0 };

  int detid = precluster[0].getDetID();
  const o2::mch::mapping::Segmentation& segment = o2::mch::mapping::segmentation(detid);

  for ( size_t i = 0; i < precluster.size(); ++i ) {
    const o2::mch::Digit& digit = precluster[i];
    int padid = digit.getPadID();

    // position and size of current pad
    double padPosition[2] = {segment.padPositionX(padid), segment.padPositionY(padid)};
    double padSize[2] = {segment.padSizeX(padid), segment.padSizeY(padid)};

    // update of xmin/max et ymin/max
    xmin = std::min(padPosition[0]-0.5*padSize[0],xmin);
    xmax = std::max(padPosition[0]+0.5*padSize[0],xmax);
    ymin = std::min(padPosition[1]-0.5*padSize[1],ymin);
    ymax = std::max(padPosition[1]+0.5*padSize[1],ymax);

    // cathode index
    int cathode = segment.isBendingPad(padid) ? 0 : 1;

    // update of the cluster position, size, charge and multiplicity
    x[cathode] += padPosition[0] * digit.getADC();
    y[cathode] += padPosition[1] * digit.getADC();
    xsize[cathode] += padSize[0];
    ysize[cathode] += padSize[1];
    charge[cathode] += digit.getADC();
      
      if(multiplicity[cathode] == 0){
          if(cathode == 0){
              padXPos[0] = padPosition[0];
              padYPos[0] = padPosition[1];
          }
          if(cathode == 1){
              padXPos[1] = padPosition[0];
              padYPos[1] = padPosition[1];
          }
      }
      else if(multiplicity[cathode] > 0){
          if((cathode == 0) && (padXPos[0] != padPosition[0])){
              isWide[0] = true;
          }
          if((cathode == 1) && (padYPos[1] != padPosition[1])){
              isWide[1] = true;
          }
      }
      
    multiplicity[cathode] += 1;
  }

  // Computation of the CoG coordinates for the two cathodes
  for ( int cathode = 0; cathode < 2; ++cathode ) {
    if ( charge[cathode] != 0 ) {
      x[cathode] /= charge[cathode];
      y[cathode] /= charge[cathode];
    }
    if ( multiplicity[cathode] != 0 ) {
      double sqrtCharge = sqrt(charge[cathode]);
      xsize[cathode] /= (multiplicity[cathode] * sqrtCharge);
      ysize[cathode] /= (multiplicity[cathode] * sqrtCharge);
    } else {
      xsize[cathode] = 1E9;
      ysize[cathode] = 1E9;
    }
  }

  // each CoG coordinate is taken from the cathode with the best precision
  Xcog = ( xsize[0] < xsize[1] ) ? x[0] : x[1];
  Ycog = ( ysize[0] < ysize[1] ) ? y[0] : y[1];
}

void PhysicsTaskPreclusters::printPreclusters(gsl::span<const o2::mch::PreCluster> preClusters, gsl::span<const o2::mch::Digit> digits)
{
  std::cout<<"\n\n============\n";
  for(auto& preCluster : preClusters) {
    // get the digits of this precluster
    auto preClusterDigits = digits.subspan(preCluster.firstDigit, preCluster.nDigits);

    bool cathode[2] = {false, false};
    float chargeSum[2] = {0, 0};
    float chargeMax[2] = {0, 0};

    int detid = preClusterDigits[0].getDetID();
    const o2::mch::mapping::Segmentation& segment = o2::mch::mapping::segmentation(detid);

    for ( size_t i = 0; i < preClusterDigits.size(); ++i ) {
      const o2::mch::Digit& digit = preClusterDigits[i];
      int padid = digit.getPadID();

      // cathode index
      int cid = segment.isBendingPad(padid) ? 0 : 1;
      cathode[cid] = true;
      chargeSum[cid] += digit.getADC();

      if (digit.getADC() > chargeMax[cid]) {
        chargeMax[cid] = digit.getADC();
      }
    }

    double Xcog, Ycog;
      bool isWide[2];
    CoG(preClusterDigits, Xcog, Ycog, isWide);

    std::cout<<"[pre-cluster] charge = "<<chargeSum[0]<<" "<<chargeSum[1]<<"   CoG = "<<Xcog<<" "<<Ycog<<std::endl;
    for (auto& d : preClusterDigits) {
      float X = segment.padPositionX(d.getPadID());
      float Y = segment.padPositionY(d.getPadID());
      bool bend = !segment.isBendingPad(d.getPadID());
      std::cout << fmt::format("  DE {:4d}  PAD {:5d}  ADC {:6d}  TIME ({} {} {:4d})",
          d.getDetID(), d.getPadID(), d.getADC(), d.getTime().orbit, d.getTime().bunchCrossing, d.getTime().sampaTime);
      std::cout << fmt::format("  CATHODE {}  PAD_XY {:+2.2f} , {:+2.2f}", (int)bend, X, Y);
      std::cout << std::endl;
    }
  }
}


bool PhysicsTaskPreclusters::plotPrecluster(const o2::mch::PreCluster& preCluster, gsl::span<const o2::mch::Digit> digits)
{
  // get the digits of this precluster
  auto preClusterDigits = digits.subspan(preCluster.firstDigit, preCluster.nDigits);

  bool cathode[2] = {false, false};
  float chargeSum[2] = {0, 0};
  float chargeMax[2] = {0, 0};
    int padhit[2]= {0,0};

  int detid = preClusterDigits[0].getDetID();
  const o2::mch::mapping::Segmentation& segment = o2::mch::mapping::segmentation(detid);

  for ( size_t i = 0; i < preClusterDigits.size(); ++i ) {
    const o2::mch::Digit& digit = preClusterDigits[i];
    int padid = digit.getPadID();

    // cathode index
    int cid = segment.isBendingPad(padid) ? 0 : 1;
    cathode[cid] = true;
      padhit[cid] += 1;
    chargeSum[cid] += digit.getADC();

    if (digit.getADC() > chargeMax[cid]) {
      chargeMax[cid] = digit.getADC();
    }
  }

  // filter out single-pad clusters
  if (preCluster.nDigits < 2) {
    return true;
  }


  float chargeTot = chargeSum[0] + chargeSum[1];
  auto hCharge = mHistogramClchgDE.find(detid);
  if ((hCharge != mHistogramClchgDE.end()) && (hCharge->second != NULL)) {
    hCharge->second->Fill(chargeTot);
  }
  auto hChargeOnCycle = mHistogramClchgDEOnCycle.find(detid);
  if ((hChargeOnCycle != mHistogramClchgDEOnCycle.end()) && (hChargeOnCycle->second != NULL)) {
        hChargeOnCycle->second->Fill(chargeTot);
  }

  // filter out clusters with small charge, which are likely to be noise
  if ((chargeSum[0]+chargeSum[1]) < 100) {
    return true;
  }
  if ((chargeMax[0] < 100) && (chargeMax[1] < 100)) {
    return true;
  }

  double Xcog, Ycog;
    bool isWide[2];
  CoG(preClusterDigits, Xcog, Ycog, isWide);

    if((cathode[0] && isWide[0]) || (cathode[1] && isWide[1]) || (cathode[0] && cathode[1])){
      auto hXY0 = mHistogramPreclustersXY[0].find(detid);
      if ((hXY0 != mHistogramPreclustersXY[0].end()) && (hXY0->second != NULL)) {
        hXY0->second->Fill(Xcog, Ycog);
      }
    }

  if(cathode[0] && isWide[0]) {
    auto hXY1 = mHistogramPreclustersXY[1].find(detid);
    if ((hXY1 != mHistogramPreclustersXY[1].end()) && (hXY1->second != NULL)) {
      hXY1->second->Fill(Xcog, Ycog);
    }
  }
  if(cathode[1] && isWide[1]) {
    auto hXY1 = mHistogramPreclustersXY[2].find(detid);
    if ((hXY1 != mHistogramPreclustersXY[2].end()) && (hXY1->second != NULL)) {
      hXY1->second->Fill(Xcog, Ycog);
    }
  }
  if(cathode[0] && cathode[1]) {
    auto hXY1 = mHistogramPreclustersXY[3].find(detid);
    if ((hXY1 != mHistogramPreclustersXY[3].end()) && (hXY1->second != NULL)) {
      hXY1->second->Fill(Xcog, Ycog);
    }
  }

  return (cathode[0] && cathode[1]);
}


void PhysicsTaskPreclusters::endOfCycle()
{
  QcInfoLogger::GetInstance() << "endOfCycle" << AliceO2::InfoLogger::InfoLogger::endm;
  for(int de = 100; de <= 1030; de++) {
    for(int i = 0; i < 3; i++) {
      //std::cout<<"DE "<<de<<"  i "<<i<<std::endl;
      auto ih = mHistogramPreclustersXY[i+1].find(de);
      if (ih == mHistogramPreclustersXY[i+1].end()) {
        continue;
      }
      TH2F* hB = ih->second;
      if (!hB) {
        continue;
      }

      ih = mHistogramPreclustersXY[0].find(de);
      if (ih == mHistogramPreclustersXY[0].end()) {
        continue;
      }
      TH2F* hAll = ih->second;
      if (!hAll) {
        continue;
      }

      ih = mHistogramPseudoeffXY[i].find(de);
      if (ih == mHistogramPseudoeffXY[i].end()) {
        continue;
      }
      TH2F* hEff = ih->second;
      if (!hEff) {
        continue;
      }

      hEff->Reset();
      hEff->Add(hB);
      hEff->Divide(hAll);
    }
  }

  mHistogramPseudoeff[0]->add(mHistogramPreclustersXY[0], mHistogramPreclustersXY[0]);
  mHistogramPseudoeff[1]->add(mHistogramPreclustersXY[1], mHistogramPreclustersXY[2]);
  mHistogramPseudoeff[1]->Divide(mHistogramPseudoeff[0]);
  mHistogramPseudoeff[2]->add(mHistogramPreclustersXY[3], mHistogramPreclustersXY[3]);
  mHistogramPseudoeff[2]->Divide(mHistogramPseudoeff[0]);

#ifdef QC_MCH_SAVE_TEMP_ROOTFILE
    TFile f("/tmp/qc.root", "RECREATE");
    
    {
        // Using PseudoeffXY to get the mean pseudoeff per DE on last cycle
         auto hMean = mMeanPseudoeffPerDE;
         auto hMeanCycle = mMeanPseudoeffPerDECycle;
        
    for(int de=0; de<1100; de++){
        auto hnum = mHistogramPreclustersXY[0].find(de);
        auto hBNB = mHistogramPreclustersXY[3].find(de);
        if ((hBNB != mHistogramPreclustersXY[3].end()) && (hBNB->second != NULL) && (hnum != mHistogramPreclustersXY[0].end()) && (hnum->second != NULL)) {
              NewPreclBNBDE[de] = 0;
              NewPreclNumDE[de] = 0;
        //  std::cout << "On va entrer dans la boucle pour lire Elec last cycle" << std::endl;
          for(int binx=1; binx<hBNB->second->GetXaxis()->GetNbins()+1; binx++){
              for(int biny=1; biny<hBNB->second->GetYaxis()->GetNbins()+1; biny++){
                          NewPreclBNBDE[de] += hBNB->second->GetBinContent(binx, biny);
                      }
                  }
          for(int binx=1; binx<hnum->second->GetXaxis()->GetNbins()+1; binx++){
              for(int biny=1; biny<hnum->second->GetYaxis()->GetNbins()+1; biny++){
                      NewPreclNumDE[de] += hnum->second->GetBinContent(binx, biny);
                    }
                }
              
          }
    }
          for(int i=0; i<1100; i++){
              MeanPseudoeffDE[i] = 0;
              MeanPseudoeffDECycle[i] = 0;
              if(NewPreclNumDE[i]>0){
                  MeanPseudoeffDE[i] = NewPreclBNBDE[i]/NewPreclNumDE[i];
              }
              if((NewPreclNumDE[i]-LastPreclNumDE[i]) > 0){
                  MeanPseudoeffDECycle[i] = (NewPreclBNBDE[i]-LastPreclBNBDE[i])/(NewPreclNumDE[i]-LastPreclNumDE[i]);
              }
              hMean->SetBinContent(i+1, MeanPseudoeffDE[i]);
              hMeanCycle->SetBinContent(i+1, MeanPseudoeffDECycle[i]);
              LastPreclBNBDE[i] = NewPreclBNBDE[i];
              LastPreclNumDE[i] = NewPreclNumDE[i];
          }
          hMean->Write();
          hMeanCycle->Write();
        std::cout << "MeanPseudoeff of DE819 since start is: " << MeanPseudoeffDE[819] << std::endl;
          std::cout << "MeanPseudoeff of DE819 in last cycle is: " << MeanPseudoeffDECycle[819] << std::endl;
      }
    
    {
      for(int i = 0; i < 4; i++) {
        for(auto& h2 : mHistogramPreclustersXY[i]) {
          if (h2.second != nullptr) {
            h2.second->Write();
          }
        }
      }
      for(auto& h : mHistogramClchgDE) {
        if (h.second != nullptr) {
          h.second->Write();
        }
      }
      for(auto& h : mHistogramClchgDEOnCycle) {
        if (h.second != nullptr) {
          h.second->Write();
                h.second->Reset();
        }
      }
    }
    
    mHistogramPseudoeff[0]->Write();
    mHistogramPseudoeff[1]->Write();
    mHistogramPseudoeff[2]->Write();

    //f.ls();
    f.Close();
    
#endif
}


void PhysicsTaskPreclusters::endOfActivity(Activity& /*activity*/)
{
  QcInfoLogger::GetInstance() << "endOfActivity" << AliceO2::InfoLogger::InfoLogger::endm;
}

void PhysicsTaskPreclusters::reset()
{
  // clean all the monitor objects here

  QcInfoLogger::GetInstance() << "Reseting the histogram" << AliceO2::InfoLogger::InfoLogger::endm;
}

} // namespace muonchambers
} // namespace quality_control_modules
} // namespace o2

