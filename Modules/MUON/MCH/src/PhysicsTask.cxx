///
/// \file   PhysicsTask.cxx
/// \author Barthelemy von Haller
/// \author Piotr Konopka
/// \author Andrea Ferrero
///

#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <algorithm>

#include "Headers/RAWDataHeader.h"
#include "DPLUtils/DPLRawParser.h"
#include "MCHRawDecoder/PageDecoder.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHMappingInterface/CathodeSegmentation.h"
#ifdef MCH_HAS_MAPPING_FACTORY
#include "MCHMappingFactory/CreateSegmentation.h"
#endif
//#include "MCHPreClustering/PreClusterFinder.h"
#include "QualityControl/QcInfoLogger.h"
#include "MCH/PhysicsTask.h"

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
  //uint16_t cru_id;
  //uint8_t dummy1;
  //uint64_t dummy2;
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
PhysicsTask::PhysicsTask() : TaskInterface(), count(1) {}

PhysicsTask::~PhysicsTask() { fclose(flog); }

void PhysicsTask::initialize(o2::framework::InitContext& /*ctx*/)
{
  QcInfoLogger::GetInstance() << "initialize PhysicsTask" << AliceO2::InfoLogger::InfoLogger::endm;
  fprintf(stdout, "initialize PhysicsTask\n");

  mDecoder.initialize();

  mPrintLevel = 0;

  flog = stdout; //fopen("/root/qc.log", "w");
  fprintf(stdout, "PhysicsTask initialization finished\n");

  mHistogramADCamplitudeVsSize = new TH2F(TString::Format("QcMuonChambers_ADCamplitude_vs_Size"),
      TString::Format("QcMuonChambers - ADC amplitude vs. size"), 500, 0, 500, 5000, 0, 5000);


  uint32_t dsid;
  for (int cruid = 0; cruid < 3; cruid++) {
    QcInfoLogger::GetInstance() << "JE SUIS ENTRÉ DANS LA BOUCLE CRUID " << cruid << AliceO2::InfoLogger::InfoLogger::endm;
    for (int linkid = 0; linkid < 24; linkid++) {
      QcInfoLogger::GetInstance() << "JE SUIS ENTRÉ DANS LA BOUCLE LINKID " << linkid << AliceO2::InfoLogger::InfoLogger::endm;

      {
        int index = 24 * cruid + linkid;
        mHistogramNhits[index] = new TH2F(TString::Format("QcMuonChambers_NHits_CRU%01d_LINK%02d", cruid, linkid),
            TString::Format("QcMuonChambers - Number of hits (CRU link %02d)", index), 40, 0, 40, 64, 0, 64);
        //mHistogramPedestals->SetDrawOption("col");
        //getObjectsManager()->startPublishing(mHistogramNhits[index]);

        mHistogramADCamplitude[index] = new TH1F(TString::Format("QcMuonChambers_ADC_Amplitude_CRU%01d_LINK%02d", cruid, linkid),
            TString::Format("QcMuonChambers - ADC amplitude (CRU link %02d)", index), 5000, 0, 5000);
        //mHistogramPedestals->SetDrawOption("col");
        //getObjectsManager()->startPublishing(mHistogramADCamplitude[index]);
      }

      int32_t link_id = mDecoder.getMapCRU(cruid, linkid);
      QcInfoLogger::GetInstance() << "  LINK_ID " << link_id << AliceO2::InfoLogger::InfoLogger::endm;
      if (link_id == -1)
        continue;
      for (int ds_addr = 0; ds_addr < 40; ds_addr++) {
        QcInfoLogger::GetInstance() << "JE SUIS ENTRÉ DANS LA BOUCLE DS_ADDR " << ds_addr << AliceO2::InfoLogger::InfoLogger::endm;
        uint32_t de;
        int32_t result = mDecoder.getMapFEC(link_id, ds_addr, de, dsid);
        QcInfoLogger::GetInstance() << "C'EST LA LIGNE APRÈS LE GETMAPFEC, DE " << de << "  RESULT " << result << AliceO2::InfoLogger::InfoLogger::endm;
        if(result < 0) continue;

        if (std::find(DEs.begin(), DEs.end(), de) == DEs.end()) {
          DEs.push_back(de);

          TH1F* h = new TH1F(TString::Format("QcMuonChambers_ADCamplitude_DE%03d", de),
              TString::Format("QcMuonChambers - ADC amplitude (DE%03d)", de), 5000, 0, 5000);
          mHistogramADCamplitudeDE.insert(make_pair(de, h));
          //getObjectsManager()->startPublishing(h);

          float Xsize = 40 * 5;
          float Xsize2 = Xsize / 2;
          float Ysize = 50;
          float Ysize2 = Ysize / 2;

          TH2F* h2 = new TH2F(TString::Format("QcMuonChambers_Nhits_DE%03d", de),
              TString::Format("QcMuonChambers - Number of hits (DE%03d)", de), Xsize * 2, -Xsize2, Xsize2, Ysize * 2, -Ysize2, Ysize2);
          mHistogramNhitsDE.insert(make_pair(de, h2));
          getObjectsManager()->startPublishing(h2);
          h2 = new TH2F(TString::Format("QcMuonChambers_Nhits_HighAmpl_DE%03d", de),
              TString::Format("QcMuonChambers - Number of hits for Csum>500 (DE%03d)", de), Xsize * 2, -Xsize2, Xsize2, Ysize * 2, -Ysize2, Ysize2);
          mHistogramNhitsHighAmplDE.insert(make_pair(de, h2));
          //getObjectsManager()->startPublishing(h2);
        }
      }
    }
  }

  for(int de = 1; de <= 1030; de++) {
    const o2::mch::mapping::Segmentation* segment = &(o2::mch::mapping::segmentation(de));
    if (segment == nullptr) continue;

    TH1F* h = new TH1F(TString::Format("QcMuonChambers_Cluster_Charge_DE%03d", de),
        TString::Format("QcMuonChambers - cluster charge (DE%03d)", de), 1000, 0, 50000);
    mHistogramClchgDE.insert(make_pair(de, h));

    float Xsize = 40 * 5;
    float Xsize2 = Xsize / 2;
    float Ysize = 50;
    float Ysize2 = Ysize / 2;
    float scale = 0.5;
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

void PhysicsTask::startOfActivity(Activity& /*activity*/)
{
  QcInfoLogger::GetInstance() << "startOfActivity" << AliceO2::InfoLogger::InfoLogger::endm;
}

void PhysicsTask::startOfCycle()
{
  QcInfoLogger::GetInstance() << "startOfCycle" << AliceO2::InfoLogger::InfoLogger::endm;
}

void PhysicsTask::monitorDataReadout(o2::framework::ProcessingContext& ctx)
{
  // Reset the containers
  mDecoder.reset();

  // For some reason the input selection doesn't work, to be investigated...
  o2::framework::DPLRawParser parser(ctx.inputs() /*, o2::framework::select("readout:MCH/RAWDATA")*/);

  for (auto it = parser.begin(), end = parser.end(); it != end; ++it) {
    // retrieving RDH v4
    auto const* rdh = it.get_if<o2::header::RAWDataHeaderV4>();
    if (!rdh)
      continue;
    // retrieving the raw pointer of the page
    auto const* raw = it.raw();
    // size of payload
    size_t payloadSize = it.size();
    if (payloadSize == 0)
      continue;

    if (mPrintLevel >= 1) {
      std::cout<<"payloadSize: "<<payloadSize<<std::endl;
      //std::cout<<"raw:     "<<(void*)raw<<std::endl;
      //std::cout<<"payload: "<<(void*)payload<<std::endl;
    }

    // Run the decoder on the CRU buffer
    mDecoder.processData((const char*)raw, (size_t)(payloadSize + sizeof(o2::header::RAWDataHeaderV4)));
  }

  std::vector<SampaHit>& hits = mDecoder.getHits();
  if (mPrintLevel >= 1)
    fprintf(flog, "hits.size()=%d\n", (int)hits.size());
  for (uint32_t i = 0; i < hits.size(); i++) {
    //continue;
    SampaHit& hit = hits[i];
    if (mPrintLevel >= 1)
      fprintf(stdout, "hit[%d]: link_id=%d, ds_addr=%d, chan_addr=%d\n",
          i, hit.link_id, hit.ds_addr, hit.chan_addr);
    if (hit.link_id >= 24 || hit.ds_addr >= 40 || hit.chan_addr >= 64) {
      fprintf(stdout, "hit[%d]: link_id=%d, ds_addr=%d, chan_addr=%d\n",
          i, hit.link_id, hit.ds_addr, hit.chan_addr);
      continue;
    }
    //if( hit.csum > 500 ) {
    mHistogramNhits[hit.link_id]->Fill(hit.ds_addr, hit.chan_addr);
    mHistogramADCamplitude[hit.link_id]->Fill(hit.csum);
    //}
  }

  std::vector<o2::mch::Digit>& digits = mDecoder.getDigits();
  if (mPrintLevel >= 1)
    fprintf(flog, "digits.size()=%d\n", (int)digits.size());
  for (uint32_t i = 0; i < digits.size(); i++) {
    o2::mch::Digit& digit = digits[i];
    plotDigit(digit);
  }
}


void PhysicsTask::monitorDataDigits(const o2::framework::DataRef& input)
{
  if (input.spec->binding != "digits")
    return;

  const auto* header = o2::header::get<header::DataHeader*>(input.header);
  if (mPrintLevel >= 1)
    fprintf(flog, "Header: %p\n", (void*)header);
  if (!header)
    return;
  //QcInfoLogger::GetInstance() << "payloadSize: " << header->payloadSize << AliceO2::InfoLogger::InfoLogger::endm;
  if (mPrintLevel >= 1)
    fprintf(flog, "payloadSize: %d\n", (int)header->payloadSize);
  if (mPrintLevel >= 1)
    fprintf(flog, "payload: %p\n", input.payload);

  std::vector<o2::mch::Digit> digits{ 0 };
  o2::mch::Digit* digitsBuffer = NULL;
  digitsBuffer = (o2::mch::Digit*)input.payload;
  size_t ndigits = ((size_t)header->payloadSize / sizeof(o2::mch::Digit));

  if (mPrintLevel >= 1)
    std::cout << "There are " << ndigits << " digits in the payload" << std::endl;

  o2::mch::Digit* ptr = (o2::mch::Digit*)digitsBuffer;
  for (size_t di = 0; di < ndigits; di++) {
    digits.push_back(*ptr);
    ptr += 1;
  }

  for (auto& d : digits) {
    if (mPrintLevel >= 1) {
      std::cout << fmt::format("  DE {:4d}  PAD {:5d}  ADC {:6d}  TIME ({} {} {:4d})",
          d.getDetID(), d.getPadID(), d.getADC(), d.getTime().orbit, d.getTime().bunchCrossing, d.getTime().sampaTime);
      std::cout << std::endl;
    }
    plotDigit(d);
  }
}


void PhysicsTask::monitorDataPreclusters(o2::framework::ProcessingContext& ctx)
{
  // get the input preclusters and associated digits
  auto preClusters = ctx.inputs().get<gsl::span<o2::mch::PreCluster>>("preclusters");
  auto digits = ctx.inputs().get<gsl::span<o2::mch::Digit>>("preclusterdigits");

  if (mPrintLevel >= 1) {
  std::cout<<"preClusters.size()="<<preClusters.size()<<std::endl;
  }

  //checkPreclusters(preClusters, digits);

  bool print = false;
  for(auto& p : preClusters) {
    if (!plotPrecluster(p, digits)) {
      print = true;
    }
  }

  if (print) {
    printPreclusters(preClusters, digits);
  }
  //for (uint32_t i = 0; i < digits.size(); i++) {
  //  o2::mch::Digit& digit = digits[i];
  //  plotDigit(digit);
  //}
}


void PhysicsTask::monitorData(o2::framework::ProcessingContext& ctx)
{
#ifdef QC_MCH_SAVE_TEMP_ROOTFILE_
  if ((count % 100) == 0) {

    TFile f("/tmp/qc.root", "RECREATE");
    for (int i = 0; i < 3 * 24; i++) {
      //mHistogramNhits[i]->Write();
      //mHistogramADCamplitude[i]->Write();
    }
    mHistogramADCamplitudeVsSize->Write();
    //std::cout<<"mHistogramADCamplitudeDE.size() = "<<mHistogramADCamplitudeDE.size()<<"  DEs.size()="<<DEs.size()<<std::endl;
    int nbDEs = DEs.size();
    for (int elem = 0; elem < nbDEs; elem++) {
      int de = DEs[elem];
      //std::cout<<"  de="<<de<<std::endl;
      {
        auto h = mHistogramADCamplitudeDE.find(de);
        if ((h != mHistogramADCamplitudeDE.end()) && (h->second != NULL)) {
          h->second->Write();
        }
      }
      {
        auto h = mHistogramNhitsDE.find(de);
        if ((h != mHistogramNhitsDE.end()) && (h->second != NULL)) {
          h->second->Write();
        }
      }
      {
        auto h = mHistogramNhitsHighAmplDE.find(de);
        if ((h != mHistogramNhitsHighAmplDE.end()) && (h->second != NULL)) {
          h->second->Write();
        }
      }
    }
    {
      for(int i = 0; i < 4; i++) {
        for(auto& h2 : mHistogramPreclustersXY[i]) {
          if (h2.second != nullptr) {
            h2.second->Write();
          }
        }
        //auto h2 = mHistogramPreclustersXY[i].find(de);
        //if ((h2 != mHistogramNhitsDE.end()) && (h2->second != NULL)) {
        //  h2->second->Write();
        //}
      }
    }

    f.ls();
    f.Close();
    if (mPrintLevel == 0) {
      printf("count: %d\n", count);
    }
  }
  if (mPrintLevel >= 1) {
    printf("count: %d\n", count);
  }
  count += 1;
#endif

  if (mPrintLevel >= 1) {
    QcInfoLogger::GetInstance() << "monitorData" << AliceO2::InfoLogger::InfoLogger::endm;
    fprintf(flog, "\n================\nmonitorData\n================\n");
  }
  monitorDataReadout(ctx);
  bool preclustersFound = false;
  bool preclusterDigitsFound = false;
  for (auto&& input : ctx.inputs()) {
    if (mPrintLevel >= 1) {
      QcInfoLogger::GetInstance() << "run PhysicsTask: input " << input.spec->binding << AliceO2::InfoLogger::InfoLogger::endm;
    }
    if (input.spec->binding == "digits") {
      monitorDataDigits(input);
    }
    if (input.spec->binding == "preclusters") {
      preclustersFound = true;
    }
    if (input.spec->binding == "preclusterdigits") {
      preclusterDigitsFound = true;
    }
  }
  //monitorDataReadout(ctx);
  if(preclustersFound && preclusterDigitsFound) {
    monitorDataPreclusters(ctx);
  }
}


void PhysicsTask::plotDigit(const o2::mch::Digit& digit)
{
  int ADC = digit.getADC();
  //int size = digit.getSize();
  int de = digit.getDetID();
  int padid = digit.getPadID();

  //fprintf(stdout, "digit[%d]: ADC=%d, DetId=%d, PadId=%d\n",
  //        i, ADC, de, padid);
  if (ADC < 0 || de < 0 || padid < 0) {
    return;
  }

  //mHistogramADCamplitudeVsSize->Fill(size, ADC);

  try {
    const o2::mch::mapping::Segmentation& segment = o2::mch::mapping::segmentation(de);

    double padX = segment.padPositionX(padid);
    double padY = segment.padPositionY(padid);
    float padSizeX = segment.padSizeX(padid);
    float padSizeY = segment.padSizeY(padid);
    int cathode = segment.isBendingPad(padid) ? 0 : 1;


    if (mPrintLevel >= 1)
      fprintf(flog, "de=%d pad=%d x=%f y=%f\n", de, padid, padX, padY);
    //if(pad.fX>=32 && pad.fX<=34 && pad.fY>=1.1 && pad.fY<=1.4)
    //  fprintf(flog, "mapping: link_id=%d ds_addr=%d chan_addr=%d  ==>  de=%d x=%f y=%f A=%d\n",
    //    hit.link_id, hit.ds_addr, hit.chan_addr, pad.fDE, pad.fX, pad.fY, hit.csum);

    auto h = mHistogramADCamplitudeDE.find(de);
    if ((h != mHistogramADCamplitudeDE.end()) && (h->second != NULL)) {
      h->second->Fill(ADC);
    }

    if (cathode == 0 && ADC > 0) {
      auto h2 = mHistogramNhitsDE.find(de);
      if ((h2 != mHistogramNhitsDE.end()) && (h2->second != NULL)) {
        int binx_min = h2->second->GetXaxis()->FindBin(padX - padSizeX / 2 + 0.1);
        int binx_max = h2->second->GetXaxis()->FindBin(padX + padSizeX / 2 - 0.1);
        int biny_min = h2->second->GetYaxis()->FindBin(padY - padSizeY / 2 + 0.1);
        int biny_max = h2->second->GetYaxis()->FindBin(padY + padSizeY / 2 - 0.1);
        for (int by = biny_min; by <= biny_max; by++) {
          float y = h2->second->GetYaxis()->GetBinCenter(by);
          for (int bx = binx_min; bx <= binx_max; bx++) {
            float x = h2->second->GetXaxis()->GetBinCenter(bx);
            h2->second->Fill(x, y);
          }
        }
      }
    }
    if (cathode == 0 && ADC > 500) {
      auto h2 = mHistogramNhitsHighAmplDE.find(de);
      if ((h2 != mHistogramNhitsHighAmplDE.end()) && (h2->second != NULL)) {
        int binx_min = h2->second->GetXaxis()->FindBin(padX - padSizeX / 2 + 0.1);
        int binx_max = h2->second->GetXaxis()->FindBin(padX + padSizeX / 2 - 0.1);
        int biny_min = h2->second->GetYaxis()->FindBin(padY - padSizeY / 2 + 0.1);
        int biny_max = h2->second->GetYaxis()->FindBin(padY + padSizeY / 2 - 0.1);
        for (int by = biny_min; by <= biny_max; by++) {
          float y = h2->second->GetYaxis()->GetBinCenter(by);
          for (int bx = binx_min; bx <= binx_max; bx++) {
            float x = h2->second->GetXaxis()->GetBinCenter(bx);
            h2->second->Fill(x, y);
          }
        }
      }
    }
  } catch (const std::exception& e) {
    QcInfoLogger::GetInstance() << "[MCH] Detection Element " << de << " not found in mapping." << AliceO2::InfoLogger::InfoLogger::endm;
    return;
  }
}




//_________________________________________________________________________________________________
static void CoG(gsl::span<const o2::mch::Digit> precluster, double& Xcog, double& Ycog)
{
  double xmin = 1E9;
  double ymin = 1E9;
  double xmax = -1E9;
  double ymax = -1E9;
  double charge[] = { 0.0, 0.0 };
  int multiplicity[] = { 0, 0 };

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



void PhysicsTask::checkPreclusters(gsl::span<const o2::mch::PreCluster> preClusters, gsl::span<const o2::mch::Digit> digits)
{
  bool doPrint = false;

  for(int pass = 0; pass < 2; pass++) {

    //std::cout<<"pass="<<pass<<"  doPrint="<<doPrint<<std::endl;

    if(doPrint) std::cout<<"\n\n============\n";
    for(auto& preCluster : preClusters) {
      // filter out single-pad clusters
      if (preCluster.nDigits < 2) {
        continue;
      }

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

      // filter out clusters with small charge, which are likely to be noise
      if ((chargeSum[0]+chargeSum[1]) < 100) {
        continue;
      }
      if ((chargeMax[0] < 100) && (chargeMax[1] < 100)) {
        continue;
      }

      double Xcog, Ycog;
      CoG(preClusterDigits, Xcog, Ycog);
      if(pass == 0 && cathode[0] && !cathode[1]) {
        if(Xcog > -30 && Xcog < -10 && Ycog > 6 && Ycog < 14) {
          doPrint = true;
          break;
        }
      }

      if (doPrint) {
        if (pass == 0) break;

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
  }
}


void PhysicsTask::printPreclusters(gsl::span<const o2::mch::PreCluster> preClusters, gsl::span<const o2::mch::Digit> digits)
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
    CoG(preClusterDigits, Xcog, Ycog);

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


bool PhysicsTask::plotPrecluster(const o2::mch::PreCluster& preCluster, gsl::span<const o2::mch::Digit> digits)
{
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

  // filter out single-pad clusters
  if (preCluster.nDigits < 2) {
    return true;
  }


  if (mPrintLevel >= 1) {
    std::cout<<"[pre-cluster] charge = "<<chargeSum[0]<<" "<<chargeSum[1]<<std::endl;
    for (auto& d : preClusterDigits) {
      std::cout << fmt::format("  DE {:4d}  PAD {:5d}  ADC {:6d}  TIME {:4d}", d.getDetID(), d.getPadID(), d.getADC(), d.getTime().sampaTime);
      std::cout << std::endl;
    }
  }

  //if(rnd.Rndm() < 0.8) cathode[0] = false;
  //if(rnd.Rndm() < 0.99) cathode[1] = false;


  float chargeTot = chargeSum[0] + chargeSum[1];
  /*if(cathode[0] && cathode[1]) {
    //chargeTot = (chargeSum[0] + chargeSum[1]) / 2;
    chargeTot = chargeSum[0] + chargeSum[1];
  } else if(cathode[0]) {
    chargeTot = chargeSum[0];
  } else if(cathode[1]) {
    chargeTot = chargeSum[1];
  }*/
  auto hCharge = mHistogramClchgDE.find(detid);
  if ((hCharge != mHistogramClchgDE.end()) && (hCharge->second != NULL)) {
    hCharge->second->Fill(chargeTot);
  }

  // filter out clusters with small charge, which are likely to be noise
  if ((chargeSum[0]+chargeSum[1]) < 100) {
    return true;
  }
  if ((chargeMax[0] < 100) && (chargeMax[1] < 100)) {
    return true;
  }

  double Xcog, Ycog;
  CoG(preClusterDigits, Xcog, Ycog);

  if (Ycog < 0) {
    return true;
  }

  auto hXY0 = mHistogramPreclustersXY[0].find(detid);
  if ((hXY0 != mHistogramPreclustersXY[0].end()) && (hXY0->second != NULL)) {
    hXY0->second->Fill(Xcog, Ycog);
  }

  int hid = 1;
  if(cathode[0]) {
    auto hXY1 = mHistogramPreclustersXY[1].find(detid);
    if ((hXY1 != mHistogramPreclustersXY[1].end()) && (hXY1->second != NULL)) {
      hXY1->second->Fill(Xcog, Ycog);
    }
  }
  if(cathode[1]) {
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


void PhysicsTask::endOfCycle()
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
    for (int i = 0; i < 3 * 24; i++) {
      //mHistogramNhits[i]->Write();
      //mHistogramADCamplitude[i]->Write();
    }
    mHistogramADCamplitudeVsSize->Write();
    //std::cout<<"mHistogramADCamplitudeDE.size() = "<<mHistogramADCamplitudeDE.size()<<"  DEs.size()="<<DEs.size()<<std::endl;
    int nbDEs = DEs.size();
    for (int elem = 0; elem < nbDEs; elem++) {
      int de = DEs[elem];
      //std::cout<<"  de="<<de<<std::endl;
      {
        auto h = mHistogramADCamplitudeDE.find(de);
        if ((h != mHistogramADCamplitudeDE.end()) && (h->second != NULL)) {
          h->second->Write();
        }
      }
      {
        auto h = mHistogramNhitsDE.find(de);
        if ((h != mHistogramNhitsDE.end()) && (h->second != NULL)) {
          h->second->Write();
        }
      }
      {
        auto h = mHistogramNhitsHighAmplDE.find(de);
        if ((h != mHistogramNhitsHighAmplDE.end()) && (h->second != NULL)) {
          h->second->Write();
        }
      }
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
    }
    mHistogramPseudoeff[0]->Write();
    mHistogramPseudoeff[1]->Write();
    mHistogramPseudoeff[2]->Write();

    //f.ls();
    f.Close();
#endif
}

void PhysicsTask::endOfActivity(Activity& /*activity*/)
{
  QcInfoLogger::GetInstance() << "endOfActivity" << AliceO2::InfoLogger::InfoLogger::endm;
}

void PhysicsTask::reset()
{
  // clean all the monitor objects here

  QcInfoLogger::GetInstance() << "Reseting the histogram" << AliceO2::InfoLogger::InfoLogger::endm;
}

} // namespace muonchambers
} // namespace quality_control_modules
} // namespace o2
