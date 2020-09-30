///
/// \file   PhysicsTaskDigits.cxx
/// \author Barthelemy von Haller
/// \author Piotr Konopka
/// \author Andrea Ferrero
/// \author Sebastien Perrin
///

#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <algorithm>

#include "DPLUtils/DPLRawParser.h"
#include "Headers/RAWDataHeader.h"
#include "MCH/PhysicsTaskDigits.h"
#ifdef MCH_HAS_MAPPING_FACTORY
#include "MCHMappingFactory/CreateSegmentation.h"
#endif
#include "MCHMappingInterface/Segmentation.h"
#include "MCHMappingInterface/CathodeSegmentation.h"
#include "MCHMappingSegContour/CathodeSegmentationContours.h"
#include "MCHRawDecoder/PageDecoder.h"
#include "QualityControl/QcInfoLogger.h"

using namespace std;

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
PhysicsTaskDigits::PhysicsTaskDigits() : TaskInterface(), count(1) {}

PhysicsTaskDigits::~PhysicsTaskDigits() { fclose(flog); }

void PhysicsTaskDigits::initialize(o2::framework::InitContext& /*ctx*/)
{
  QcInfoLogger::GetInstance() << "initialize PhysicsTaskDigits" << AliceO2::InfoLogger::InfoLogger::endm;
  fprintf(stdout, "initialize PhysicsTaskDigits\n");

  mDecoder.initialize();

  uint32_t dsid;
  for (int cruid = 0; cruid < 32; cruid++) {
      for (int linkid = 0; linkid < 24; linkid++) {
      {
          int index = 24 * cruid + linkid;
          mHistogramNhits[index] = new TH2F(TString::Format("QcMuonChambers_NHits_CRU%01d_LINK%02d", cruid, linkid),
                                        TString::Format("QcMuonChambers - Number of hits (CRU link %02d)", index), 40, 0, 40, 64, 0, 64);
          mHistogramADCamplitude[index] = new TH1F(TString::Format("QcMuonChambers_ADC_Amplitude_CRU%01d_LINK%02d", cruid, linkid),
                                        TString::Format("QcMuonChambers - ADC amplitude (CRU link %02d)", index), 5000, 0, 5000);
      }
      int32_t solar_id = mDecoder.getMapCRU(cruid, linkid);
      if (solar_id == -1)
          continue;
      for (int ds_addr = 0; ds_addr < 40; ds_addr++) {
          uint32_t de;
          int32_t result = mDecoder.getMapFEC(solar_id, ds_addr, de, dsid);
      if(result < 0) continue;
          if (std::find(DEs.begin(), DEs.end(), de) == DEs.end()) {
              DEs.push_back(de);
              TH1F* h = new TH1F(TString::Format("QcMuonChambers_ADCamplitude_DE%03d", de),
                            TString::Format("QcMuonChambers - ADC amplitude (DE%03d)", de), 5000, 0, 5000);
              mHistogramADCamplitudeDE.insert(make_pair(de, h));
              
              float Xsize = 40 * 5;
              float Xsize2 = Xsize / 2;
              float Ysize = 50;
              float Ysize2 = Ysize / 2;

              TH2F* h2 = new TH2F(TString::Format("QcMuonChambers_Nhits_DE%03d_B", de),
                    TString::Format("QcMuonChambers - Number of hits (DE%03d B)", de), Xsize * 2, -Xsize2, Xsize2, Ysize * 2, -Ysize2, Ysize2);
              mHistogramNhitsDE[0].insert(make_pair(de, h2));
              getObjectsManager()->startPublishing(h2);
              h2 = new TH2F(TString::Format("QcMuonChambers_Nhits_DE%03d_NB", de),
                    TString::Format("QcMuonChambers - Number of hits (DE%03d NB)", de), Xsize * 2, -Xsize2, Xsize2, Ysize * 2, -Ysize2, Ysize2);
              mHistogramNhitsDE[1].insert(make_pair(de, h2));
              getObjectsManager()->startPublishing(h2);
              h2 = new TH2F(TString::Format("QcMuonChambers_Nhits_HighAmpl_DE%03d_B", de),
                    TString::Format("QcMuonChambers - Number of hits for Csum>500 (DE%03d B)", de), Xsize * 2, -Xsize2, Xsize2, Ysize * 2, -Ysize2, Ysize2);
              mHistogramNhitsHighAmplDE[0].insert(make_pair(de, h2));
              h2 = new TH2F(TString::Format("QcMuonChambers_Nhits_HighAmpl_DE%03d_NB", de),
                    TString::Format("QcMuonChambers - Number of hits for Csum>500 (DE%03d NB)", de), Xsize * 2, -Xsize2, Xsize2, Ysize * 2, -Ysize2, Ysize2);
              mHistogramNhitsHighAmplDE[1].insert(make_pair(de, h2));
              h2 = new TH2F(TString::Format("QcMuonChambers_Norbits_DE%03d_B", de),
                    TString::Format("QcMuonChambers - Number of orbits (DE%03d B)", de), Xsize * 2, -Xsize2, Xsize2, Ysize * 2, -Ysize2, Ysize2);
              mHistogramNorbitsDE[0].insert(make_pair(de, h2));
              h2 = new TH2F(TString::Format("QcMuonChambers_Norbits_DE%03d_NB", de),
                    TString::Format("QcMuonChambers - Number of orbits (DE%03d NB)", de), Xsize * 2, -Xsize2, Xsize2, Ysize * 2, -Ysize2, Ysize2);
              mHistogramNorbitsDE[1].insert(make_pair(de, h2));
        }
      }
    }
  }
    
  for(int fee = 0; fee <= MCH_FFEID_MAX; fee++) {
      for(int link = 0; link < 12; link++) {
          norbits[fee][link] = lastorbitseen[fee][link] = 0;
    }
  }
  for(int de = 0; de < 1100; de++) {
    MeanOccupancyDE[de] = MeanOccupancyDECycle[de] = LastMeanNhitsDE[de] = LastMeanNorbitsDE[de] = NewMeanNhitsDE[de] = NewMeanNorbitsDE[de] = NbinsDE[de] = 0;
  }

    // Histograms using the Electronic Mapping

    mHistogramNorbitsElec = new TH2F("QcMuonChambers_Norbits_Elec", "QcMuonChambers - Norbits",
        (MCH_FFEID_MAX+1)*12*40, 0, (MCH_FFEID_MAX+1)*12*40, 64, 0, 64);
    getObjectsManager()->startPublishing(mHistogramNorbitsElec);
    mHistogramNHitsElec = new TH2F("QcMuonChambers_NHits_Elec", "QcMuonChambers - NHits",
        (MCH_FFEID_MAX+1)*12*40, 0, (MCH_FFEID_MAX+1)*12*40, 64, 0, 64);
    getObjectsManager()->startPublishing(mHistogramNHitsElec);
    mHistogramOccupancyElec = new TH2F("QcMuonChambers_Occupancy_Elec", "QcMuonChambers - Occupancy (MHz)",
        (MCH_FFEID_MAX+1)*12*40, 0, (MCH_FFEID_MAX+1)*12*40, 64, 0, 64);
    getObjectsManager()->startPublishing(mHistogramOccupancyElec);

    // 1D histograms for mean occupancy per DE (integrated or per elapsed cycle)
    
    mMeanOccupancyPerDE = new TH1F("QcMuonChambers_MeanOccupancy", "Mean Occupancy of each DE (MHz)", 1100, -0.5, 1099.5);
    getObjectsManager()->startPublishing(mMeanOccupancyPerDE);
    mMeanOccupancyPerDECycle = new TH1F("QcMuonChambers_MeanOccupancy_OnCycle", "Mean Occupancy of each DE during the cycle (MHz)", 1100, -0.5, 1099.5);
    getObjectsManager()->startPublishing(mMeanOccupancyPerDECycle);

    for (int de = 0; de < 1030; de++){
      const o2::mch::mapping::Segmentation* segment = &(o2::mch::mapping::segmentation(de));
      if (segment == nullptr) continue;

      float Xsize = 40 * 5;
      float Xsize2 = Xsize / 2;
      float Ysize = 50;
      float Ysize2 = Ysize / 2;
      float scale = 0.5;

      // Histograms using the XY Mapping

      {
        TH2F* hXY = new TH2F(TString::Format("QcMuonChambers_Occupancy_B_XY_%03d", de),
            TString::Format("QcMuonChambers - Occupancy XY (DE%03d B) (MHz)", de), Xsize / scale, -Xsize2, Xsize2, Ysize / scale, -Ysize2, Ysize2);
        mHistogramOccupancyXY[0].insert(make_pair(de, hXY));
        getObjectsManager()->startPublishing(hXY);
        hXY = new TH2F(TString::Format("QcMuonChambers_Occupancy_NB_XY_%03d", de),
            TString::Format("QcMuonChambers - Occupancy XY (DE%03d NB) (MHz)", de), Xsize / scale, -Xsize2, Xsize2, Ysize / scale, -Ysize2, Ysize2);
        mHistogramOccupancyXY[1].insert(make_pair(de, hXY));
        getObjectsManager()->startPublishing(hXY);
      }
    }

    mHistogramOccupancy[0] = new GlobalHistogram("QcMuonChambers_Occupancy_den", "Occupancy (MHz)");
    mHistogramOccupancy[0]->init();

    mHistogramOrbits[0] = new GlobalHistogram("QcMuonChambers_Orbits_den", "Orbits");
    mHistogramOrbits[0]->init();

  flog = stdout;
  fprintf(stdout, "PhysicsTaskDigits initialization finished\n");
}

void PhysicsTaskDigits::startOfActivity(Activity& /*activity*/)
{
  QcInfoLogger::GetInstance() << "startOfActivity" << AliceO2::InfoLogger::InfoLogger::endm;
}

void PhysicsTaskDigits::startOfCycle()
{
  QcInfoLogger::GetInstance() << "startOfCycle" << AliceO2::InfoLogger::InfoLogger::endm;
}

void PhysicsTaskDigits::monitorDataDigits(o2::framework::ProcessingContext& ctx)
{
    fprintf(flog, "\n================\nmonitorDataDigits\n================\n");
    
  // get the input preclusters and associated digits with the orbit information
  auto digits = ctx.inputs().get<gsl::span<o2::mch::Digit>>("digits");
  auto orbits = ctx.inputs().get<gsl::span<uint64_t>>("orbits");
    
    for (auto& orb : orbits){ //Normalement une seule fois
        uint32_t orbitnumber = (orb & 0xFFFFFFFF);
        uint32_t link = (orb >> 32) & 0xFF;
        uint32_t fee = (orb >> 40) & 0xFF;
        if(link != 15){
            if(orbitnumber != lastorbitseen[fee][link]) {
              norbits[fee][link] += 1;
            }
            lastorbitseen[fee][link] = orbitnumber;
        }
        else if(link == 15){
            for(int li=0; li<12; li++){
                if(orbitnumber != lastorbitseen[fee][li]) {
                  norbits[fee][li] += 1;
                }
                lastorbitseen[fee][li] = orbitnumber;
            }
        }
    }
        
    for (auto& d : digits) {
        plotDigit(d);
    }
}

void PhysicsTaskDigits::monitorData(o2::framework::ProcessingContext& ctx)
{
  bool digitsFound = false;
  bool orbitsFound = false;
  for (auto&& input : ctx.inputs()) {
    if (input.spec->binding == "digits") {
      digitsFound = true;
    }
    if (input.spec->binding == "orbits") {
      orbitsFound = true;
    }
  }
  if (digitsFound) {
    monitorDataDigits(ctx);
  }
}

void PhysicsTaskDigits::plotDigit(const o2::mch::Digit& digit)
{
  int ADC = digit.getADC();
  int de = digit.getDetID();
  int padid = digit.getPadID();

  if (ADC < 0 || de <= 0 || padid < 0) {
    return;
  }

  try {
    const o2::mch::mapping::Segmentation& segment = o2::mch::mapping::segmentation(de);

    double padX = segment.padPositionX(padid);
    double padY = segment.padPositionY(padid);
    float padSizeX = segment.padSizeX(padid);
    float padSizeY = segment.padSizeY(padid);
    int cathode = segment.isBendingPad(padid) ? 0 : 1;
    int dsid = segment.padDualSampaId(padid);
    int chan_addr = segment.padDualSampaChannel(padid);

    uint32_t solar_id = 0;
    uint32_t ds_addr = 0;
    int32_t cruid = 0;
    int32_t linkid = 0;
    int32_t fee_id = 0;

    // get the unique solar ID and the DS address associated to this digit
    solar_id = mDecoder.getMapFECinv(de, dsid, solar_id, ds_addr);

    // get the CRU ID and CRU link ID from the inverse CRU mapping, and compute the FEE ID
    if(mDecoder.getMapCRUInv(solar_id, cruid, linkid)){
      fee_id = cruid * 2 + (linkid / 12);
    }

    int xbin = fee_id * 12 * 40 + (linkid % 12) * 40 + ds_addr + 1;
    int ybin = chan_addr + 1;


    mHistogramNHitsElec->Fill(xbin-0.5, ybin-0.5);

    auto h = mHistogramADCamplitudeDE.find(de);
    if ((h != mHistogramADCamplitudeDE.end()) && (h->second != NULL)) {
      h->second->Fill(ADC);
    }


      for(int cat=0; cat<2; cat++){
        if (cathode == cat && ADC > 0) {
          auto h2 = mHistogramNhitsDE[cat].find(de);
          if ((h2 != mHistogramNhitsDE[cat].end()) && (h2->second != NULL)) {
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
          h2 = mHistogramNorbitsDE[cat].find(de);
          if ((h2 != mHistogramNorbitsDE[cat].end()) && (h2->second != NULL)) {
            int NYbins = h2->second->GetYaxis()->GetNbins();
            int NXbins = h2->second->GetXaxis()->GetNbins();
            for (int by = 0; by < NYbins; by++) {
              float y = h2->second->GetYaxis()->GetBinCenter(by);
              for (int bx = 0; bx < NXbins; bx++) {
                float x = h2->second->GetXaxis()->GetBinCenter(bx);

                int bpad = 0;
                int nbpad = 0;
                uint32_t solar_id_boucle = 0;
                uint32_t ds_addr_boucle = 0;
                int32_t cruid_boucle = 0;
                int32_t linkid_boucle = 0;
                int32_t fee_id_boucle = 0;

                if(segment.findPadPairByPosition(x, y, bpad, nbpad)){
                    int dsid_boucle = 0;
                    int chan_addr_boucle = 0;
                    if(cat == 0){
                      dsid_boucle = segment.padDualSampaId(bpad);
                      chan_addr_boucle = segment.padDualSampaChannel(bpad);
                    }
                    else if(cat == 1){
                        dsid_boucle = segment.padDualSampaId(nbpad);
                        chan_addr_boucle = segment.padDualSampaChannel(nbpad);
                    }
                  solar_id_boucle = mDecoder.getMapFECinv(de, dsid_boucle, solar_id_boucle, ds_addr_boucle);
                  if(solar_id_boucle == solar_id){
                    h2->second->SetBinContent(bx, by, norbits[fee_id][linkid]);
                    if(!mDecoder.getMapCRUInv(solar_id_boucle, cruid_boucle, linkid_boucle)){
                      std::cout << "Problem getMapCRUInv !!!!!!!!!!!!!" << std::endl;
                    }
                    if(mDecoder.getMapCRUInv(solar_id_boucle, cruid_boucle, linkid_boucle)){
                      fee_id_boucle = cruid_boucle * 2 + (linkid_boucle / 12);
                      int xbin = fee_id_boucle * 12 * 40 + (linkid_boucle % 12) * 40 + ds_addr_boucle + 1;
                      int ybin = chan_addr_boucle + 1;
                    }
                  }
                }

              }
            }
          }
        }
        if (cathode == cat && ADC > 500) {
          auto h2 = mHistogramNhitsHighAmplDE[cat].find(de);
          if ((h2 != mHistogramNhitsHighAmplDE[cat].end()) && (h2->second != NULL)) {
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
      }
      
  } catch (const std::exception& e) {
    QcInfoLogger::GetInstance() << "[MCH] Detection Element " << de << " not found in mapping." << AliceO2::InfoLogger::InfoLogger::endm;
    return;
  }
}

void PhysicsTaskDigits::endOfCycle()
{
  QcInfoLogger::GetInstance() << "endOfCycle" << AliceO2::InfoLogger::InfoLogger::endm;
  for(int de = 100; de <= 1030; de++) {
    for(int i = 0; i < 2; i++) {
        auto ih = mHistogramOccupancyXY[i+1].find(de);
        ih = mHistogramOccupancyXY[i].find(de);
        if (ih == mHistogramOccupancyXY[i].end()) {
          continue;
        }
    }
  }

  mHistogramOrbits[0]->set(mHistogramNorbitsDE[0], mHistogramNorbitsDE[1]);
  mHistogramOccupancy[0]->set(mHistogramNhitsDE[0], mHistogramNhitsDE[1]);
  mHistogramOccupancy[0]->Divide(mHistogramOrbits[0]);
  mHistogramOccupancy[0]->Scale(1/87.5); // Converting Occupancy from NbHits/NbOrbits to MHz

  // Fill norbits values for electronics channels associated to readout pads
  for(int fee_id = 0; fee_id <= MCH_FFEID_MAX; fee_id++) {
    int cruid = fee_id / 2;
    
    // loop on CRU link and check if it corresponds to an existing SOLAR board
    for(int linkid = 0; linkid < 12; linkid++) {
      int crulinkid = linkid + 12 * (fee_id % 2);

      int32_t solar_id = mDecoder.getMapCRU(cruid, crulinkid);
      if (solar_id == -1) {
        continue;
      }

      // loop on DS boards and check if it exists in the mapping
      for (int ds_addr = 0; ds_addr < 40; ds_addr++) {
        uint32_t de, dsid;
        int32_t result = mDecoder.getMapFEC(solar_id, ds_addr, de, dsid);
        if(result < 0) continue;

        int xbin = fee_id * 12 * 40 + (linkid % 12) * 40 + ds_addr + 1;

        // loop on DS channels and check if it is associated to a readout pad
        for(int chan_addr = 0; chan_addr < 64; chan_addr++) {

          const o2::mch::mapping::Segmentation& segment = o2::mch::mapping::segmentation(de);
          int padid = segment.findPadByFEE(dsid, chan_addr);
          if(padid < 0) {
            continue;
          }

          int ybin = chan_addr + 1;
          mHistogramNorbitsElec->SetBinContent(xbin, ybin, norbits[fee_id][linkid]);
        }
      }
    }
  }


  mHistogramOccupancyElec->Reset();
  mHistogramOccupancyElec->Add(mHistogramNHitsElec);
  mHistogramOccupancyElec->Divide(mHistogramNorbitsElec);
  mHistogramOccupancyElec->Scale(1/87.5); // Converting Occupancy from NbHits/NbOrbits to MHz

#ifdef QC_MCH_SAVE_TEMP_ROOTFILE
    TFile f("/tmp/qc.root", "RECREATE");
    
    mHistogramNorbitsElec->Write();
    mHistogramNHitsElec->Write();
    mHistogramOccupancyElec->Write();
    
    int nbDEs = DEs.size();
    for (int elem = 0; elem < nbDEs; elem++) {
      int de = DEs[elem];
      {
        auto h = mHistogramADCamplitudeDE.find(de);
        if ((h != mHistogramADCamplitudeDE.end()) && (h->second != NULL)) {
          h->second->Write();
        }
      }
      for(int i = 0; i < 2; i++) {
          {
            auto h = mHistogramNhitsDE[i].find(de);
            if ((h != mHistogramNhitsDE[i].end()) && (h->second != NULL)) {
              h->second->Write();
            }
          }
          {
            auto h = mHistogramNorbitsDE[i].find(de);
            if ((h != mHistogramNorbitsDE[i].end()) && (h->second != NULL)) {
              h->second->Write();
            }
          }
          {
              auto h = mHistogramOccupancyXY[i].find(de);
              if ((h != mHistogramOccupancyXY[i].end()) && (h->second != NULL)) {
                  auto hhits = mHistogramNhitsDE[i].find(de);
                  auto horbits = mHistogramNorbitsDE[i].find(de);
                  if ((hhits != mHistogramNhitsDE[i].end()) && (horbits != mHistogramNorbitsDE[i].end()) && (hhits->second != NULL) && (horbits->second != NULL)) {
                      h->second->Divide(hhits->second,horbits->second);
                      h->second->Scale(1/87.5); // Converting Occupancy from NbHits/NbOrbits to MHz
                      h->second->Write();
              }
            }
          }
      }
    }

    {
      // Using OccupancyElec to get the mean occupancy per DE
      auto h = mHistogramOccupancyElec;
      auto horbits = mHistogramNorbitsElec;
      auto h1 = mMeanOccupancyPerDE;
      if (h && h1 && horbits) {
        for(int de = 0; de < 1100; de++) {
          MeanOccupancyDE[de] = 0;
          NbinsDE[de] = 0;
        }
        for(int binx=1; binx<h->GetXaxis()->GetNbins()+1; binx++){
          for(int biny=1; biny<h->GetYaxis()->GetNbins()+1; biny++){

            int norbits = horbits->GetBinContent(binx, biny);
            if(norbits <= 0) {
              // no orbits detected for this channel, skip it
              continue;
            }

            uint32_t ds_addr =  (binx-1) % 40;
            uint32_t linkid = ( (binx-1-ds_addr) / 40 ) % 12;
            uint32_t fee_id = (binx-1-ds_addr-40*linkid) / (12*40);
            uint32_t chan_addr = biny-1;
            uint32_t de;
            uint32_t dsid;
            uint32_t cru_id = fee_id / 2;
            uint32_t crulinkid = linkid + 12 * (fee_id % 2);
            int32_t solar_id = mDecoder.getMapCRU(cru_id, crulinkid);
            if(solar_id > -1){
              int32_t result = mDecoder.getMapFEC(solar_id, ds_addr, de, dsid);
              if(result > -1){
                MeanOccupancyDE[de] += h->GetBinContent(binx, biny);
                NbinsDE[de] += 1;
              }
            }
          }
        }
        for(int i=0; i<1100; i++){
          if(NbinsDE[i]>0){
            MeanOccupancyDE[i] /= NbinsDE[i];
          }
          h1->SetBinContent(i+1, MeanOccupancyDE[i]);
        }
        h1->Write();
        std::cout << "MeanOccupancy of DE502 is: " << MeanOccupancyDE[502] << std::endl;
      }
    }

    {
        // Using NHitsElec and Norbits to get the mean occupancy per DE on last cycle
      auto hhits = mHistogramNHitsElec;
      auto horbits = mHistogramNorbitsElec;
      auto h1 = mMeanOccupancyPerDECycle;
      if (hhits && horbits && h1) {
          for(int de = 0; de < 1100; de++) {
              NewMeanNhitsDE[de] = 0;
              NewMeanNorbitsDE[de] = 0;
          }
          for(int binx=1; binx<hhits->GetXaxis()->GetNbins()+1; binx++){
              for(int biny=1; biny<hhits->GetYaxis()->GetNbins()+1; biny++){
                  uint32_t ds_addr =  (binx-1) % 40;
                  uint32_t linkid = ( (binx-1-ds_addr) / 40 ) % 12;
                  uint32_t fee_id = (binx-1-ds_addr-40*linkid) / (12*40);
                  uint32_t chan_addr = biny-1;
                  uint32_t de;
                  uint32_t dsid;
                  uint32_t cru_id = fee_id / 2;
                  uint32_t crulinkid = linkid + 12 * (fee_id % 2);
                  int32_t solar_id = mDecoder.getMapCRU(cru_id, crulinkid);
                  if(solar_id > -1){
                      int32_t result = mDecoder.getMapFEC(solar_id, ds_addr, de, dsid);
                      if(result > -1){
                          NewMeanNhitsDE[de] += hhits->GetBinContent(binx, biny);
                          NewMeanNorbitsDE[de] += horbits->GetBinContent(binx, biny);
                          NbinsDE[de] += 1;
                      }
                  }
              }
          }
          for(int i=0; i<1100; i++){
              MeanOccupancyDECycle[i] = 0;
              if(NbinsDE[i]>0){
                  NewMeanNhitsDE[i] /= NbinsDE[i];
                  NewMeanNorbitsDE[i] /= NbinsDE[i];
              }
              if((NewMeanNorbitsDE[i]-LastMeanNorbitsDE[i]) > 0){
                  MeanOccupancyDECycle[i] = (NewMeanNhitsDE[i]-LastMeanNhitsDE[i])/(NewMeanNorbitsDE[i]-LastMeanNorbitsDE[i])/87.5; //Scaling to MHz
              }
              h1->SetBinContent(i+1, MeanOccupancyDECycle[i]);
              LastMeanNhitsDE[i] = NewMeanNhitsDE[i];
              LastMeanNorbitsDE[i] = NewMeanNorbitsDE[i];
          }
          h1->Write();
          std::cout << "MeanOccupancy of DE502 in last cycle is: " << MeanOccupancyDECycle[502] << std::endl;
      }
    }
    
    mHistogramOrbits[0]->Write();
    mHistogramOccupancy[0]->Write();

    //f.ls();
    f.Close();
#endif
}


void PhysicsTaskDigits::endOfActivity(Activity& /*activity*/)
{
  QcInfoLogger::GetInstance() << "endOfActivity" << AliceO2::InfoLogger::InfoLogger::endm;
}

void PhysicsTaskDigits::reset()
{
  // clean all the monitor objects here

  QcInfoLogger::GetInstance() << "Reseting the histogram" << AliceO2::InfoLogger::InfoLogger::endm;
}

} // namespace muonchambers
} // namespace quality_control_modules
} // namespace o2
