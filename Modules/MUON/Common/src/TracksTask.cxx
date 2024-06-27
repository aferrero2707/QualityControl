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

#include "MUONCommon/TracksTask.h"

#include "MUONCommon/Helpers.h"
#include "QualityControl/ObjectsManager.h"
#include <DataFormatsMCH/Cluster.h>
#include <DataFormatsMCH/Digit.h>
#include <DataFormatsMCH/ROFRecord.h>
#include <DataFormatsMCH/TrackMCH.h>
#include <DetectorsBase/GeometryManager.h>
#include "QualityControl/QcInfoLogger.h"
#include <Framework/DataRefUtils.h>
#include <Framework/InputRecord.h>
#include <Framework/TimingInfo.h>
#include <MCHGeometryTransformer/Transformations.h>
#include <ReconstructionDataFormats/TrackMCHMID.h>
#include <ReconstructionDataFormats/GlobalFwdTrack.h>
#include <DataFormatsZDC/RecEventFlat.h>
#include <CommonConstants/LHCConstants.h>
#include <gsl/span>

#include <fstream>

namespace o2::quality_control_modules::muon
{

TracksTask::TracksTask()
{
}

TracksTask::~TracksTask() = default;

GID::mask_t adaptSource(GID::mask_t src)
{
  if (src[GID::Source::MFTMCHMID] == 1) {
    src.reset(GID::Source::MFTMCHMID); // does not exist
    src.set(GID::Source::MFTMCH);
    // ensure we request the individual tracks as we use their information in the plotter
    src.set(GID::Source::MFT);
    src.set(GID::Source::MCH);
    src.set(GID::Source::MID);
  }
  if (src[GID::Source::MFTMCH] == 1) {
    // ensure we request the individual tracks as we use their information in the plotter
    src.set(GID::Source::MFT);
    src.set(GID::Source::MCH);
  }
  if (src[GID::Source::MCHMID] == 1) {
    // ensure we request the individual tracks as we use their information in the plotter
    src.set(GID::Source::MCH);
    src.set(GID::Source::MID);
  }
  return src;
}

void TracksTask::initialize(o2::framework::InitContext& /*ic*/)
{
  ILOG(Debug, Devel) << "initialize TracksTask" << ENDM; // QcInfoLogger is used. FairMQ logs will go to there as well.

  mBackgroundITS.clear();
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 108626425));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 108657381));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 108679092));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 108679092));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 108713194));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 108713194));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 108722471));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 108734899));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 108761049));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 108761049));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 108799971));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 108830949));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 108830949));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 108852658));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 108869515));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 108886778));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 108896055));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 109047915));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 109082036));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 109086490));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 109233912));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109255604));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 109255604));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 109303425));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 109303425));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109303442));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 109303442));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109329978));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 109351653));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 109416763));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109472564));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109477005));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109498699));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109503540));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 109503540));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109525240));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 109525240));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109537642));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 109537642));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109542092));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 109559336));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109559349));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 109581028));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 109581028));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 109607195));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 109607195));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 109624427));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109628875));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 109628875));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109689510));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 109689510));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 109693965));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 109693965));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109715653));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 109715653));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 109732919));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 109732919));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 109737348));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 109742197));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 109824138));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 109824138));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109824151));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 109863101));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 109863101));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 109889224));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 109910913));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109915773));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109932620));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 109932620));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109932621));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 109932633));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 109937443));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 109937443));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 109937452));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 109949872));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 109976003));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 109976003));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 110036659));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 110036659));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 110045927));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 110045927));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 110067636));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 110149568));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 110149568));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 110210215));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 110241190));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 110241190));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 110253610));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 110275303));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 110275303));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 110297009));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 110297009));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 110318708));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 110318708));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 110327980));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 110340408));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 110349678));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 110362087));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 110362097));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 110383778));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 110393056));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 110393074));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 110393074));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 110427188));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 110448867));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 110448867));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 110448879));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 110448879));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 110448883));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 110453322));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 110492269));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 110557368));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 110561792));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 110579059));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 110631714));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 110691981));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 110691981));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 110740212));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 110757062));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 110796007));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 110843858));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 110908933));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 110913781));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 110913781));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 110952333));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 111000556));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 111022252));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111034684));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111039115));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 111039115));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 111056366));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 111060825));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 111060825));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 111121445));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 111147609));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 111164837));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 111164860));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111164861));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 111164861));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111195817));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 111195817));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 111273318));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111277770));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111295036));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 111295036));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111304312));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 111342868));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111360098));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111364549));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111364561));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 111386253));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 111391072));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111391090));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 111391090));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111391097));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 111407957));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 111407957));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 111429632));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 111446895));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 111468582));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111473053));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111490298));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 111490298));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111538129));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 111542946));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 111581527));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 111624926));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 111663850));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 111668292));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 111668294));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 111668294));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 111750648));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 111755085));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 111815725));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 111820186));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 111825000));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 111825000));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 111825009));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 111880809));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 111880809));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111885266));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 111902521));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111906959));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 111906959));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111924193));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 111924193));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 111955184));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 111972060));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 111972060));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 111976869));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 111976869));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 111976889));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 112010985));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 112015425));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112037123));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112037126));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 112037149));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 112037149));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112058838));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112058838));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112063669));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112063669));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112102209));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112102209));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112107040));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 112107041));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 112128744));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 112172128));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112227954));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112227954));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112232411));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 112237216));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 112293028));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 112297494));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112379814));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 112379815));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 112379815));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112405962));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 112405962));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 112405962));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112405962));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 112423219));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 112475892));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112497582));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112514457));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112531688));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112536158));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112596771));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112601230));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112622914));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112649456));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112649456));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112683552));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 112688005));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112736240));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 112774815));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112801318));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 112839880));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 112839880));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 112839885));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 112878819));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 112878819));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112900518));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112900518));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 112909802));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112922212));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 112926657));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 112926657));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 112926658));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 112926658));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 112953186));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 112953214));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 112953214));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 112965625));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 113035146));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 113039986));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 113039986));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 113039989));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 113083381));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 113121922));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 113121922));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 113126758));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 113139193));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 113143635));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 113182562));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 113182562));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 113204271));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 113204271));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 113235240));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 113256952));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 113273793));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 113291054));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 113291054));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 113291058));
  mBackgroundITS.emplace_back(o2::InteractionRecord(0, 113300338));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 113360586));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 113360586));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 113365430));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 113365430));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2376, 113377831));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 113377831));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 113377836));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 113377849));
  mBackgroundITS.emplace_back(o2::InteractionRecord(594, 113377849));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 113382300));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 113403972));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 113403978));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 113408801));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 113408803));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 113408805));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1782, 113430517));
  mBackgroundITS.emplace_back(o2::InteractionRecord(2970, 113452198));
  mBackgroundITS.emplace_back(o2::InteractionRecord(1188, 113473915));

  ILOG(Info, Support) << "loading sources" << ENDM; // QcInfoLogger is used. FairMQ logs will go to there as well.

  auto srcFixed = mSrc;
  // For track type selection
  if (auto param = mCustomParameters.find("GID"); param != mCustomParameters.end()) {
    ILOG(Info, Devel) << "Custom parameter - GID (= sources by user): " << param->second << ENDM;
    ILOG(Info, Devel) << "Allowed sources           = " << mAllowedSources << " " << GID::getSourcesNames(mAllowedSources) << ENDM;
    auto requested = GID::getSourcesMask(param->second);
    ILOG(Info, Devel) << "Requested Sources         = " << requested << " " << GID::getSourcesNames(requested) << ENDM;
    mSrc = mAllowedSources & requested;
    srcFixed = adaptSource(mSrc);
    ILOG(Info, Devel) << "Allowed requested sources = " << mSrc << " " << GID::getSourcesNames(mSrc) << ENDM;
    ILOG(Info, Devel) << "Sources for data request  = " << srcFixed << " " << GID::getSourcesNames(srcFixed) << ENDM;
  }

  ILOG(Info, Support) << "Will do DataRequest for " << GID::getSourcesNames(srcFixed) << ENDM;
  if (srcFixed[GID::Source::MFTMCHMID] == 1) {
    srcFixed.reset(GID::Source::MFTMCHMID);
    srcFixed.set(GID::Source::MFTMCH);
  }
  mDataRequest = std::make_shared<o2::globaltracking::DataRequest>();
  mDataRequest->requestTracks(srcFixed, false);
}

void TracksTask::createTrackHistos(const Activity& activity)
{
  bool fullHistos = getConfigurationParameter<bool>(mCustomParameters, "fullHistos", false, activity);

  double maxTracksPerTF = getConfigurationParameter<double>(mCustomParameters, "maxTracksPerTF", 400, activity);
  double cutRAbsMin = getConfigurationParameter<double>(mCustomParameters, "cutRAbsMin", 17.6, activity);
  double cutRAbsMax = getConfigurationParameter<double>(mCustomParameters, "cutRAbsMax", 89.5, activity);
  double cutEtaMin = getConfigurationParameter<double>(mCustomParameters, "cutEtaMin", -4.0, activity);
  double cutEtaMax = getConfigurationParameter<double>(mCustomParameters, "cutEtaMax", -2.5, activity);
  double cutPtMin = getConfigurationParameter<double>(mCustomParameters, "cutPtMin", 0.5, activity);
  double cutChi2Min = getConfigurationParameter<double>(mCustomParameters, "cutChi2Min", 0, activity);
  double cutChi2Max = getConfigurationParameter<double>(mCustomParameters, "cutChi2Max", 1000, activity);
  double nSigmaPDCA = getConfigurationParameter<double>(mCustomParameters, "nSigmaPDCA", 6, activity);
  double matchScoreMaxMFT = getConfigurationParameter<double>(mCustomParameters, "matchScoreMaxMFT", 1000, activity);
  double diMuonTimeCut = getConfigurationParameter<double>(mCustomParameters, "diMuonTimeCut", 100, activity) / 1000;

  int etaBins = getConfigurationParameter<int>(mCustomParameters, "etaBins", 200, activity);
  int phiBins = getConfigurationParameter<int>(mCustomParameters, "phiBins", 180, activity);
  int ptBins = getConfigurationParameter<int>(mCustomParameters, "ptBins", 300, activity);

  //======================================
  // Track plotters without cuts

  auto createPlotter = [&](GID::Source source, std::string path) {
    if (mSrc[source] == 1) {
      ILOG(Info, Devel) << "Creating plotter for path " << path << ENDM;
      mTrackPlotters[source] = std::make_unique<muon::TrackPlotter>(maxTracksPerTF, etaBins, phiBins, ptBins, source, path, fullHistos);
      mTrackPlotters[source]->publish(getObjectsManager());
    }
  };

  createPlotter(GID::Source::MCH, "");
  createPlotter(GID::Source::MCHMID, "MCH-MID/");
  createPlotter(GID::Source::MFTMCH, "MFT-MCH/");
  createPlotter(GID::Source::MFTMCHMID, "MFT-MCH-MID/");

  //======================================
  // Track plotters with cuts

  std::vector<MuonCutFunc> muonCuts{
    // Rabs cut
    [cutRAbsMin, cutRAbsMax](const MuonTrack& t) { return ((t.getRAbs() >= cutRAbsMin) && (t.getRAbs() <= cutRAbsMax)); },
    // Eta cut
    [cutEtaMin, cutEtaMax](const MuonTrack& t) { return ((t.getMuonMomentumAtVertexMCH().eta() >= cutEtaMin) && (t.getMuonMomentumAtVertexMCH().eta() <= cutEtaMax)); },
    // Pt cut
    [cutPtMin](const MuonTrack& t) { return ((t.getMuonMomentumAtVertexMCH().Pt() >= cutPtMin)); },
    // pDCA cut
    [nSigmaPDCA](const MuonTrack& t) {
      static const double sigmaPDCA23 = 80.;
      static const double sigmaPDCA310 = 54.;
      static const double relPRes = 0.0004;
      static const double slopeRes = 0.0005;

      double thetaAbs = TMath::ATan(t.getRAbs() / 505.) * TMath::RadToDeg();

      double pUncorr = t.getTrackParamMCH().p();
      double p = t.getMuonMomentumAtVertexMCH().P();

      double pDCA = pUncorr * t.getDCAMCH();
      double sigmaPDCA = (thetaAbs < 3) ? sigmaPDCA23 : sigmaPDCA310;
      double nrp = nSigmaPDCA * relPRes * p;
      double pResEffect = sigmaPDCA / (1. - nrp / (1. + nrp));
      double slopeResEffect = 535. * slopeRes * p;
      double sigmaPDCAWithRes = TMath::Sqrt(pResEffect * pResEffect + slopeResEffect * slopeResEffect);
      if (pDCA > nSigmaPDCA * sigmaPDCAWithRes) {
        return false;
      }

      return true;
    },
    // MFT-MCH match score
    [matchScoreMaxMFT](const MuonTrack& t) {
      if (t.hasMFT() && t.hasMCH() && t.getMatchInfoFwd().getMFTMCHMatchingScore() > matchScoreMaxMFT)
        return false;
      return true;
    },
    // MCH chi2 cut
    [cutChi2Min, cutChi2Max](const MuonTrack& t) { return ((t.getChi2OverNDFMCH() >= cutChi2Min) && (t.getChi2OverNDFMCH() <= cutChi2Max)); }
  };

  // ZDC background selection/rejection
  std::vector<MuonCutFunc> zdcBdgSelection {
    [&](const MuonTrack& t) {
      auto muonIR = t.getIRMCH();

      for (const auto& zdcIR : mBackgroundZDC) {
        auto diffIR = muonIR.toLong() - zdcIR.toLong() - 31;
        if (std::abs(diffIR) <= 5) {
          std::cout << "[PIPPO] ZDC background found at " << zdcIR << " (" << muonIR << ")" << std::endl;
          return true;
        }
      }
      return false;
    }
  };

  std::vector<MuonCutFunc> zdcBdgRejection {
    [&](const MuonTrack& t) {
      auto muonIR = t.getIRMCH();

      for (const auto& zdcIR : mBackgroundZDC) {
        auto diffIR = muonIR.toLong() - zdcIR.toLong() - 31;
        if (std::abs(diffIR) <= 5) {
          return false;
        }
      }
      return true;
    }
  };

  // ITS background selection/rejection
  std::vector<MuonCutFunc> itsBdgSelection {
    [&](const MuonTrack& t) {
      auto muonIR = t.getIRMCH();

      for (const auto& itsIR : mBackgroundITS) {
        auto diffIR = muonIR.toLong() - itsIR.toLong();
        if (diffIR >= 0 && diffIR <= 594) {
          std::cout << "[PIPPO] ITS background found at " << itsIR << " (" << muonIR << ")" << std::endl;
          return true;
        }
      }
      return false;
    }
  };

  std::vector<MuonCutFunc> itsBdgRejection {
    [&](const MuonTrack& t) {
      auto muonIR = t.getIRMCH();

      for (const auto& itsIR : mBackgroundITS) {
        auto diffIR = muonIR.toLong() - itsIR.toLong();
        if (diffIR >= 0 && diffIR <= 594) {
          std::cout << "[PIPPO] ITS background rejected " << std::endl;
          return false;
        }
      }
      return true;
    }
  };

  std::vector<MuonCutFunc> bdgSelection{ zdcBdgSelection[0], itsBdgSelection[0] };
  std::vector<MuonCutFunc> bdgRejection{ zdcBdgRejection[0], itsBdgRejection[0] };

  std::vector<MuonCutFunc> zdcOnlyBdgSelection{ zdcBdgSelection[0], itsBdgRejection[0] };
  std::vector<MuonCutFunc> itsOnlyBdgSelection{ zdcBdgRejection[0], itsBdgSelection[0] };

  std::vector<DiMuonCutFunc> diMuonCuts{
    // cut on time difference between the two muon tracks
    [diMuonTimeCut](const MuonTrack& t1, const MuonTrack& t2) { return (std::abs(t1.getTime().getTimeStamp() - t2.getTime().getTimeStamp()) < diMuonTimeCut); }
  };

  auto createPlotterWithCuts = [&](GID::Source source, std::string path, std::vector<MuonCutFunc>& cuts) {
    if (mSrc[source] == 1) {
      ILOG(Info, Devel) << "Creating plotter for path " << path << ENDM;
      mTrackPlottersWithCuts[source] = std::make_unique<muon::TrackPlotter>(maxTracksPerTF, etaBins, phiBins, ptBins, source, path, fullHistos);
      mTrackPlottersWithCuts[source]->setMuonCuts(cuts);
      mTrackPlottersWithCuts[source]->setDiMuonCuts(diMuonCuts);
      mTrackPlottersWithCuts[source]->publish(getObjectsManager());
    }
  };

  createPlotterWithCuts(GID::Source::MCH, "WithCuts/", muonCuts);
  createPlotterWithCuts(GID::Source::MCHMID, "MCH-MID/WithCuts/", muonCuts);
  createPlotterWithCuts(GID::Source::MFTMCH, "MFT-MCH/WithCuts/", muonCuts);
  createPlotterWithCuts(GID::Source::MFTMCHMID, "MFT-MCH-MID/WithCuts/", muonCuts);

  mTrackPlottersBgdZDC[0] = std::make_unique<muon::TrackPlotter>(maxTracksPerTF, etaBins, phiBins, ptBins, GID::Source::MCH, "BgdZDC/", fullHistos);
  mTrackPlottersBgdZDC[0]->setMuonCuts(zdcBdgSelection);
  mTrackPlottersBgdZDC[0]->setDiMuonCuts(diMuonCuts);
  mTrackPlottersBgdZDC[0]->publish(getObjectsManager());

  mTrackPlottersBgdZDC[1] = std::make_unique<muon::TrackPlotter>(maxTracksPerTF, etaBins, phiBins, ptBins, GID::Source::MCH, "NoBgdZDC/", fullHistos);
  mTrackPlottersBgdZDC[1]->setMuonCuts(zdcBdgRejection);
  mTrackPlottersBgdZDC[1]->setDiMuonCuts(diMuonCuts);
  mTrackPlottersBgdZDC[1]->publish(getObjectsManager());

  mTrackPlottersBgdZDC[2] = std::make_unique<muon::TrackPlotter>(maxTracksPerTF, etaBins, phiBins, ptBins, GID::Source::MCH, "BgdITS/", fullHistos);
  mTrackPlottersBgdZDC[2]->setMuonCuts(itsBdgSelection);
  mTrackPlottersBgdZDC[2]->setDiMuonCuts(diMuonCuts);
  mTrackPlottersBgdZDC[2]->publish(getObjectsManager());

  mTrackPlottersBgdZDC[3] = std::make_unique<muon::TrackPlotter>(maxTracksPerTF, etaBins, phiBins, ptBins, GID::Source::MCH, "NoBgdITS/", fullHistos);
  mTrackPlottersBgdZDC[3]->setMuonCuts(itsBdgRejection);
  mTrackPlottersBgdZDC[3]->setDiMuonCuts(diMuonCuts);
  mTrackPlottersBgdZDC[3]->publish(getObjectsManager());

  mTrackPlottersBgdZDC[4] = std::make_unique<muon::TrackPlotter>(maxTracksPerTF, etaBins, phiBins, ptBins, GID::Source::MCH, "BgdITSZDC/", fullHistos);
  mTrackPlottersBgdZDC[4]->setMuonCuts(bdgSelection);
  mTrackPlottersBgdZDC[4]->setDiMuonCuts(diMuonCuts);
  mTrackPlottersBgdZDC[4]->publish(getObjectsManager());

  mTrackPlottersBgdZDC[5] = std::make_unique<muon::TrackPlotter>(maxTracksPerTF, etaBins, phiBins, ptBins, GID::Source::MCH, "NoBgdITSZDC/", fullHistos);
  mTrackPlottersBgdZDC[5]->setMuonCuts(bdgRejection);
  mTrackPlottersBgdZDC[5]->setDiMuonCuts(diMuonCuts);
  mTrackPlottersBgdZDC[5]->publish(getObjectsManager());

  mTrackPlottersBgdZDC[6] = std::make_unique<muon::TrackPlotter>(maxTracksPerTF, etaBins, phiBins, ptBins, GID::Source::MCH, "BgdZDC-NoBgdITS/", fullHistos);
  mTrackPlottersBgdZDC[6]->setMuonCuts(zdcOnlyBdgSelection);
  mTrackPlottersBgdZDC[6]->setDiMuonCuts(diMuonCuts);
  mTrackPlottersBgdZDC[6]->publish(getObjectsManager());

  mTrackPlottersBgdZDC[7] = std::make_unique<muon::TrackPlotter>(maxTracksPerTF, etaBins, phiBins, ptBins, GID::Source::MCH, "BgdITS-NoBgdZDC/", fullHistos);
  mTrackPlottersBgdZDC[7]->setMuonCuts(itsOnlyBdgSelection);
  mTrackPlottersBgdZDC[7]->setDiMuonCuts(diMuonCuts);
  mTrackPlottersBgdZDC[7]->publish(getObjectsManager());

  mBcZDC = std::make_unique<TH1F>("BcZDC", "BcZDC;bc", o2::constants::lhc::LHCMaxBunches, 0, o2::constants::lhc::LHCMaxBunches);
  getObjectsManager()->startPublishing(mBcZDC.get());

  mDCAvsBcZDC = std::make_unique<TH2F>("DCAvsBcZDC", "DCAvsBcZDC;bc;DCA (cm)", 200, -100, 100, 100, 0, 100);
  getObjectsManager()->startPublishing(mDCAvsBcZDC.get());
}

void TracksTask::removeTrackHistos()
{
  ILOG(Debug, Devel) << "Un-publishing objects" << ENDM;
  for (auto& p : mTrackPlotters) {
    p.second->unpublish(getObjectsManager());
  }
  for (auto& p : mTrackPlottersWithCuts) {
    p.second->unpublish(getObjectsManager());
  }

  ILOG(Debug, Devel) << "Destroying objects" << ENDM;
  mTrackPlotters.clear();
  mTrackPlottersWithCuts.clear();
}

void TracksTask::startOfActivity(const Activity& activity)
{
  ILOG(Debug, Devel) << "startOfActivity : " << activity << ENDM;
  createTrackHistos(activity);
}

void TracksTask::startOfCycle()
{
  ILOG(Debug, Devel) << "startOfCycle" << ENDM;
}

bool TracksTask::assertInputs(o2::framework::ProcessingContext& ctx)
{
  if (!ctx.inputs().isValid("trackMCH")) {
    ILOG(Info, Support) << "no mch tracks available on input" << ENDM;
    return false;
  }
  if (!ctx.inputs().isValid("trackMCHROF")) {
    ILOG(Info, Support) << "no mch track rofs available on input" << ENDM;
    return false;
  }
  if (!ctx.inputs().isValid("trackMCHTRACKCLUSTERS")) {
    ILOG(Info, Support) << "no mch track clusters available on input" << ENDM;
    return false;
  }
  if (mSrc[GID::Source::MCHMID] == 1) {
    if (!ctx.inputs().isValid("matchMCHMID")) {
      ILOG(Info, Support) << "no muon (mch+mid) track available on input" << ENDM;
      return false;
    }
    if (!ctx.inputs().isValid("trackMID")) {
      ILOG(Info, Support) << "no mid track available on input" << ENDM;
      return false;
    }
  }
  if (mSrc[GID::Source::MFTMCH] == 1) {
    if (!ctx.inputs().isValid("fwdtracks")) {
      ILOG(Info, Support) << "no muon (mch+mft) track available on input" << ENDM;
      return false;
    }
  }
  if (mSrc[GID::Source::MFTMCHMID] == 1) {
    if (!ctx.inputs().isValid("fwdtracks")) {
      ILOG(Info, Support) << "no muon (mch+mft+mid) track available on input" << ENDM;
      return false;
    }
  }
  return true;
}

void TracksTask::monitorData(o2::framework::ProcessingContext& ctx)
{
  ILOG(Debug, Devel) << "Debug: MonitorData" << ENDM;

  int firstTForbit = ctx.services().get<o2::framework::TimingInfo>().firstTForbit;
  ILOG(Debug, Devel) << "Debug: firstTForbit=" << firstTForbit << ENDM;

  if (!assertInputs(ctx)) {
    return;
  }

  ILOG(Debug, Devel) << "Debug: Asserted inputs" << ENDM;

  mRecoCont.collectData(ctx, *mDataRequest.get());

  // ZDC data
  auto bcrecZDC = ctx.inputs().get<gsl::span<o2::zdc::BCRecData>>("zdc-bcrec");
  auto energyZDC = ctx.inputs().get<gsl::span<o2::zdc::ZDCEnergy>>("zdc-energyrec");
  auto tdcZDC = ctx.inputs().get<gsl::span<o2::zdc::ZDCTDCData>>("zdc-tdcrec");
  auto infoZDC = ctx.inputs().get<gsl::span<uint16_t>>("zdc-inforec");

  o2::zdc::RecEventFlat ev;
  ev.init(bcrecZDC, energyZDC, tdcZDC, infoZDC);

  mBackgroundZDC.clear();
  while (ev.next()) {
    int32_t itdc = o2::zdc::TDCZNAC;
    int nhit = ev.NtdcV(itdc);

    for (int32_t ipos = 0; ipos < nhit; ipos++) {
      double mytdc = o2::zdc::FTDCVal * ev.TDCVal[itdc][ipos];

      if (mytdc > 5.7 && mytdc < 8.7) {
        // Backgroud event found here!
        std::cout << "[TOTO] ZDC background event found: " << ev.ir << std::endl;
        mBackgroundZDC.emplace_back(ev.ir);
        mBcZDC->Fill(ev.ir.bc);
        for (const auto& itsIR : mBackgroundITS) {
          if (itsIR.orbit == ev.ir.orbit) {
            std::cout << "ITS+ZDC background event in same orbit: ZDC [" << ev.ir << "]  ITS [ " << itsIR << "]  BC difference: " << (ev.ir.toLong() - itsIR.toLong()) << std::endl;
          }
        }
      }
    }
  }

  ILOG(Debug, Devel) << "Debug: Collected data" << ENDM;

  for (auto& p : mTrackPlotters) {
    if (p.second) {
      p.second->setFirstTForbit(firstTForbit);
    }
  }
  for (auto& p : mTrackPlottersWithCuts) {
    if (p.second) {
      p.second->setFirstTForbit(firstTForbit);
    }
  }
  for (auto& p : mTrackPlottersBgdZDC) {
    if (p) {
      p->setFirstTForbit(firstTForbit);
    }
  }

  if (mSrc[GID::MCH] == 1) {
    ILOG(Debug, Devel) << "Debug: MCH requested" << ENDM;
    if (mRecoCont.isTrackSourceLoaded(GID::MCH)) {
      ILOG(Debug, Devel) << "Debug: MCH source loaded" << ENDM;
      mTrackPlotters[GID::MCH]->fillHistograms(mRecoCont);
      mTrackPlottersWithCuts[GID::MCH]->fillHistograms(mRecoCont);
      for (auto& p : mTrackPlottersBgdZDC) {
        p->fillHistograms(mRecoCont);
      }
      for (const auto& zdcIR : mBackgroundZDC) {
        for (const auto& track : mTrackPlotters[GID::MCH]->getMuonTracks()) {
          auto muonIR = track.first.getIRMCH();

          auto diffIR = muonIR.toLong() - zdcIR.toLong();
          if (std::abs(diffIR) < 100) {
            mDCAvsBcZDC->Fill(diffIR, track.first.getDCAMCH());
          }
        }
      }
    }
  }
  if (mSrc[GID::MCHMID] == 1) {
    ILOG(Debug, Devel) << "Debug: MCHMID requested" << ENDM;
    if (mRecoCont.isMatchSourceLoaded(GID::MCHMID)) {
      ILOG(Debug, Devel) << "Debug: MCHMID source loaded" << ENDM;
      mTrackPlotters[GID::MCHMID]->fillHistograms(mRecoCont);
      mTrackPlottersWithCuts[GID::MCHMID]->fillHistograms(mRecoCont);
    }
  }
  if (mSrc[GID::MFTMCH] == 1) {
    ILOG(Debug, Devel) << "Debug: MFTMCH requested" << ENDM;
    if (mRecoCont.isTrackSourceLoaded(GID::MFTMCH)) {
      mTrackPlotters[GID::MFTMCH]->fillHistograms(mRecoCont);
      mTrackPlottersWithCuts[GID::MFTMCH]->fillHistograms(mRecoCont);
    }
  }
  if (mSrc[GID::MFTMCHMID] == 1) {
    ILOG(Debug, Devel) << "Debug: MFTMCHMID requested" << ENDM;
    if (mRecoCont.isTrackSourceLoaded(GID::MFTMCH)) {
      mTrackPlotters[GID::MFTMCHMID]->fillHistograms(mRecoCont);
      mTrackPlottersWithCuts[GID::MFTMCHMID]->fillHistograms(mRecoCont);
    }
  }
}

void TracksTask::endOfCycle()
{
  ILOG(Debug, Devel) << "endOfCycle" << ENDM;
  for (auto& p : mTrackPlotters) {
    p.second->endOfCycle();
  }
  for (auto& p : mTrackPlottersWithCuts) {
    p.second->endOfCycle();
  }
  for (auto& p : mTrackPlottersBgdZDC) {
    p->endOfCycle();
  }
}

void TracksTask::endOfActivity(const Activity& /*activity*/)
{
  ILOG(Debug, Devel) << "endOfActivity" << ENDM;
  removeTrackHistos();
}

void TracksTask::reset()
{
  ILOG(Debug, Devel) << "reset" << ENDM;
  for (auto& p : mTrackPlotters) {
    p.second->reset();
  }
  for (auto& p : mTrackPlottersWithCuts) {
    p.second->reset();
  }
  for (auto& p : mTrackPlottersBgdZDC) {
    p->reset();
  }
}

} // namespace o2::quality_control_modules::muon
