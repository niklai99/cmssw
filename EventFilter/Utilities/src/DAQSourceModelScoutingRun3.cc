#include "EventFilter//Utilities/interface/DAQSourceModelsScoutingRun3.h"

void DataModeScoutingRun3::makeDirectoryEntries(std::vector<std::string> const& baseDirs, std::string const& runDir) {
  std::filesystem::path runDirP(runDir);
  for (auto& baseDir : baseDirs) {
    std::filesystem::path baseDirP(baseDir);
    buPaths_.emplace_back(baseDirP / runDirP);
  }
}

std::pair<bool, std::vector<std::string>> DataModeScoutingRun3::defineAdditionalFiles(std::string const& primaryName,
                                                                                  bool fileListMode) const {
  std::vector<std::string> additionalFiles;

  if (fileListMode) {
    //for the unit test
    additionalFiles.push_back(primaryName + "_1");
    return std::make_pair(true, additionalFiles);
  }

  auto fullpath = std::filesystem::path(primaryName);
  auto fullname = fullpath.filename();

  for (size_t i = 1; i < buPaths_.size(); i++) {
    std::filesystem::path newPath = buPaths_[i] / fullname;
    additionalFiles.push_back(newPath.generic_string());
  }
  return std::make_pair(true, additionalFiles);
}

void DataModeScoutingRun3::readEvent(edm::EventPrincipal& eventPrincipal) {
  assert(!events_.empty());

  edm::TimeValue_t time;
  timeval stv;
  gettimeofday(&stv, nullptr);
  time = stv.tv_sec;
  time = (time << 32) + stv.tv_usec;
  edm::Timestamp tstamp(time);

  // set provenance helpers
  uint32_t hdrEventID = currOrbit;
  edm::EventID eventID = edm::EventID(daqSource_->eventRunNumber(), daqSource_->currentLumiSection(), hdrEventID);
  edm::EventAuxiliary aux(eventID, daqSource_->processGUID(), tstamp, events_[0]->isRealData(), edm::EventAuxiliary::PhysicsTrigger);

  aux.setProcessHistoryID(daqSource_->processHistoryID());
  daqSource_->makeEventWrapper(eventPrincipal, aux);

  // create scouting raw data collection
  std::unique_ptr<SRDCollection> rawData(new SRDCollection);

  // loop over selected sources
  for(auto const& event: events_) {
    fillSRDCollection(*rawData, (char*)event->payload(), event->eventSize());
  }  

  std::unique_ptr<edm::WrapperBase> edp(new edm::Wrapper<SRDCollection>(std::move(rawData)));
  eventPrincipal.put(daqProvenanceHelpers_[0]->branchDescription(), std::move(edp), daqProvenanceHelpers_[0]->dummyProvenance());
 
  eventCached_ = false;
}

void DataModeScoutingRun3::fillSRDCollection(
    SRDCollection& rawData, char* buff, size_t len
  ){
  
  size_t pos = 0;

  // get the source ID
  // TODO: it would be better to include it in the header
  int sourceId  = *((uint32_t*)(buff + pos));
  pos += 4;

  // size of the orbit paylod
  size_t orbitSize = len-pos;

  // set the size (=orbit size) in the SRDColletion
  // of the current source. Currently is not possible to use
  // the FRD size because it needs to be resized, since it 
  // is expecting to receive payloads which are multiple of 8 bytes.
  // (not true for scouting payloads)
  FEDRawData& fedData = rawData.FEDData(sourceId);
  fedData.resize(orbitSize, 4);
  //rawData.setSourceSize(sourceId, orbitSize);

  // FEDRawData is expecting multiples of 8 bytes,
  // add a padding.
  //fedData.resize(orbitSize + orbitSize%8);
  memcpy(fedData.data(), buff+pos, orbitSize);

  return;
}

std::vector<std::shared_ptr<const edm::DaqProvenanceHelper>>& DataModeScoutingRun3::makeDaqProvenanceHelpers() {
  //set FRD data collection
  daqProvenanceHelpers_.clear();
  daqProvenanceHelpers_.emplace_back(std::make_shared<const edm::DaqProvenanceHelper>(
      edm::TypeID(typeid(SRDCollection)), "SRDCollection", "SRDCollection", "DAQSource"));
  return daqProvenanceHelpers_;
}

bool DataModeScoutingRun3::nextEventView() {
  blockCompleted_ = false;
  if (eventCached_)
    return true;

  int evt_idx = 0;
  for (int i = 0; i < numFiles_; i++){
    // add last event length if it was a valid orbit
    // i.e. it has been processed
    if (validOrbits_[i]){
      dataBlockAddrs_[i] += events_[evt_idx]->size();
      evt_idx++;
    }
  }

  return makeEvents();
}

bool DataModeScoutingRun3::makeEvents() {
  // clear events and reset current orbit
  events_.clear();
  currOrbit = 0xFFFFFFFF; // max uint
  assert(!blockCompleted_);

  // create current "events" (= orbits) list from each data source
  // Check if one dataBlock terminated earlier than others.
  for (int i = 0; i < numFiles_; i++) {
    if (dataBlockAddrs_[i] >= dataBlockMaxAddrs_[i]) {
      completedBlocks_[i] = true;
      continue;
    } 

    events_.emplace_back(std::make_unique<FRDEventMsgView>(dataBlockAddrs_[i]));
    if (dataBlockAddrs_[i] + events_.back()->size() > dataBlockMaxAddrs_[i])
      throw cms::Exception("DAQSource::getNextEvent")
          << " event id:" << events_.back()->event() << " lumi:" << events_.back()->lumi() << " run:" << events_.back()->run()
          << " of size:" << events_.back()->size() << " bytes does not fit into the buffer or has corrupted header";
    
    // find the min orbit for the current event
    if ((events_.back()->event()<currOrbit) && (!completedBlocks_[i])){
      currOrbit = events_.back()->event();
    }
  }

  // mark valid orbits from each data source
  // = find when orbit is missing from one source
  bool allBlocksCompleted = true;
  int evt_idx = 0;
  for (int i=0; i < numFiles_; i++){
    if (completedBlocks_[i]) {
      validOrbits_[i] = false;
      continue;
    }

    if(events_[evt_idx]->event() != currOrbit){
      validOrbits_[i] = false;
    } else {
      validOrbits_[i] = true;
      allBlocksCompleted = false;
    }

    evt_idx++;
  }

  if (allBlocksCompleted){
    blockCompleted_ = true;
  }
  return !allBlocksCompleted;
}

bool DataModeScoutingRun3::checksumValid() { return true; }

std::string DataModeScoutingRun3::getChecksumError() const { return std::string(); }
