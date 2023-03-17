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
    std::cout << "Added new file: " << primaryName + "_1" << std::endl;
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
  std::cout << "Read new event!" << std::endl;

  edm::TimeValue_t time;
  timeval stv;
  gettimeofday(&stv, nullptr);
  time = stv.tv_sec;
  time = (time << 32) + stv.tv_usec;
  edm::Timestamp tstamp(time);

  // set provenance helpers
  uint32_t hdrEventID = events_[0]->event();//currOrbit;
  edm::EventID eventID = edm::EventID(daqSource_->eventRunNumber(), daqSource_->currentLumiSection(), hdrEventID);
  edm::EventAuxiliary aux(eventID, daqSource_->processGUID(), tstamp, events_[0]->isRealData(), edm::EventAuxiliary::PhysicsTrigger);

  aux.setProcessHistoryID(daqSource_->processHistoryID());
  daqSource_->makeEventWrapper(eventPrincipal, aux);

  // create scouting raw data collection
  std::unique_ptr<SDSRawDataCollection> rawData(new SDSRawDataCollection);

  std::unique_ptr<edm::WrapperBase> edp(new edm::Wrapper<SDSRawDataCollection>(std::move(rawData)));
  eventPrincipal.put(daqProvenanceHelpers_[0]->branchDescription(), std::move(edp), daqProvenanceHelpers_[0]->dummyProvenance());
 
  eventCached_ = false;
}

std::vector<std::shared_ptr<const edm::DaqProvenanceHelper>>& DataModeScoutingRun3::makeDaqProvenanceHelpers() {
  //set FRD data collection
  daqProvenanceHelpers_.clear();
  daqProvenanceHelpers_.emplace_back(std::make_shared<const edm::DaqProvenanceHelper>(
      edm::TypeID(typeid(SDSRawDataCollection)), "SDSRawDataCollection", "SDSRawDataCollection", "DAQSource"));
  return daqProvenanceHelpers_;
}

bool DataModeScoutingRun3::nextEventView() {
  blockCompleted_ = false;
  if (eventCached_)
    return true;
  for (unsigned int i = 0; i < events_.size(); i++) {
    //add last event length to each stripe
    dataBlockAddrs_[i] += events_[i]->size();
  }
  return makeEvents();
}

bool DataModeScoutingRun3::makeEvents() {
  events_.clear();
  assert(!blockCompleted_);
  for (int i = 0; i < numFiles_; i++) {
    if (dataBlockAddrs_[i] >= dataBlockMaxAddrs_[i]) {
      //must be exact
      assert(dataBlockAddrs_[i] == dataBlockMaxAddrs_[i]);
      blockCompleted_ = true;
      return false;
    } else {
      if (blockCompleted_)
        throw cms::Exception("DataModeScoutingRun3::makeEvents")
            << "not all blocks were completed at the same time";
    }
    if (blockCompleted_)
      continue;
    events_.emplace_back(std::make_unique<FRDEventMsgView>(dataBlockAddrs_[i]));
    if (dataBlockAddrs_[i] + events_[i]->size() > dataBlockMaxAddrs_[i])
      throw cms::Exception("DAQSource::getNextEvent")
          << " event id:" << events_[i]->event() << " lumi:" << events_[i]->lumi() << " run:" << events_[i]->run()
          << " of size:" << events_[i]->size() << " bytes does not fit into the buffer or has corrupted header";
  }
  return !blockCompleted_;
}

bool DataModeScoutingRun3::checksumValid() { return true; }

std::string DataModeScoutingRun3::getChecksumError() const { return std::string(); }
