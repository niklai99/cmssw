#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/L1Scouting/interface/SRawDataCollection.h"
#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/L1Scouting/interface/SDSNumbering.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"

#include "DataFormats/L1TMuon/interface/L1MuKBMTCombinedStub.h"

#include "EventFilter/L1ScoutingRawToDigi/interface/shifts.h"
#include "EventFilter/L1ScoutingRawToDigi/interface/scales.h"
#include "EventFilter/L1ScoutingRawToDigi/interface/masks.h"
#include "EventFilter/L1ScoutingRawToDigi/interface/blocks.h"

class ScBMTFRawToDigi : public edm::stream::EDProducer<> {
public:
  explicit ScBMTFRawToDigi(const edm::ParameterSet&);
  ~ScBMTFRawToDigi() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  //void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  //void endStream() override;

  void unpackOrbit(
    scoutingRun3::BmtfStubOrbitCollection* stubs,
    const unsigned char* buf, size_t len,
    int SDSID
  );

  int calculateEta(uint i, int wheel, uint sector, uint station);

  L1MuKBMTCombinedStub buildStub(int wheel, int sector, int station,
                                 int phi, int phiB, bool tag,
                                 int eta, int qeta, int bx,
                                 int quality);

  L1MuKBMTCombinedStub buildStubNoEta(int wheel, int sector, int station,
                                      int phi, int phiB, bool tag,
                                      int bx, int quality);

  std::vector<L1MuKBMTCombinedStub> bx_stubs;
  std::unique_ptr<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>> dummyLVec_;

  bool debug = false;

  edm::InputTag srcInputTag;
  edm::EDGetToken rawToken;

  std::vector<int> eta1_;
  std::vector<int> eta2_;
  std::vector<int> eta3_;
};
