#include <DataFormats/Common/interface/RefProd.h>
#include "DataFormats/Common/interface/Wrapper.h"

#include <DataFormats/L1Scouting/interface/SRawDataCollection.h>
#include "DataFormats/L1Scouting/interface/OrbitCollection.h"

namespace scoutingRun3 {
  edm::Wrapper<scoutingRun3::OrbitCollection<l1t::Muon>>            MuonOrbitCollectionWrapper;
  edm::Wrapper<scoutingRun3::OrbitCollection<l1t::Jet>>             JetOrbitCollectionWrapper;
  edm::Wrapper<scoutingRun3::OrbitCollection<l1t::EGamma>>          EGammaOrbitCollectionWrapper;
  edm::Wrapper<scoutingRun3::OrbitCollection<l1t::Tau>>             TauOrbitCollectionWrapper;
  edm::Wrapper<scoutingRun3::OrbitCollection<l1t::EtSum>>           EtSumOrbitCollectionWrapper;
  edm::Wrapper<scoutingRun3::OrbitCollection<L1MuKBMTCombinedStub>> BmtfStubOrbitCollectionWrapper;
}