#include "EventFilter/L1ScoutingRawToDigi/plugins/ScBMTFRawToDigi.h"



ScBMTFRawToDigi::ScBMTFRawToDigi(const edm::ParameterSet& iConfig) {
  using namespace edm;
  srcInputTag  = iConfig.getParameter<InputTag>( "srcInputTag" );
  debug = iConfig.getUntrackedParameter<bool>("debug", false);

  // produces<l1t::MuonBxCollection>().setBranchAlias( "MuonBxCollection" );
  produces<scoutingRun3::BmtfStubOrbitCollection>().setBranchAlias( "BmtfStubOrbitCollection" );
  rawToken = consumes<SRDCollection>(srcInputTag);

  dummyLVec_.reset( new ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>() );

  eta1_ = iConfig.getParameter<std::vector<int>>("cotTheta_1");
  eta2_ = iConfig.getParameter<std::vector<int>>("cotTheta_2");
  eta3_ = iConfig.getParameter<std::vector<int>>("cotTheta_3");
}



ScBMTFRawToDigi::~ScBMTFRawToDigi() {};



void ScBMTFRawToDigi::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  Handle<SRDCollection> ScoutingRawDataCollection;
  iEvent.getByToken( rawToken, ScoutingRawDataCollection );

  // create pointer to bmtf stubs collection
  std::unique_ptr<scoutingRun3::BmtfStubOrbitCollection> unpackedStubs(new scoutingRun3::BmtfStubOrbitCollection);

  // loop over all Scouting Data Source IDs
  for (unsigned int i=SDSNumbering::BmtfMinSDSID; i<SDSNumbering::BmtfMaxSDSID; i++) {
    // get data and orbit size from i^th source
    const FEDRawData& sourceRawData = ScoutingRawDataCollection->FEDData(i);
    size_t orbitSize = sourceRawData.size();

    if((sourceRawData.size()==0) && debug){
      std::cout << "No raw data for BMTF FED " << i << std::endl;
    }

    unpackOrbit(unpackedStubs.get(), sourceRawData.data(), orbitSize, i);
  }

  unpackedStubs.get()->flatten();

  // store collection in the event
  iEvent.put( std::move(unpackedStubs) );
}



void ScBMTFRawToDigi::unpackOrbit(
  //l1t::MuonBxCollection* muons,
  scoutingRun3::BmtfStubOrbitCollection* stubs,
  const unsigned char* buf, size_t len,
  int SDSID
  ){

  using namespace scoutingRun3;
  size_t pos = 0;

  //muons->setBXRange(0,3565);

  while (pos < len) {
    assert(pos+4 <= len);

    // get BX header
    uint32_t header = *((uint32_t*)(buf + pos));
    pos += 4;
    // decode header
    uint32_t bx     = (header & 0xffff0000) >> 16;
    uint32_t sCount = (header & 0x000000ff);
    // decode orbit number
    uint32_t orbit = *((uint32_t*)(buf + pos));
    pos += 4;
    orbit &= 0x7FFFFFFF;

    // declare block to read
    bmtf::block *bl = (bmtf::block *)(buf + pos);
    pos += sCount*8;
    assert(pos <= len);

    if (debug){
      std::cout  << " BMTF #" << SDSID << " Orbit " << orbit << ", BX -> "<< bx << ", nStubs -> " << sCount << std::endl;
    }

    // Unpack stubs for this BX
    // TODO: apply filters
    // int32_t valid, phi, phiB, tag, qual, eta, qeta, station, wheel, reserved, sector;
    int32_t phi, phiB, tag, qual, eta, qeta, station, wheel, sector;
    // map for station and wheel, to find chambers with 2 stubs
    std::vector<std::vector<bool>> stwh_matrix(4, std::vector<bool>(5,false));
    for (unsigned int i=0; i<sCount; i++) {

      // shifts and masks
      // valid    = ((bl->stub[i] >> bmtf::shiftsStubs::valid   ) & bmtf::masksStubs::valid   );
      phi      = ((bl->stub[i] >> bmtf::shiftsStubs::phi     ) & bmtf::masksStubs::phi     );
      phiB     = ((bl->stub[i] >> bmtf::shiftsStubs::phiB    ) & bmtf::masksStubs::phiB    );
      qual     = ((bl->stub[i] >> bmtf::shiftsStubs::qual    ) & bmtf::masksStubs::qual    );
      eta      = ((bl->stub[i] >> bmtf::shiftsStubs::eta     ) & bmtf::masksStubs::eta     );
      qeta     = ((bl->stub[i] >> bmtf::shiftsStubs::qeta    ) & bmtf::masksStubs::qeta    );
      station  = ((bl->stub[i] >> bmtf::shiftsStubs::station ) & bmtf::masksStubs::station ) + 1;
      wheel    = ((bl->stub[i] >> bmtf::shiftsStubs::wheel   ) & bmtf::masksStubs::wheel   );
      // reserved = ((bl->stub[i] >> bmtf::shiftsStubs::reserved) & bmtf::masksStubs::reserved);
      sector   = static_cast<uint16_t>(SDSID - SDSNumbering::BmtfMinSDSID);


      // trick to set Ts2tag for L1MuKBMTFCombinedStub
      if (stwh_matrix[station-1][wheel+2]==false) {
        tag = 1;
      } else {
        tag = 0;
      }
      stwh_matrix[station-1][wheel+2] = true;

      phi      = phi   >= 2048 ? phi   - 4096 : phi;
      phiB     = phiB  >=  512 ? phiB  - 1024 : phiB;
      wheel    = wheel >=    4 ? wheel -    8 : wheel;


      if (eta==0) {
        L1MuKBMTCombinedStub comb_stub = buildStubNoEta(wheel, sector, station, phi, phiB, tag, bx, qual);
        stubs->push_back(bx, comb_stub);
      } else {
        L1MuKBMTCombinedStub comb_stub = buildStub(wheel, sector, station, phi, phiB, tag, eta, qeta, bx, qual);
        stubs->push_back(bx, comb_stub);
      }

    } // end of bx

  } // end orbit while loop

}



int ScBMTFRawToDigi::calculateEta(uint i, int wheel, uint sector, uint station) {
  int eta = 0;
  if (wheel > 0) {
    eta = 7 * wheel + 3 - i;
  } else if (wheel < 0) {
    eta = 7 * wheel + i - 3;
  } else {
    if (sector == 0 || sector == 3 || sector == 4 || sector == 7 || sector == 8 || sector == 11)
      eta = i - 3;
    else
      eta = 3 - i;
  }

  if (station == 1)
    eta = -eta1_[eta + 17];
  else if (station == 2)
    eta = -eta2_[eta + 17];
  else
    eta = -eta3_[eta + 17];

  return eta;
}



L1MuKBMTCombinedStub ScBMTFRawToDigi::buildStub(int wheel, int sector, int station,
                                                int phi, int phiB, bool tag,
                                                int eta, int qeta, int bx,
                                                int quality) {
  // convert eta hw values to global units
  int qeta1 = 0;
  int qeta2 = 0;
  int eta1 = 255;
  int eta2 = 255;
  int mask = 0;

  bool hasEta = false;
  for (uint i = 0; i < 7; ++i) {
    mask = (i << i);
    if ((eta & mask) == 0)
      continue;
    if (!hasEta) {
      hasEta = true;
      eta1 = calculateEta(i, wheel, sector, station);
      if ((qeta & mask) == 1)
        qeta1 = 2;
      else
        qeta1 = 1;
    } else {
      eta2 = calculateEta(i, wheel, sector, station);
      if ((qeta & mask) == 1)
        qeta2 = 2;
      else
        qeta2 = 1;
    }
  }

  // TODO: tag = true is hardcode for now
  L1MuKBMTCombinedStub stub(wheel, sector, station, phi, phiB, tag, bx, quality, eta1, eta2, qeta1, qeta2);

  return stub;
}



L1MuKBMTCombinedStub ScBMTFRawToDigi::buildStubNoEta(int wheel, int sector, int station,
                                                     int phi, int phiB, bool tag,
                                                     int bx, int quality) {

  int qeta1 = 0;
  int qeta2 = 0;
  int eta1 = 7;
  int eta2 = 7;
  L1MuKBMTCombinedStub stub(wheel, sector, station, phi, phiB, tag, bx, quality, eta1, eta2, qeta1, qeta2);

  return stub;
}



void ScBMTFRawToDigi::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScBMTFRawToDigi);
