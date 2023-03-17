#include "EventFilter/L1ScoutingRawToDigi/plugins/ScGMTRawToDigi.h"

ScGMTRawToDigi::ScGMTRawToDigi(const edm::ParameterSet& iConfig) {
  using namespace edm;
  srcInputTag  = iConfig.getParameter<InputTag>( "srcInputTag" );
  debug = iConfig.getUntrackedParameter<bool>("debug", false);

  produces<l1t::MuonBxCollection>().setBranchAlias( "MuonBxCollection" );
  rawToken = consumes<SRDCollection>(srcInputTag);
  
  //  bx_muons.reserve(8);
  dummyLVec_.reset( new ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>() );
}

ScGMTRawToDigi::~ScGMTRawToDigi() {};

void ScGMTRawToDigi::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  Handle<SRDCollection> ScoutingRawDataCollection;
  iEvent.getByToken( rawToken, ScoutingRawDataCollection );

  const FEDRawData& sourceRawData = ScoutingRawDataCollection->FEDData(SDSNumbering::GmtSDSID);
  size_t orbitSize = ScoutingRawDataCollection->getSourceSize(SDSNumbering::GmtSDSID);

  std::unique_ptr<l1t::MuonBxCollection> unpackedMuons(new l1t::MuonBxCollection);
  
  if((sourceRawData.size()==0) && debug){
    std::cout << "No raw data for GMT FED\n";  
  }

  unpackOrbit(unpackedMuons.get(), sourceRawData.data(), orbitSize); 

  // store collection in the event
  iEvent.put( std::move(unpackedMuons) );
}

void ScGMTRawToDigi::unpackOrbit(
  l1t::MuonBxCollection* muons, 
  const unsigned char* buf, size_t len
  ){
  
  using namespace scoutingRun3;
  size_t pos = 0;

  muons->setBXRange(0,3565);
  
  while (pos < len) {
    assert(pos+4 <= len);
    
    // get BX header
    uint32_t header = *((uint32_t*)(buf + pos));
    pos += 4;
    // count mA and mB
    uint32_t mAcount = (header & header_masks::mAcount) >> header_shifts::mAcount;
    uint32_t mBcount = (header & header_masks::mBcount) >> header_shifts::mBcount;
    
    // declare block to read
    ugmt::block *bl = (ugmt::block *)(buf + pos);
    pos += 4 + 4 + (mAcount+mBcount)*12; 
    assert(pos <= len);

    uint32_t orbit = bl->orbit & 0x7FFFFFFF;  
    uint32_t bx = bl->bx;
   
    if (debug){
      std::cout  << " GMT Orbit " << orbit << ", BX -> "<< bx << ", nMuons -> " << mAcount+mBcount << std::endl;
    }
    
    // Unpack muons for this BX
    
    // bx_muons.clear();
    
    // cuts should be applied
    bool excludeIntermediate=true;
    int ptcut=0;
    unsigned int qualcut=0;

    for (unsigned int i=0; i<mAcount+mBcount; i++) {

      uint32_t interm = (bl->mu[i].extra >> ugmt::shiftsMuon::interm) & ugmt::masksMuon::interm;
      if (excludeIntermediate && (interm == 1)) continue;

      uint32_t index    = (bl->mu[i].s >> ugmt::shiftsMuon::index)  & ugmt::masksMuon::index;
      uint32_t ietaextu = (bl->mu[i].f >> ugmt::shiftsMuon::etaext) & ugmt::masksMuon::etaextv;
      int32_t ietaext;
      if (((bl->mu[i].f >> ugmt::shiftsMuon::etaext) & ugmt::masksMuon::etaexts)!=0) {
          ietaext = ietaextu -= 256;
      } else {
          ietaext = ietaextu;
      }
      
      // extract pt and quality and apply cut if required
      int32_t iptuncon = (bl->mu[i].s >> ugmt::shiftsMuon::ptuncon) & ugmt::masksMuon::ptuncon;
      int32_t ipt      = (bl->mu[i].f >> ugmt::shiftsMuon::pt)      & ugmt::masksMuon::pt;
      if ((ipt-1) < ptcut) {
          continue;
      }
      uint32_t qual = (bl->mu[i].f >> ugmt::shiftsMuon::qual) & ugmt::masksMuon::qual;
      if (qual < qualcut) {
          continue;
      }
      
      // extract integer value for extrapolated phi
      int32_t iphiext = ((bl->mu[i].f >> ugmt::shiftsMuon::phiext) & ugmt::masksMuon::phiext);

      // extract integer value for extrapolated phi
      int32_t idxy = ((bl->mu[i].s >> ugmt::shiftsMuon::dxy) & ugmt::masksMuon::dxy);

      // extract iso bits and charge
      uint32_t iso = (bl->mu[i].s >> ugmt::shiftsMuon::iso) & ugmt::masksMuon::iso;
      int32_t chrg = 0;
      if (((bl->mu[i].s >> ugmt::shiftsMuon::chrgv) & ugmt::masksMuon::chrgv)==1)
          chrg=((bl->mu[i].s >> ugmt::shiftsMuon::chrg) & ugmt::masksMuon::chrg)==1 ? -1 : 1 ;

      // extract eta and phi at muon station
      int32_t  iphi  = (bl->mu[i].s >> ugmt::shiftsMuon::phi)      & ugmt::masksMuon::phi;
      uint32_t ieta1 = (bl->mu[i].extra >> ugmt::shiftsMuon::eta1) & ugmt::masksMuon::eta;
      uint32_t ieta2 = (bl->mu[i].extra >> ugmt::shiftsMuon::eta2) & ugmt::masksMuon::eta;


      uint32_t ieta_u;
      int32_t ieta;
      // checking if raw eta should be taken from muon 1 or muon 2
      if ( (bl->mu[i].extra & 0x1) == 0 ) {
          ieta_u = ieta1;
      } else {
          ieta_u = ieta2;
      }

      // two's complement
      if ( ieta_u > 256 ) {
          ieta = ieta_u - 512;
      } else {
          ieta = ieta_u;
      }

      // convert to physical units using scales
      //float fpt      = (ipt     -1) * ugmt::scales::pt_scale;                 // -1 since bin 0 is for invalid muons
      float fptuncon = (iptuncon-1) * ugmt::scales::ptunconstrained_scale;    // -1 since bin 0 is for invalid muons
      float fphi     = iphi         * ugmt::scales::phi_scale;
      float fphiext  = iphiext      * ugmt::scales::phi_scale;
      //float feta     = ieta         * ugmt::scales::eta_scale;
      float fetaext  = ietaext      * ugmt::scales::eta_scale;

      if (fphiext>M_PI) fphiext = fphiext - 2.*M_PI;
      if (fphi   >M_PI) fphi    = fphi    - 2.*M_PI;


      l1t::Muon muon(
          *dummyLVec_, ipt, ieta, iphi, qual, chrg, chrg != 0, iso,
          index, 0, false, 0, 0, 0, 0,
          ietaext, iphiext, fetaext, fphiext,
          iptuncon, fptuncon, idxy
      ); 

      
      //bx_muons.push_back(muon);
      muons->push_back(bx, muon);

    } // end of bx
    
    // add muons to the collection
    //muons->push_back(bx, bx_muons);

  } // end orbit while loop

}

void ScGMTRawToDigi::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScGMTRawToDigi);
