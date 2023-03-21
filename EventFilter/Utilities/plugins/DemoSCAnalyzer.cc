
// system include files
#include <fstream>
#include <iomanip>
#include <memory>

#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/MessageLogger/interface/MessageDrop.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"

using namespace edm;
using namespace std;

class DemoSCAnalyzer: public edm::one::EDAnalyzer<edm::one::SharedResources> {
  
  public:
    explicit DemoSCAnalyzer(const edm::ParameterSet&);
    ~DemoSCAnalyzer() override{};

    void analyze(const edm::Event&, const edm::EventSetup&) override;

  private:
    // EDM tokens:
    edm::EDGetTokenT<l1t::MuonBxCollection> muonToken_;

    int minBx=0;
    int maxBx=3565;

    edm::Service<TFileService> fs;
    TH1D* muonBxOccupancy;
};

DemoSCAnalyzer::DemoSCAnalyzer(const edm::ParameterSet& iConfig){
  muonToken_  = consumes<l1t::MuonBxCollection>(iConfig.getParameter<InputTag>("muInputTag"));

  minBx = iConfig.getParameter<int>("minBx");
  maxBx = iConfig.getParameter<int>("maxBx");

  muonBxOccupancy = fs->make<TH1D>("nMuonsPerBx" , "nMuonsPerBx" , maxBx-minBx , minBx , maxBx );

  //usesResource();
}

void DemoSCAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& evSetup) {
  // handle inputs
  Handle<BXVector<l1t::Muon>> muons;
  iEvent.getByToken(muonToken_, muons);

  for (int bxi = minBx; bxi <= maxBx; ++bxi) {
    //cout << " ========== BX = " << std::dec << i << " =============================" << endl; 

    if (muons->size(bxi)>0){
      muonBxOccupancy->Fill(bxi, muons->size(bxi));  
    }

    //Loop over Muons
    int nObj = 0;
    cout << " ------ Muons --------" << endl;
    if (muons.isValid()) {
      if (bxi >= muons->getFirstBX() && bxi <= muons->getLastBX()) {
        for (std::vector<l1t::Muon>::const_iterator mu = muons->begin(bxi); mu != muons->end(bxi); ++mu) {
          cout << "  " << std::dec << std::setw(2) << std::setfill(' ') << nObj << std::setfill('0') << ")";
          cout << "   Pt " << std::dec << std::setw(3) << mu->hwPt() << " (0x" << std::hex << std::setw(3)
               << std::setfill('0') << mu->hwPt() << ")";
          cout << "   EtaAtVtx " << std::dec << std::setw(3) << mu->hwEtaAtVtx() << " (0x" << std::hex << std::setw(3)
               << std::setfill('0') << (mu->hwEtaAtVtx() & 0x1ff) << ")";
          cout << "   Eta " << std::dec << std::setw(3) << mu->hwEta() << " (0x" << std::hex << std::setw(3)
               << std::setfill('0') << (mu->hwEta() & 0x1ff) << ")";
          cout << "   PhiAtVtx " << std::dec << std::setw(3) << mu->hwPhiAtVtx() << " (0x" << std::hex << std::setw(3)
               << std::setfill('0') << mu->hwPhiAtVtx() << ")";
          cout << "   Phi " << std::dec << std::setw(3) << mu->hwPhi() << " (0x" << std::hex << std::setw(3)
               << std::setfill('0') << mu->hwPhi() << ")";
          cout << "   Iso " << std::dec << std::setw(1) << mu->hwIso();
          cout << "   Qual " << std::dec << std::setw(1) << mu->hwQual();
          cout << endl;
          nObj++;
        }
      } else {
        cout << "No Muons stored for this bx " << bxi << endl;
      }
    } else {
      cout << "No Muon Data in this event " << endl;
    }

  } // end BX loop
  
}

DEFINE_FWK_MODULE(DemoSCAnalyzer);
