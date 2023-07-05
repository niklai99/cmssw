#include "FWCore/Framework/interface/MakerMacros.h"

// system include files
#include <fstream>
#include <iomanip>
#include <memory>

// user include files
//   base class
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "DataFormats/L1Trigger/interface/Muon.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/MessageLogger/interface/MessageDrop.h"

#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/L1Scouting/interface/SRawDataCollection.h"
#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/L1Scouting/interface/SDSNumbering.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"

#include "DataFormats/L1TMuon/interface/L1MuKBMTCombinedStub.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"



// class declaration
class KBmtfMuonAnalysis : public edm::stream::EDAnalyzer<> {
  public:
    explicit KBmtfMuonAnalysis(const edm::ParameterSet&);
    ~KBmtfMuonAnalysis() override{};
    void analyze(const edm::Event&, const edm::EventSetup&) override;

    unsigned int calcGlobalPhi(const l1t::RegionalMuonCand*);
    double calcDr(const l1t::RegionalMuonCand*, const l1t::Muon*);

    void updateTotGmtM() { ++totGmtM_; };
    void updateTotMatches() { ++totMatches_; };
    int getTotGmtM() { return totGmtM_; };
    int getTotMatches() { return totMatches_; };

  private:
    edm::EDGetTokenT<scoutingRun3::MuonOrbitCollection> gmtMuonToken_;
    edm::EDGetTokenT<l1t::RegionalMuonCandBxCollection> bmtfMuonToken_;

    double drCut_;
    double phiMult_;
    double etaMult_;

    int minBx_;
    int maxBx_;
    bool debug_;

    int totGmtM_;
    int totMatches_;

    string filenameResults_;
    std::ofstream fileResults_;

};



KBmtfMuonAnalysis::KBmtfMuonAnalysis(const edm::ParameterSet& iConfig)
    : gmtMuonToken_(consumes<scoutingRun3::MuonOrbitCollection>(iConfig.getParameter<edm::InputTag>("gmtMuonInputTag"))),
      bmtfMuonToken_(consumes<l1t::RegionalMuonCandBxCollection>(iConfig.getParameter<edm::InputTag>("bmtfMuonInputTag"))),
      drCut_(iConfig.getParameter<double>("drCut")),
      phiMult_(iConfig.getParameter<double>("phiMult")),
      etaMult_(iConfig.getParameter<double>("etaMult")),
      minBx_(iConfig.getParameter<int>("minBx")),
      maxBx_(iConfig.getParameter<int>("maxBx")),
      debug_(iConfig.getParameter<bool>("debug")),
      totGmtM_(0),
      totMatches_(0),
      filenameResults_(iConfig.getParameter<string>("filenameResults")),
      fileResults_(filenameResults_)
{

}



// loop over events
void KBmtfMuonAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& evSetup) {
  using namespace edm;

  //input
  Handle<scoutingRun3::MuonOrbitCollection> gmtMuons;
  Handle<l1t::RegionalMuonCandBxCollection> bmtfMuons;
  iEvent.getByToken(gmtMuonToken_, gmtMuons);
  iEvent.getByToken(bmtfMuonToken_, bmtfMuons);

  if (debug_) {
    std::cout << " -----------------------------------------------------  " << std::endl;
    std::cout << " *********** Run " << std::dec << iEvent.id().run() << " Event " << iEvent.id().event()
              << " **************  " << std::endl;
    std::cout << " ----------------------------------------------------- " << std::endl;
  }

  const std::vector<int>* gmtMuonsIndex = gmtMuons->getIndex();
  int bx = 0;
  double l1dr = 0.0, l1dr_min = 0.0;
  int l1_match_i = -1, l1_i = -1;
  int n_matches = 0;
  int n_gmt_m = 0;
  double l1_reg_m_physPhi, l1_reg_m_physEta;
  double l1_gmt_m_physPhi, l1_gmt_m_physEta;
  for (int i=0; i < gmtMuons->sizeFlatData(); ++i) {
    const l1t::Muon *gmt_m = gmtMuons->getFlatData(i);

    // barrel gmt muons
    if ((gmt_m->tfMuonIndex()>=36) && (gmt_m->tfMuonIndex()<=70)) {
      ++n_gmt_m;
      updateTotGmtM();

      // loop over BMTF muons in same BX
      l1dr_min = 100.0;
      l1_match_i = -1;
      l1_i = -1;
      l1_reg_m_physPhi = -100.0;
      l1_gmt_m_physPhi = gmt_m->hwPhi()/phiMult_;
      l1_reg_m_physEta = -100.0;
      l1_gmt_m_physEta = gmt_m->hwEta()/etaMult_;
      // for (std::vector<l1t::RegionalMuonCand>::const_iterator bmtf_m=bmtfMuons->begin(bx-1); bmtf_m!=bmtfMuons->end(bx-1); ++bmtf_m) {
      for (int bx=0; bx<maxBx_; ++bx) {
        if (bmtfMuons->size(bx)==0) continue;

        for (size_t k=0; k<bmtfMuons->size(bx); ++k) {
          const l1t::RegionalMuonCand *bmtf_m = &(bmtfMuons->at(bx, k));
          ++l1_i;
          l1dr = calcDr(bmtf_m, gmt_m);
          if (l1dr < l1dr_min) {
            l1_match_i = l1_i;
            l1dr_min = l1dr;
            l1_reg_m_physPhi = calcGlobalPhi(bmtf_m)/phiMult_;
            l1_reg_m_physEta = bmtf_m->hwEta()/etaMult_;
          }
        }
      }

      if (l1dr_min < drCut_) {
        ++n_matches;
        updateTotMatches();
        std::cout << "Match: " << std::endl;
        std::cout << "    dr = " << l1dr_min << "    l1_match " << l1_match_i
                  << "    #gmt muons: " << n_gmt_m << "    #matches " << n_matches
                  << "    #Tot gmt muons: " << getTotGmtM() << "    #Tot matches " << getTotMatches()
                  << std::endl;

        fileResults_ << l1dr_min << ","
                     << l1_reg_m_physPhi << ","
                     << l1_gmt_m_physPhi << ","
                     << l1_reg_m_physEta << ","
                     << l1_gmt_m_physEta << std::endl;

      }
    }

  }

  // loop over BX
  // std::cout << "Check 1 " << gmtMuonsIndex->size() << " " << gmtMuons->sizeFlatData() << std::endl;
  // if (gmtMuonsIndex->size()!=0) {
  //   for (size_t i=0; i < gmtMuonsIndex->size()-1; ++i) {
  //     // get BX number
  //     bx = i;
  //     std::cout << "Check 2-1 " << gmtMuons->getIndex(i) << std::endl;
  //     std::cout << "Check 2-1 " << gmtMuons->getIndex(i+1) << std::endl;
  //     // loop over gmt muons in BX
  //     for (int j=gmtMuons->getIndex(i); j<gmtMuons->getIndex(i+1); ++j) {
  //       std::cout << "Check 2-1-1" << std::endl;
  //       const l1t::Muon *gmt_m = gmtMuons->getFlatData(j);
  //       std::cout << "Check 2-1-2" << std::endl;

  //       // barrel gmt muons
  //       if ((gmt_m->tfMuonIndex()>=36) && (gmt_m->tfMuonIndex()<=70)) {
  //         std::cout << "Check 2-1-2-1" << std::endl;
  //         ++n_gmt_m;

  //         // loop over BMTF muons in same BX
  //         l1dr_min = 100.0;
  //         l1_match_i = -1;
  //         l1_i = -1;
  //         // for (std::vector<l1t::RegionalMuonCand>::const_iterator bmtf_m=bmtfMuons->begin(bx-1); bmtf_m!=bmtfMuons->end(bx-1); ++bmtf_m) {
  //         for (size_t k=0; k<=bmtfMuons->size(bx-1); ++k) {
  //           std::cout << "Check 2-1-2-1-1" << std::endl;
  //           const l1t::RegionalMuonCand *bmtf_m = &(bmtfMuons->at(bx-1, k));
  //           std::cout << "Check 2-1-2-1-2" << std::endl;
  //           ++l1_i;
  //           l1dr = calcDr(bmtf_m, gmt_m);
  //           std::cout << "Check 2-1-2-1-3" << std::endl;
  //           if (l1dr < l1dr_min) {
  //             l1_match_i = l1_i;
  //             l1dr_min = l1dr;
  //           }
  //         }
  //         std::cout << "Check 2-1-2-2" << std::endl;

  //         if (l1dr_min < drCut_) {
  //           ++n_matches;
  //           std::cout << "Match: " << std::endl;
  //           std::cout << "    dr = " << l1dr_min << "    l1_match " << l1_match_i
  //                     << "    #gmt muons: " << n_gmt_m << "    #matches " << n_matches
  //                     << std::endl;
  //         }
  //       }
  //       std::cout << "Check 2-1-3" << std::endl;

  //     }

  //     std::cout << "Check 2-2" << std::endl;
  //   }

  //   std::cout << "Check 3" << std::endl;
  // }

}



unsigned int KBmtfMuonAnalysis::calcGlobalPhi(const l1t::RegionalMuonCand *l1_reg_m) {

  unsigned int globalPhi = l1_reg_m->processor()*48 + l1_reg_m->hwPhi();
  globalPhi += 576 - 24;      // first processor starts at -15degrees in cms phi
  globalPhi = globalPhi%576;  // wrap around

  return globalPhi;
}



double KBmtfMuonAnalysis::calcDr(const l1t::RegionalMuonCand *l1_reg_m, const l1t::Muon *l1_gmt_m) {

  double l1_reg_m_physPhi = calcGlobalPhi(l1_reg_m)/phiMult_;
  if (l1_reg_m_physPhi > M_PI) { l1_reg_m_physPhi -= 2*M_PI; }

  double l1_gmt_m_physPhi = l1_gmt_m->hwPhi()/phiMult_;
  if (l1_gmt_m_physPhi > M_PI) { l1_gmt_m_physPhi -= 2*M_PI; }

  double dphi = abs(l1_reg_m_physPhi - l1_gmt_m_physPhi);
  if (dphi > M_PI) { dphi = std::abs(dphi - 2*M_PI); }

  double deta = l1_reg_m->hwEta()/etaMult_ - l1_gmt_m->hwEta()/etaMult_;

  double dr = std::sqrt(dphi*dphi + deta*deta);

  if (debug_) {
    std::cout << "****** dr calc debug ******" << std::endl;
    std::cout << "l1_reg_m_phi = " << l1_reg_m_physPhi << std::endl;
    std::cout << "l1_gmt_m_phi = " << l1_gmt_m_physPhi << std::endl;
    std::cout << "l1_reg_m_eta = " << l1_reg_m->hwEta()/etaMult_ << std::endl;
    std::cout << "l1_gmt_m_eta = " << l1_gmt_m->hwEta()/etaMult_ << std::endl;
    std::cout << "dphi         = " << dphi << std::endl;
    std::cout << "deta         = " << deta << std::endl;
    std::cout << "dr           = " << dr   << std::endl;
  }

  return dr;
}

DEFINE_FWK_MODULE(KBmtfMuonAnalysis);
