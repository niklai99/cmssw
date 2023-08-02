#include "FWCore/Framework/interface/MakerMacros.h"

// system include files
#include <fstream>
#include <iomanip>
#include <memory>
#include <string>
#include <cmath>

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
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"


// root include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
// #include "TDirectory.h"



// ----------------------------- CLASS DECLARATION  ----------------------------
class kBmtfMatchGmt : public edm::stream::EDAnalyzer<> {


    public:

        // constructor and destructor
        explicit kBmtfMatchGmt(const edm::ParameterSet&); 
        ~kBmtfMatchGmt() override{};

        // method for analyzing the events
        void analyze(const edm::Event&, const edm::EventSetup&) override;


    private:
        
        // structure to store the BMTF muon variables
        struct BmtfMuon {
            int bx, idx, dxy, qual, sign, signValid, processor;
            double pt, ptUnconstrained, phi, eta;
        };

        // structure to store the GMT muon variables
        struct GmtMuon {
            int bx, idx, tfidx, dxy, charge, chargeValid, quality;
            double pt, ptUnconstrained, phi, eta, phiAtVtx, etaAtVtx;
        };


        // methods to initialize the analysis type
        void initMatchMuons(); 


        // conversion methods from hardware to physical coordinates
        int calcGlobalPhi(const l1t::RegionalMuonCand&) const;
        double calcPhysPhi(int) const;
        double calcPhysEta(int) const;
        double calcPtGeV(unsigned int) const;

        // methods to match BMTF muons to GMT muons
        void matchBmtfOnGmt(const edm::Handle<l1t::RegionalMuonCandBxCollection>&, const edm::Handle<scoutingRun3::MuonOrbitCollection>&, int);
        void matchGmtOnBmtf(const edm::Handle<l1t::RegionalMuonCandBxCollection>&, const edm::Handle<scoutingRun3::MuonOrbitCollection>&, int);
        double calcDr(const l1t::RegionalMuonCand&, const l1t::Muon*) const; 
        bool isBarrel(const l1t::Muon*) const;

        // method to save the two matched muons if they pass the deltaR threshold
        void saveBmtfMuons(const l1t::RegionalMuonCand&, int, int);
        void saveGmtMuons(const l1t::Muon*, int, int);



        // the tokens to access the data
        edm::EDGetTokenT<l1t::RegionalMuonCandBxCollection>     bmtfMuonToken_;
        edm::EDGetTokenT<scoutingRun3::MuonOrbitCollection>     gmtMuonToken_;

        // the min and max BX to be analyzed
        int minBx_;
        int maxBx_;

        // phi and eta conversion factors
        double phiMult_;
        double etaMult_;

        // deltaR threshold for matching BMTF muons to GMT muons
        double drThreshold_;

        // deltaR between the two muons
        double dr_;

        // debug and verbose flags
        bool debug_;
        bool vverbose_;

        // the root file service to handle the output file
        edm::Service<TFileService> fs;

        // the muon data structures
        BmtfMuon bmtfMuonData_;
        GmtMuon gmtMuonData_;

        // the tree for the matched muons
        TTree* bmtfOnGmtBmtfMuonsTree_;
        TTree* bmtfOnGmtGmtMuonsTree_;
        TTree* gmtOnBmtfBmtfMuonsTree_;
        TTree* gmtOnBmtfGmtMuonsTree_;

        // the tree for unmatched muons
        TTree* bmtfOnGmtUnmatchedBmtfMuonsTree_;
        TTree* gmtOnBmtfUnmatchedGmtMuonsTree_;

        // the orbit number
        unsigned int orbitNum_;

};
// -----------------------------------------------------------------------------



// -------------------------------- constructor  -------------------------------

kBmtfMatchGmt::kBmtfMatchGmt(const edm::ParameterSet& iConfig):
    bmtfMuonToken_(consumes<l1t::RegionalMuonCandBxCollection>(iConfig.getParameter<edm::InputTag>("bmtfMuonTag"))),                // the token to access the muons from the BMTF
    gmtMuonToken_(consumes<scoutingRun3::MuonOrbitCollection>(iConfig.getParameter<edm::InputTag>("gmtMuonTag"))),                  // the token to access the muons from the GMT
    minBx_(iConfig.getParameter<int>("minBx")),                                                                                     // the min BX to be analyzed [0]
    maxBx_(iConfig.getParameter<int>("maxBx")),                                                                                     // the max BX to be analyzed [3564]
    phiMult_(iConfig.getParameter<double>("phiMult")),                                                                              // the phi conversion factor from hardware to physical units 
    etaMult_(iConfig.getParameter<double>("etaMult")),                                                                              // the eta conversion factor from hardware to physical units
    drThreshold_(iConfig.getParameter<double>("drThreshold")),                                                                      // the deltaR threshold for matching BMTF muons to GMT muons
    debug_(iConfig.getParameter<bool>("debug")),                                                                                    // the debug flag (currently not used)
    vverbose_(iConfig.getParameter<bool>("vverbose"))                                                                               // the very verbose flag (currently not used)
{   

    // initialize the analysis type
    initMatchMuons();

}

// -----------------------------------------------------------------------------





// ----------------------- method called for each orbit  -----------------------
void kBmtfMatchGmt::analyze(const edm::Event& iEvent, const edm::EventSetup& evSetup) {

    // get the muons
    edm::Handle<l1t::RegionalMuonCandBxCollection> bmtfMuonHandle; // a BXVector of RegionalMuonCand
    iEvent.getByToken(bmtfMuonToken_, bmtfMuonHandle);

    // get the GMT muons
    edm::Handle<scoutingRun3::MuonOrbitCollection> gmtMuonHandle; 
    iEvent.getByToken(gmtMuonToken_, gmtMuonHandle);

    // update the orbit number
    orbitNum_++;

    // loop over the BXs
    for (int bx = minBx_; bx < maxBx_; ++bx) {

        // match the kBMTF muons to the GMT muons
        matchBmtfOnGmt(bmtfMuonHandle, gmtMuonHandle, bx); 

        // match the GMT muons to the kBMTF muons
        matchGmtOnBmtf(bmtfMuonHandle, gmtMuonHandle, bx);
    
    }

}
// -----------------------------------------------------------------------------



// ------------------------ methods called for each BX  ------------------------

void kBmtfMatchGmt::matchBmtfOnGmt(const edm::Handle<l1t::RegionalMuonCandBxCollection>& bmtfMuonHandle, const edm::Handle<scoutingRun3::MuonOrbitCollection>& gmtMuonHandle, int bx) {

    int nBmtfMuons = bmtfMuonHandle->size(bx);
    if (nBmtfMuons == 0) return;

    const std::vector<int>* gmtMuonIndex = gmtMuonHandle->getIndex();
    int current_bx_index = gmtMuonIndex->at(bx);
    int next_bx_index    = gmtMuonIndex->at(bx+1);
    int nGmtMuons        = next_bx_index - current_bx_index;
    if (nGmtMuons == 0) return;

    std::vector<std::tuple<double, int, int>> drValues; // <deltaR, BMTF index, GMT index>

    // Get all valid pairings and their deltaR
    for (int iBmtfMuon = 0; iBmtfMuon < nBmtfMuons; ++iBmtfMuon) {
        const l1t::RegionalMuonCand& bmtfMuon = bmtfMuonHandle->at(bx, iBmtfMuon);  
        for (int iGmtMuon = current_bx_index; iGmtMuon < next_bx_index; ++iGmtMuon) {
            const l1t::Muon* gmtMuon = gmtMuonHandle->getFlatData(iGmtMuon);
            if (!isBarrel(gmtMuon)) continue;
            double currentDr = calcDr(bmtfMuon, gmtMuon);
            drValues.push_back({currentDr, iBmtfMuon, iGmtMuon});
        }
    }

    // Sort pairs by deltaR
    std::sort(drValues.begin(), drValues.end(), [](const auto &a, const auto &b) {
        return std::get<0>(a) < std::get<0>(b);
    });

    std::set<int> matchedBmtfIndices;
    std::set<int> matchedGmtIndices;

    for (auto &[dr, iBmtf, iGmt] : drValues) {
        if (dr >= drThreshold_) continue; // Skip if deltaR is too large
        if (matchedBmtfIndices.find(iBmtf) == matchedBmtfIndices.end() && matchedGmtIndices.find(iGmt) == matchedGmtIndices.end()) {
            matchedBmtfIndices.insert(iBmtf);
            matchedGmtIndices.insert(iGmt);

            // update the dr_ variable for the tree
            dr_ = dr;

            const l1t::RegionalMuonCand& matchedBmtfMuon = bmtfMuonHandle->at(bx, iBmtf);
            const l1t::Muon* matchedGmtMuon = gmtMuonHandle->getFlatData(iGmt);
            
            saveBmtfMuons(matchedBmtfMuon, bx, iBmtf);
            saveGmtMuons(matchedGmtMuon, bx, iGmt - current_bx_index);
            bmtfOnGmtBmtfMuonsTree_->Fill();
            bmtfOnGmtGmtMuonsTree_->Fill();
        }
    }

    // Handle unmatched BMTF muons
    for (int iBmtfMuon = 0; iBmtfMuon < nBmtfMuons; ++iBmtfMuon) {
        if (matchedBmtfIndices.find(iBmtfMuon) == matchedBmtfIndices.end()) {
            const l1t::RegionalMuonCand& unmatchedBmtfMuon = bmtfMuonHandle->at(bx, iBmtfMuon);
            saveBmtfMuons(unmatchedBmtfMuon, bx, iBmtfMuon);
            bmtfOnGmtUnmatchedBmtfMuonsTree_->Fill();
        }
    }
}




void kBmtfMatchGmt::matchGmtOnBmtf(const edm::Handle<l1t::RegionalMuonCandBxCollection>& bmtfMuonHandle, const edm::Handle<scoutingRun3::MuonOrbitCollection>& gmtMuonHandle, int bx) {

    int nBmtfMuons = bmtfMuonHandle->size(bx);
    if (nBmtfMuons == 0) return;

    const std::vector<int>* gmtMuonIndex = gmtMuonHandle->getIndex();
    int current_bx_index = gmtMuonIndex->at(bx);
    int next_bx_index    = gmtMuonIndex->at(bx+1);
    int nGmtMuons        = next_bx_index - current_bx_index;
    if (nGmtMuons == 0) return;

    std::vector<std::tuple<double, int, int>> drValues; // <deltaR, GMT index, BMTF index>

    // Get all valid pairings and their deltaR
    for (int iGmtMuon = current_bx_index; iGmtMuon < next_bx_index; ++iGmtMuon) {
        const l1t::Muon* gmtMuon = gmtMuonHandle->getFlatData(iGmtMuon);
        if (!isBarrel(gmtMuon)) continue;
        for (int iBmtfMuon = 0; iBmtfMuon < nBmtfMuons; ++iBmtfMuon) {
            const l1t::RegionalMuonCand& bmtfMuon = bmtfMuonHandle->at(bx, iBmtfMuon);  
            double currentDr = calcDr(bmtfMuon, gmtMuon);
            drValues.push_back({currentDr, iGmtMuon, iBmtfMuon});
        }
    }

    // Sort pairs by deltaR
    std::sort(drValues.begin(), drValues.end(), [](const auto &a, const auto &b) {
        return std::get<0>(a) < std::get<0>(b);
    });

    std::set<int> matchedBmtfIndices;
    std::set<int> matchedGmtIndices;

    for (auto &[dr, iGmt, iBmtf] : drValues) {
        if (dr >= drThreshold_) continue; // Skip if deltaR is too large
        if (matchedGmtIndices.find(iGmt) == matchedGmtIndices.end() && matchedBmtfIndices.find(iBmtf) == matchedBmtfIndices.end()) {
            matchedBmtfIndices.insert(iBmtf);
            matchedGmtIndices.insert(iGmt);

            // update the dr_ variable for the tree
            dr_ = dr;

            const l1t::RegionalMuonCand& matchedBmtfMuon = bmtfMuonHandle->at(bx, iBmtf);
            const l1t::Muon* matchedGmtMuon = gmtMuonHandle->getFlatData(iGmt);
            
            saveGmtMuons(matchedGmtMuon, bx, iGmt - current_bx_index);
            saveBmtfMuons(matchedBmtfMuon, bx, iBmtf);
            gmtOnBmtfGmtMuonsTree_->Fill();
            gmtOnBmtfBmtfMuonsTree_->Fill();
        }
    }

    // Handle unmatched GMT muons
    for (int iGmtMuon = current_bx_index; iGmtMuon < next_bx_index; ++iGmtMuon) {
        if (matchedGmtIndices.find(iGmtMuon) == matchedGmtIndices.end()) {
            const l1t::Muon* unmatchedGmtMuon = gmtMuonHandle->getFlatData(iGmtMuon);
            if (!isBarrel(unmatchedGmtMuon)) continue;
            saveGmtMuons(unmatchedGmtMuon, bx, iGmtMuon - current_bx_index);
            gmtOnBmtfUnmatchedGmtMuonsTree_->Fill();
        }
    }
}



// -----------------------------------------------------------------------------


// --------------------------- analysis initializers  --------------------------

// method to initialize the muon matching
void kBmtfMatchGmt::initMatchMuons() {

    // initialize the orbit number
    orbitNum_ = 0;


    // BMTF ON GMT

    bmtfOnGmtBmtfMuonsTree_ = fs->make<TTree>("bmtfOnGmtBmtfMuons", "bmtfOnGmtBmtfMuons");
    bmtfOnGmtGmtMuonsTree_  = fs->make<TTree>("bmtfOnGmtGmtMuons",  "bmtfOnGmtGmtMuons");
    bmtfOnGmtUnmatchedBmtfMuonsTree_ = fs->make<TTree>("bmtfOnGmtUnmatchedBmtfMuons", "bmtfOnGmtUnmatchedBmtfMuons");

    // initialize the tree branches for the BMTF muons (bx, idx, pt, ptUnconstrained, dxy, phi, eta, sign, signValid)
    bmtfOnGmtBmtfMuonsTree_->Branch("bmtfMuonBx",              &bmtfMuonData_.bx,              "bmtfMuonBx/I");
    bmtfOnGmtBmtfMuonsTree_->Branch("bmtfMuonIdx",             &bmtfMuonData_.idx,             "bmtfMuonIdx/I");
    bmtfOnGmtBmtfMuonsTree_->Branch("bmtfMuonPt",              &bmtfMuonData_.pt,              "bmtfMuonPt/D");
    bmtfOnGmtBmtfMuonsTree_->Branch("bmtfMuonPtUnconstrained", &bmtfMuonData_.ptUnconstrained, "bmtfMuonPtUnconstrained/D");
    bmtfOnGmtBmtfMuonsTree_->Branch("bmtfMuonDxy",             &bmtfMuonData_.dxy,             "bmtfMuonDxy/I");
    bmtfOnGmtBmtfMuonsTree_->Branch("bmtfMuonPhi",             &bmtfMuonData_.phi,             "bmtfMuonPhi/D");
    bmtfOnGmtBmtfMuonsTree_->Branch("bmtfMuonEta",             &bmtfMuonData_.eta,             "bmtfMuonEta/D");
    bmtfOnGmtBmtfMuonsTree_->Branch("bmtfMuonSign",            &bmtfMuonData_.sign,            "bmtfMuonSign/I");
    bmtfOnGmtBmtfMuonsTree_->Branch("bmtfMuonSignValid",       &bmtfMuonData_.signValid,       "bmtfMuonSignValid/I");
    bmtfOnGmtBmtfMuonsTree_->Branch("bmtfMuonDr",              &dr_,                           "bmtfMuonDr/D");
    bmtfOnGmtBmtfMuonsTree_->Branch("bmtfMuonQual",            &bmtfMuonData_.qual,            "bmtfMuonQual/I");
    bmtfOnGmtBmtfMuonsTree_->Branch("bmtfMuonProcessor",       &bmtfMuonData_.processor,       "bmtfMuonProcessor/I");


    bmtfOnGmtBmtfMuonsTree_->Branch("orbitNum", &orbitNum_, "orbitNum/I");

    // initialize the tree branches for the GMT muon (bx, idx, pt, ptUnconstrained, dxy, phi, eta, charge, chargeValid)
    bmtfOnGmtGmtMuonsTree_->Branch("gmtMuonBx",              &gmtMuonData_.bx,              "gmtMuonBx/I");
    bmtfOnGmtGmtMuonsTree_->Branch("gmtMuonIdx",             &gmtMuonData_.idx,             "gmtMuonIdx/I");
    bmtfOnGmtGmtMuonsTree_->Branch("gmtMuonPt",              &gmtMuonData_.pt,              "gmtMuonPt/D");
    bmtfOnGmtGmtMuonsTree_->Branch("gmtMuonPtUnconstrained", &gmtMuonData_.ptUnconstrained, "gmtMuonPtUnconstrained/D");
    bmtfOnGmtGmtMuonsTree_->Branch("gmtMuonDxy",             &gmtMuonData_.dxy,             "gmtMuonDxy/I");
    bmtfOnGmtGmtMuonsTree_->Branch("gmtMuonPhi",             &gmtMuonData_.phi,             "gmtMuonPhi/D");
    bmtfOnGmtGmtMuonsTree_->Branch("gmtMuonEta",             &gmtMuonData_.eta,             "gmtMuonEta/D");
    bmtfOnGmtGmtMuonsTree_->Branch("gmtMuonCharge",          &gmtMuonData_.charge,          "gmtMuonCharge/I");
    bmtfOnGmtGmtMuonsTree_->Branch("gmtMuonChargeValid",     &gmtMuonData_.chargeValid,     "gmtMuonChargeValid/I");
    bmtfOnGmtGmtMuonsTree_->Branch("gmtMuonDr",              &dr_,                          "gmtMuonDr/D");
    bmtfOnGmtGmtMuonsTree_->Branch("gmtMuonQuality",         &gmtMuonData_.quality,         "gmtMuonQuality/I");

    bmtfOnGmtGmtMuonsTree_->Branch("orbitNum", &orbitNum_, "orbitNum/I");

    // initialize the tree branches for the unmatched BMTF muons (bx, idx, pt, ptUnconstrained, dxy, phi, eta, sign, signValid)
    bmtfOnGmtUnmatchedBmtfMuonsTree_->Branch("bmtfMuonBx",              &bmtfMuonData_.bx,              "bmtfMuonBx/I");
    bmtfOnGmtUnmatchedBmtfMuonsTree_->Branch("bmtfMuonIdx",             &bmtfMuonData_.idx,             "bmtfMuonIdx/I");
    bmtfOnGmtUnmatchedBmtfMuonsTree_->Branch("bmtfMuonPt",              &bmtfMuonData_.pt,              "bmtfMuonPt/D");
    bmtfOnGmtUnmatchedBmtfMuonsTree_->Branch("bmtfMuonPtUnconstrained", &bmtfMuonData_.ptUnconstrained, "bmtfMuonPtUnconstrained/D");
    bmtfOnGmtUnmatchedBmtfMuonsTree_->Branch("bmtfMuonDxy",             &bmtfMuonData_.dxy,             "bmtfMuonDxy/I");
    bmtfOnGmtUnmatchedBmtfMuonsTree_->Branch("bmtfMuonPhi",             &bmtfMuonData_.phi,             "bmtfMuonPhi/D");
    bmtfOnGmtUnmatchedBmtfMuonsTree_->Branch("bmtfMuonEta",             &bmtfMuonData_.eta,             "bmtfMuonEta/D");
    bmtfOnGmtUnmatchedBmtfMuonsTree_->Branch("bmtfMuonSign",            &bmtfMuonData_.sign,            "bmtfMuonSign/I");
    bmtfOnGmtUnmatchedBmtfMuonsTree_->Branch("bmtfMuonSignValid",       &bmtfMuonData_.signValid,       "bmtfMuonSignValid/I");
    bmtfOnGmtUnmatchedBmtfMuonsTree_->Branch("bmtfMuonDr",              &dr_,                           "bmtfMuonDr/D");
    bmtfOnGmtUnmatchedBmtfMuonsTree_->Branch("bmtfMuonQual",            &bmtfMuonData_.qual,            "bmtfMuonQual/I");
    bmtfOnGmtUnmatchedBmtfMuonsTree_->Branch("bmtfMuonProcessor",       &bmtfMuonData_.processor,       "bmtfMuonProcessor/I");

    bmtfOnGmtUnmatchedBmtfMuonsTree_->Branch("orbitNum", &orbitNum_, "orbitNum/I");


    // GMT ON BMTF

    gmtOnBmtfBmtfMuonsTree_ = fs->make<TTree>("gmtOnBmtfBmtfMuons", "gmtOnBmtfBmtfMuons");
    gmtOnBmtfGmtMuonsTree_  = fs->make<TTree>("gmtOnBmtfGmtMuons",  "gmtOnBmtfGmtMuons");
    gmtOnBmtfUnmatchedGmtMuonsTree_ = fs->make<TTree>("gmtOnBmtfUnmatchedGmtMuons", "gmtOnBmtfUnmatchedGmtMuons");

    // initialize the tree branches for the BMTF muons (bx, idx, pt, ptUnconstrained, dxy, phi, eta, sign, signValid)
    gmtOnBmtfBmtfMuonsTree_->Branch("bmtfMuonBx",              &bmtfMuonData_.bx,              "bmtfMuonBx/I");
    gmtOnBmtfBmtfMuonsTree_->Branch("bmtfMuonIdx",             &bmtfMuonData_.idx,             "bmtfMuonIdx/I");
    gmtOnBmtfBmtfMuonsTree_->Branch("bmtfMuonPt",              &bmtfMuonData_.pt,              "bmtfMuonPt/D");
    gmtOnBmtfBmtfMuonsTree_->Branch("bmtfMuonPtUnconstrained", &bmtfMuonData_.ptUnconstrained, "bmtfMuonPtUnconstrained/D");
    gmtOnBmtfBmtfMuonsTree_->Branch("bmtfMuonDxy",             &bmtfMuonData_.dxy,             "bmtfMuonDxy/I");
    gmtOnBmtfBmtfMuonsTree_->Branch("bmtfMuonPhi",             &bmtfMuonData_.phi,             "bmtfMuonPhi/D");
    gmtOnBmtfBmtfMuonsTree_->Branch("bmtfMuonEta",             &bmtfMuonData_.eta,             "bmtfMuonEta/D");
    gmtOnBmtfBmtfMuonsTree_->Branch("bmtfMuonSign",            &bmtfMuonData_.sign,            "bmtfMuonSign/I");
    gmtOnBmtfBmtfMuonsTree_->Branch("bmtfMuonSignValid",       &bmtfMuonData_.signValid,       "bmtfMuonSignValid/I");
    gmtOnBmtfBmtfMuonsTree_->Branch("bmtfMuonDr",              &dr_,                           "bmtfMuonDr/D");
    gmtOnBmtfBmtfMuonsTree_->Branch("bmtfMuonQual",            &bmtfMuonData_.qual,            "bmtfMuonQual/I");
    gmtOnBmtfBmtfMuonsTree_->Branch("bmtfMuonProcessor",       &bmtfMuonData_.processor,       "bmtfMuonProcessor/I");

    gmtOnBmtfBmtfMuonsTree_->Branch("orbitNum", &orbitNum_, "orbitNum/I");


    // initialize the tree branches for the GMT muon (bx, idx, pt, ptUnconstrained, dxy, phi, eta, charge, chargeValid)
    gmtOnBmtfGmtMuonsTree_->Branch("gmtMuonBx",              &gmtMuonData_.bx,              "gmtMuonBx/I");
    gmtOnBmtfGmtMuonsTree_->Branch("gmtMuonIdx",             &gmtMuonData_.idx,             "gmtMuonIdx/I");
    gmtOnBmtfGmtMuonsTree_->Branch("gmtMuonPt",              &gmtMuonData_.pt,              "gmtMuonPt/D");
    gmtOnBmtfGmtMuonsTree_->Branch("gmtMuonPtUnconstrained", &gmtMuonData_.ptUnconstrained, "gmtMuonPtUnconstrained/D");
    gmtOnBmtfGmtMuonsTree_->Branch("gmtMuonDxy",             &gmtMuonData_.dxy,             "gmtMuonDxy/I");
    gmtOnBmtfGmtMuonsTree_->Branch("gmtMuonPhi",             &gmtMuonData_.phi,             "gmtMuonPhi/D");
    gmtOnBmtfGmtMuonsTree_->Branch("gmtMuonEta",             &gmtMuonData_.eta,             "gmtMuonEta/D");
    gmtOnBmtfGmtMuonsTree_->Branch("gmtMuonCharge",          &gmtMuonData_.charge,          "gmtMuonCharge/I");
    gmtOnBmtfGmtMuonsTree_->Branch("gmtMuonChargeValid",     &gmtMuonData_.chargeValid,     "gmtMuonChargeValid/I");
    gmtOnBmtfGmtMuonsTree_->Branch("gmtMuonDr",              &dr_,                          "gmtMuonDr/D");
    gmtOnBmtfGmtMuonsTree_->Branch("gmtMuonQuality",         &gmtMuonData_.quality,         "gmtMuonQuality/I");

    gmtOnBmtfGmtMuonsTree_->Branch("orbitNum", &orbitNum_, "orbitNum/I");

    // initialize the tree branches for the unmatched BMTF muons (bx, idx, pt, ptUnconstrained, dxy, phi, eta, sign, signValid)
    gmtOnBmtfUnmatchedGmtMuonsTree_->Branch("gmtMuonBx",              &gmtMuonData_.bx,              "gmtMuonBx/I");
    gmtOnBmtfUnmatchedGmtMuonsTree_->Branch("gmtMuonIdx",             &gmtMuonData_.idx,             "gmtMuonIdx/I");
    gmtOnBmtfUnmatchedGmtMuonsTree_->Branch("gmtMuonPt",              &gmtMuonData_.pt,              "gmtMuonPt/D");
    gmtOnBmtfUnmatchedGmtMuonsTree_->Branch("gmtMuonPtUnconstrained", &gmtMuonData_.ptUnconstrained, "gmtMuonPtUnconstrained/D");
    gmtOnBmtfUnmatchedGmtMuonsTree_->Branch("gmtMuonDxy",             &gmtMuonData_.dxy,             "gmtMuonDxy/I");
    gmtOnBmtfUnmatchedGmtMuonsTree_->Branch("gmtMuonPhi",             &gmtMuonData_.phi,             "gmtMuonPhi/D");
    gmtOnBmtfUnmatchedGmtMuonsTree_->Branch("gmtMuonEta",             &gmtMuonData_.eta,             "gmtMuonEta/D");
    gmtOnBmtfUnmatchedGmtMuonsTree_->Branch("gmtMuonCharge",          &gmtMuonData_.charge,          "gmtMuonCharge/I");
    gmtOnBmtfUnmatchedGmtMuonsTree_->Branch("gmtMuonChargeValid",     &gmtMuonData_.chargeValid,     "gmtMuonChargeValid/I");
    gmtOnBmtfUnmatchedGmtMuonsTree_->Branch("gmtMuonDr",              &dr_,                          "gmtMuonDr/D");
    gmtOnBmtfUnmatchedGmtMuonsTree_->Branch("gmtMuonQuality",         &gmtMuonData_.quality,         "gmtMuonQuality/I");

    gmtOnBmtfUnmatchedGmtMuonsTree_->Branch("orbitNum", &orbitNum_, "orbitNum/I");


}


// -----------------------------------------------------------------------------




// ------------------------------ utility methods  -----------------------------

// method that computes the globalPhi coordinate of the muon in the hardware units
int kBmtfMatchGmt::calcGlobalPhi(const l1t::RegionalMuonCand& muon) const {
    int globalPhi = muon.processor()*48 + muon.hwPhi();
    globalPhi += 576 - 24;
    globalPhi = globalPhi % 576;
    return globalPhi;
}

// method that computes the physical phi coordinate of the muon
double kBmtfMatchGmt::calcPhysPhi(int globalPhi) const {
    double physPhi = 1.0 * globalPhi / phiMult_;
    if (physPhi > M_PI) {
        physPhi = physPhi - 2*M_PI;
    }
    return physPhi;
}

// method that computes the physical eta coordinate of the muon
double kBmtfMatchGmt::calcPhysEta(int hwEta) const {
    double physEta = 1.0 * hwEta / etaMult_;
    return physEta;
}

// method that computes the physical pt of the muon
double kBmtfMatchGmt::calcPtGeV(unsigned int pt) const {
    return pt * 0.5;
}

// method that computes the deltaR between a BMTF muon and a GMT muon
double kBmtfMatchGmt::calcDr(const l1t::RegionalMuonCand& bmtfMuon, const l1t::Muon* gmtMuon) const {
    double dPhi = fabs(calcPhysPhi(calcGlobalPhi(bmtfMuon)) - calcPhysPhi(gmtMuon->hwPhi()));
    double dEta = fabs(calcPhysEta(bmtfMuon.hwEta()) - calcPhysEta(gmtMuon->hwEta()));
    return sqrt(dPhi*dPhi + dEta*dEta);
}

// method to check if the GMT muon is in the barrel to match it to a BMTF muon
bool kBmtfMatchGmt::isBarrel(const l1t::Muon* gmtMuon) const {
    // the GMT muon is in the barrel it has tfIndex between 36 and 70 (included)
    return (gmtMuon->tfMuonIndex() >= 36 && gmtMuon->tfMuonIndex() <= 70);
}

// method to save the BMTF muon data in the tree
void kBmtfMatchGmt::saveBmtfMuons(const l1t::RegionalMuonCand& muon, int bx, int idx) {

    // fill the tree branches
    bmtfMuonData_.bx              = bx;
    bmtfMuonData_.pt              = calcPtGeV(muon.hwPt());
    bmtfMuonData_.ptUnconstrained = calcPtGeV(muon.hwPtUnconstrained());
    bmtfMuonData_.dxy             = muon.hwDXY();
    bmtfMuonData_.phi             = calcPhysPhi(calcGlobalPhi(muon));
    bmtfMuonData_.eta             = calcPhysEta(muon.hwEta());
    bmtfMuonData_.sign            = muon.hwSign();
    bmtfMuonData_.signValid       = muon.hwSignValid();
    bmtfMuonData_.idx             = idx;
    bmtfMuonData_.qual            = muon.hwQual();
    bmtfMuonData_.processor       = muon.processor();

}

// method to save the GMT muon data in the tree
void kBmtfMatchGmt::saveGmtMuons(const l1t::Muon* gmtMuon, int bx, int idx) {

    // fill the tree branches
    gmtMuonData_.bx              = bx;
    gmtMuonData_.pt              = calcPtGeV(gmtMuon->hwPt());
    gmtMuonData_.ptUnconstrained = calcPtGeV(gmtMuon->hwPtUnconstrained());
    gmtMuonData_.dxy             = gmtMuon->hwDXY();
    gmtMuonData_.phi             = calcPhysPhi(gmtMuon->hwPhi());
    gmtMuonData_.eta             = calcPhysEta(gmtMuon->hwEta());
    gmtMuonData_.charge          = gmtMuon->hwCharge();
    gmtMuonData_.chargeValid     = gmtMuon->hwChargeValid();
    gmtMuonData_.idx             = idx;
    gmtMuonData_.quality         = gmtMuon->hwQual();

}



// -----------------------------------------------------------------------------


DEFINE_FWK_MODULE(kBmtfMatchGmt);


