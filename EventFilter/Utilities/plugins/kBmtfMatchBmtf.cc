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
class kBmtfMatchBmtf : public edm::stream::EDAnalyzer<> {


    public:

        // constructor and destructor
        explicit kBmtfMatchBmtf(const edm::ParameterSet&); 
        ~kBmtfMatchBmtf() override{};

        // method for analyzing the events
        void analyze(const edm::Event&, const edm::EventSetup&) override;


    private:
        
        // structure to store the BMTF muon variables
        struct BmtfMuon {
            int bx, idx, dxy, qual, sign, signValid, processor;
            double pt, ptUnconstrained, phi, eta;
        };


        // methods to initialize the analysis type
        void initMatchMuons(); 

        // conversion methods from hardware to physical coordinates
        int calcGlobalPhi(const l1t::RegionalMuonCand&) const;
        double calcPhysPhi(int) const;
        double calcPhysEta(int) const;
        double calcPtGeV(unsigned int) const;

        // methods to match BMTF muons
        void matchBmtfMuons(const edm::Handle<l1t::RegionalMuonCandBxCollection>&, int);
        double calcDr(const l1t::RegionalMuonCand&, const l1t::RegionalMuonCand&) const; 

        // method to save the two matched muons if they pass the deltaR threshold
        void saveBmtfMuons(const l1t::RegionalMuonCand&, int, int);

        edm::EDGetTokenT<l1t::RegionalMuonCandBxCollection>     bmtfMuonToken_;

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

        // the structure to store the BMTF muon variables
        BmtfMuon bmtfMuonData_;

        // tree for matched muons 
        TTree* bmtfMatchedMuonsTree_;

        // tree for unmatched muons
        TTree* bmtfUnmatchedMuonsTree_;

        // the orbit number
        unsigned int orbitNum_;

};
// -----------------------------------------------------------------------------



// -------------------------------- constructor  -------------------------------

kBmtfMatchBmtf::kBmtfMatchBmtf(const edm::ParameterSet& iConfig):
    bmtfMuonToken_(consumes<l1t::RegionalMuonCandBxCollection>(iConfig.getParameter<edm::InputTag>("bmtfMuonTag"))),                // the token to access the muons from the BMTF
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
void kBmtfMatchBmtf::analyze(const edm::Event& iEvent, const edm::EventSetup& evSetup) {

    // get the muons
    edm::Handle<l1t::RegionalMuonCandBxCollection> bmtfMuonHandle; // a BXVector of RegionalMuonCand
    iEvent.getByToken(bmtfMuonToken_, bmtfMuonHandle);

    // update the orbit number
    orbitNum_ = iEvent.id().event();

    // loop over the BXs
    for (int bx = minBx_; bx < maxBx_; ++bx) {

        // match BMTF muons 
        matchBmtfMuons(bmtfMuonHandle, bx);
    
    }

}
// -----------------------------------------------------------------------------



// ------------------------ methods called for each BX  ------------------------

void kBmtfMatchBmtf::matchBmtfMuons(const edm::Handle<l1t::RegionalMuonCandBxCollection>& bmtfMuonHandle, int bx) {

    int nBmtfMuons = bmtfMuonHandle->size(bx);

    // if there are less than 2 muons in the BX then we cannot have a match
    if (nBmtfMuons < 2) return;

    std::vector<std::tuple<double, int, int>> drValues; // <deltaR, BMTF index 1, BMTF index 2>

    // Get all valid pairings and their deltaR
    for (int iBmtfMuon = 0; iBmtfMuon < nBmtfMuons; ++iBmtfMuon) {
        const l1t::RegionalMuonCand& bmtfMuon1 = bmtfMuonHandle->at(bx, iBmtfMuon);  
        for (int jBmtfMuon = 0; jBmtfMuon < nBmtfMuons; ++jBmtfMuon) {
            if (iBmtfMuon == jBmtfMuon) continue;  // Skip the comparison with itself
            const l1t::RegionalMuonCand& bmtfMuon2 = bmtfMuonHandle->at(bx, jBmtfMuon);  

            int processorDiff = abs(bmtfMuon1.processor() - bmtfMuon2.processor());

            // Ensure the muons are from neighboring boards
            if (processorDiff > 1 && !(bmtfMuon1.processor() == 0 && bmtfMuon2.processor() == 11) && !(bmtfMuon1.processor() == 11 && bmtfMuon2.processor() == 0)) {
                continue;
            }

            double currentDr = calcDr(bmtfMuon1, bmtfMuon2);
            drValues.push_back({currentDr, iBmtfMuon, jBmtfMuon});
        }
    }



    // Sort pairs by deltaR
    std::sort(drValues.begin(), drValues.end(), [](const auto &a, const auto &b) {
        return std::get<0>(a) < std::get<0>(b);
    });

    std::set<int> matchedBmtfIndices;

    for (auto &[dr, iBmtf, jBmtf] : drValues) {
        if (dr >= drThreshold_) continue; // Skip if deltaR is too large
        if (matchedBmtfIndices.find(iBmtf) == matchedBmtfIndices.end() && matchedBmtfIndices.find(jBmtf) == matchedBmtfIndices.end()) {
            matchedBmtfIndices.insert(iBmtf);
            matchedBmtfIndices.insert(jBmtf);

            // update the dr_ variable for the tree
            dr_ = dr;

            const l1t::RegionalMuonCand& matchedBmtfMuon1 = bmtfMuonHandle->at(bx, iBmtf);
            const l1t::RegionalMuonCand& matchedBmtfMuon2 = bmtfMuonHandle->at(bx, jBmtf);
                
            saveBmtfMuons(matchedBmtfMuon1, bx, iBmtf);
            bmtfMatchedMuonsTree_->Fill();

            saveBmtfMuons(matchedBmtfMuon2, bx, jBmtf);
            bmtfMatchedMuonsTree_->Fill();
        }
    }

    // Handle unmatched BMTF muons
    for (int iBmtfMuon = 0; iBmtfMuon < nBmtfMuons; ++iBmtfMuon) {
        if (matchedBmtfIndices.find(iBmtfMuon) == matchedBmtfIndices.end()) {
            const l1t::RegionalMuonCand& unmatchedBmtfMuon = bmtfMuonHandle->at(bx, iBmtfMuon);
            saveBmtfMuons(unmatchedBmtfMuon, bx, iBmtfMuon);
            bmtfUnmatchedMuonsTree_->Fill();
        }
    }

}






// -----------------------------------------------------------------------------


// --------------------------- analysis initializers  --------------------------

// method to initialize the muon matching
void kBmtfMatchBmtf::initMatchMuons() {


    bmtfMatchedMuonsTree_   = fs->make<TTree>("bmtfMatchedMuons", "bmtfMatchedMuons");
    bmtfUnmatchedMuonsTree_ = fs->make<TTree>("bmtfUnmatchedMuons", "bmtfUnmatchedMuons");

    // initialize the tree branches for the BMTF muons (bx, idx, pt, ptUnconstrained, dxy, phi, eta, sign, signValid)
    bmtfMatchedMuonsTree_->Branch("bmtfMuonBx",              &bmtfMuonData_.bx,              "bmtfMuonBx/I");
    bmtfMatchedMuonsTree_->Branch("bmtfMuonIdx",             &bmtfMuonData_.idx,             "bmtfMuonIdx/I");
    bmtfMatchedMuonsTree_->Branch("bmtfMuonPt",              &bmtfMuonData_.pt,              "bmtfMuonPt/D");
    bmtfMatchedMuonsTree_->Branch("bmtfMuonPtUnconstrained", &bmtfMuonData_.ptUnconstrained, "bmtfMuonPtUnconstrained/D");
    bmtfMatchedMuonsTree_->Branch("bmtfMuonDxy",             &bmtfMuonData_.dxy,             "bmtfMuonDxy/I");
    bmtfMatchedMuonsTree_->Branch("bmtfMuonPhi",             &bmtfMuonData_.phi,             "bmtfMuonPhi/D");
    bmtfMatchedMuonsTree_->Branch("bmtfMuonEta",             &bmtfMuonData_.eta,             "bmtfMuonEta/D");
    bmtfMatchedMuonsTree_->Branch("bmtfMuonSign",            &bmtfMuonData_.sign,            "bmtfMuonSign/I");
    bmtfMatchedMuonsTree_->Branch("bmtfMuonSignValid",       &bmtfMuonData_.signValid,       "bmtfMuonSignValid/I");
    bmtfMatchedMuonsTree_->Branch("bmtfMuonDr",              &dr_,                           "bmtfMuonDr/D");
    bmtfMatchedMuonsTree_->Branch("bmtfMuonQual",            &bmtfMuonData_.qual,            "bmtfMuonQual/I");
    bmtfMatchedMuonsTree_->Branch("bmtfMuonProcessor",       &bmtfMuonData_.processor,       "bmtfMuonProcessor/I");
    bmtfMatchedMuonsTree_->Branch("orbitNum",                &orbitNum_,                     "orbitNum/I");


    // initialize the tree branches for the unmatched BMTF muons (bx, idx, pt, ptUnconstrained, dxy, phi, eta, sign, signValid)
    bmtfUnmatchedMuonsTree_->Branch("bmtfMuonBx",              &bmtfMuonData_.bx,              "bmtfMuonBx/I");
    bmtfUnmatchedMuonsTree_->Branch("bmtfMuonIdx",             &bmtfMuonData_.idx,             "bmtfMuonIdx/I");
    bmtfUnmatchedMuonsTree_->Branch("bmtfMuonPt",              &bmtfMuonData_.pt,              "bmtfMuonPt/D");
    bmtfUnmatchedMuonsTree_->Branch("bmtfMuonPtUnconstrained", &bmtfMuonData_.ptUnconstrained, "bmtfMuonPtUnconstrained/D");
    bmtfUnmatchedMuonsTree_->Branch("bmtfMuonDxy",             &bmtfMuonData_.dxy,             "bmtfMuonDxy/I");
    bmtfUnmatchedMuonsTree_->Branch("bmtfMuonPhi",             &bmtfMuonData_.phi,             "bmtfMuonPhi/D");
    bmtfUnmatchedMuonsTree_->Branch("bmtfMuonEta",             &bmtfMuonData_.eta,             "bmtfMuonEta/D");
    bmtfUnmatchedMuonsTree_->Branch("bmtfMuonSign",            &bmtfMuonData_.sign,            "bmtfMuonSign/I");
    bmtfUnmatchedMuonsTree_->Branch("bmtfMuonSignValid",       &bmtfMuonData_.signValid,       "bmtfMuonSignValid/I");
    bmtfUnmatchedMuonsTree_->Branch("bmtfMuonDr",              &dr_,                           "bmtfMuonDr/D");
    bmtfUnmatchedMuonsTree_->Branch("bmtfMuonQual",            &bmtfMuonData_.qual,            "bmtfMuonQual/I");
    bmtfUnmatchedMuonsTree_->Branch("bmtfMuonProcessor",       &bmtfMuonData_.processor,       "bmtfMuonProcessor/I");
    bmtfUnmatchedMuonsTree_->Branch("orbitNum",                &orbitNum_,                     "orbitNum/I");

}


// -----------------------------------------------------------------------------




// ------------------------------ utility methods  -----------------------------

// method that computes the globalPhi coordinate of the muon in the hardware units
int kBmtfMatchBmtf::calcGlobalPhi(const l1t::RegionalMuonCand& muon) const {
    int globalPhi = muon.processor()*48 + muon.hwPhi();
    globalPhi += 576 - 24;
    globalPhi = globalPhi % 576;
    return globalPhi;
}

// method that computes the physical phi coordinate of the muon
double kBmtfMatchBmtf::calcPhysPhi(int globalPhi) const {
    double physPhi = 1.0 * globalPhi / phiMult_;
    if (physPhi > M_PI) {
        physPhi = physPhi - 2*M_PI;
    }
    return physPhi;
}

// method that computes the physical eta coordinate of the muon
double kBmtfMatchBmtf::calcPhysEta(int hwEta) const {
    double physEta = 1.0 * hwEta / etaMult_;
    return physEta;
}

// method that computes the physical pt of the muon
double kBmtfMatchBmtf::calcPtGeV(unsigned int pt) const {
    return pt * 0.5;
}

// method that computes the deltaR between a BMTF muon and a GMT muon
double kBmtfMatchBmtf::calcDr(const l1t::RegionalMuonCand& bmtfMuon1, const l1t::RegionalMuonCand& bmtfMuon2) const {
    double dPhi = fabs(calcPhysPhi(calcGlobalPhi(bmtfMuon1)) - calcPhysPhi(calcGlobalPhi(bmtfMuon2)));
    double dEta = fabs(calcPhysEta(bmtfMuon1.hwEta()) - calcPhysEta(bmtfMuon2.hwEta()));
    return sqrt(dPhi*dPhi + dEta*dEta);
}


// method to save the BMTF muon data in the tree
void kBmtfMatchBmtf::saveBmtfMuons(const l1t::RegionalMuonCand& muon, int bx, int idx) {

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


// -----------------------------------------------------------------------------


DEFINE_FWK_MODULE(kBmtfMatchBmtf);


