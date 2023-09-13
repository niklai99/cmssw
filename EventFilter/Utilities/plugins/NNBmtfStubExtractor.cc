#include "FWCore/Framework/interface/MakerMacros.h"

// system include files
#include <fstream>
#include <iomanip>
#include <memory>
#include <string>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/L1Trigger/interface/Muon.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/MessageLogger/interface/MessageDrop.h"

#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/L1Scouting/interface/SRawDataCollection.h"
#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/L1Scouting/interface/SDSNumbering.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"

#include "DataFormats/L1TMuon/interface/L1MuKBMTrack.h"
#include "DataFormats/L1TMuon/interface/L1MuKBMTCombinedStub.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"

#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuonSetup.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"

#include "L1Trigger/L1TMuonBarrel/interface/L1TMuonBarrelKalmanStubProcessor.h"
#include "CondFormats/L1TObjects/interface/L1MuDTTFParameters.h"
#include "CondFormats/DataRecord/interface/L1MuDTTFParametersRcd.h"
#include "CondFormats/L1TObjects/interface/L1MuDTTFMasks.h"
#include "CondFormats/DataRecord/interface/L1MuDTTFMasksRcd.h"

#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhDigi.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThDigi.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"
#include "DataFormats/L1TMuon/interface/L1MuKBMTCombinedStub.h"
#include "CondFormats/L1TObjects/interface/L1TMuonBarrelParams.h"
#include "CondFormats/DataRecord/interface/L1TMuonBarrelParamsRcd.h"



// Using statements
using reco::MuonCollection;
using l1t::MuonBxCollection;
using l1t::RegionalMuonCandBxCollection;



// Class declaration
class NNBmtfStubExtractor : public edm::stream::EDAnalyzer<> {
public:

    // Constructor
    explicit NNBmtfStubExtractor(const edm::ParameterSet&);
    // Destructor
    ~NNBmtfStubExtractor() override = default;
    // Analyzer
    void analyze(const edm::Event&, const edm::EventSetup&) override;

private:

    // Helper methods to compute simple quantities
    int calcGlobalPhi(const l1t::RegionalMuonCand&) const;
    double calcPhysPhi(int) const;
    double calcPhysEta(int) const;
    double calcPtGeV(unsigned int) const;
    double calcDr(const l1t::RegionalMuonCand&, const l1t::Muon*) const;
    double calcDr(const l1t::Muon*, const TrajectoryStateOnSurface*) const;
    bool isBarrel(const l1t::Muon*) const;

    // Helper methods to initialize the analyzer
    void initCsvHeader();
    void initializeMuPropagators(const edm::EventSetup&);
    FreeTrajectoryState createFTSForRecoMuon(const reco::Muon&, edm::ESHandle<MagneticField>);
    std::tuple<TrajectoryStateOnSurface, TrajectoryStateOnSurface> propagateToStations(const FreeTrajectoryState&);

    // Helper methods to check validity of the muon
    bool recoMuonCheck(const reco::Muon&) const;

    // Matching methods
    std::tuple<double, int, int, int> findClosestL1Match(const edm::Event&, const TrajectoryStateOnSurface&, const BXVector<l1t::Muon>&);
    bool doL1Matching(const l1t::RegionalMuonCand&, const L1MuKBMTrack&, const l1t::Muon&);

    // Helper methods to write data to file
    void writeEventInfoToFile(const edm::Event&);
    void writeRecoToFile(const reco::Muon&, const TrajectoryStateOnSurface&);
    void writeL1ToFile(const l1t::RegionalMuonCand&);
    void writeKBmtfToFile(const L1MuKBMTrack&);
    void writeGmtToFile(const l1t::Muon&);
    void writeStubToFile(const L1MuKBMTCombinedStub&);
    void writeStubsToFile(const L1MuKBMTCombinedStubRefVector&);


    // CMSSW tokens
    edm::EDGetTokenT<MuonCollection>                                recoMuonsToken_;
    edm::EDGetTokenT<MuonBxCollection>                              simGmtStage2DigisToken_;
    edm::EDGetTokenT<RegionalMuonCandBxCollection>                  simKBmtfDigisToken_;
    edm::EDGetTokenT<L1MuKBMTrackBxCollection>                      simKBmtfTracksToken_;
    const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord>  magfieldToken_;

    // Processors
    L1TMuonBarrelKalmanStubProcessor stubProcessor_;

    // Propagators
    PropagateToMuonSetup const muPropagator1stSetup_;
    PropagateToMuonSetup const muPropagator2ndSetup_;
    PropagateToMuon muPropagator1st_;
    PropagateToMuon muPropagator2nd_;
    
    // Input parameters
    double drCut_;  
    double phiMult_;
    double etaMult_;
    std::vector<int> cotTheta_1_;
    std::vector<int> cotTheta_2_;
    std::vector<int> cotTheta_3_;
    int minBx_;
    int maxBx_;
    bool useDummyValues_;

    // File names
    std::string filenameNNBmtfStubs_;

    // Output file
    std::ofstream outputFile;

    // Constants
    static const double PT_CONVERSION_FACTOR;

};


// Initialize constants
const double NNBmtfStubExtractor::PT_CONVERSION_FACTOR = 0.5;



// Constructor
NNBmtfStubExtractor::NNBmtfStubExtractor(const edm::ParameterSet& iConfig):
    recoMuonsToken_(consumes<MuonCollection>(iConfig.getParameter<edm::InputTag>("recoMuons"))),
    simGmtStage2DigisToken_(consumes<MuonBxCollection>(iConfig.getParameter<edm::InputTag>("l1tMuonBXVector"))),
    simKBmtfDigisToken_(consumes<RegionalMuonCandBxCollection>(iConfig.getParameter<edm::InputTag>("l1tRegionalMuonCandBXVector"))),
    simKBmtfTracksToken_(consumes<L1MuKBMTrackBxCollection>(iConfig.getParameter<edm::InputTag>("L1MuKBMTrackBXVector"))),  
    magfieldToken_(consumesCollector().esConsumes()),
    muPropagator1stSetup_(iConfig.getParameter<edm::ParameterSet>("muProp1st"), consumesCollector()),
    muPropagator2ndSetup_(iConfig.getParameter<edm::ParameterSet>("muProp2nd"), consumesCollector()),
    drCut_(iConfig.getParameter<double>("drCut")),
    phiMult_(iConfig.getParameter<double>("phiMult")),
    etaMult_(iConfig.getParameter<double>("etaMult")),
    minBx_(iConfig.getParameter<int>("minBx")),
    maxBx_(iConfig.getParameter<int>("maxBx")),
    useDummyValues_(iConfig.getParameter<bool>("useDummyValues")),
    filenameNNBmtfStubs_(iConfig.getParameter<std::string>("filenameNNBmtfStubs")),
    outputFile(filenameNNBmtfStubs_)
{
    initCsvHeader();
}



// ----------------------- method called for each orbit  -----------------------
void NNBmtfStubExtractor::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    using namespace edm;

    // Initialization
    initializeMuPropagators(iSetup);

    // Fetch data from the event
    BXVector<l1t::Muon> vGmtMuons = iEvent.get(simGmtStage2DigisToken_);
    vector<reco::Muon> vRecoMuons = iEvent.get(recoMuonsToken_);
    auto const theMagField        = iSetup.getHandle(magfieldToken_);

    // Iterate over reconstructed muons
    int recoMuonIndex = -1;
    for (const auto& recoMuon : vRecoMuons) {
        ++recoMuonIndex;

        // Skip if not a global muon (currently not implemented)
        if (!recoMuonCheck(recoMuon)) continue;

        // TODO: tight muon selection
        // isRecoTight = ?
        bool isRecoTight = true;
        if (!isRecoTight) continue;

        // Setup for propagation and matching
        auto ftsRecoMuon = createFTSForRecoMuon(recoMuon, theMagField);
        auto [stateAtStation1, stateAtStation2] = propagateToStations(ftsRecoMuon);

        // Validate propagation results
        if (!stateAtStation1.isValid() || !stateAtStation2.isValid()) continue;

        // Match reco muon to L1 muon
        auto [minDeltaR, matchedRegionalMuonIndex, matchedKBmtfTrackIndex, matchedGmtMuonIndex] = findClosestL1Match(iEvent, stateAtStation2, vGmtMuons);

        // Skip if no L1 muon was matched
        if ( !(minDeltaR < drCut_) ) continue;
        
        // Get RECO matches (L1, KBMTF, GMT)
        const auto& l1tRegMuon    = iEvent.get(simKBmtfDigisToken_)[matchedRegionalMuonIndex];
        const auto& l1tKBmtfTrack = iEvent.get(simKBmtfTracksToken_)[matchedKBmtfTrackIndex];
        const auto& l1tGmtMuon    = vGmtMuons[matchedGmtMuonIndex];

        // Write data to file 
        writeEventInfoToFile(iEvent);
        writeRecoToFile(recoMuon, stateAtStation2);
        writeL1ToFile(l1tRegMuon);
        writeKBmtfToFile(l1tKBmtfTrack);
        writeGmtToFile(l1tGmtMuon);

        unsigned int nStubs = l1tKBmtfTrack.stubs().size();

        // write the number of stubs
        outputFile << nStubs << ",";

        // Sort the stubs by their station numbers
        L1MuKBMTCombinedStubRefVector sortedStubs = l1tKBmtfTrack.stubs(); // Assuming this returns a copy or a const reference

        // Sorting lambda
        std::sort(sortedStubs.begin(), sortedStubs.end(), [](const auto& a, const auto& b) {
            return a->stNum() < b->stNum();
        });

        // Write stubs to file
        writeStubsToFile(sortedStubs);

        outputFile << std::endl;
        
    }
    //--------------------------------------------------------------------------

}
// -----------------------------------------------------------------------------


// ------------------------------ Helper methods  -----------------------------

// Helper function to initialize the CSV header
void NNBmtfStubExtractor::initCsvHeader() {
    outputFile    << "run,ls,event,"
                        // RECO quantities
                        << "ptReco,etaExtRecoSt2,phiExtRecoSt2,chargeReco,"
                        << "etaVtxReco,phiVtxReco,dXYReco,isTight,"
                        // L1 quantities
                        << "ptL1,pt2L1,hwDXYL1,phiL1,etaL1,hwSignL1,hwSignValidL1,hwQualityL1,"
                        // KBMTF quantities
                        << "kbmtfRVtx,kbmtfPhiVtx,kbmtfRMu,kbmtfPhiMu,kbmtfPhiBendMu,kbmtfFloatPTUnconstrained,kbmtfDXY,"
                        << "kbmtfApproxChi2,kbmtfCoarseEta,kbmtfRank,kbmtfFineEta,kbmtfHasFineEta,"
                        << "kbmtfSector,kbmtfHitPattern,kbmtfQuality,kbmtfWheel,kbmtfStep,"
                        // GMT quantities
                        << "GmtHwEtaAtVtx,GmtHwPhiAtVtx,GmtEtaAtVtx,GmtPhiAtVtx,GmtPt,GmtPtU,GmtEta,GmtPhi,"
                        // STUBS quantities
                        << "n_stubs,"
                        << "s1_stNum,s1_scNum,s1_whNum,s1_eta,s1_qeta,s1_phi,s1_phiB,s1_quality,"
                        << "s2_stNum,s2_scNum,s2_whNum,s2_eta,s2_qeta,s2_phi,s2_phiB,s2_quality,"
                        << "s3_stNum,s3_scNum,s3_whNum,s3_eta,s3_qeta,s3_phi,s3_phiB,s3_quality,"
                        << "s4_stNum,s4_scNum,s4_whNum,s4_eta,s4_qeta,s4_phi,s4_phiB,s4_quality,"
                        << std::endl;
}



// Helper function to initialize the propagators
void NNBmtfStubExtractor::initializeMuPropagators(const edm::EventSetup& iSetup) {
    muPropagator1st_ = muPropagator1stSetup_.init(iSetup);
    muPropagator2nd_ = muPropagator2ndSetup_.init(iSetup);
}



// Helper function to create a FreeTrajectoryState for a reco muon
FreeTrajectoryState NNBmtfStubExtractor::createFTSForRecoMuon(const reco::Muon& recoMuon, edm::ESHandle<MagneticField> magneticField) {
    return FreeTrajectoryState(
        GlobalPoint(recoMuon.vx(), recoMuon.vy(), recoMuon.vz()),
        GlobalVector(recoMuon.px(), recoMuon.py(), recoMuon.pz()),
        recoMuon.charge(),
        magneticField.product()
    );
}



// Helper function to propagate the reco muon to the muon stations
std::tuple<TrajectoryStateOnSurface, TrajectoryStateOnSurface> NNBmtfStubExtractor::propagateToStations(const FreeTrajectoryState& ftsRecoMuon) {
    auto stateAtStation1 = muPropagator1st_.extrapolate(ftsRecoMuon);
    auto stateAtStation2 = muPropagator2nd_.extrapolate(ftsRecoMuon);
    return {stateAtStation1, stateAtStation2};
}



// Check if the reco muon is valid
bool NNBmtfStubExtractor::recoMuonCheck(const reco::Muon& recoMuon) const {
    return true; 
}



// Find the closest L1 match to the reco muon
std::tuple<double, int, int, int> NNBmtfStubExtractor::findClosestL1Match(
    const edm::Event& iEvent, 
    const TrajectoryStateOnSurface& stateAtStation2, 
    const BXVector<l1t::Muon>& detectedGmtMuons) {

    double minDeltaR = 100.0;
    int matchedRegionalMuonIndex = -1, matchedKBmtfTrackIndex = -1, matchedGmtMuonIndex = -1;

    // Loop over KBMTF l1tRegionalMuonCands
    int regionalMuonIndex = -1;
    for (const auto& regionalMuon : iEvent.get(simKBmtfDigisToken_)) {
        ++regionalMuonIndex;

        // Loop over KBMTF L1MuKBMTrack
        int KBmtfTrackIndex = -1;
        for (const auto& KBmtfTrack : iEvent.get(simKBmtfTracksToken_)) {
            ++KBmtfTrackIndex;

            // Loop over uGMT muons
            int gmtMuonIndex = -1;
            for(const auto& gmtMuon : detectedGmtMuons) {
                ++gmtMuonIndex;

                if (doL1Matching(regionalMuon, KBmtfTrack, gmtMuon)) {
                    double deltaR = calcDr(&gmtMuon, &stateAtStation2);
                    if (deltaR < minDeltaR) {
                        minDeltaR                = deltaR;
                        matchedRegionalMuonIndex = regionalMuonIndex;
                        matchedKBmtfTrackIndex   = KBmtfTrackIndex;
                        matchedGmtMuonIndex      = gmtMuonIndex;
                    }
                }
            }
        }
    }
    return {minDeltaR, matchedRegionalMuonIndex, matchedKBmtfTrackIndex, matchedGmtMuonIndex};
}



// Check if the L1 muon matches the reco muon
bool NNBmtfStubExtractor::doL1Matching(const l1t::RegionalMuonCand& regionalMuon, const L1MuKBMTrack& KBmtfTrack, const l1t::Muon& gmtMuon) {
    bool isCoarseEtaMatch = !KBmtfTrack.hasFineEta() && (regionalMuon.hwEta() == KBmtfTrack.coarseEta());
    bool isFineEtaMatch   =  KBmtfTrack.hasFineEta() && (regionalMuon.hwEta() == KBmtfTrack.fineEta());

    bool isEtaMatch             = isFineEtaMatch || isCoarseEtaMatch;
    bool isPtUnconstrainedMatch = regionalMuon.hwPtUnconstrained() == gmtMuon.hwPtUnconstrained();
    bool isDxyMatch             = regionalMuon.hwDXY() == gmtMuon.hwDXY();

    return isEtaMatch && isPtUnconstrainedMatch && isDxyMatch;
}



// method that computes the globalPhi coordinate of the muon in the hardware units
int NNBmtfStubExtractor::calcGlobalPhi(const l1t::RegionalMuonCand& muon) const {
    int globalPhi = muon.processor()*48 + muon.hwPhi();
    globalPhi += 576 - 24;
    globalPhi = globalPhi % 576;
    return globalPhi;
}



// method that computes the physical phi coordinate of the muon
double NNBmtfStubExtractor::calcPhysPhi(int globalPhi) const {
    double physPhi = 1.0 * globalPhi / phiMult_;
    if (physPhi > M_PI) {
        physPhi = physPhi - 2*M_PI;
    }
    return physPhi;
}



// method that computes the physical eta coordinate of the muon
double NNBmtfStubExtractor::calcPhysEta(int hwEta) const {
    double physEta = 1.0 * hwEta / etaMult_;
    return physEta;
}



// method that computes the physical pt of the muon
double NNBmtfStubExtractor::calcPtGeV(unsigned int pt) const {
    return pt * PT_CONVERSION_FACTOR;
}



// method that computes the deltaR between a BMTF muon and a GMT muon
double NNBmtfStubExtractor::calcDr(const l1t::RegionalMuonCand& bmtfMuon, const l1t::Muon* gmtMuon) const {
    double dPhi = fabs(calcPhysPhi(calcGlobalPhi(bmtfMuon)) - calcPhysPhi(gmtMuon->hwPhi()));
    double dEta = fabs(calcPhysEta(bmtfMuon.hwEta()) - calcPhysEta(gmtMuon->hwEta()));
    return sqrt(dPhi*dPhi + dEta*dEta);
}



// method that computes the deltaR between a GMT muon and a reco track propagated to St2
double NNBmtfStubExtractor::calcDr(const l1t::Muon* l1tGmtMuon, const TrajectoryStateOnSurface* track) const {
    double l1tPhysPhi = l1tGmtMuon->hwPhi() / phiMult_;
    if (l1tPhysPhi > M_PI) {
        l1tPhysPhi = l1tPhysPhi - 2*M_PI;
    }

    double dphi = abs(l1tPhysPhi - track->globalPosition().phi());
    if (dphi > M_PI) {
        dphi = std::abs(dphi - 2*M_PI);
    }

    double deta = l1tGmtMuon->hwEta() / etaMult_ - track->globalPosition().eta();
    double dr = std::sqrt(dphi*dphi + deta*deta);

    return dr;
}



// method to check if the GMT muon is in the barrel to match it to a BMTF muon
bool NNBmtfStubExtractor::isBarrel(const l1t::Muon* gmtMuon) const {
    // the GMT muon is in the barrel it has tfIndex between 36 and 70 (included)
    return (gmtMuon->tfMuonIndex() >= 36) && (gmtMuon->tfMuonIndex() <= 70);
}



void NNBmtfStubExtractor::writeEventInfoToFile(const edm::Event& iEvent) {
    outputFile << iEvent.id().run() << "," << iEvent.luminosityBlock() << "," << iEvent.id().event() << ",";
}



void NNBmtfStubExtractor::writeRecoToFile(const reco::Muon& recoMuon, const TrajectoryStateOnSurface& stateAtMuSt2) {
    outputFile  << recoMuon.pt() << ","                                            // "ptReco,"
                << stateAtMuSt2.globalPosition().eta() << ","                      // "etaExtRecoSt2,"
                << stateAtMuSt2.globalPosition().phi() << ","                      // "phiExtRecoSt2,"
                << recoMuon.charge() << ","                                        // "chargeReco,"
                << recoMuon.eta() << ","                                           // "etaVtxReco,"
                << recoMuon.phi() << ","                                           // "phiVtxReco,"
                << "-1" << ","                                                     // "dXYReco,"
                << "1"  << ",";                                                    // "isTight,"
}



void NNBmtfStubExtractor::writeL1ToFile(const l1t::RegionalMuonCand& l1tRegMuon) {
    outputFile  << calcPtGeV(l1tRegMuon.hwPt()) << ","                             // "ptL1,"
                << calcPtGeV(l1tRegMuon.hwPtUnconstrained()) << ","                // "pt2L1,"
                << l1tRegMuon.hwDXY() << ","                                       // "hwDXYL1,"
                << calcPhysPhi(calcGlobalPhi(l1tRegMuon)) << ","                   // "phiL1,"
                << calcPhysEta(l1tRegMuon.hwEta()) << ","                          // "etaL1,"
                << l1tRegMuon.hwSign() << ","                                      // "hwSignL1,"
                << l1tRegMuon.hwSignValid() << ","                                 // "hwSignValidL1,"
                << l1tRegMuon.hwQual() << ",";                                     // "hwQualityL1,"
}



void NNBmtfStubExtractor::writeKBmtfToFile(const L1MuKBMTrack& l1tKBmtfTrack) {
    outputFile  << l1tKBmtfTrack.curvatureAtVertex() << ","                        // "kbmtfRVtx,"
                << l1tKBmtfTrack.phiAtVertex() << ","                              // "kbmtfPhiVtx,"
                << l1tKBmtfTrack.curvature() << ","                                // "kbmtfRMu,"
                << l1tKBmtfTrack.positionAngle() << ","                            // "kbmtfPhiMu,"
                << l1tKBmtfTrack.bendingAngle() << ","                             // "kbmtfPhiBendMu,"
                << l1tKBmtfTrack.ptUnconstrained() << ","                          // "kbmtfFloatPTUnconstrained,"
                << l1tKBmtfTrack.dxy() << ","                                      // "kbmtfDXY,"
                << l1tKBmtfTrack.approxChi2() << ","                               // "kbmtfApproxChi2,"
                << l1tKBmtfTrack.coarseEta() << ","                                // "kbmtfCoarseEta,"
                << l1tKBmtfTrack.rank() << ","                                     // "kbmtfRank,"
                << l1tKBmtfTrack.fineEta() << ","                                  // "kbmtfFineEta,"
                << l1tKBmtfTrack.hasFineEta() << ","                               // "kbmtfHasFineEta,"
                << l1tKBmtfTrack.sector() << ","                                   // "kbmtfSector,"
                << l1tKBmtfTrack.hitPattern() << ","                               // "kbmtfHitPattern,"
                << l1tKBmtfTrack.quality() << ","                                  // "kbmtfQuality,"
                << l1tKBmtfTrack.wheel() << ","                                    // "kbmtfWheel,"
                << l1tKBmtfTrack.step() << ",";                                    // "kbmtfStep,"
}

void NNBmtfStubExtractor::writeGmtToFile(const l1t::Muon& l1tGmtMuon) {
    outputFile  << l1tGmtMuon.hwEtaAtVtx() << ","                                  // "GmtHwEtaAtVtx,"
                << l1tGmtMuon.hwPhiAtVtx() << ","                                  // "GmtHwPhiAtVtx,"
                << l1tGmtMuon.etaAtVtx() << ","                                    // "GmtEtaAtVtx,"
                << l1tGmtMuon.phiAtVtx() << ","                                    // "GmtPhiAtVtx,"
                << calcPtGeV(l1tGmtMuon.hwPt()) << ","                             // "GmtPt,"
                << calcPtGeV(l1tGmtMuon.hwPtUnconstrained()) << ","                // "GmtPtU,"
                << calcPhysEta(l1tGmtMuon.hwEta()) << ","                          // "GmtEta,"
                << calcPhysPhi(l1tGmtMuon.hwPhi()) << ",";                         // "GmtPhi,"
}



void NNBmtfStubExtractor::writeStubsToFile(const L1MuKBMTCombinedStubRefVector& sortedStubs) {
    
    // Map to store stubs according to their station number
    std::map<int, L1MuKBMTCombinedStubRef> stationToStubMap;  

    // Fill the map
    for (const auto& stub : sortedStubs) {
        stationToStubMap[stub->stNum()] = stub;
    }

    // Initialize the iterator for cycling through the stubs
    auto cyclicIt = sortedStubs.begin();

    // Write stubs or dummy/duplicated stubs to file in station order
    for (int station = 1; station <= 4; ++station) {
        auto it = stationToStubMap.find(station);
        
        if (it != stationToStubMap.end()) {
            // Write stub information
            const auto& stubRef = it->second;
            writeStubToFile(*(stubRef));
        } else {        
            if (useDummyValues_) {
                // Write dummy stub information
                outputFile << "-1,-1,-3,-1,0,0,0,0,";
            } else {
                // Write duplicated stub information
                if (cyclicIt == sortedStubs.end()) {
                    cyclicIt = sortedStubs.begin();
                }
                const auto& stubRef = *cyclicIt;
                writeStubToFile(*(stubRef));

                // Move to the next stub for potential further duplication
                ++cyclicIt;
            }
        }
    }
}



void NNBmtfStubExtractor::writeStubToFile(const L1MuKBMTCombinedStub& stub) {
    outputFile  << stub.stNum() << ","
                << stub.scNum() << ","
                << stub.whNum() << ",";

    if (stub.tag()==0) {
        outputFile << stub.eta1() << ",";
        outputFile << stub.qeta1() << ",";
    } else {
        outputFile << stub.eta2() << ",";
        outputFile << stub.qeta2() << ",";
    }

    outputFile  << stub.phi() << "," 
                << stub.phiB() << ","
                << stub.quality() << ",";
}



// define this as a plug-in
DEFINE_FWK_MODULE(NNBmtfStubExtractor);
