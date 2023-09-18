#-------------------------------------------------------------------------------
#
# Packages
#
#-------------------------------------------------------------------------------
import os
import math
import FWCore.ParameterSet.Config as cms
# import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.Eras.Era_Run3_cff import Run3
from Configuration.Eras.Modifier_pf_badHcalMitigationOff_cff import pf_badHcalMitigationOff
#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
#
# - Process and Configuration
#
#-------------------------------------------------------------------------------
# process
# process = cms.Process("Extractor", Run2_2018, pf_badHcalMitigation)
process = cms.Process("Extractor", Run3, pf_badHcalMitigationOff)

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")
process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# KBMTF emulator
process.load('L1Trigger.L1TMuonBarrel.simKBmtfStubs_cfi')
process.simKBmtfStubs.srcPhi = cms.InputTag('bmtfDigis')
process.simKBmtfStubs.srcTheta = cms.InputTag('bmtfDigis')
process.load('L1Trigger.L1TMuonBarrel.simKBmtfDigis_cfi')
process.load('L1Trigger.L1TMuon.simGmtStage2Digis_cfi')
process.simGmtStage2Digis.barrelTFInput  = cms.InputTag("simKBmtfDigis", "BMTF")
process.simGmtStage2Digis.overlapTFInput = cms.InputTag("omtfStage2Digis", "", "RECO")
process.simGmtStage2Digis.forwardTFInput = cms.InputTag("emtfStage2Digis", "", "RECO")
process.simGmtStage2Digis.autoBxRange    = cms.bool(True) # if True the output BX range is calculated from the inputs and 'bxMin' and 'bxMax' are ignored
process.simGmtStage2Digis.bxMin          = cms.int32(-2)
process.simGmtStage2Digis.bxMax          = cms.int32(2)
process.simGmtStage2Digis.autoCancelMode = cms.bool(False) # if True the cancel out methods are configured depending on the FW version number and 'emtfCancelMode' is ignored
process.simGmtStage2Digis.bmtfCancelMode = cms.string("kftracks")
process.simGmtStage2Digis.emtfCancelMode = cms.string("coordinate") # 'tracks' or 'coordinate'

# maximum number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# report every n events
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


# parser
options = VarParsing.VarParsing('analysis')
options.parseArguments()

if len(options.inputFiles)==0:
    print("[ERROR] No input file given")
    exit()
    
    
# input files
# options.inputFiles is the filename that contains for each line a file name
# we need to read the lines and put them in a list
fileList = []
with open(f"/home/nilai/CMSSW_13_1_0_pre4/src/_RECO_STUDIES/file_lists/{options.inputFiles[0]}") as f:
    for line in f:
        fileList.append("file:" + line.strip())

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(*fileList),
    secondaryFileNames = cms.untracked.vstring()
)

process.dumpED = cms.EDAnalyzer("EventContentAnalyzer")


# other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_Prompt_v3', '')
# process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v32', '')
#-------------------------------------------------------------------------------





#-------------------------------------------------------------------------------
#
# EDAnalyzer with all args for module constructor
#
#-------------------------------------------------------------------------------
process.extractor = cms.EDAnalyzer('NNBmtfStubExtractor',
    recoMuons       = cms.InputTag("muons", "", "RECO"),
    # L1MuonBXVector = cms.InputTag("gmtStage2Digis", "Muon", "RECO"),
    l1tMuonBXVector             = cms.InputTag("simGmtStage2Digis", "", "Extractor"),
    l1tRegionalMuonCandBXVector = cms.InputTag("simKBmtfDigis","BMTF"),
    L1MuKBMTrackBXVector        = cms.InputTag("simKBmtfDigis", "", "Extractor"),

    muProp1st = cms.PSet(
        useTrack          = cms.string("global"),   # 'none' to use Candidate P4; or 'tracker', 'muon', 'global'
        useState          = cms.string("atVertex"), # 'innermost' and 'outermost' require the TrackExtra
        useSimpleGeometry = cms.bool(True),
        useStation2       = cms.bool(False),
        fallbackToME1     = cms.bool(False),
        # new
        cosmicPropagationHypothesis = cms.bool(False),
        useMB2InOverlap             = cms.bool(False),
        propagatorAlong             = cms.ESInputTag("", "SteppingHelixPropagatorAlong"),
        propagatorAny               = cms.ESInputTag("", "SteppingHelixPropagatorAny"),
        propagatorOpposite          = cms.ESInputTag("", "SteppingHelixPropagatorOpposite")
    ),

    muProp2nd = cms.PSet(
        useTrack          = cms.string("global"),   # 'none' to use Candidate P4; or 'tracker', 'muon', 'global'
        useState          = cms.string("atVertex"), # 'innermost' and 'outermost' require the TrackExtra
        useSimpleGeometry = cms.bool(True),
        useStation2       = cms.bool(True),
        fallbackToME1     = cms.bool(False),
        # new
        cosmicPropagationHypothesis = cms.bool(False),
        useMB2InOverlap             = cms.bool(False),
        propagatorAlong             = cms.ESInputTag("", "SteppingHelixPropagatorAlong"),
        propagatorAny               = cms.ESInputTag("", "SteppingHelixPropagatorAny"),
        propagatorOpposite          = cms.ESInputTag("", "SteppingHelixPropagatorOpposite")
    ),

    drCut   = cms.double(0.1),
    phiMult = cms.double(576./(2*math.pi)),
    etaMult = cms.double(1./0.010875),

    minPhiQuality   = cms.int32(0),
    minThetaQuality = cms.int32(0),
    minBx           = cms.int32(-2),
    maxBx           = cms.int32(2),
    
    useDummyValues  = cms.bool(False),
    
    filenameNNBmtfStubs = cms.string("/home/nilai/CMSSW_13_1_0_pre4/src/_RECO_STUDIES/dataset/" + options.outputFile.replace(".root","") + ".csv")
)

process.ps = cms.Sequence(process.simKBmtfStubs + process.simKBmtfDigis + process.simGmtStage2Digis + process.extractor)
process.p  = cms.Path(process.ps)

process.schedule = cms.Schedule(process.p)
#-------------------------------------------------------------------------------
