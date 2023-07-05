from __future__ import print_function
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

options = VarParsing.VarParsing ('analysis')

options.register ('runNumber',
                  368318, # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Run Number")

options.register ('daqSourceMode',
                  'ScoutingRun3', # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "DAQ source data mode")

options.register ('buBaseDir',
                  '.', # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "BU base directory")

options.register ('fuBaseDir',
                  '.', # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "BU base directory")

options.register ('fffBaseDir',
                  '.', # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "FFF base directory")

options.register ('numThreads',
                  2, # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Number of CMSSW threads")

options.register ('numFwkStreams',
                  1, # default value
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Number of CMSSW streams")


options.parseArguments()

cmsswbase = os.path.expandvars("$CMSSW_BASE/")

process = cms.Process("SCPU")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(options.numThreads),
    numberOfStreams = cms.untracked.uint32(options.numFwkStreams),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1) # ShmStreamConsumer requires synchronization at LuminosityBlock boundaries
)
process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(threshold = cms.untracked.string( "DEBUG" )),
    destinations = cms.untracked.vstring( 'cout' )
)

process.FastMonitoringService = cms.Service("FastMonitoringService",
    sleepTime = cms.untracked.int32(1)
)

process.Timing = cms.Service("Timing",
  summaryOnly = cms.untracked.bool(False),
  useJobReport = cms.untracked.bool(True)
)

process.EvFDaqDirector = cms.Service("EvFDaqDirector",
    useFileBroker = cms.untracked.bool(False),
    fileBrokerHostFromCfg = cms.untracked.bool(True),
    fileBrokerHost = cms.untracked.string("htcp40.cern.ch"),
    runNumber = cms.untracked.uint32(options.runNumber),
    baseDir = cms.untracked.string(options.fffBaseDir+"/"+options.fuBaseDir),
    buBaseDir = cms.untracked.string(options.fffBaseDir+"/"+options.buBaseDir),
    directorIsBU = cms.untracked.bool(False),
)

try:
  os.makedirs(options.fffBaseDir+"/"+options.fuBaseDir+"/run"+str(options.runNumber).zfill(6))
except Exception as ex:
  print(str(ex))
  pass

ram_dir_path=options.buBaseDir+"/run"+str(options.runNumber).zfill(6)+"/"

process.source = cms.Source("DAQSource",
    testing = cms.untracked.bool(True),
    dataMode = cms.untracked.string(options.daqSourceMode),
    verifyChecksum = cms.untracked.bool(False),
    useL1EventID = cms.untracked.bool(False),
    eventChunkBlock = cms.untracked.uint32(2 * 1024),
    eventChunkSize = cms.untracked.uint32(3 * 1024),
    maxChunkSize = cms.untracked.uint32(10 * 1024),
    numBuffers = cms.untracked.uint32(2),
    maxBufferedFiles = cms.untracked.uint32(2),
    fileListMode = cms.untracked.bool(True),
    fileNames = cms.untracked.vstring(
        ram_dir_path + "run" + str(options.runNumber) + "_ls0150_index000000.raw"
    )

)

## test pluging
process.GmtUnpacker = cms.EDProducer('ScGMTRawToDigi',
  srcInputTag = cms.InputTag('rawDataCollector'),
  debug=cms.untracked.bool(True)
)
process.CaloUnpacker = cms.EDProducer('ScCaloRawToDigi',
  srcInputTag = cms.InputTag('rawDataCollector'),
  debug=cms.untracked.bool(False)
)
process.BmtfUnpacker = cms.EDProducer('ScBMTFRawToDigi',
  srcInputTag = cms.InputTag('rawDataCollector'),
  debug=cms.untracked.bool(True),
  cotTheta_1=cms.vint32(105,101,97,93,88,84,79,69,64,58,52,46,40,34,21,14,7,0,-7,-14,-21,-34,-40,-46,-52,-58,-64,-69,-79,-84,-88,-93,-97,-101,-105),
  cotTheta_2=cms.vint32(93,89,85,81,77,73,68,60,55,50,45,40,34,29,17,12,6,0,-6,-12,-17,-29,-34,-40,-45,-50,-55,-60,-68,-73,-77,-81,-85,-89,-93),
  cotTheta_3=cms.vint32(81,77,74,70,66,62,58,51,46,42,38,33,29,24,15,10,5,0,-5,-10,-15,-24,-29,-33,-38,-42,-46,-51,-58,-62,-66,-70,-74,-77,-81),
)

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('file:./PoolOutputTest.root'),
    outputCommands = cms.untracked.vstring(
        "keep *",
        "keep *_rawDataCollector_*_*",
        # "keep *_GmtUnpacker_*_*",
        # "keep *_CaloUnpacker_*_*",
        "keep *_BmtfUnpacker_*_*"
        ),
    #compressionLevel = cms.untracked.int32(1)
    )

rawToDigiTask = cms.Task(
  # process.GmtUnpacker,
  # process.CaloUnpacker,
  process.BmtfUnpacker
)

"""
bmtfKalmanTrackingSettings = cms.PSet(
    verbose = cms.bool(False),  #
    lutFile = cms.string("L1Trigger/L1TMuon/data/bmtf_luts/kalmanLUTs_v302.root"),
    initialK = cms.vdouble(-1.196,-1.581,-2.133,-2.263),
    initialK2 = cms.vdouble(-3.26e-4,-7.165e-4,2.305e-3,-5.63e-3),
#    eLoss = cms.vdouble(-2.85e-4,-6.21e-5,-1.26e-4,-1.23e-4),
    eLoss = cms.vdouble(+0.000765,0,0,0),

    aPhi = cms.vdouble(1.942, .01511, .01476, .009799),
    aPhiB = cms.vdouble(-1.508,-0.1237,-0.1496,-0.1333),
    aPhiBNLO = cms.vdouble(0.000331,0,0,0),

    bPhi = cms.vdouble(-1,.18245,.20898,.17286),
    bPhiB = cms.vdouble(-1,1.18245,1.20898,1.17286),
    phiAt2 = cms.double(0.15918),
    etaLUT0 = cms.vdouble(8.946,7.508,6.279,6.399),
    etaLUT1 = cms.vdouble(0.159,0.116,0.088,0.128),
    #generic cuts
    chiSquare = cms.vdouble(0.0,0.109375,0.234375,0.359375),
    chiSquareCutPattern = cms.vint32(7,11,13,14,15),
    chiSquareCutCurvMax = cms.vint32(2500,2500,2500,2500,2500),
    chiSquareCut = cms.vint32(126,126,126,126,126),


    #vertex cuts
    trackComp = cms.vdouble(1.75,1.25,0.625,0.250),
    trackCompErr1 = cms.vdouble(2.0,2.0,2.0,2.0),
    trackCompErr2 = cms.vdouble(0.218750,0.218750,0.218750,0.3125),
    trackCompCutPattern = cms.vint32(3,5,6,9,10,12),
    trackCompCutCurvMax = cms.vint32(34,34,34,34,34,34),   #this is shifted<<4
    trackCompCut        = cms.vint32(15,15,15,15,15,15),
    chiSquareCutTight   = cms.vint32(40,126,60,126,126,126),

    combos4=cms.vint32(9,10,11,12,13,14,15),
    combos3=cms.vint32(5,6,7),
    combos2=cms.vint32(3),
    combos1=cms.vint32(), #for future possible usage

    useOfflineAlgo = cms.bool(False),
    ###Only for the offline algo -not in firmware --------------------
    mScatteringPhi = cms.vdouble(2.49e-3,5.47e-5,3.49e-5,1.37e-5),
    mScatteringPhiB = cms.vdouble(7.22e-3,3.461e-3,4.447e-3,4.12e-3),
    pointResolutionPhi = cms.double(1.),
    pointResolutionPhiB = cms.double(500.),
    pointResolutionPhiBH = cms.vdouble(151., 173., 155., 153.),
    pointResolutionPhiBL = cms.vdouble(17866., 19306., 23984., 23746.),
    pointResolutionVertex = cms.double(1.)
)



process.simKBmtfDigis = cms.EDProducer("L1TMuonBarrelKalmanTrackProducer",
    src = cms.InputTag("BmtfUnpacker"),
    bx = cms.vint32(-2,-1,0,1,2),
#    bx = cms.vint32(0),
    algoSettings = bmtfKalmanTrackingSettings,
    trackFinderSettings = cms.PSet(
        sectorsToProcess = cms.vint32(0,1,2,3),
        verbose = cms.int32(0),
        sectorSettings = cms.PSet(
#            verbose = cms.int32(1),
            verbose = cms.int32(0),
            wheelsToProcess = cms.vint32(-2,-1,0,1,2),
            regionSettings = cms.PSet(
                verbose=cms.int32(0)
            )
        )

    )
)
"""


process.p1 = cms.Path(rawToDigiTask)
#process.p = cms.Path(process.GmtUnpacker*process.CaloUnpacker)
# process.p2 = cms.Path(process.simKBmtfDigis)

process.ep = cms.EndPath(
    process.output
)