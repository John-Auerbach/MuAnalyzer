import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('analysis')

options.register( 'isMC',
                                  False,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.bool,
				  "True if is MC dataset")

options.register( 'isSig',
                                  False,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "True if is signal dataset")


options.register( 'runLocally',
                                  False,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "True if running locally")

options.parseArguments()

from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process("TestAnalyzer", Run3)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
if(options.isMC):
        process.GlobalTag=GlobalTag(process.GlobalTag, '124X_mcRun3_2022_realistic_v12', '')
        #process.GlobalTag.globaltag = cms.string('106X_upgrade2018_realistic_v15_L1v1')
else:
        process.GlobalTag.globaltag = cms.string('106X_dataRun2_v32')

process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(0)
)

if(options.runLocally):
	process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
	import FWCore.Utilities.FileUtils as FileUtils
	if(options.isMC):
	    mylist = FileUtils.loadListFromFile ('datafiles/bkg_gen_data.txt')#DYJets_mumuFilter_combined.txt
	else:
	   mylist = FileUtils.loadListFromFile('datafiles/Filtered_Files_Doublemuon.txt')
	readFiles = cms.untracked.vstring( *mylist)

process.options = cms.untracked.PSet(
 SkipEvent = cms.untracked.vstring( "Error: uninitialized ProxyBase used" ),
 IgnoreCompletely = cms.untracked.vstring( "ProductNotFound" )
 #SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    fileNames = readFiles
)

process.demo = cms.EDAnalyzer('MuAnalyzer',
    recoMuons = cms.InputTag("muons"),
    tracks = cms.InputTag("generalTracks"),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    genParticles = cms.InputTag("genParticles"),
    DTSegmentLabel = cms.InputTag("dt4DSegments"),
    RPCRecHitTag = cms.InputTag("rpcRecHits"),
    edmPileupInfo = cms.InputTag("mixData"),
    trigResults = cms.InputTag("TriggerResults","","HLT"),
    muonPathsToPass = cms.vstring("HLT_IsoMu24_v","HLT_IsoMu27_v"),
    HBHERecHits = cms.InputTag("hbhereco"),
    CaloJetSource = cms.InputTag("ak4CaloJets"),
    StandAloneTracks = cms.InputTag("standAloneMuons"),
    EERecHits = cms.InputTag("reducedEcalRecHitsEE"),
    EBRecHits = cms.InputTag("reducedEcalRecHitsEB"),
    dbremweighttag = cms.InputTag("g4SimHits","DBremEventWeight"),
    isMC = cms.untracked.bool(options.isMC),
    isSig = cms.untracked.bool(options.isSig),
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("bkg_DiMuon_Histos.root"),
				   closeFileFast = cms.untracked.bool(False)
)

process.p = cms.Path(process.demo)
