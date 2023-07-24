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

process = cms.Process("TEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

if(options.isMC):
        process.GlobalTag.globaltag = cms.string('106X_upgrade2018_realistic_v15_L1v1')
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
	    mylist = FileUtils.loadListFromFile ('datafiles/1p0sig_gen_data.txt')#DYJets_mumuFilter_combined.txt
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
                                   fileName = cms.string("DiMuon_Histos.root"),
				   closeFileFast = cms.untracked.bool(False)
)

process.p = cms.Path(process.demo)
