#ifndef EventInfo_h
#define EventInfo_h

#include <vector>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TH1F.h"

class EventInfo {
public:
  EventInfo();
  bool passTriggers(const edm::Event& iEvent,
		    edm::EDGetToken m_trigResultsToken,
		    edm::EDGetToken m_trigEventToken,
		    std::vector<std::string> m_muonPathsToPass);
  bool matchTagAndTrigger(const reco::Muon* selectedMuon);
  double eventWeight;
  bool filledWeights;
  std::vector<double> PuUpDownWeights;
  double muRocWeight;
  std::vector<double> muIdWeights;
  std::vector<double> muIsoWeights;
  std::vector<double> muTrigWeights;
  int nTagMuons;
  double probeTrackIso;
  double probeEcalIso;
  double probeReducedChi;
  double tagMuonEta;
  double tagMuonPt;
  double probeTrackEta;
  double probeTrackPhi;
  double probeTrackPt;
  double probeTrackP;
  double DtHitPhiByDepth[4];
  double DtHitZByDepth[4];
  double pairVtxChi;
  double diMuonMass;
  //Transient Track Studies
  double closestApproach;
  double staMinDr;
  double standaloneDEoverE;
  double dxy;
  double dsz;
  double dca;
  double dcal;
  //Calo Jets
  double minCaloJetDr;
  double caloJetHcalE;
  double caloJetEcalE;
  double caloJetTotalE;
  double caloConeE;
  int caloConeHits;
  //Gen Matching
  double minGenMuDr;
  double minGenMuDE;
  //Standalone muon
  double probeMuondE;
  double probeMuondR;
  double standaloneE;
  double globalMuonE;
  int probeCharge;
  double stadEta;
  double stadPhi;
  double staChi2;
  int staNHits;
  double staTransDCA;
  double staLongDCA;
  //pileup
  double pileupWeight;
  int nPUmean;
  //Trigger objects
  std::vector<trigger::TriggerObject> passingMuons;
  //HCAL
  int nCellsFound;
  int nNeighborCellsFound;
  bool foundDepths[7];
  int iEta;
  double cellCenterDEta;
  double cellCenterDPhi;
  double hitEnergies[7];
  double neighborHitEnergies[7];
  int bremDepth;
  int hitsOverThresh;
  int coneHits;
  double coneEnergy;
  int expectedHits;
  double HOMuonHitEnergy;
  double HOMuonHitDr;
};
#endif
