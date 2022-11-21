#ifndef Histograms_h
#define Histograms_h

#include <vector>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HBHEChannelInfo.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/Common/interface/SortedCollection.h"
#include "DarkPhoton/MuAnalyzer/interface/EventInfo.h"

class TTree;
class TH1F;
class TH2F;

class Histograms {
public:
  Histograms();

  void book(TFileDirectory histFolder,bool MC);
  void FillHists(EventInfo info);
  void Normalize();
  void IncCutFlow(double weight);
  void FailLocation(double ieta, double iphi, double weight);
  void ResetCutFlow();

public:
  bool isMC;
  //Histograms
  TH1F* m_eventCount;
  TH1F* m_eventWeight;
  TH1F* m_muRocWeight;
  TH1F* m_dbremWeight;
  //Scale factors
  TH1F* m_IDSF;
  TH1F* m_IDSFup;
  TH1F* m_IDSFdown;
  TH1F* m_ISOSF;
  TH1F* m_ISOSFup;
  TH1F* m_ISOSFdown;
  TH1F* m_TrigSF;
  TH1F* m_TrigSFup;
  TH1F* m_TrigSFdown;
  TH1F* m_eventCount_PUup;
  TH1F* m_eventCount_PUdown;
  TH1F* m_eventCount_IDup;
  TH1F* m_eventCount_IDdown;
  TH1F* m_eventCount_ISOup;
  TH1F* m_eventCount_ISOdown;
  TH1F* m_eventCount_Trigup;
  TH1F* m_eventCount_Trigdown;
  TH1F* m_cutProgress;
  TH1F* m_ProbeTrackPt;
  TH1F* m_ProbeTrackEta;
  TH2F* m_ProbeTrackEtaPhi;
  TH1F* m_ProbeTrackP;
  TH1F* m_ProbeMuonEta;
  TH1F* m_ProbeMuonPt;
  TH2F* m_ProbeTrackPtEta;
  TH1F* m_minGenMuDr;
  TH1F* m_minGenMuDE;
  TH1F* m_initialMuE;
  TH1F* m_ProbeTrackIsolation;
  TH1F* m_ProbeEcalIsolation;
  int cutProgress;
  TH1F* m_StandaloneMuonE;
  TH1F* m_GlobalMuonE;
  TH1F* m_StandalonedE;
  TH1F* m_StandalonedE_PUup;
  TH1F* m_StandalonedE_PUdown;
  TH1F* m_StandalonedE_IDup;
  TH1F* m_StandalonedE_IDdown;
  TH1F* m_StandalonedE_ISOup;
  TH1F* m_StandalonedE_ISOdown;
  TH1F* m_StandalonedE_Trigup;
  TH1F* m_StandalonedE_Trigdown;
  TH1F* m_AllVertexOffset;
  TH1F* m_EventWeights;
  TH1F* m_DrtoStandalone;
  TH1F* m_VtxZ;
  TH1F* m_CaloJetMinDr;
  TH1F* m_CaloJetEcalE;
  TH1F* m_CaloJetHcalE;
  TH1F* m_CaloJetTotalE;
  TH1F* m_DiMuonMass;
  TH1F* m_TagMuonEta;
  TH1F* m_TagMuonPt;
  TH1F* m_nSelectedMuons;
  TH1F* m_PairVtxChi;
  TH1F* m_TrackIsolation;
  TH1F* m_TrackReducedChi;
  TH1F* m_Probedxy;
  TH1F* m_Probedsz;
  TH1F* m_ProbeDCA;
  TH1F* m_ProbeLongDCA;
  TH2F* m_ProbeDCAphi;
  TH1F* m_standalonedPhi;
  TH1F* m_standalonedPhiPosQ;
  TH1F* m_standalonedPhiNegQ;
  TH1F* m_standaloneChi2;
  TH1F* m_standaloneNHits;
  TH1F* m_staTransDCA;
  TH1F* m_staLongDCA;
  TH1F* m_PileupWeights;
  TH1F* m_PileupWeightUp;
  TH1F* m_PileupWeightDown;
  TH1F* m_PUmean;
};
#endif
