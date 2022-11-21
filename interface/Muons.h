#ifndef Muons_h
#define Muons_h

#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DarkPhoton/MuAnalyzer/interface/RoccoR.h"
#include "TH2D.h"
#include "TRandom3.h"

class Muons {
public:
  Muons(bool xrootd = false);

  int SelectMuons(const edm::Event& iEvent, edm::EDGetToken m_recoMuonToken);
  double MuonRoccorWeight(const reco::Muon* muon, edm::Handle<reco::GenParticleCollection> genParticles);
  std::vector<double> MuonTightIsoWeight(double MuonPt, double MuonEta);
  std::vector<double> MuonTightIDWeight(double MuonPt, double MuonEta);
  std::vector<double> MuonTrigWeight(double MuonPt, double MuonEta);

  RoccoR rc;
  TRandom3 randEng;
  TH2D* Muon_TightID;
  TH2D* Muon_TightISO;
  TH2D* Muon_Trig;
  std::vector<const reco::Muon*> selectedMuons;
  std::vector<const reco::Muon*> selectedEndcapMuons;
  const reco::Muon* highPtSelectedMuon;
  const reco::Muon* highPtSelectedEndcapMuon;
  const reco::Muon* secondhighPtSelectedMuon;

  double highmuonpt;
  double secondhighmuonpt;
  double highendcappt;
  double muonisolation;

  double muidweight;
  double muidweightUp;
  double muidweightDown;
  double muisoweight;
  double muisoweightUp;
  double muisoweightDown;
  double mutrigweight;
  double mutrigweightUp;
  double mutrigweightDown;
};

#endif
