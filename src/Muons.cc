#include "DarkPhoton/MuAnalyzer/interface/Muons.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <TMath.h>
#include <TRandom3.h>

Muons::Muons(bool xrootd) {
  highmuonpt = 0;
  secondhighmuonpt = 0;
  rc.init(edm::FileInPath("DarkPhoton/MuAnalyzer/text/RoccoR2018UL.txt").fullPath());
  randEng = TRandom3();
  std::string iMuon_ID_eff, iMuon_ISO_eff, iMuon_Trig_eff;
  if (xrootd) {
    iMuon_ID_eff =
        "root://cmseos.fnal.gov//store/user/revering/muonefficiencies/"
        "Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root";
    iMuon_ISO_eff =
        "root://cmseos.fnal.gov//store/user/revering/muonefficiencies/"
        "Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root";
    iMuon_Trig_eff =
        "root://cmseos.fnal.gov//store/user/revering/muonefficiencies/"
        "Efficiencies_muon_generalTracks_Z_Run2018_UL_SingleMuonTriggers.root";
  } else {
    iMuon_ID_eff =
        "${CMSSW_BASE}/src/DarkPhoton/MuAnalyzer/muonefficiencies/Run2/UL/2018/2018_Z/"
        "Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root";
    iMuon_ISO_eff =
        "${CMSSW_BASE}/src/DarkPhoton/MuAnalyzer/muonefficiencies/Run2/UL/2018/2018_Z/"
        "Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root";
    iMuon_Trig_eff =
        "${CMSSW_BASE}/src/DarkPhoton/MuAnalyzer/muonefficiencies/Run2/UL/2018/2018_trigger/"
        "Efficiencies_muon_generalTracks_Z_Run2018_UL_SingleMuonTriggers.root";
  }

  TFile* lFile_ID = TFile::Open(iMuon_ID_eff.c_str(), "READ");
  Muon_TightID = (TH2D*)lFile_ID->Get("NUM_TightID_DEN_TrackerMuons_abseta_pt");
  Muon_TightID->SetDirectory(0);

  lFile_ID->Close();

  TFile* lFile_ISO = TFile::Open(iMuon_ISO_eff.c_str(), "READ");
  Muon_TightISO = (TH2D*)lFile_ISO->Get("NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt");
  Muon_TightISO->SetDirectory(0);

  lFile_ISO->Close();
  TFile* lFile_Trig = TFile::Open(iMuon_Trig_eff.c_str(), "READ");
  Muon_Trig = (TH2D*)lFile_Trig->Get("NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt");
  Muon_Trig->SetDirectory(0);
  lFile_Trig->Close();
}

double Muons::MuonRoccorWeight(const reco::Muon* muon, edm::Handle<reco::GenParticleCollection> genParticles) {
  double genPt = -1;
  for (const auto& particle : *(genParticles.product())) {
    if (fabs(particle.pdgId()) != 13) {
      continue;
    }
    if (deltaR(particle.eta(), particle.phi(), muon->eta(), muon->phi()) > 0.1) {
      continue;
    }
    if (fabs(muon->pt() - genPt) > fabs(muon->pt() - particle.pt())) {
      genPt = particle.pt();
    }
  }

  int Q = muon->charge();
  double pt = muon->pt();
  double eta = muon->eta();
  double phi = muon->phi();
  if (genPt > 0) {
    return rc.kSpreadMC(Q, pt, eta, phi, genPt, 0, 0);
  }
  return 1;
  //   return rc.kSmearMC(Q,pt,eta,phi,nl,u,0,0);
}

std::vector<double> Muons::MuonTightIDWeight(double MuonPt, double MuonEta) {
  double muPtForId = 0.;
  double muEtaForId = 0.;

  if (MuonPt > 100) {
    muPtForId = 100.;
  } else if (MuonPt < 20.) {
    muPtForId = 21.;
  } else {
    muPtForId = MuonPt;
  }

  if (std::abs(MuonEta) > 2.3) {
    muEtaForId = 2.29;
  } else {
    muEtaForId = std::abs(MuonEta);
  }

  muidweight = Muon_TightID->GetBinContent(Muon_TightID->FindBin(muEtaForId, muPtForId));
  muidweightUp = muidweight +
                 sqrt(pow(Muon_TightID->GetBinError(Muon_TightID->FindBin(muEtaForId, muPtForId)), 2) + pow((.01), 2));
  muidweightDown = muidweight - sqrt(pow(Muon_TightID->GetBinError(Muon_TightID->FindBin(muEtaForId, muPtForId)), 2) +
                                     pow((.01), 2));

  std::vector<double> muWeights;
  muWeights.push_back(muidweight);
  muWeights.push_back(muidweightUp);
  muWeights.push_back(muidweightDown);

  return muWeights;
}

std::vector<double> Muons::MuonTightIsoWeight(double MuonPt, double MuonEta) {
  double muPtForId = 0.;
  double muEtaForId = 0.;

  if (MuonPt > 100) {
    muPtForId = 100.;
  } else if (MuonPt < 20.) {
    muPtForId = 21.;
  } else {
    muPtForId = MuonPt;
  }

  if (std::abs(MuonEta) > 2.3) {
    muEtaForId = 2.29;
  } else {
    muEtaForId = std::abs(MuonEta);
  }

  muisoweight = Muon_TightISO->GetBinContent(Muon_TightISO->FindBin(muEtaForId, muPtForId));
  muisoweightUp = muisoweight + sqrt(pow(Muon_TightISO->GetBinError(Muon_TightISO->FindBin(muEtaForId, muPtForId)), 2) +
                                     pow((.01), 2));
  muisoweightDown =
      muisoweight -
      sqrt(pow(Muon_TightISO->GetBinError(Muon_TightISO->FindBin(muEtaForId, muPtForId)), 2) + pow((.01), 2));

  std::vector<double> muWeights;
  muWeights.push_back(muisoweight);
  muWeights.push_back(muisoweightUp);
  muWeights.push_back(muisoweightDown);

  return muWeights;
}

std::vector<double> Muons::MuonTrigWeight(double MuonPt, double MuonEta) {
  double muPtForId = 0.;
  double muEtaForId = 0.;

  if (MuonPt > 100) {
    muPtForId = 100.;
  } else if (MuonPt < 20.) {
    muPtForId = 21.;
  } else {
    muPtForId = MuonPt;
  }

  if (std::abs(MuonEta) > 2.3) {
    muEtaForId = 2.29;
  } else {
    muEtaForId = std::abs(MuonEta);
  }
  mutrigweight = Muon_Trig->GetBinContent(Muon_Trig->FindBin(muEtaForId, muPtForId));
  mutrigweightUp =
      mutrigweight + sqrt(pow(Muon_Trig->GetBinError(Muon_Trig->FindBin(muEtaForId, muPtForId)), 2) + pow((.01), 2));
  mutrigweightDown =
      mutrigweight - sqrt(pow(Muon_Trig->GetBinError(Muon_Trig->FindBin(muEtaForId, muPtForId)), 2) + pow((.01), 2));

  std::vector<double> muWeights;
  muWeights.push_back(mutrigweight);
  muWeights.push_back(mutrigweightUp);
  muWeights.push_back(mutrigweightDown);
  return muWeights;
}

int Muons::SelectMuons(const edm::Event& iEvent, edm::EDGetToken m_recoMuonToken) {
  int cutPro = 0;
  edm::Handle<std::vector<reco::Muon>> recoMuons;
  iEvent.getByToken(m_recoMuonToken, recoMuons);

  for (std::vector<reco::Muon>::const_iterator iMuon = recoMuons->begin(); iMuon != recoMuons->end(); iMuon++) {
    //Tight ID
    if (!(iMuon->isPFMuon() && iMuon->isGlobalMuon()))
      continue;
    if (!(iMuon->passed(reco::Muon::CutBasedIdTight)))
      continue;
    if (cutPro == 0) {
      cutPro = 1;
    }
    //PFIso Loose requirement
    if (!(iMuon->passed(reco::Muon::PFIsoTight)))
      continue;
    if (cutPro == 1) {
      cutPro = 2;
    }
    //Track iso
    //    if(!(iMuon->passed(reco::Muon::TkIsoTight))) continue;
    //    if ((iMuon->pfIsolationR04().sumChargedHadronPt + TMath::Max(0., iMuon->pfIsolationR04().sumNeutralHadronEt + iMuon->pfIsolationR04().sumPhotonEt - 0.5*iMuon->pfIsolationR04().sumPUPt))/iMuon->pt() > 0.25) continue;

    //pT above trigger turnon
    if (iMuon->pt() < 26 || fabs(iMuon->eta()) > 2.4)
      continue;
    if (cutPro == 2) {
      cutPro = 3;
    }
    if (!iMuon->isGlobalMuon())
      continue;
    if (cutPro == 3) {
      cutPro = 4;
    }
    selectedMuons.push_back(&(*iMuon));

    if (selectedMuons.size() == 1) {
      highPtSelectedMuon = selectedMuons[0];
      highmuonpt = iMuon->pt();
      secondhighPtSelectedMuon = selectedMuons[0];
      secondhighmuonpt = 0;
    } else {
      if (iMuon->pt() > highmuonpt) {
        secondhighPtSelectedMuon = highPtSelectedMuon;
        secondhighmuonpt = iMuon->pt();
        highPtSelectedMuon = selectedMuons.back();
        highmuonpt = iMuon->pt();
      } else if (iMuon->pt() > secondhighmuonpt) {
        secondhighPtSelectedMuon = selectedMuons.back();
        secondhighmuonpt = iMuon->pt();
      }
    }
    if (selectedEndcapMuons.size() == 0) {
      highendcappt = 0;
    }
    if (fabs(iMuon->eta()) > 1.653) {
      selectedEndcapMuons.push_back(&(*iMuon));
      if (iMuon->pt() > highendcappt) {
        highendcappt = iMuon->pt();
        highPtSelectedEndcapMuon = selectedMuons.back();
      }
    }
  }
  return cutPro;
}
