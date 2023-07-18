#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "DarkPhoton/MuAnalyzer/interface/Histograms.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include "Geometry/CSCGeometry/interface/CSCLayerGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "TH2F.h"
#include "TTree.h"
#include <iostream>

Histograms::Histograms() {
 cutProgress = 0; 
}

void Histograms::book(TFileDirectory histFolder, bool MC) {
  isMC=MC;
  m_missingHits = histFolder.make<TH1F>("missingCount", "; ;Events",8,-0.5,7.5);
  m_eventCount = histFolder.make<TH1F>("eventCount", "; ;Events", 1, 0, 1);
  m_eventCount_PUup = histFolder.make<TH1F>("eventCount_PUup", "; ;Events", 1, 0, 1);
  m_eventCount_PUdown = histFolder.make<TH1F>("eventCount_PUdown", "; ;Events", 1, 0, 1);
  m_eventCount_IDup = histFolder.make<TH1F>("eventCount_IDup", "; ;Events", 1, 0, 1);
  m_eventCount_IDdown = histFolder.make<TH1F>("eventCount_IDdown", "; ;Events", 1, 0, 1);
  m_eventCount_ISOup = histFolder.make<TH1F>("eventCount_ISOup", "; ;Events", 1, 0, 1);
  m_eventCount_ISOdown = histFolder.make<TH1F>("eventCount_ISOdown", "; ;Events", 1, 0, 1);
  m_eventCount_Trigup = histFolder.make<TH1F>("eventCount_Trigup", "; ;Events", 1, 0, 1);
  m_eventCount_Trigdown = histFolder.make<TH1F>("eventCount_Trigdown", "; ;Events", 1, 0, 1);
  m_muRocWeight = histFolder.make<TH1F>("MuRocWeight", ";Muon Rochester Correction Weight; Events", 100, 0.9, 1.1);
  m_IDSF = histFolder.make<TH1F>("IDSF", ";Tight muon ID scale factor; Events", 100, 0.9, 1.1);
  m_IDSFup = histFolder.make<TH1F>("IDSFup", ";Tight muon ID scale factor up; Events", 100, 0.9, 1.1);
  m_IDSFdown = histFolder.make<TH1F>("IDSFdown", ";Tight muon ID scale factor down; Events", 100, 0.9, 1.1);
  m_ISOSF = histFolder.make<TH1F>("ISOSF", ";Tight muon ISO scale factor; Events", 100, 0.9, 1.1);
  m_ISOSFup = histFolder.make<TH1F>("ISOSFup", ";Tight muon ISO scale factor up; Events", 100, 0.9, 1.1);
  m_ISOSFdown = histFolder.make<TH1F>("ISOSFdown", ";Tight muon ISO scale factor down; Events", 100, 0.9, 1.1);
  m_TrigSF = histFolder.make<TH1F>("TrigSF", ";Tight muon Trigger scale factor; Events", 100, 0.9, 1.1);
  m_TrigSFup = histFolder.make<TH1F>("TrigSFup", ";Tight muon Trigger scale factor up; Events", 100, 0.9, 1.1);
  m_TrigSFdown = histFolder.make<TH1F>("TrigSFdown", ";Tight muon Trigger scale factor down; Events", 100, 0.9, 1.1);

  m_cutProgress = histFolder.make<TH1F>("cutProgress", ";# Cut Progress; Events passing cut level", 40, -.5, 39.5);
  std::vector<std::string> cutLabels = {"AllEvents",
                                        "Tight ID Muon",
                                        "Tight ISO Muon",
                                        "Muon Pt/Eta Req.",
                                        "Tag is global muon",
                                        "High Pt Track",
                                        "Track |eta| < 2.4",
                                        "High Purity Track",
                                        "Track Belongs to PV",
                                        "Track Long DCA",
                                        "Track Trans DCA",
                                        "Tag/Probe Opposite Charge",
                                        "Tag/Probe Fit Vtx",
                                        "Tag/Probe Vtx Chi",
                                        "Tag/Probe Inv Mass",
                                        "Probe Tracker Iso",
                                        "Probe Ecal Iso",
                                        "Tag Passed Trigger",
                                        "No Other Probes Pair"};
  for (uint8_t i = 0; i < cutLabels.size(); i++) {
    m_cutProgress->GetXaxis()->SetBinLabel(i + 1, cutLabels[i].c_str());
  }
  m_DiMuonMass = histFolder.make<TH1F>("DiMuonInvariantMass", "; Dimuon Invariant Mass (GeV); Events", 100, 50, 150);
  m_TagMuonEta = histFolder.make<TH1F>("TaggingMuonEta", "; Tagging Muon #eta; Events", 100, -2.6, 2.6);
  m_TagMuonPt = histFolder.make<TH1F>("TaggingMuonPt", "; Tagging Muon pt (GeV); Events", 200, 0, 120);
  m_nSelectedMuons = histFolder.make<TH1F>("NumberOfSelectedMuons", "; # of Selected Muons; Events", 20, -0.5, 19.5);
  m_ProbeTrackPt = histFolder.make<TH1F>("ProbeTrackPt", ";#Pt (GeV); Events", 200, 0, 100);
  m_ProbeTrackP = histFolder.make<TH1F>("ProbeTrackP",";Momentum (GeV); Events", 200,0,400);
  m_ProbeTrackEta = histFolder.make<TH1F>("ProbeTrackEta", ";#eta; Events", 200, -4, 4);
  m_ProbeTrackPtEta = histFolder.make<TH2F>("ProbeTrackPtEta", ";#Pt (GeV);#eta; Events", 50, 20, 70, 80, -4, 4);
  m_ProbeTrackEtaPhi = histFolder.make<TH2F>("ProbeTrackEtaPhi", ";#eta;#phi;Events", 100, -2.6, 2.6, 100, -3.15, 3.15);
  m_ProbeEcalIsolation = histFolder.make<TH1F>("ProbeEcalIsolation", ";Ecal Isolation; Events", 100, 0, 100);
  m_minGenMuDr = histFolder.make<TH1F>("minGenMuDr", "; #Delta R; Events", 200, 0, 4);
  m_minGenMuDE = histFolder.make<TH1F>("minGenMuDE", "; #Delta E over E; Events", 200, -1, 1);
  m_StandaloneMuonE = histFolder.make<TH1F>("StandaloneMuonE", "; Energy (GeV); Events", 200, 0, 200);
  m_GlobalMuonE = histFolder.make<TH1F>("GlobaldMuonE", "; #Delta E over E; Events", 200, 0, 200);
  m_StandalonedE =
      histFolder.make<TH1F>("StandaloneMuondEoverE", "; (Standalone E - Probe E)/Probe E; Events", 200, -1, 1);
  m_StandalonedE_PUup =
      histFolder.make<TH1F>("StandaloneMuondEoverE_PUup", "; (Standalone E - Probe E)/Probe E; Events", 200, -1, 1);
  m_StandalonedE_PUdown =
      histFolder.make<TH1F>("StandaloneMuondEoverE_PUdown", "; (Standalone E - Probe E)/Probe E; Events", 200, -1, 1);
  m_StandalonedE_IDup =
      histFolder.make<TH1F>("StandaloneMuondEoverE_IDup", "; (Standalone E - Probe E)/Probe E; Events", 200, -1, 1);
  m_StandalonedE_IDdown =
      histFolder.make<TH1F>("StandaloneMuondEoverE_IDdown", "; (Standalone E - Probe E)/Probe E; Events", 200, -1, 1);
  m_StandalonedE_ISOup =
      histFolder.make<TH1F>("StandaloneMuondEoverE_ISOup", "; (Standalone E - Probe E)/Probe E; Events", 200, -1, 1);
  m_StandalonedE_ISOdown =
      histFolder.make<TH1F>("StandaloneMuondEoverE_ISOdown", "; (Standalone E - Probe E)/Probe E; Events", 200, -1, 1);
  m_StandalonedE_Trigup =
      histFolder.make<TH1F>("StandaloneMuondEoverE_Trigup", "; (Standalone E - Probe E)/Probe E; Events", 200, -1, 1);
  m_StandalonedE_Trigdown =
      histFolder.make<TH1F>("StandaloneMuondEoverE_Trigdown", "; (Standalone E - Probe E)/Probe E; Events", 200, -1, 1);
  m_AllVertexOffset = histFolder.make<TH1F>("AllClosestApproach", "; Distance of Closest Approach; Events", 200, 0, 10);
  m_EventWeights = histFolder.make<TH1F>("EventWeights", "; Weight; Number of Events", 200, 0, 10);
  m_DrtoStandalone = histFolder.make<TH1F>(
      "DrtoNearestStandalone", "; #Delta R to Nearest Standalone muon; Number of Events", 140, 0, 7);
  m_CaloJetMinDr = histFolder.make<TH1F>("MinCaloJetDr", ";#Delta R from Probe to Nearest Calo Jet; Events", 140, 0, 7);
  m_CaloJetEcalE = histFolder.make<TH1F>("CaloJetEcalE", ";Calo Jet Ecal E (GeV); Events", 100, 0, 100);
  m_CaloJetHcalE = histFolder.make<TH1F>("CaloJetHcalE", ";Calo Jet Ecal E (GeV); Events", 100, 0, 100);
  m_CaloJetTotalE = histFolder.make<TH1F>("CaloJetTotalE", ";Calo Jet Total E (GeV); Events", 100, 0, 100);
  m_standalonedPhi = histFolder.make<TH1F>("StandalonedPhi", ";#Delta#phi;Events", 100, -1.0, 1.0);
  m_standalonedPhiPosQ = histFolder.make<TH1F>("StandalonedPhiPosQ", ";#Delta#phi;Events", 100, -1.0, 1.0);
  m_standalonedPhiNegQ = histFolder.make<TH1F>("StandalonedPhiNegQ", ";#Delta#phi;Events", 100, -1.0, 1.0);
  m_PairVtxChi =
      histFolder.make<TH1F>("PairedVtxReducedChiSq", "; Paired Vertex Reduced #Chi Squared; Events", 100, 0, 10);
  m_TrackIsolation = histFolder.make<TH1F>("TrackIsolation", "; Track Based Isolation; Events", 200, 0, 5);
  m_TrackReducedChi = histFolder.make<TH1F>("TrackReducedChi", "; Track Reduced #chi^{2}; Events", 100, 0, 30);
  m_Probedxy = histFolder.make<TH1F>("ProbeTrackdxy", "; Probe Track dxy; Events", 200, -0.3, 0.3);
  m_Probedsz = histFolder.make<TH1F>("ProbeTrackdsz", "; Probe Track dsz; Events", 200, -7, 7);
  m_ProbeDCA = histFolder.make<TH1F>("ProbeTrackDCA", "; Probe Track Transverse DCA; Events", 200, 0, 0.05);
  m_ProbeLongDCA = histFolder.make<TH1F>("ProbeTrackLongDCA", "; Probe Track Long DCA; Events", 200, 0, 0.3);
  m_ProbeDCAphi = histFolder.make<TH2F>(
      "ProbeTrackDCAphi", "; Probe Track DCA; Probe Track Phi; Events", 200, 0, 0.4, 100, -3.14, 3.14);
  m_standaloneChi2 = histFolder.make<TH1F>("standaloneChi2", ";Standalone #Chi^2;Events", 200, 0, 100);
  m_standaloneNHits = histFolder.make<TH1F>("standaloneNHits", ";Standalone N Hits;Events", 100, 0, 20);
  m_staTransDCA = histFolder.make<TH1F>("standaloneTDCA", ";Standalone transverse DCA;Events", 100, 0, 5);
  m_staLongDCA = histFolder.make<TH1F>("standaloneLongDCA", ";Standalone longitudinal DCA;Events", 100, 0, 5);
  m_PileupWeights = histFolder.make<TH1F>("PileupWeights", ";N Pileup Interactions; Weight", 100, 0, 4);
  m_PUmean = histFolder.make<TH1F>("NPUMean", ";MC N Pileup Mean; Events", 100, 0, 100);
}

void Histograms::FillHists(EventInfo info) {

  int missingHits = 0;
  for  (int depth = 0; depth < 7; depth++) {
    if (info.foundDepths[depth] && info.hitEnergies[depth] < 0.1) {
      missingHits++;
    }
  }
  m_missingHits->Fill(missingHits);

  double eventWeightPUup, eventWeightPUdown, eventWeightIDup, eventWeightIDdown, eventWeightISOup, eventWeightISOdown,
      eventWeightTrigup, eventWeightTrigdown;
  if (info.filledWeights) {
    eventWeightPUup = info.eventWeight * info.PuUpDownWeights[0] / info.pileupWeight;
    eventWeightPUdown = info.eventWeight * info.PuUpDownWeights[1] / info.pileupWeight;
    eventWeightIDup = info.eventWeight * info.muIdWeights[1] / info.muIdWeights[0];
    eventWeightIDdown = info.eventWeight * info.muIdWeights[2] / info.muIdWeights[0];
    eventWeightISOup = info.eventWeight * info.muIsoWeights[1] / info.muIsoWeights[0];
    eventWeightISOdown = info.eventWeight * info.muIsoWeights[2] / info.muIsoWeights[0];
    eventWeightTrigup = info.eventWeight * info.muTrigWeights[1] / info.muTrigWeights[0];
    eventWeightTrigdown = info.eventWeight * info.muTrigWeights[2] / info.muTrigWeights[0];

  }
  if(isMC)
  {
    m_PileupWeights->Fill(info.pileupWeight);
    m_eventCount_PUup->Fill(0.5, eventWeightPUup);
    m_eventCount_PUdown->Fill(0.5, eventWeightPUdown);
    m_eventCount_IDup->Fill(0.5, eventWeightIDup);
    m_eventCount_IDdown->Fill(0.5, eventWeightIDdown);
    m_eventCount_ISOup->Fill(0.5, eventWeightISOup);
    m_eventCount_ISOdown->Fill(0.5, eventWeightISOdown);
    m_eventCount_Trigup->Fill(0.5, eventWeightTrigup);
    m_eventCount_Trigdown->Fill(0.5, eventWeightTrigdown);
    m_muRocWeight->Fill(info.muRocWeight, info.eventWeight);
  }
  m_PUmean->Fill(info.nPUmean);
  m_eventCount->Fill(0.5, info.eventWeight);
  if (info.filledWeights) {
    m_IDSF->Fill(info.muIdWeights[0], info.eventWeight);
    m_IDSFup->Fill(info.muIdWeights[1], eventWeightIDup);
    m_IDSFdown->Fill(info.muIdWeights[2], eventWeightIDdown);
    m_ISOSF->Fill(info.muIsoWeights[0], info.eventWeight);
    m_ISOSFup->Fill(info.muIsoWeights[1], eventWeightISOup);
    m_ISOSFdown->Fill(info.muIsoWeights[2], eventWeightISOdown);
    m_TrigSF->Fill(info.muTrigWeights[0], info.eventWeight);
    m_TrigSFup->Fill(info.muTrigWeights[1], eventWeightTrigup);
    m_TrigSFdown->Fill(info.muTrigWeights[2], eventWeightTrigdown);
  }
  m_EventWeights->Fill(info.eventWeight);
  m_nSelectedMuons->Fill(info.nTagMuons, info.eventWeight);
  m_TrackIsolation->Fill(info.probeTrackIso, info.eventWeight);
  m_Probedxy->Fill(info.dxy, info.eventWeight);
  m_Probedsz->Fill(info.dsz, info.eventWeight);
  m_ProbeDCA->Fill(info.dca, info.eventWeight);
  m_ProbeLongDCA->Fill(info.dcal, info.eventWeight);
  m_ProbeDCAphi->Fill(info.dca, info.probeTrackPhi, info.eventWeight);
  m_TrackReducedChi->Fill(info.probeReducedChi, info.eventWeight);
  m_TagMuonEta->Fill(info.tagMuonEta, info.eventWeight);
  m_PairVtxChi->Fill(info.pairVtxChi, info.eventWeight);
  m_ProbeEcalIsolation->Fill(info.probeEcalIso, info.eventWeight);
  m_TagMuonPt->Fill(info.tagMuonPt, info.eventWeight);
  m_ProbeTrackPt->Fill(info.probeTrackPt, info.eventWeight);
  m_ProbeTrackP->Fill(info.probeTrackP, info.eventWeight);
  m_ProbeTrackEta->Fill(info.probeTrackEta, info.eventWeight);
  m_ProbeTrackPtEta->Fill(info.probeTrackPt, info.probeTrackEta, info.eventWeight);
  m_ProbeTrackEtaPhi->Fill(info.probeTrackEta, info.probeTrackPhi, info.eventWeight);
  m_DiMuonMass->Fill(info.diMuonMass, info.eventWeight);
  if(isMC)
  {
    m_minGenMuDr->Fill(info.minGenMuDr, info.eventWeight);
    m_minGenMuDE->Fill(info.minGenMuDE, info.eventWeight);
  }
  m_GlobalMuonE->Fill(info.globalMuonE, info.eventWeight);
  m_standalonedPhi->Fill(info.stadPhi, info.eventWeight);
  if (info.probeCharge > 0) 
  {
    m_standalonedPhiPosQ->Fill(info.stadPhi, info.eventWeight);
  } 
  else 
  {
    m_standalonedPhiNegQ->Fill(info.stadPhi, info.eventWeight);
  }

  if (info.standaloneE > 0) {
    m_StandaloneMuonE->Fill(info.standaloneE, info.eventWeight);
    m_StandalonedE->Fill(info.standaloneDEoverE, info.eventWeight);
    if (info.filledWeights) {
      m_StandalonedE_PUup->Fill(info.standaloneDEoverE, eventWeightPUup);
      m_StandalonedE_PUdown->Fill(info.standaloneDEoverE, eventWeightPUdown);
      m_StandalonedE_IDup->Fill(info.standaloneDEoverE, eventWeightIDup);
      m_StandalonedE_IDdown->Fill(info.standaloneDEoverE, eventWeightIDdown);
      m_StandalonedE_ISOup->Fill(info.standaloneDEoverE, eventWeightISOup);
      m_StandalonedE_ISOdown->Fill(info.standaloneDEoverE, eventWeightISOdown);
      m_StandalonedE_Trigup->Fill(info.standaloneDEoverE, eventWeightTrigup);
      m_StandalonedE_Trigdown->Fill(info.standaloneDEoverE, eventWeightTrigdown);
    }
  }
  m_AllVertexOffset->Fill(info.closestApproach, info.eventWeight);
  m_DrtoStandalone->Fill(info.staMinDr, info.eventWeight);
  m_CaloJetMinDr->Fill(info.minCaloJetDr, info.eventWeight);
  m_CaloJetEcalE->Fill(info.caloJetEcalE, info.eventWeight);
  m_CaloJetHcalE->Fill(info.caloJetHcalE, info.eventWeight);
  m_CaloJetTotalE->Fill(info.caloJetTotalE, info.eventWeight);
  m_standaloneChi2->Fill(info.staChi2, info.eventWeight);
  m_standaloneNHits->Fill(info.staNHits, info.eventWeight);
  m_staTransDCA->Fill(info.staTransDCA, info.eventWeight);
  m_staLongDCA->Fill(info.staLongDCA, info.eventWeight);

}

void Histograms::IncCutFlow(double weight) {
  m_cutProgress->Fill(cutProgress, weight);
  cutProgress++;
  return;
}

void Histograms::ResetCutFlow() { cutProgress = 0; }
