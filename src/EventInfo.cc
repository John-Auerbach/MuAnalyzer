#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DarkPhoton/MuAnalyzer/interface/EventInfo.h"
#include <iostream>

EventInfo::EventInfo() {
  eventWeight = 1;
  muRocWeight = 1;
  filledWeights = false;
  nTagMuons = -1;
  probeTrackIso = -1;
  probeEcalIso = -1;
  minGenMuDr = -1;
  minGenMuDE = -10;
  probeReducedChi = -1;
  tagMuonEta = 0;
  tagMuonPt = -1;
  probeTrackEta = 0;
  probeTrackPhi = 0;
  probeTrackPt = -1.;
  pairVtxChi = -1;
  diMuonMass = -1;
  closestApproach = -1;
  standaloneE = -1;
  staMinDr = -1;
  minCaloJetDr = -1;
  caloJetHcalE = -1;
  caloJetEcalE = -1;
  caloJetTotalE = -1;
  probeMuondE = -1;
  globalMuonE = -1;
  probeCharge = 0;
  stadEta = -10;
  stadPhi = -10;
  staNHits = 0;
  staChi2 = -1;
  staTransDCA = -1;
  staLongDCA = -1;
  nPUmean = -1.;
  pileupWeight = -1.;
  expectedHits = 0;
}

bool EventInfo::passTriggers(const edm::Event& iEvent, 
		                edm::EDGetToken m_trigResultsToken,
				edm::EDGetToken m_trigEventToken,
				std::vector<std::string> m_muonPathsToPass) 
{ 
   bool passTriggers = false;
   edm::Handle<edm::TriggerResults> triggerResults;
   edm::Handle<trigger::TriggerEvent> trigEvent;
   iEvent.getByToken(m_trigEventToken,trigEvent);
   iEvent.getByToken(m_trigResultsToken, triggerResults);
   const edm::TriggerNames& trigNames = iEvent.triggerNames(*triggerResults);
   for (size_t i = 0; i < trigNames.size(); ++i) {
      const std::string& name = trigNames.triggerName(i);
      for (auto& pathName : m_muonPathsToPass) {
         if ((name.find(pathName) != std::string::npos)) {
            if (triggerResults->accept(i)) { 
               passTriggers = true;
               trigger::size_type filterIndex = trigEvent->filterIndex(edm::InputTag("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07","","HLT"));
	       if(filterIndex<trigEvent->sizeFilters())
	       {
	          const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex);
		  const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
		  for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt)
		  {
		     const trigger::TriggerObject& obj = trigObjColl[*keyIt];
                     passingMuons.push_back(obj);
		  }
	       }
            }
         }
      }
   }

   return passTriggers;
}   

bool EventInfo::matchTagAndTrigger(const reco::Muon* selectedMuon)
{
   for(std::vector<trigger::TriggerObject>::const_iterator trigObject = passingMuons.begin(); trigObject!=passingMuons.end();++trigObject)
   {
      double dR=deltaR(trigObject->eta(),trigObject->phi(),selectedMuon->eta(),selectedMuon->phi());
      if(dR<0.1) return true;
   }
   return false;
}
