#include "DarkPhoton/MuAnalyzer/interface/CSC.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DarkPhoton/MuAnalyzer/interface/Histograms.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

CSC::CSC(bool loadSFs) {
  minDR = 10;
  recMinDR = 100;
  minSimCSCDr = -1;
  nSimCSCHits = 0;
  CSCdEta = -10;
  CSCdPhi = -10;
  minTotalImpactParameter = 10000;
  for (int i = 0; i < 4; i++) {
    minDrByDepth[i] = -1;
  }
  if(loadSFs)
  {
      printf("Trying to open the CSC scale factors.\n");
      std::string sfName = "/data/cmszfs1/user/revering/dphoton/slc7/CMSSW_10_6_17_patch1/src/DarkPhoton/MuAnalyzer/histograms/missingCSCscaleFactors.root";
    TFile* lFile_SF = TFile::Open(sfName.c_str(),"READ");
    printf("Opened the file.\n");
    for (int station=0;station<4;station++)
    {
      std::string name="missingCSCSf_station"+std::to_string(station+1);
      missingCSCScaleFactors[station] = (TH2F*)lFile_SF->Get(name.c_str());
      missingCSCScaleFactors[station]->SetDirectory(0);
      printf("Filled a station thingy.\n");
    }
    lFile_SF->Close();
  }
}

void CSC::ExtrapolateTrackToDT(const edm::Event& iEvent,
                                const edm::EventSetup& iSetup,
                                edm::EDGetTokenT<DTRecSegment4DCollection> DTSegment_Label,
                                edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitToken,
                                edm::ESGetToken<DTGeometry, MuonGeometryRecord> dtGeomToken,
                                const reco::Track* iTrack,
                                reco::TransientTrack tracksToVertex,
                                GlobalPoint VertexPosition) {

  edm::Handle<DTRecSegment4DCollection> TheDTSegments = iEvent.getHandle(DTSegment_Label);
  iEvent.getByToken(DTSegment_Label, TheDTSegments);
  minTotalImpactParameter = 10000;
  int nStations = 4;
  for (int i = 0; i < nStations; i++) {
    minDtDrByDepth[i] = -1;
    minDtDEtaByDepth[i] = -1;
    minDtDPhiByDepth[i] = -1;
    minDtDZByDepth[i] = -1;
    DtHitPhiByDepth[i] = -100;
    DtHitZByDepth[i] = -100;
  }

  const DTGeometry* dtGeometry = &iSetup.getData(dtGeomToken);
  if (TheDTSegments.isValid()) {
    DTRecSegment4DCollection::id_iterator chamberIdIt;
    for (chamberIdIt = TheDTSegments->id_begin(); chamberIdIt != TheDTSegments->id_end(); ++chamberIdIt)
    {
       const DTChamber* chamber = dtGeometry->chamber(*chamberIdIt);

       DTRecSegment4DCollection::range range = TheDTSegments->get((*chamberIdIt));
       for (DTRecSegment4DCollection::const_iterator iSegment = range.first; iSegment != range.second; ++iSegment){

         GlobalPoint TheGlobalPosition = chamber->toGlobal((*iSegment).localPosition());
         TrajectoryStateClosestToPoint traj = tracksToVertex.trajectoryStateClosestToPoint(TheGlobalPosition);
         double thisMinDr =
             deltaR(traj.position().eta(), traj.position().phi(), TheGlobalPosition.eta(), TheGlobalPosition.phi());

         int stationId = chamber->id().station()-1;
         if (thisMinDr < minDtDrByDepth[stationId] || minDtDrByDepth[stationId] < 0) {
           minDtDrByDepth[stationId] = thisMinDr;
           minDtDPhiByDepth[stationId] = traj.position().phi()-TheGlobalPosition.phi();
           minDtDEtaByDepth[stationId] = traj.position().eta()-TheGlobalPosition.eta();
           minDtDZByDepth[stationId] = traj.position().z()-TheGlobalPosition.z();
           DtHitPhiByDepth[stationId] = traj.position().phi();
           DtHitZByDepth[stationId] = traj.position().z();
         }
       }
     }
  }
}

void CSC::ExtrapolateTrackToCSC(const edm::Event& iEvent,
                                const edm::EventSetup& iSetup,
                                edm::EDGetTokenT<CSCSegmentCollection> CSCSegment_Label,
                                std::vector<const reco::Track*>::const_iterator& iTrack,
                                GlobalVector one_momentum,
                                std::vector<reco::TransientTrack> tracksToVertex,
                                GlobalPoint VertexPosition) {
  edm::Handle<CSCSegmentCollection> TheCSCSegments;
  iEvent.getByToken(CSCSegment_Label, TheCSCSegments);

  edm::ESHandle<CSCGeometry> TheCSCGeometry;
  iSetup.get<MuonGeometryRecord>().get(TheCSCGeometry);
  if (TheCSCSegments.isValid()) {
    for (CSCSegmentCollection::const_iterator iSegment = TheCSCSegments->begin(); iSegment != TheCSCSegments->end();
         iSegment++) {
      CSCDetId iDetId = (CSCDetId)(*iSegment).cscDetId();
      if ((*iTrack)->eta() < 0 && iDetId.endcap() == 1)
        continue;
      if ((*iTrack)->eta() > 0 && iDetId.endcap() == 2)
        continue;

      //Only using 1st layer of CSCs
      //if(iDetId.station() != 1) continue;

      DetId TheDetUnitId(iSegment->cscDetId());
      const GeomDetUnit* TheUnit = (*TheCSCGeometry).idToDetUnit(TheDetUnitId);

      double dPhi = fabs(one_momentum.phi() - TheUnit->toGlobal(iSegment->localPosition()).phi());
      //double dPhi = fabs(one_momentum.phi() - TheUnit->position().phi());
      if (dPhi > ROOT::Math::Pi())
        dPhi -= 2 * ROOT::Math::Pi();

      LocalPoint TheLocalPosition = iSegment->localPosition();
      const BoundPlane& TheSurface = TheUnit->surface();
      GlobalPoint TheGlobalPosition = TheSurface.toGlobal(TheLocalPosition);

      if (minDR > sqrt((pow((one_momentum.eta() - TheGlobalPosition.eta()), 2.0) + pow(dPhi, 2.0)))) {
        minDR = sqrt((pow((one_momentum.eta() - TheGlobalPosition.eta()), 2.0) + pow(dPhi, 2.0)));
        //       if(minDR > sqrt(( pow((one_momentum.eta() - TheUnit->position().eta()),2.0) + pow(dPhi, 2.0)))){
        //         minDR = sqrt(( pow((one_momentum.eta() - TheUnit->position().eta()),2.0) + pow(dPhi, 2.0)));
        TrackEta_dR = one_momentum.eta();
        TrackPhi_dR = one_momentum.phi();
        CSCdEta = one_momentum.eta() - TheGlobalPosition.eta();
        CSCdPhi = one_momentum.phi() - TheGlobalPosition.phi();
        TrackGlobalPoint =
            GlobalPoint(GlobalPoint::Polar(one_momentum.theta(), one_momentum.phi(), TheGlobalPosition.mag()));
        TrackP_dR = sqrt(pow(one_momentum.x(), 2) + pow(one_momentum.y(), 2));
      }

      TrajectoryStateClosestToPoint traj = tracksToVertex[0].trajectoryStateClosestToPoint(TheGlobalPosition);

      if (minTotalImpactParameter > sqrt(traj.perigeeParameters().transverseImpactParameter() *
                                             traj.perigeeParameters().transverseImpactParameter() +
                                         traj.perigeeParameters().longitudinalImpactParameter() *
                                             traj.perigeeParameters().longitudinalImpactParameter())) {
        minTotalImpactParameter = sqrt(traj.perigeeParameters().transverseImpactParameter() *
                                           traj.perigeeParameters().transverseImpactParameter() +
                                       traj.perigeeParameters().longitudinalImpactParameter() *
                                           traj.perigeeParameters().longitudinalImpactParameter());
        TrackEta = one_momentum.eta();
        TrackPhi = one_momentum.phi();
        TrackP = sqrt(pow(one_momentum.x(), 2) + pow(one_momentum.y(), 2) + pow(one_momentum.z(), 2));
      }
    }
  }
}

void CSC::ExtrapolateTrackToCSC(const edm::Event& iEvent,
                                const edm::EventSetup& iSetup,
                                edm::EDGetTokenT<CSCSegmentCollection> CSCSegment_Label,
                                edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitToken,
                                const reco::Track* iTrack,
                                reco::TransientTrack tracksToVertex,
                                GlobalPoint VertexPosition) {
  ExtrapolateTrackToCSC(iEvent, iSetup, CSCSegment_Label, iTrack, tracksToVertex, VertexPosition);

  edm::Handle<RPCRecHitCollection> rpcRecHits;
  iEvent.getByToken(rpcRecHitToken, rpcRecHits);
  edm::ESHandle<RPCGeometry> TheRPCGeometry;
  iSetup.get<MuonGeometryRecord>().get(TheRPCGeometry);

  RPCRecHitCollection::const_iterator recIt;
  for (recIt = rpcRecHits->begin(); recIt != rpcRecHits->end(); recIt++) {
    DetId TheDetUnitId(recIt->rpcId());
    const GeomDetUnit* TheUnit = (*TheRPCGeometry).idToDetUnit(TheDetUnitId);

    LocalPoint TheLocalPosition = recIt->localPosition();
    const BoundPlane& TheSurface = TheUnit->surface();
    GlobalPoint TheGlobalPosition = TheSurface.toGlobal(TheLocalPosition);
    TrajectoryStateClosestToPoint traj = tracksToVertex.trajectoryStateClosestToPoint(TheGlobalPosition);
    double RHdEta = traj.position().eta() - TheGlobalPosition.eta();
    double RHdPhi = traj.position().phi() - TheGlobalPosition.phi();
    if (recMinDR > sqrt((pow(RHdEta, 2.0) + pow(RHdPhi, 2.0)))) {
      recMinDR = sqrt((pow(RHdEta, 2.0) + pow(RHdPhi, 2.0)));
    }
  }
}
void CSC::ExtrapolateTrackToCSC(const edm::Event& iEvent,
                                const edm::EventSetup& iSetup,
                                edm::EDGetTokenT<CSCSegmentCollection> CSCSegment_Label,
                                edm::EDGetTokenT<CSCRecHit2DCollection> CSCRecHitToken,
                                const reco::Track* iTrack,
                                reco::TransientTrack tracksToVertex,
                                GlobalPoint VertexPosition) {
  ExtrapolateTrackToCSC(iEvent, iSetup, CSCSegment_Label, iTrack, tracksToVertex, VertexPosition);

  edm::Handle<CSCRecHit2DCollection> cscRecHits;
  iEvent.getByToken(CSCRecHitToken, cscRecHits);
  edm::ESHandle<CSCGeometry> TheCSCGeometry;
  iSetup.get<MuonGeometryRecord>().get(TheCSCGeometry);

  CSCRecHit2DCollection::const_iterator recIt;
  for (recIt = cscRecHits->begin(); recIt != cscRecHits->end(); recIt++) {
    CSCDetId iDetId = (CSCDetId)(*recIt).cscDetId();
    if (iTrack->eta() < 0 && iDetId.endcap() == 1)
      continue;
    if (iTrack->eta() > 0 && iDetId.endcap() == 2)
      continue;
    DetId TheDetUnitId(recIt->cscDetId());
    const GeomDetUnit* TheUnit = (*TheCSCGeometry).idToDetUnit(TheDetUnitId);

    LocalPoint TheLocalPosition = recIt->localPosition();
    const BoundPlane& TheSurface = TheUnit->surface();
    GlobalPoint TheGlobalPosition = TheSurface.toGlobal(TheLocalPosition);
    TrajectoryStateClosestToPoint traj = tracksToVertex.trajectoryStateClosestToPoint(TheGlobalPosition);
    double RHdEta = traj.position().eta() - TheGlobalPosition.eta();
    double RHdPhi = traj.position().phi() - TheGlobalPosition.phi();
    if (recMinDR > sqrt((pow(RHdEta, 2.0) + pow(RHdPhi, 2.0)))) {
      recMinDR = sqrt((pow(RHdEta, 2.0) + pow(RHdPhi, 2.0)));
    }
  }
}

void CSC::ExtrapolateTrackToCSC(const edm::Event& iEvent,
                                const edm::EventSetup& iSetup,
                                edm::EDGetTokenT<CSCSegmentCollection> CSCSegment_Label,
                                const reco::Track* iTrack,
                                reco::TransientTrack tracksToVertex,
                                GlobalPoint VertexPosition) {
  edm::Handle<CSCSegmentCollection> TheCSCSegments;
  iEvent.getByToken(CSCSegment_Label, TheCSCSegments);
  edm::ESHandle<CSCGeometry> TheCSCGeometry;
  iSetup.get<MuonGeometryRecord>().get(TheCSCGeometry);
  minDR = 10;
  minTotalImpactParameter = 10000;
  int nStations = 4;
  for (int i = 0; i < nStations; i++) {
    minDrByDepth[i] = -1;
  }

  if (TheCSCSegments.isValid()) {
    for (CSCSegmentCollection::const_iterator iSegment = TheCSCSegments->begin(); iSegment != TheCSCSegments->end();
         iSegment++) {
      CSCDetId iDetId = (CSCDetId)(*iSegment).cscDetId();
      if (iTrack->eta() < 0 && iDetId.endcap() == 1)
        continue;
      if (iTrack->eta() > 0 && iDetId.endcap() == 2)
        continue;

      //Only using 1st layer of CSCs
      //if(iDetId.station() != 1) continue;

      DetId TheDetUnitId(iSegment->cscDetId());
      const GeomDetUnit* TheUnit = (*TheCSCGeometry).idToDetUnit(TheDetUnitId);

      LocalPoint TheLocalPosition = iSegment->localPosition();
      const BoundPlane& TheSurface = TheUnit->surface();
      GlobalPoint TheGlobalPosition = TheSurface.toGlobal(TheLocalPosition);
      TrajectoryStateClosestToPoint traj = tracksToVertex.trajectoryStateClosestToPoint(TheGlobalPosition);
      double thisMinDr =
          deltaR(traj.position().eta(), traj.position().phi(), TheGlobalPosition.eta(), TheGlobalPosition.phi());

      if (minDR > thisMinDr) {
        minDR = thisMinDr;

        CSCdEta = traj.position().eta() - TheGlobalPosition.eta();
        CSCdPhi = traj.position().phi() - TheGlobalPosition.phi();

        TrackEta_dR = iTrack->momentum().eta();
        TrackPhi_dR = iTrack->momentum().phi();
        TrackGlobalPoint = GlobalPoint(
            GlobalPoint::Polar(iTrack->momentum().theta(), iTrack->momentum().phi(), TheGlobalPosition.mag()));
        TrackP_dR = sqrt(pow(iTrack->momentum().x(), 2) + pow(iTrack->momentum().y(), 2));
        CSCTraj = TheUnit->toGlobal(iSegment->localDirection());
        CSCChiSq = iSegment->chi2();
      }
      int stationId = iDetId.station() - 1;
      if (thisMinDr < minDrByDepth[stationId] || minDrByDepth[stationId] < 0) {
        minDrByDepth[stationId] = thisMinDr;
      }

      if (minTotalImpactParameter > sqrt(traj.perigeeParameters().transverseImpactParameter() *
                                             traj.perigeeParameters().transverseImpactParameter() +
                                         traj.perigeeParameters().longitudinalImpactParameter() *
                                             traj.perigeeParameters().longitudinalImpactParameter())) {
        minTotalImpactParameter = sqrt(traj.perigeeParameters().transverseImpactParameter() *
                                           traj.perigeeParameters().transverseImpactParameter() +
                                       traj.perigeeParameters().longitudinalImpactParameter() *
                                           traj.perigeeParameters().longitudinalImpactParameter());
        TrackEta = iTrack->momentum().eta();
        TrackPhi = iTrack->momentum().phi();
        TrackP = sqrt(pow(iTrack->momentum().x(), 2) + pow(iTrack->momentum().y(), 2) + pow(iTrack->momentum().z(), 2));
      }
    }
  }
}

void CSC::ExtrapolateMuonToCSC(const edm::Event& iEvent,
                               const edm::EventSetup& iSetup,
                               edm::EDGetTokenT<CSCSegmentCollection> CSCSegment_Label,
                               const reco::Muon* iMuon,
                               GlobalVector two_momentum,
                               std::vector<reco::TransientTrack> tracksToVertex) {
  edm::Handle<CSCSegmentCollection> TheCSCSegments;
  iEvent.getByToken(CSCSegment_Label, TheCSCSegments);

  edm::ESHandle<CSCGeometry> TheCSCGeometry;
  iSetup.get<MuonGeometryRecord>().get(TheCSCGeometry);

  if (TheCSCSegments.isValid()) {
    for (CSCSegmentCollection::const_iterator iSegment = TheCSCSegments->begin(); iSegment != TheCSCSegments->end();
         iSegment++) {
      CSCDetId iDetId = (CSCDetId)(*iSegment).cscDetId();

      //Only using 1st layer of CSCs
      //       if(iDetId.station() != 1) continue;

      if (iMuon->eta() < 0 && iDetId.endcap() == 1)
        continue;
      if (iMuon->eta() > 0 && iDetId.endcap() == 2)
        continue;

      DetId TheDetUnitId(iSegment->cscDetId());
      const GeomDetUnit* TheUnit = (*TheCSCGeometry).idToDetUnit(TheDetUnitId);

      LocalPoint TheLocalPosition = iSegment->localPosition();
      //       const BoundPlane& TheSurface = TheUnit->surface();
      GlobalPoint TheGlobalPosition = TheUnit->toGlobal(TheLocalPosition);

      if (fabs(iMuon->eta()) > 1.653) {
        TrajectoryStateClosestToPoint MuonTraj = tracksToVertex[1].trajectoryStateClosestToPoint(TheGlobalPosition);
        if (minTotalImpactParameter_Muon > sqrt(MuonTraj.perigeeParameters().transverseImpactParameter() *
                                                    MuonTraj.perigeeParameters().transverseImpactParameter() +
                                                MuonTraj.perigeeParameters().longitudinalImpactParameter() *
                                                    MuonTraj.perigeeParameters().longitudinalImpactParameter())) {
          minTotalImpactParameter_Muon = sqrt(MuonTraj.perigeeParameters().transverseImpactParameter() *
                                                  MuonTraj.perigeeParameters().transverseImpactParameter() +
                                              MuonTraj.perigeeParameters().longitudinalImpactParameter() *
                                                  MuonTraj.perigeeParameters().longitudinalImpactParameter());
          MuonEta = iMuon->eta();
          MuonPhi = iMuon->phi();
          MuonP = sqrt(pow(two_momentum.x(), 2) + pow(two_momentum.y(), 2) + pow(two_momentum.z(), 2));
        }
        double dPhi_Muon = fabs(two_momentum.phi() - TheGlobalPosition.phi());
        if (dPhi_Muon > ROOT::Math::Pi())
          dPhi_Muon -= 2 * ROOT::Math::Pi();
        if (minDR_Muon > sqrt((pow((two_momentum.eta() - TheGlobalPosition.eta()), 2.0) + pow(dPhi_Muon, 2.0)))) {
          minDR_Muon = sqrt((pow((two_momentum.eta() - TheGlobalPosition.eta()), 2.0) + pow(dPhi_Muon, 2.0)));
          MuonEta_dR = iMuon->eta();
          MuonPhi_dR = iMuon->phi();
          MuonGlobalPoint =
              GlobalPoint(GlobalPoint::Polar(two_momentum.theta(), two_momentum.phi(), TheGlobalPosition.mag()));
          MuonP_dR = sqrt(pow(two_momentum.x(), 2) + pow(two_momentum.y(), 2) + pow(two_momentum.z(), 2));
        }
      }
    }
  }
}

void CSC::ExtrapolateTrackToSimCSC(const edm::Event& iEvent,
                                   const edm::EventSetup& iSetup,
                                   edm::EDGetToken cscSimHitsToken,
                                   const reco::Track* selectedTrack) {
  edm::Handle<edm::PSimHitContainer> cscSimHits;
  iEvent.getByToken(cscSimHitsToken, cscSimHits);
  edm::ESHandle<CSCGeometry> cscGeom;
  iSetup.get<MuonGeometryRecord>().get(cscGeom);
  bool goodDepths[4];
  for (int depth = 0; depth < 4; depth++) {
    goodDepths[depth] = false;
  }
  if (cscSimHits.isValid()) {
    for (unsigned int simHitCounter = 0; simHitCounter < cscSimHits->size(); simHitCounter++)
    {
      const PSimHit& simHit = (*cscSimHits)[simHitCounter];
      //if (abs(simHit.particleType()) == 13) {
      const CSCLayer* csclayer = cscGeom->layer((CSCDetId)simHit.detUnitId());
      int stationId = csclayer->id().station() - 1;
      const auto& thisPoint = csclayer->toGlobal(simHit.localPosition());
      double SimCSCDr = deltaR(thisPoint.eta(), thisPoint.phi(), selectedTrack->eta(), selectedTrack->phi());
      if (SimCSCDr < 0.05) {
        goodDepths[stationId] = true;
      }
      if (SimCSCDr < minSimCSCDr || minSimCSCDr < 0) {
        minSimCSCDr = SimCSCDr;
      }
    }
    //}
  }
  for (int station = 0; station < 4; station++) {
    if (goodDepths[station]) {
      nSimCSCHits++;
    }
  }
}

double CSC::GetMissingHitWeight(int station, double eta, double phi)
{
  if(missingCSCScaleFactors[station]->GetBinError(missingCSCScaleFactors[station]->FindBin(eta,phi))==0) return 1;
  double scaleFactor=missingCSCScaleFactors[station]->GetBinContent(missingCSCScaleFactors[station]->FindBin(eta,phi));
  if(scaleFactor>1){scaleFactor=1;}
  return scaleFactor;
}

