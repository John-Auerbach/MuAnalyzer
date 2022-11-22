#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DarkPhoton/MuAnalyzer/interface/Tracks.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Framework/interface/Event.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

Tracks::Tracks() {}

int Tracks::SelectTracks(const edm::Event& iEvent,
                         edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label,
                         edm::Handle<reco::VertexCollection> vtxHandle,
                         const reco::Muon* Tag,
                         edm::ESHandle<TransientTrackBuilder> transientTrackBuilder) {
  int cutProgress = 0;
  edm::Handle<std::vector<reco::Track>> thePATTrackHandle;
  iEvent.getByToken(trackCollection_label, thePATTrackHandle);
  highendcaptrackpt = 0;
  for (std::vector<reco::Track>::const_iterator iTrack = thePATTrackHandle->begin(); iTrack != thePATTrackHandle->end();
       ++iTrack) {
    if (deltaR(iTrack->eta(), iTrack->phi(), Tag->eta(), Tag->phi()) < 0.2)
      continue;
    if (iTrack->pt() < 20)
      continue;
    if (cutProgress == 0) {
      cutProgress = 1;
    }
    if (fabs(iTrack->eta()) > 2.4)
      continue;
    if (cutProgress == 1) {
      cutProgress = 2;
    }
    if (!iTrack->quality(reco::Track::qualityByName("highPurity")))
      continue;
    if (cutProgress == 2) {
      cutProgress = 3;
    }
    bool foundtrack = false;
    GlobalPoint tkVtx;
    for (unsigned int i = 0; i < vtxHandle->size(); i++) {
      reco::VertexRef vtx(vtxHandle, i);
      if (!vtx->isValid()) {
        continue;
      }
      for (unsigned int j = 0; j < vtx->tracksSize(); j++) {
        if (vtx->trackRefAt(j)->pt() == iTrack->pt() && vtx->trackRefAt(j)->eta() == iTrack->eta() &&
            vtx->trackRefAt(j)->phi() == iTrack->phi()) {
          foundtrack = true;
          GlobalPoint vert(vtx->x(), vtx->y(), vtx->z());
          tkVtx = vert;
        }
      }
    }
    if (!foundtrack)
      continue;
    if (cutProgress == 3) {
      cutProgress = 4;
    }
    reco::TransientTrack tk = transientTrackBuilder->build(*iTrack);
    TrajectoryStateClosestToPoint traj = tk.trajectoryStateClosestToPoint(tkVtx);
    double transDCA = traj.perigeeParameters().transverseImpactParameter();
    double longDCA = traj.perigeeParameters().longitudinalImpactParameter();
    if (longDCA > 0.05)
      continue;
    if (cutProgress == 4) {
      cutProgress = 5;
    }
    if (transDCA > 0.005)
      continue;
    if (cutProgress == 5) {
      cutProgress = 6;
    }
    selectedTracks.push_back(&(*iTrack));
    trackSelected = true;
    if (iTrack->pt() > highendcaptrackpt) {
      highptendcaptrack = &(*iTrack);
      highendcaptrackpt = iTrack->pt();
    }
    selectedEndcapTracks.push_back(&(*iTrack));
  }
  return cutProgress;
}

int Tracks::SelectUnpairedTracks(const edm::Event& iEvent,
                                 edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label,
                                 edm::Handle<reco::VertexCollection> vtxHandle) {
  int cutProgress = 0;
  edm::Handle<std::vector<reco::Track>> thePATTrackHandle;
  iEvent.getByToken(trackCollection_label, thePATTrackHandle);
  highendcaptrackpt = 0;
  for (std::vector<reco::Track>::const_iterator iTrack = thePATTrackHandle->begin(); iTrack != thePATTrackHandle->end();
       ++iTrack) {
    if (iTrack->pt() < 20)
      continue;
    if (cutProgress == 0) {
      cutProgress = 1;
    }
    selectedTracks.push_back(&(*iTrack));
    if (fabs(iTrack->eta()) > 2.4 || fabs(iTrack->eta()) < 1.45)
      continue;
    if (cutProgress == 1) {
      cutProgress = 2;
    }
    if (!iTrack->quality(reco::Track::qualityByName("highPurity")))
      continue;
    if (cutProgress == 2) {
      cutProgress = 3;
    }
    trackSelected = true;
    if (iTrack->pt() > highendcaptrackpt) {
      highptendcaptrack = &(*iTrack);
      highendcaptrackpt = iTrack->pt();
    }
    selectedEndcapTracks.push_back(&(*iTrack));
  }
  return cutProgress;
}

int Tracks::MatchMuonToTrack(const edm::Event& iEvent,
                             edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label,
                             edm::Handle<reco::VertexCollection> vtxHandle,
                             const reco::Muon* Tag,
                             edm::ESHandle<TransientTrackBuilder> transientTrackBuilder) {
  int cutProgress = 0;
  edm::Handle<std::vector<reco::Track>> thePATTrackHandle;
  iEvent.getByToken(trackCollection_label, thePATTrackHandle);
  highendcaptrackpt = 0;
  for (std::vector<reco::Track>::const_iterator iTrack = thePATTrackHandle->begin(); iTrack != thePATTrackHandle->end();
       ++iTrack) {
    if (deltaR(iTrack->eta(), iTrack->phi(), Tag->eta(), Tag->phi()) > 0.2)
      continue;
    if (iTrack->pt() < 20)
      continue;
    if (cutProgress == 0) {
      cutProgress = 1;
    }
    if (fabs(iTrack->eta()) > 2.4 || fabs(iTrack->eta()) < 1.45)
      continue;
    if (cutProgress == 1) {
      cutProgress = 2;
    }
    if (!iTrack->quality(reco::Track::qualityByName("highPurity")))
      continue;
    if (cutProgress == 2) {
      cutProgress = 3;
    }
    bool foundtrack = false;
    GlobalPoint tkVtx;
    for (unsigned int i = 0; i < vtxHandle->size(); i++) {
      reco::VertexRef vtx(vtxHandle, i);
      if (!vtx->isValid()) {
        continue;
      }
      for (unsigned int j = 0; j < vtx->tracksSize(); j++) {
        if (vtx->trackRefAt(j)->pt() == iTrack->pt() && vtx->trackRefAt(j)->eta() == iTrack->eta() &&
            vtx->trackRefAt(j)->phi() == iTrack->phi()) {
          foundtrack = true;
          GlobalPoint vert(vtx->x(), vtx->y(), vtx->z());
          tkVtx = vert;
        }
      }
    }
    if (!foundtrack)
      continue;
    if (cutProgress == 3) {
      cutProgress = 4;
    }
    reco::TransientTrack tk = transientTrackBuilder->build(*iTrack);
    TrajectoryStateClosestToPoint traj = tk.trajectoryStateClosestToPoint(tkVtx);
    double transDCA = traj.perigeeParameters().transverseImpactParameter();
    double longDCA = traj.perigeeParameters().longitudinalImpactParameter();
    if (longDCA > 0.05)
      continue;
    if (cutProgress == 4) {
      cutProgress = 5;
    }
    if (transDCA > 0.005)
      continue;
    if (cutProgress == 5) {
      cutProgress = 6;
    }
    MuonMatchedTrack = &(*iTrack);
  }
  return cutProgress;
}
double Tracks::GetIsolation(const edm::Event& iEvent,
                            edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label,
                            double conesize,
                            edm::Handle<reco::VertexCollection> vtxHandle,
                            const reco::Track* MainTrack) {
  bool foundtrack = false;
  int vtxindex = -1;
  int trackindex = -1;
  double Isolation = 0;
  for (unsigned int i = 0; i < vtxHandle->size(); i++) {
    reco::VertexRef vtx(vtxHandle, i);
    if (!vtx->isValid()) {
      continue;
    }
    for (unsigned int j = 0; j < vtx->tracksSize(); j++) {
      if (vtx->trackRefAt(j)->pt() == MainTrack->pt() && vtx->trackRefAt(j)->eta() == MainTrack->eta() &&
          vtx->trackRefAt(j)->phi() == MainTrack->phi()) {
        vtxindex = (int)i;
        trackindex = (int)j;
        foundtrack = true;
        GlobalPoint vert(vtx->x(), vtx->y(), vtx->z());
        probeVertPoint = vert;
      }
    }
  }
  probeVtx = vtxindex;
  if (!foundtrack) {
    return -1;
  }

  reco::VertexRef primaryVtx(vtxHandle, vtxindex);

  for (unsigned int i = 0; i < primaryVtx->tracksSize(); i++) {
    if ((int)i == trackindex) {
      continue;
    }
    reco::TrackBaseRef secondarytrack = primaryVtx->trackRefAt(i);
    double dphi = fabs(MainTrack->phi() - secondarytrack->phi());
    if (dphi > ROOT::Math::Pi())
      dphi -= 2 * ROOT::Math::Pi();
    double Dr = sqrt(pow(MainTrack->eta() - secondarytrack->eta(), 2.0) + pow(dphi, 2.0));
    if (Dr > conesize || Dr < 0.01) {
      continue;
    }
    Isolation += secondarytrack->pt();
  }
  return Isolation;
}

double Tracks::GetIsolation(const edm::Event& iEvent,
                            edm::EDGetTokenT<std::vector<reco::Track>> trackCollection_label,
                            double conesize,
                            edm::Handle<reco::VertexCollection> vtxHandle,
                            const reco::TrackRef MainTrack) {
  bool foundtrack = false;
  int vtxindex = -1;
  int trackindex = -1;
  double Isolation = 0;
  //double Isolation = 0;
  for (unsigned int i = 0; i < vtxHandle->size(); i++) {
    reco::VertexRef vtx(vtxHandle, i);
    if (!vtx->isValid()) {
      continue;
    }
    for (unsigned int j = 0; j < vtx->tracksSize(); j++) {
      if (vtx->trackRefAt(j)->pt() == MainTrack->pt() && vtx->trackRefAt(j)->eta() == MainTrack->eta() &&
          vtx->trackRefAt(j)->phi() == MainTrack->phi())

      // if((*vtx->trackRefAt(j))==(*MainTrack))
      {
        vtxindex = (int)i;
        trackindex = (int)j;
        foundtrack = true;
        GlobalPoint vert(vtx->x(), vtx->y(), vtx->z());
        probeVertPoint = vert;
      }
    }
  }
  tagVtx = vtxindex;
  NVertices = (int)vtxHandle->size();
  if (!foundtrack || vtxindex != 0) {
    return -1;
  }

  reco::VertexRef primaryVtx(vtxHandle, vtxindex);

  for (unsigned int i = 0; i < primaryVtx->tracksSize(); i++) {
    if ((int)i == trackindex) {
      continue;
    }
    reco::TrackBaseRef secondarytrack = primaryVtx->trackRefAt(i);
    double dphi = fabs(MainTrack->phi() - secondarytrack->phi());
    if (dphi > ROOT::Math::Pi())
      dphi -= 2 * ROOT::Math::Pi();
    double Dr = sqrt(pow(MainTrack->eta() - secondarytrack->eta(), 2.0) + pow(dphi, 2.0));
    if (Dr > conesize || Dr < 0.01) {
      continue;
    }
    Isolation += secondarytrack->pt();
  }
  return Isolation;
}

int Tracks::PairTracks(std::vector<const reco::Track*>::const_iterator& Track,
                       const reco::TrackRef MuonTrack,
                       edm::ESHandle<TransientTrackBuilder> transientTrackBuilder) {
  tracksToVertex.clear();
  tracksToVertex.push_back(transientTrackBuilder->build(*Track));
  tracksToVertex.push_back(transientTrackBuilder->build(*MuonTrack));

  // try to fit these two tracks to a common vertex
  KalmanVertexFitter vertexFitter;
  fittedVertex = vertexFitter.vertex(tracksToVertex);

  // some poor fits will simply fail to find a common vertex
  if (fittedVertex.isValid() && fittedVertex.totalChiSquared() >= 0. && fittedVertex.degreesOfFreedom() > 0) {
    // others we can exclude by their poor fit
    if (fittedVertex.totalChiSquared() / fittedVertex.degreesOfFreedom() < 3.) {
      // important! evaluate momentum vectors AT THE VERTEX
      TrajectoryStateClosestToPoint one_TSCP = tracksToVertex[0].trajectoryStateClosestToPoint(fittedVertex.position());
      TrajectoryStateClosestToPoint two_TSCP = tracksToVertex[1].trajectoryStateClosestToPoint(fittedVertex.position());
      one_momentum = one_TSCP.momentum();
      //if(sqrt(pow(one_momentum.x(),2)+pow(one_momentum.y(),2))<15){return false;}
      two_momentum = two_TSCP.momentum();

      double total_energy = sqrt(one_momentum.mag2() + 0.106 * 0.106) + sqrt(two_momentum.mag2() + 0.106 * 0.106);
      double total_px = one_momentum.x() + two_momentum.x();
      double total_py = one_momentum.y() + two_momentum.y();
      double total_pz = one_momentum.z() + two_momentum.z();
      MuonTrackMass = sqrt(pow(total_energy, 2) - pow(total_px, 2) - pow(total_py, 2) - pow(total_pz, 2));
    } else {
      return 1;
    }
    if (MuonTrackMass < 50 || MuonTrackMass > 150)
      return 2;
  } else {
    return 0;
  }
  pairvertexchi = fittedVertex.totalChiSquared() / fittedVertex.degreesOfFreedom();
  return 3;
}

int Tracks::PairTracks(const reco::Track* Track,
                       const reco::TrackRef MuonTrack,
                       edm::ESHandle<TransientTrackBuilder> transientTrackBuilder) {
  tracksToVertex.clear();
  tracksToVertex.push_back(transientTrackBuilder->build(*Track));
  tracksToVertex.push_back(transientTrackBuilder->build(*MuonTrack));

  // try to fit these two tracks to a common vertex
  KalmanVertexFitter vertexFitter;
  fittedVertex = vertexFitter.vertex(tracksToVertex);

  // some poor fits will simply fail to find a common vertex
  if (fittedVertex.isValid() && fittedVertex.totalChiSquared() >= 0. && fittedVertex.degreesOfFreedom() > 0) {
    // others we can exclude by their poor fit
    if (fittedVertex.totalChiSquared() / fittedVertex.degreesOfFreedom() < 3.) {
      // important! evaluate momentum vectors AT THE VERTEX
      TrajectoryStateClosestToPoint one_TSCP = tracksToVertex[0].trajectoryStateClosestToPoint(fittedVertex.position());
      TrajectoryStateClosestToPoint two_TSCP = tracksToVertex[1].trajectoryStateClosestToPoint(fittedVertex.position());
      one_momentum = one_TSCP.momentum();
      //if(sqrt(pow(one_momentum.x(),2)+pow(one_momentum.y(),2))<15){return false;}
      two_momentum = two_TSCP.momentum();

      double total_energy = sqrt(one_momentum.mag2() + 0.106 * 0.106) + sqrt(two_momentum.mag2() + 0.106 * 0.106);
      double total_px = one_momentum.x() + two_momentum.x();
      double total_py = one_momentum.y() + two_momentum.y();
      double total_pz = one_momentum.z() + two_momentum.z();
      MuonTrackMass = sqrt(pow(total_energy, 2) - pow(total_px, 2) - pow(total_py, 2) - pow(total_pz, 2));
    } else {
      return 1;
    }
    if (MuonTrackMass < 80 || MuonTrackMass > 100)
      return 2;
  } else {
    return 0;
  }

  return 3;
}

void Tracks::PairSimTracks(const reco::Track* Track,
                           const reco::TrackRef MuonTrack,
                           edm::ESHandle<TransientTrackBuilder> transientTrackBuilder) {
  tracksToVertex.push_back(transientTrackBuilder->build(*Track));
  tracksToVertex.push_back(transientTrackBuilder->build(*MuonTrack));
  // try to fit these two tracks to a common vertex
  KalmanVertexFitter vertexFitter;
  fittedVertex = vertexFitter.vertex(tracksToVertex);

  // some poor fits will simply fail to find a common vertex
  if (fittedVertex.isValid() && fittedVertex.totalChiSquared() >= 0. && fittedVertex.degreesOfFreedom() > 0) {
    // others we can exclude by their poor fit
    simVtxChi = fittedVertex.totalChiSquared() / fittedVertex.degreesOfFreedom();
  } else {
    simVtxChi = -1;
  }
  // important! evaluate momentum vectors AT THE VERTEX
  TrajectoryStateClosestToPoint one_TSCP = tracksToVertex[0].trajectoryStateClosestToPoint(fittedVertex.position());
  TrajectoryStateClosestToPoint two_TSCP = tracksToVertex[1].trajectoryStateClosestToPoint(fittedVertex.position());
  one_momentum = one_TSCP.momentum();
  two_momentum = two_TSCP.momentum();

  double total_energy = sqrt(one_momentum.mag2() + 0.106 * 0.106) + sqrt(two_momentum.mag2() + 0.106 * 0.106);
  double total_px = one_momentum.x() + two_momentum.x();
  double total_py = one_momentum.y() + two_momentum.y();
  double total_pz = one_momentum.z() + two_momentum.z();
  SimMuonTrackMass = sqrt(pow(total_energy, 2) - pow(total_px, 2) - pow(total_py, 2) - pow(total_pz, 2));
  SimProbeTrackPt = one_momentum.perp();
  SimProbeTrackEta = one_momentum.eta();
}

bool Tracks::PairTrackerTracks(std::vector<const reco::Track*>::const_iterator& Track,
                               std::vector<const reco::Track*>::const_iterator& Track_2nd,
                               edm::ESHandle<TransientTrackBuilder> transientTrackBuilder) {
  tracksToVertex.push_back(transientTrackBuilder->build(*Track));
  tracksToVertex.push_back(transientTrackBuilder->build(*Track_2nd));

  // try to fit these two tracks to a common vertex
  KalmanVertexFitter vertexFitter;
  fittedVertex = vertexFitter.vertex(tracksToVertex);

  // some poor fits will simply fail to find a common vertex
  if (fittedVertex.isValid() && fittedVertex.totalChiSquared() >= 0. && fittedVertex.degreesOfFreedom() > 0) {
    // others we can exclude by their poor fit
    if (fittedVertex.totalChiSquared() / fittedVertex.degreesOfFreedom() < 3.) {
      TrajectoryStateClosestToPoint one_TSCP = tracksToVertex[0].trajectoryStateClosestToPoint(fittedVertex.position());
      TrajectoryStateClosestToPoint two_TSCP = tracksToVertex[1].trajectoryStateClosestToPoint(fittedVertex.position());
      one_momentum = one_TSCP.momentum();
      two_momentum = two_TSCP.momentum();
      return true;
    } else {
      return false;
    }

  } else {
    return false;
  }
}
