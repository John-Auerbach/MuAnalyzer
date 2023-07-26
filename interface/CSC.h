#ifndef CSC_H
#define CSC_H

#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include "Geometry/CSCGeometry/interface/CSCLayerGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCRoll.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TH2F.h"

class CSC {
public:
  CSC(bool loadSFs=false);

  void ExtrapolateTrackToCSC(const edm::Event& iEvent,
                             const edm::EventSetup& iSetup,
                             edm::EDGetTokenT<CSCSegmentCollection> CSCSegment_Label,
                             std::vector<const reco::Track*>::const_iterator& iTrack,
                             GlobalVector one_momentum,
                             std::vector<reco::TransientTrack> tracksToVertex,
                             GlobalPoint VertexPosition);
  void ExtrapolateTrackToCSC(const edm::Event& iEvent,
                             const edm::EventSetup& iSetup,
                             edm::EDGetTokenT<CSCSegmentCollection> CSCSegment_Label,
                             edm::EDGetTokenT<CSCRecHit2DCollection> CSCRecHitToken,
                             const reco::Track* iTrack,
                             reco::TransientTrack tracksToVertex,
                             GlobalPoint VertexPosition);
  void ExtrapolateTrackToCSC(const edm::Event& iEvent,
                             const edm::EventSetup& iSetup,
                             edm::EDGetTokenT<CSCSegmentCollection> CSCSegment_Label,
                             const reco::Track* iTrack,
                             reco::TransientTrack tracksToVertex,
                             GlobalPoint VertexPosition);
  void ExtrapolateTrackToCSC(const edm::Event& iEvent,
                             const edm::EventSetup& iSetup,
                             edm::EDGetTokenT<CSCSegmentCollection> CSCSegment_Label,
                             edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitToken,
                             const reco::Track* iTrack,
                             reco::TransientTrack tracksToVertex,
                             GlobalPoint VertexPosition);
  void ExtrapolateTrackToDTs(const edm::Event& iEvent,
                             const edm::EventSetup& iSetup,
                             edm::EDGetTokenT<DTRecSegment4DCollection> DTSegment_Label,
                             edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitToken,
                             const reco::Track* iTrack,
                             reco::TransientTrack tracksToVertex,
                             GlobalPoint VertexPosition);
  void ExtrapolateMuonToCSC(const edm::Event& iEvent,
                            const edm::EventSetup& iSetup,
                            edm::EDGetTokenT<CSCSegmentCollection> CSCSegment_Label,
                            const reco::Muon* iMuon,
                            GlobalVector two_momentum,
                            std::vector<reco::TransientTrack> tracksToVertex);
  void ExtrapolateTrackToSimCSC(const edm::Event& iEvent,
                                const edm::EventSetup& iSetup,
                                edm::EDGetToken cscSimHitsToken,
                                const reco::Track* selectedTrack);
 
  double GetMissingHitWeight(int station, double eta, double phi);

  GlobalPoint MuonGlobalPoint;
  GlobalPoint TrackGlobalPoint;
  GlobalPoint TrackVertex;
  GlobalVector CSCTraj;
  TH2F* missingCSCScaleFactors[4];
  bool loadSFs;
  double CSCChiSq;
  double MuonEta;
  double MuonPhi;
  double MuonP;
  double MuonEta_dR;
  double MuonPhi_dR;
  double MuonP_dR;
  double minDR;
  int nSimCSCHits;
  double CSCdEta;
  double minDrByDepth[4];
  double minDPhiByDepth[4];
  double minDEtaByDepth[4];
  double minDZByDepth[4];
  double minDtDrByDepth[4];
  double minDtDPhiByDepth[4];
  double minDtDEtaByDepth[4];
  double minDtDZByDepth[4];
  double DtHitPhiByDepth[4];
  double DtHitZByDepth[4];
  double minSimCSCDr;
  double CSCdPhi;
  double recMinDR;
  double minDR_Muon;
  double minTotalImpactParameter;
  double minTotalImpactParameter_Muon;
  double TrackP;
  double TrackEta;
  double TrackPhi;
  double TrackP_dR;
  double TrackEta_dR;
  double TrackPhi_dR;
};

#endif
