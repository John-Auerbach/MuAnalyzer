#ifndef HCAL_H
#define HCAL_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HBHEChannelInfo.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/HcalRecNumberingRecord.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DarkPhoton/MuAnalyzer/interface/Histograms.h"
#include "TH1F.h"
#include "TH2F.h"

class HCAL {
public:
  HCAL();
  void CheckHCAL(
      const edm::Event& iEvent,
      const edm::EventSetup& iSetup,
      edm::EDGetTokenT<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> HBHERecHit_Label, edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometryToken);
  std::pair<double, int> EnergyInCone(
      const edm::Event& iEvent,
      const edm::EventSetup& iSetup,
      edm::EDGetTokenT<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> HBHERecHit_Label,
      edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometryToken,
      GlobalPoint track);
  double MuonMindR(
      const edm::Event& iEvent,
      const edm::EventSetup& iSetup,
      edm::EDGetTokenT<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> HBHERecHit_Label,
      edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometryToken,
      GlobalPoint MuonGlobalPoint);
  bool HitsPlots(
      const edm::Event& iEvent,
      const edm::EventSetup& iSetup,
      edm::EDGetTokenT<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> HBHERecHit_Label,
      edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometryToken,
      edm::ESGetToken<HcalTopology, HcalRecNumberingRecord> topoToken_,
      GlobalPoint TrackGlobalPoint,
      GlobalPoint RandGlobalPoint,
      bool GoodRand,
      Histograms myHistograms,
      reco::TransientTrack track);
  void SetCenterCellDistance(const edm::EventSetup& iSetup, edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometryToken, GlobalPoint TrackGlobalPoint);
  bool FindMuonHits(
      const edm::Event& iEvent,
      const edm::EventSetup& iSetup,
      edm::EDGetTokenT<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> HBHERecHit_Label,
      edm::EDGetTokenT<HORecHitCollection> horecoToken_,
      edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometryToken,
      edm::ESGetToken<HcalTopology, HcalRecNumberingRecord> topoToken_,
      GlobalPoint TrackGlobalPoint,
      double charge,
      reco::TransientTrack track);
  double GetIsolation(const edm::Event&,
                      const edm::EventSetup&,
                      edm::EDGetTokenT<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>>,
                      edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometryToken,
                      const reco::TransientTrack,
		      double coneSize=0.2);
  double HOMuonHitEnergy;
  double HOMuonHitDr;
  double MuonHitEnergy;
  int MuonHitDepth;
  double MuonMinDr;
  HcalDetId minHCALDetId;
  double m_coneEnergy;
  int m_coneHits;
  int m_bremDepth;
  int TrackiEta;
  int centerInCone;
  int centerOutOfCone;
  double Tdeta;
  double Tdphi;
  double cellDetaByDepth[7];
  double cellDphiByDepth[7];
  double cellEtaEdgeDistance[7];
  double cellPhiEdgeDistance[7];
  double BigConeEnergy;
  int m_HitsOverThresh;
  double m_hitEnergies[7];
  double m_neighborHitEnergies[7];
  double m_neighborHitEnergiesPhi[7];
  bool m_failAdjacent;
  bool m_crossDepths;
  bool foundDepths[7];
  bool onEdge;
  double etaEdgeMin;
  double phiEdgeMin;
  int CellsFound;
  int NeighborCellsFound;
  int NeighborPhiCellsFound;
  int ieta;
  int iphi;

private:
  const double Hit_Thresholds[7] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
  void GetConeIDs(const HcalTopology* theHBHETopology,
                  HcalDetId* MuonAlignedCells,
                  HcalDetId ClosestCell,
                  const int Ndepths,
                  const int CellsPerDepth);
  void GetCornerIDs(const HcalTopology* theHBHETopology,
                    HcalDetId* MuonAlignedCells,
                    HcalDetId ClosestCell,
                    const int Ndepths);
  void GetCenterCells(const HcalTopology* theHBHETopology,
                      HcalDetId* MuonAlignedCells,
                      HcalDetId ClosestCell,
                      const int Ndepths,
                      const int CellsPerDepth);
  void GetAdjacentCells(const HcalTopology* theHBHETopology,
                        HcalDetId* MuonAlignedCells,
                        HcalDetId ClosestCell,
                        const int Ndepths,
                        int ieta,
                        double deta,
                        double dphi,
                        const int CellsPerDepth);
  int GetProjectedCells(const CaloSubdetectorGeometry* HEGeom,
                        HcalDetId* TrackAlignedCells,
                        double vtxz,
                        GlobalPoint direction);
  int GetTransientProjectedCellsNeighbors(const HcalTopology* theHBHETopology,
                                          const CaloSubdetectorGeometry* HEGeom,
                                          HcalDetId* TrackAlignedCells,
                                          reco::TransientTrack muTrack,
                                          bool etaPlus);
  int GetTransientProjectedCellsNeighborsPhi(const HcalTopology* theHBHETopology,
                                             const CaloSubdetectorGeometry* HEGeom,
                                             HcalDetId* TrackAlignedCells,
                                             reco::TransientTrack muTrack,
                                             bool phiPlus);
  int GetTransientProjectedCells(const CaloSubdetectorGeometry* HEGeom,
                                 HcalDetId* TrackAlignedCells,
                                 reco::TransientTrack muTrack);
};

#endif
