#include "DarkPhoton/MuAnalyzer/interface/HCAL.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/HcalRecNumberingRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DarkPhoton/MuAnalyzer/interface/Histograms.h"
#include "TH1F.h"

#include "Math/VectorUtil.h"
#include <algorithm>
#include <iostream>

HCAL::HCAL() {
  CellsFound = 0;
  centerInCone = 0;
  centerOutOfCone = 0;
  m_crossDepths = false;
  onEdge = false;
  etaEdgeMin = -1;
  phiEdgeMin = -1;
  HOMuonHitDr = -1;
  HOMuonHitEnergy = -1;
  for (int depth = 0; depth < 7; depth++) {
    cellDetaByDepth[depth] = 10.;
    cellDphiByDepth[depth] = 10.;
    cellEtaEdgeDistance[depth] = 10.;
    cellPhiEdgeDistance[depth] = 10.;
    m_hitEnergies[depth] = 0;
    m_hitDrs[depth] = 0;
    foundDepths[depth] = false;
  }
}

void HCAL::CheckHCAL(
    const edm::Event& iEvent,
    const edm::EventSetup& iSetup,
    edm::EDGetTokenT<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> HBHERecHit_Label,
    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometryToken) {
  std::cout << "Inside CheckHCAL" << std::endl;

  edm::Handle<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> hcalRecHits;
  iEvent.getByToken(HBHERecHit_Label, hcalRecHits);

  const CaloGeometry* caloGeom = &iSetup.getData(geometryToken);

  if (!hcalRecHits.isValid()) {
    //    std::cout << "Could not find HCAL RecHits" << std::endl;
  } else {
    const HBHERecHitCollection* hbhe = hcalRecHits.product();
    for (HBHERecHitCollection::const_iterator hbherechit = hbhe->begin(); hbherechit != hbhe->end(); hbherechit++) {
      HcalDetId id(hbherechit->detid());

      std::shared_ptr<const CaloCellGeometry> hbhe_cell = caloGeom->getGeometry(hbherechit->id());
      //     Global3DPoint hbhe_position = hbhe_cell->getPosition();

      /*      std::cout << "hbherechit->energy(): " << hbherechit->energy() << std::endl;
       std::cout << "hbherechit->time(): " << hbherechit->time() << std::endl;

       std::cout << "id.depth(): " << id.depth() << std::endl;
       std::cout << "id.subdet(): " << id.subdet() <<  std::endl;
       std::cout << "hbhe_position.eta(): " << hbhe_position.eta() << std::endl;
       std::cout << "hbhe_position.phi(): " << hbhe_position.phi() << std::endl;
*/
    }
  }
}

std::pair<double, int> HCAL::EnergyInCone(
    const edm::Event& iEvent,
    const edm::EventSetup& iSetup,
    edm::EDGetTokenT<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> HBHERecHit_Label,
    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometryToken,
    GlobalPoint track) {
  edm::Handle<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> hcalRecHits;
  iEvent.getByToken(HBHERecHit_Label, hcalRecHits);

  const CaloGeometry* caloGeom = &iSetup.getData(geometryToken);
  double energyInCone = 0;
  int nhits = 0;

  if (!hcalRecHits.isValid()) {
    //    std::cout << "Could not find HCAL RecHits" << std::endl;
  } else {
    const HBHERecHitCollection* hbhe = hcalRecHits.product();
    for (HBHERecHitCollection::const_iterator hbherechit = hbhe->begin(); hbherechit != hbhe->end(); hbherechit++) {
      HcalDetId id(hbherechit->detid());

      std::shared_ptr<const CaloCellGeometry> hbhe_cell = caloGeom->getGeometry(hbherechit->id());
      Global3DPoint hbhe_position = hbhe_cell->getPosition();
      double trackDr = deltaR(hbhe_position.eta(), hbhe_position.phi(), track.eta(), track.phi());
      if (trackDr < 0.4 && hbherechit->energy() > 0.) {
        energyInCone = energyInCone + hbherechit->energy();
        nhits++;
      }
    }
  }
  return std::make_pair(energyInCone, nhits);
}

double HCAL::MuonMindR(
    const edm::Event& iEvent,
    const edm::EventSetup& iSetup,
    edm::EDGetTokenT<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> HBHERecHit_Label,
    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometryToken,
    GlobalPoint MuonGlobalPoint) {
  double MuonEta = MuonGlobalPoint.eta();
  double MuonPhi = MuonGlobalPoint.phi();
  double minHCALdR = 1000;
  std::cout << "inside MuonMindR" << std::endl;

  edm::Handle<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> hcalRecHits;
  iEvent.getByToken(HBHERecHit_Label, hcalRecHits);

  const CaloGeometry* caloGeom = &iSetup.getData(geometryToken);

  if (!hcalRecHits.isValid()) {
    std::cout << "Could not find HCAL RecHits" << std::endl;
  } else {
    const HBHERecHitCollection* hbhe = hcalRecHits.product();
    for (HBHERecHitCollection::const_iterator hbherechit = hbhe->begin(); hbherechit != hbhe->end(); hbherechit++) {
      HcalDetId id(hbherechit->detid());
      std::shared_ptr<const CaloCellGeometry> hbhe_cell = caloGeom->getGeometry(hbherechit->id());
      Global3DPoint hbhe_position = hbhe_cell->getPosition();
      double dPhi = fabs(MuonPhi - hbhe_position.phi());
      if (dPhi > ROOT::Math::Pi())
        dPhi -= 2 * ROOT::Math::Pi();
      if (fabs(id.ieta()) < 19) {
        continue;
      }

      if (minHCALdR > sqrt((pow((MuonEta - hbhe_position.eta()), 2.0) + pow(dPhi, 2.0)))) {
        if (hbherechit->energy() != 0) {
          minHCALdR = sqrt((pow((MuonEta - hbhe_position.eta()), 2.0) + pow(dPhi, 2.0)));
          MuonMinDr = minHCALdR;
          MuonHitEnergy = hbherechit->energy();
          MuonHitDepth = id.depth();
          minHCALDetId = id;
        }
      }
      //       hbhe_cell->reset();
    }
  }
  return minHCALdR;
}

void HCAL::GetCenterCells(const HcalTopology* theHBHETopology,
                          HcalDetId* TrackAlignedCells,
                          HcalDetId ClosestCell,
                          const int Ndepths,
                          const int CellsPerDepth) {
  int startdepth = ClosestCell.depth();
  HcalDetId IteratingId = ClosestCell;
  for (int i = startdepth; i > 0; i--) {
    TrackAlignedCells[(i - 1) * CellsPerDepth] = IteratingId;
    theHBHETopology->decrementDepth(IteratingId);
  }
  IteratingId = ClosestCell;
  for (int i = startdepth + 1; i <= Ndepths; i++) {
    if (!theHBHETopology->validHcal(TrackAlignedCells[CellsPerDepth * (i - 2)])) {
      continue;
    }
    theHBHETopology->incrementDepth(IteratingId);
    TrackAlignedCells[(i - 1) * CellsPerDepth] = IteratingId;
  }

  return;
}

void HCAL::GetConeIDs(const HcalTopology* theHBHETopology,
                      HcalDetId* TrackAlignedCells,
                      HcalDetId ClosestCell,
                      const int Ndepths,
                      const int CellsPerDepth) {
  GetCenterCells(theHBHETopology, TrackAlignedCells, ClosestCell, Ndepths, CellsPerDepth);

  for (int i = 1; i <= Ndepths; i++) {
    if (!theHBHETopology->validHcal(TrackAlignedCells[CellsPerDepth * (i - 1)])) {
      continue;
    }
    HcalDetId incIEta[2];
    HcalDetId decIEta[2];
    theHBHETopology->incIEta(TrackAlignedCells[CellsPerDepth * (i - 1)], incIEta);
    theHBHETopology->decIEta(TrackAlignedCells[CellsPerDepth * (i - 1)], decIEta);
    HcalDetId incIPhi;
    HcalDetId decIPhi;
    if (theHBHETopology->incIPhi(TrackAlignedCells[CellsPerDepth * (i - 1)], incIPhi)) {
      TrackAlignedCells[CellsPerDepth * (i - 1) + 5] = incIPhi;
    }
    if (theHBHETopology->decIPhi(TrackAlignedCells[CellsPerDepth * (i - 1)], decIPhi)) {
      TrackAlignedCells[CellsPerDepth * (i - 1) + 6] = decIPhi;
    }
    TrackAlignedCells[CellsPerDepth * (i - 1) + 1] = incIEta[0];
    TrackAlignedCells[CellsPerDepth * (i - 1) + 2] = incIEta[1];
    TrackAlignedCells[CellsPerDepth * (i - 1) + 3] = decIEta[0];
    TrackAlignedCells[CellsPerDepth * (i - 1) + 4] = decIEta[1];
  }
  return;
}

void HCAL::GetAdjacentCells(const HcalTopology* theHBHETopology,
                            HcalDetId* TrackAlignedCells,
                            HcalDetId ClosestCell,
                            const int Ndepths,
                            int ieta,
                            double deta,
                            double dphi,
                            const int CellsPerDepth) {
  GetCenterCells(theHBHETopology, TrackAlignedCells, ClosestCell, Ndepths, CellsPerDepth);

  for (int i = 1; i <= Ndepths; i++) {
    if (!theHBHETopology->validHcal(TrackAlignedCells[CellsPerDepth * (i - 1)])) {
      continue;
    }
    HcalDetId IEta[2];
    HcalDetId Corner;
    if (deta < ROOT::Math::Pi()) {
      if (deta > 0) {
        theHBHETopology->incIEta(TrackAlignedCells[CellsPerDepth * (i - 1)], IEta);
      } else {
        theHBHETopology->decIEta(TrackAlignedCells[CellsPerDepth * (i - 1)], IEta);
      }
      if (dphi < ROOT::Math::Pi()) {
        int high = 0;
        if (theHBHETopology->validHcal(IEta[1])) {
          high = 1;
        }
        if (dphi > 0) {
          theHBHETopology->incIPhi(IEta[high], Corner);
        } else {
          theHBHETopology->decIPhi(IEta[0], Corner);
        }
      }
    }
    HcalDetId IPhi;
    if (dphi < ROOT::Math::Pi()) {
      if (dphi > 0) {
        if (theHBHETopology->incIPhi(TrackAlignedCells[CellsPerDepth * (i - 1)], IPhi)) {
          TrackAlignedCells[CellsPerDepth * (i - 1) + 3] = IPhi;
        }
      } else {
        if (theHBHETopology->decIPhi(TrackAlignedCells[CellsPerDepth * (i - 1)], IPhi)) {
          TrackAlignedCells[CellsPerDepth * (i - 1) + 3] = IPhi;
        }
      }
    }
    TrackAlignedCells[CellsPerDepth * (i - 1) + 1] = IEta[0];
    TrackAlignedCells[CellsPerDepth * (i - 1) + 2] = IEta[1];
    TrackAlignedCells[CellsPerDepth * (i - 1) + 4] = Corner;
  }
  return;
}
void HCAL::GetCornerIDs(const HcalTopology* theHBHETopology,
                        HcalDetId* CornerCells,
                        HcalDetId ClosestCell,
                        const int Ndepths) {
  const int CellsPerDepth = 8;
  int startdepth = ClosestCell.depth();
  HcalDetId IteratingId = ClosestCell;
  HcalDetId CenterIds[Ndepths];
  for (int i = startdepth; i > 0; i--) {
    CenterIds[i - 1] = IteratingId;
    theHBHETopology->decrementDepth(IteratingId);
  }
  IteratingId = ClosestCell;
  for (int i = startdepth + 1; i <= Ndepths; i++) {
    if (!theHBHETopology->validHcal(CenterIds[i - 2])) {
      continue;
    }
    theHBHETopology->incrementDepth(IteratingId);
    CenterIds[i - 1] = IteratingId;
  }
  for (int i = 1; i <= Ndepths; i++) {
    if (!theHBHETopology->validHcal(CenterIds[i - 1])) {
      continue;
    }
    HcalDetId incIPhi;
    HcalDetId decIPhi;
    theHBHETopology->incIPhi(CenterIds[i - 1], incIPhi);
    theHBHETopology->decIPhi(CenterIds[i - 1], decIPhi);
    HcalDetId incIEtaincIPhi[2];
    HcalDetId decIEtadecIPhi[2];
    HcalDetId incIEtadecIPhi[2];
    HcalDetId decIEtaincIPhi[2];
    theHBHETopology->incIEta(incIPhi, incIEtaincIPhi);
    theHBHETopology->decIEta(incIPhi, decIEtaincIPhi);
    theHBHETopology->incIEta(decIPhi, incIEtadecIPhi);
    theHBHETopology->decIEta(decIPhi, decIEtadecIPhi);
    CornerCells[CellsPerDepth * (i - 1)] = incIEtaincIPhi[0];
    CornerCells[CellsPerDepth * (i - 1) + 1] = incIEtaincIPhi[1];
    CornerCells[CellsPerDepth * (i - 1) + 2] = decIEtaincIPhi[0];
    CornerCells[CellsPerDepth * (i - 1) + 3] = decIEtaincIPhi[1];
    CornerCells[CellsPerDepth * (i - 1) + 4] = incIEtadecIPhi[0];
    CornerCells[CellsPerDepth * (i - 1) + 5] = incIEtadecIPhi[1];
    CornerCells[CellsPerDepth * (i - 1) + 6] = decIEtadecIPhi[0];
    CornerCells[CellsPerDepth * (i - 1) + 7] = decIEtadecIPhi[1];
  }
  return;
}

bool HCAL::HitsPlots(
    const edm::Event& iEvent,
    const edm::EventSetup& iSetup,
    edm::EDGetTokenT<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> HBHERecHit_Label,
    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometryToken,
    edm::ESGetToken<HcalTopology, HcalRecNumberingRecord> topoToken_,
    GlobalPoint TrackGlobalPoint,
    GlobalPoint RandGlobalPoint,
    bool GoodRand,
    Histograms myHistograms,
    reco::TransientTrack track) {
  edm::Handle<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> hcalRecHits;
  iEvent.getByToken(HBHERecHit_Label, hcalRecHits);
  static double Hits[4];
  Hits[0] = 0;
  Hits[1] = 0;
  Hits[2] = 0;
  Hits[3] = 0;

  const CaloGeometry* caloGeom = &iSetup.getData(geometryToken);
  const HcalTopology* theHBHETopology = &iSetup.getData(topoToken_);

  const CaloSubdetectorGeometry* HEGeom = caloGeom->getSubdetectorGeometry(DetId::Hcal, 2);
  int TrackiPhi;
  HcalDetId ClosestCell = (HcalDetId)HEGeom->getClosestCell(TrackGlobalPoint);
  TrackiEta = ClosestCell.ieta();
  TrackiPhi = ClosestCell.iphi();
  if (fabs(TrackiEta) < 17) {
    return false;
  }
  if (TrackiEta < -16 && TrackiPhi > 52 && TrackiPhi < 64) {
    return false;
  }
  if (fabs(TrackiEta) > 28) {
    return false;
  }
  const int Ndepths = 10;
  HcalDetId CenterCells[Ndepths];
  HcalDetId NeighborCells[Ndepths];
  HcalDetId NeighborPhiCells[Ndepths];
  Tdeta = TrackGlobalPoint.eta() - caloGeom->getGeometry(ClosestCell)->etaPos();
  Tdphi = TrackGlobalPoint.phi() - caloGeom->getGeometry(ClosestCell)->phiPos();

  CellsFound = GetTransientProjectedCells(HEGeom, CenterCells, track);
  int lastDepth = -1;
  int depthOffset = 0;
  for (int cell = 0; cell < CellsFound; cell++) {
    if (cell >= Ndepths) {
      break;
    }
    GlobalPoint cellGlobalPoint = caloGeom->getGeometry(CenterCells[cell])->getPosition();
    int thisDepth = CenterCells[cell].depth() - 1;
    if (thisDepth < 7) {
      foundDepths[thisDepth] = true;
    }
    if (thisDepth == lastDepth) {
      m_crossDepths = true;
      depthOffset++;
      continue;
    }
    lastDepth = thisDepth;
    double etaSize = caloGeom->getGeometry(CenterCells[cell])->etaSpan();
    double phiSize = caloGeom->getGeometry(CenterCells[cell])->phiSpan();
    TrajectoryStateClosestToPoint traj = track.trajectoryStateClosestToPoint(cellGlobalPoint);
    cellDetaByDepth[cell - depthOffset] = traj.position().eta() - cellGlobalPoint.eta();
    cellDphiByDepth[cell - depthOffset] = traj.position().phi() - cellGlobalPoint.phi();
    double etaEdge = etaSize / 2. - fabs(cellDetaByDepth[cell - depthOffset]);
    double phiEdge = phiSize / 2. - fabs(cellDphiByDepth[cell - depthOffset]);
    if (etaEdgeMin < 0 || etaEdge < etaEdgeMin) {
      etaEdgeMin = etaEdge;
    }
    if (phiEdgeMin < 0 || phiEdge < phiEdgeMin) {
      phiEdgeMin = phiEdge;
    }
    if (cellDetaByDepth[cell - depthOffset] > 0) {
      cellEtaEdgeDistance[cell - depthOffset] = cellDetaByDepth[cell - depthOffset] - etaSize / 2.;
    } else {
      cellEtaEdgeDistance[cell - depthOffset] = cellDetaByDepth[cell - depthOffset] + etaSize / 2.;
    }
    if (cellDphiByDepth[cell - depthOffset] > 0) {
      cellPhiEdgeDistance[cell - depthOffset] = cellDphiByDepth[cell - depthOffset] - phiSize / 2.;
    } else {
      cellPhiEdgeDistance[cell - depthOffset] = cellDphiByDepth[cell - depthOffset] + phiSize / 2.;
    }
  }
  if (fabs(etaEdgeMin) < 0.016 || fabs(phiEdgeMin) < 0.004) {
    onEdge = true;
  }
  bool etaPlus = Tdeta > 0;
  bool phiPlus = Tdphi > 0;
  NeighborCellsFound = GetTransientProjectedCellsNeighbors(theHBHETopology, HEGeom, NeighborCells, track, etaPlus);
  NeighborPhiCellsFound =
      GetTransientProjectedCellsNeighborsPhi(theHBHETopology, HEGeom, NeighborPhiCells, track, phiPlus);

  if (!hcalRecHits.isValid()) {
    std::cout << "Could not find HCAL RecHits" << std::endl;
    return false;
  } else {
    const HBHERecHitCollection* hbhe = hcalRecHits.product();
    std::deque<std::tuple<int, int, double>> MuonHits[7];
    double layerenergies[7];
    double neighborenergies[7];
    double neighborphienergies[7];
    for (int i = 0; i < 7; i++) {
      layerenergies[i] = 0;
      neighborenergies[i] = 0;
      neighborphienergies[i] = 0;
    }
    TrajectoryStateClosestToPoint t0 = track.impactPointTSCP();
    for (HBHERecHitCollection::const_iterator hbherechit = hbhe->begin(); hbherechit != hbhe->end(); hbherechit++) {
      HcalDetId id(hbherechit->detid());
      std::shared_ptr<const CaloCellGeometry> hbhe_cell = caloGeom->getGeometry(hbherechit->id());
      HcalDetId* centermatch = std::find(std::begin(CenterCells), std::end(CenterCells), id);
      HcalDetId* neighbormatch = std::find(std::begin(NeighborCells), std::end(NeighborCells), id);
      HcalDetId* neighborphimatch = std::find(std::begin(NeighborPhiCells), std::end(NeighborPhiCells), id);

      int HitiEta = id.ieta();
      if (fabs(HitiEta) < 16) {
        continue;
      }
      const GlobalPoint hitPos = hbhe_cell->getPosition();
      TrajectoryStateClosestToPoint traj = track.trajectoryStateClosestToPoint(hitPos);
      math::XYZVector idPositionRoot(hitPos.x(), hitPos.y(), hitPos.z());
      math::XYZVector trajRoot(traj.position().x(), traj.position().y(), traj.position().z());

      if (hbherechit->energy() != 0) {
        if (centermatch != std::end(CenterCells)) {
          Hits[1] += hbherechit->energy();
          Hits[0]++;
          if (id.depth() < 8) {
            layerenergies[id.depth() - 1] += hbherechit->energy();
          }
          if ((ROOT::Math::VectorUtil::DeltaR(idPositionRoot, trajRoot) < 0.2) &&
              ((idPositionRoot.eta() * t0.momentum().eta()) > 0)) {
            centerInCone++;
          } else {
            centerOutOfCone++;
          }
        }
        if (neighbormatch != std::end(NeighborCells)) {
          if (id.depth() < 8) {
            neighborenergies[id.depth() - 1] += hbherechit->energy();
          }
        }
        if (neighborphimatch != std::end(NeighborPhiCells)) {
          if (id.depth() < 8) {
            neighborphienergies[id.depth() - 1] += hbherechit->energy();
          }
        }
      }
    }
    int hitsoverthresh = 0;
    for (int i = 0; i < 7; i++) {
      m_hitEnergies[i] = layerenergies[i];
      m_neighborHitEnergies[i] = neighborenergies[i];
      m_neighborHitEnergiesPhi[i] = neighborphienergies[i];
      if (layerenergies[i] > Hit_Thresholds[i]) {
        hitsoverthresh++;
      }
    }
    m_coneHits = Hits[0];
    m_coneEnergy = Hits[1];
    m_HitsOverThresh = hitsoverthresh;
  }
  return true;
}

void HCAL::SetCenterCellDistance(const edm::EventSetup& iSetup, edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometryToken, GlobalPoint TrackGlobalPoint) {
  const CaloGeometry* caloGeom = &iSetup.getData(geometryToken);
  const CaloSubdetectorGeometry* HEGeom = caloGeom->getSubdetectorGeometry(DetId::Hcal, 2);
  HcalDetId ClosestCell = (HcalDetId)HEGeom->getClosestCell(TrackGlobalPoint);
  TrackiEta = ClosestCell.ieta();
  int TrackiPhi = ClosestCell.iphi();
  ieta = TrackiEta;
  iphi = TrackiPhi;
  Tdphi = TrackGlobalPoint.phi() - caloGeom->getGeometry(ClosestCell)->phiPos();
  Tdeta = TrackGlobalPoint.eta() - caloGeom->getGeometry(ClosestCell)->etaPos();
  if (Tdphi > ROOT::Math::Pi())
    Tdphi -= 2 * ROOT::Math::Pi();
  if (Tdphi < -ROOT::Math::Pi())
    Tdphi += 2 * ROOT::Math::Pi();
}

bool HCAL::FindMuonHits(
    const edm::Event& iEvent,
    const edm::EventSetup& iSetup,
    edm::EDGetTokenT<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> HBHERecHit_Label,
    edm::EDGetTokenT<HORecHitCollection> horecoToken_,
    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometryToken,
    edm::ESGetToken<HcalTopology, HcalRecNumberingRecord> topoToken_,
    GlobalPoint TrackGlobalPoint,
    double charge,
    reco::TransientTrack track) {
  edm::Handle<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> hcalRecHits;
  iEvent.getByToken(HBHERecHit_Label, hcalRecHits);
  double ReducedConeE = 0;
  static double Hits[4];
  Hits[0] = 0;
  Hits[1] = 0;
  Hits[2] = 0;
  Hits[3] = 0;
  m_HitsOverThresh = 0;
  const CaloGeometry* caloGeom = &iSetup.getData(geometryToken);
  const HcalTopology* theHBHETopology = &iSetup.getData(topoToken_);

  const CaloSubdetectorGeometry* HEGeom = caloGeom->getSubdetectorGeometry(DetId::Hcal, 2);
  HcalDetId ClosestCell = (HcalDetId)HEGeom->getClosestCell(TrackGlobalPoint);
  TrackiEta = ClosestCell.ieta();
  int TrackiPhi = ClosestCell.iphi();
  int BremDepth = ClosestCell.depth();
  m_bremDepth = BremDepth;
  ieta = TrackiEta;
  iphi = TrackiPhi;
  if (fabs(TrackiEta) < 17) {
    return false;
  }
  if (TrackiEta < -16 && TrackiPhi > 52 && TrackiPhi < 64) {
    return false;
  }
  if (fabs(TrackiEta) > 28) {
    return false;
  }

  const int Ndepths = 7;
  const int CellsPerDepth = 5;
  HcalDetId AdjacentCells[CellsPerDepth * Ndepths];
  HcalDetId LowThreshAdjacentCells[CellsPerDepth * Ndepths];
  HcalDetId CenterCells[CellsPerDepth * Ndepths];
  HcalDetId NeighborCells[CellsPerDepth * Ndepths];
  //GetCenterCells(theHBHETopology,TrackAlignedCells,ClosestCell,Ndepths,CellsPerDepth);
  /*  double highetathresh, lowetathresh, phithresh;
  if(fabs(TrackiEta)>20)
  {
     highetathresh = 0.015;
     lowetathresh = 0.03;
     phithresh = 0.055;
  }
  else
  {
     lowetathresh = 0.024;
     highetathresh = 0.01;
     phithresh = 0.008;
  }*/
  /*lowetathresh = 0.0;
  highetathresh = 0.0;
  phithresh = 0.0;*/
  //GetCornerIDs(theHBHETopology,CornerAlignedCells,ClosestCell,Ndepths);

  Tdphi = TrackGlobalPoint.phi() - caloGeom->getGeometry(ClosestCell)->phiPos();
  Tdeta = TrackGlobalPoint.eta() - caloGeom->getGeometry(ClosestCell)->etaPos();
  if (Tdphi > ROOT::Math::Pi())
    Tdphi -= 2 * ROOT::Math::Pi();
  if (Tdphi < -ROOT::Math::Pi())
    Tdphi += 2 * ROOT::Math::Pi();
  /*  bool etaedge = false;
  if(TrackiEta>0)
  {
     if(Tdeta>highetathresh){etaedge=true;}
     if(Tdeta<-lowetathresh){etaedge=true;}
  }
  else
  {
     if(Tdeta<-highetathresh){etaedge=true;}
     if(Tdeta>lowetathresh){etaedge=true;}
//       if(fabs(Tdeta)>lowetathresh){etaedge=true;}
  }*/
  //  bool phiedge = (fabs(Tdphi)>phithresh);
  GetAdjacentCells(
      theHBHETopology, LowThreshAdjacentCells, ClosestCell, Ndepths, TrackiEta, Tdeta, Tdphi, CellsPerDepth);
  //  if(!etaedge){Tdeta=10;}
  //  if(!phiedge){Tdphi=10;}
  GetAdjacentCells(theHBHETopology, AdjacentCells, ClosestCell, Ndepths, TrackiEta, Tdeta, Tdphi, CellsPerDepth);
  //  GetCenterCells(theHBHETopology,CenterCells,ClosestCell,Ndepths,1);
  CellsFound = GetTransientProjectedCells(HEGeom, CenterCells, track);
  int lastDepth = -1;
  int depthOffset = 0;
  for (int cell = 0; cell < CellsFound; cell++) {
    if (cell >= (CellsPerDepth * Ndepths)) {
      break;
    }
    GlobalPoint cellGlobalPoint = caloGeom->getGeometry(CenterCells[cell])->getPosition();
    int thisDepth = CenterCells[cell].depth() - 1;
    if (thisDepth < 7) {
      foundDepths[thisDepth] = true;
    }
    if (thisDepth == lastDepth) {
      m_crossDepths = true;
      depthOffset++;
      continue;
    }
    lastDepth = thisDepth;
    double etaSize = caloGeom->getGeometry(CenterCells[cell])->etaSpan();
    double phiSize = caloGeom->getGeometry(CenterCells[cell])->phiSpan();
    TrajectoryStateClosestToPoint traj = track.trajectoryStateClosestToPoint(cellGlobalPoint);
    cellDetaByDepth[cell - depthOffset] = traj.position().eta() - cellGlobalPoint.eta();
    cellDphiByDepth[cell - depthOffset] = traj.position().phi() - cellGlobalPoint.phi();
    double etaEdge = etaSize / 2. - fabs(cellDetaByDepth[cell - depthOffset]);
    double phiEdge = phiSize / 2. - fabs(cellDphiByDepth[cell - depthOffset]);
    if (etaEdgeMin < 0 || etaEdge < etaEdgeMin) {
      etaEdgeMin = etaEdge;
    }
    if (phiEdgeMin < 0 || phiEdge < phiEdgeMin) {
      phiEdgeMin = phiEdge;
    }
    if (cellDetaByDepth[cell - depthOffset] > 0) {
      cellEtaEdgeDistance[cell - depthOffset] = cellDetaByDepth[cell - depthOffset] - etaSize / 2.;
    } else {
      cellEtaEdgeDistance[cell - depthOffset] = cellDetaByDepth[cell - depthOffset] + etaSize / 2.;
    }
    if (cellDphiByDepth[cell - depthOffset] > 0) {
      cellPhiEdgeDistance[cell - depthOffset] = cellDphiByDepth[cell - depthOffset] - phiSize / 2.;
    } else {
      cellPhiEdgeDistance[cell - depthOffset] = cellDphiByDepth[cell - depthOffset] + phiSize / 2.;
    }
  }
  if (fabs(etaEdgeMin) < 0.016 || fabs(phiEdgeMin) < 0.004) {
    onEdge = true;
  }

  Tdphi = TrackGlobalPoint.phi() - caloGeom->getGeometry(ClosestCell)->phiPos();
  Tdeta = TrackGlobalPoint.eta() - caloGeom->getGeometry(ClosestCell)->etaPos();
  bool etaPlus = Tdeta > 0;
  NeighborCellsFound = GetTransientProjectedCellsNeighbors(theHBHETopology, HEGeom, NeighborCells, track, etaPlus);

  //---------------------------------------------------------------------------------------

  auto const& hoht = iEvent.getHandle(horecoToken_);
  for (HORecHitCollection::const_iterator hohtrechit = (*hoht).begin(); hohtrechit != (*hoht).end(); hohtrechit++) {
    std::shared_ptr<const CaloCellGeometry> hoht_cell = caloGeom->getGeometry(hohtrechit->id());
    Global3DPoint hoht_position = hoht_cell->getPosition();
    const GlobalPoint hitPos = hoht_cell->getPosition();
    TrajectoryStateClosestToPoint traj = track.trajectoryStateClosestToPoint(hitPos);
    double trackDr = deltaR(hoht_position.eta(), hoht_position.phi(), traj.position().eta(), traj.position().phi());
    //std::cout << "HO HCAL.cc\n";
    //std::cout << "trackDr: "+std::to_string(trackDr)+"\n";
    if ((trackDr < HOMuonHitDr || HOMuonHitDr < 0)) { // && (trackDr < 0.2)
      HOMuonHitEnergy = hohtrechit->energy();
      HOMuonHitDr = trackDr;
    }
    //std::cout << "HOMuonHitDr: "+std::to_string(HOMuonHitDr)+"\n";
    //std::cout << "HO HCAL.cc+\n";
  }


  //---------------------------------------------------------------------------------------



  if (!hcalRecHits.isValid()) {
    printf("Could not find HCAL RecHits.\n");
  } else {
    const HBHERecHitCollection* hbhe = hcalRecHits.product();
    double layerenergies[7], layerdrs[7], lowthreshadjacent[7], neighborenergies[7];
    for (int i = 0; i < 7; i++) {
      layerenergies[i] = 0;
      layerdrs[i] = 0;
      lowthreshadjacent[i] = 0;
      neighborenergies[i] = 0;
    }

    for (HBHERecHitCollection::const_iterator hbherechit = hbhe->begin(); hbherechit != hbhe->end(); hbherechit++) {
      HcalDetId id(hbherechit->detid());
      std::shared_ptr<const CaloCellGeometry> hbhe_cell = caloGeom->getGeometry(hbherechit->id());
      Global3DPoint hbhe_position = hbhe_cell->getPosition();
      const GlobalPoint hitPos = hbhe_cell->getPosition();
      TrajectoryStateClosestToPoint traj = track.trajectoryStateClosestToPoint(hitPos);
      double trackDr = deltaR(hbhe_position.eta(), hbhe_position.phi(), traj.position().eta(), traj.position().phi());

      HcalDetId* trackmatch = std::find(std::begin(AdjacentCells), std::end(AdjacentCells), id);
      HcalDetId* centermatch = std::find(std::begin(CenterCells), std::end(CenterCells), id);
      HcalDetId* neighbormatch = std::find(std::begin(NeighborCells), std::end(NeighborCells), id);
      HcalDetId* lowthreshmatch = std::find(std::begin(LowThreshAdjacentCells), std::end(LowThreshAdjacentCells), id);
      int HitiEta = id.ieta();
      if (fabs(HitiEta) < 16) {
        continue;
      }
      if (hbherechit->energy() != 0) {
        if (centermatch != std::end(CenterCells)) {
          Hits[0]++;
          Hits[1] += hbherechit->energy();
          if (hbherechit->energy() > 0.8) {
            ReducedConeE += hbherechit->energy();
          }
          if (id.depth() < 8) { // && trackDr <
            layerenergies[id.depth() - 1] += hbherechit->energy();
            layerdrs[id.depth() - 1] +=  trackDr;
          }
        }
        if (neighbormatch != std::end(NeighborCells)) {
          if (id.depth() < 8) {
            neighborenergies[id.depth() - 1] += hbherechit->energy();
          }
        }
        if (trackmatch != std::end(AdjacentCells)) {
          if (centermatch != std::end(CenterCells)) {
            if (id.depth() < 8) {
              //if(hbherechit->energy()>layerenergies[id.depth()-1]){layerenergies[id.depth()-1]=hbherechit->energy();}
            }
          }
        }
        if (lowthreshmatch != std::end(LowThreshAdjacentCells) && centermatch == std::end(CenterCells) &&
            (hbherechit->energy() > lowthreshadjacent[id.depth()])) {
          lowthreshadjacent[id.depth()] = hbherechit->energy();
        }
      }
      //       hbhe_cell->reset();
    }
    int hitsoverthresh = 0;
    m_failAdjacent = false;
    for (int i = 0; i < 7; i++) {
      m_hitEnergies[i] = layerenergies[i];
      m_hitDrs[i] = layerdrs[i];
      m_neighborHitEnergies[i] = neighborenergies[i];
      if (fabs(TrackiEta) > 25 || i < 6) {
        if (layerenergies[i] < Hit_Thresholds[i]) {
          if (lowthreshadjacent[i] > Hit_Thresholds[i]) {
            m_failAdjacent = true;
          }
        }
      }
    }
    for (int i = 0; i < 7; i++) {
      if (layerenergies[i] > Hit_Thresholds[i]) {
        hitsoverthresh++;
      }
    }
    m_coneHits = Hits[0];
    m_coneEnergy = Hits[1];
    m_HitsOverThresh = hitsoverthresh;
  }
  return true;
}

double HCAL::GetIsolation(
    const edm::Event& iEvent,
    const edm::EventSetup& iSetup,
    edm::EDGetTokenT<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> HBHERecHit_Label,
    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometryToken,
    const reco::TransientTrack track,
    double coneSize) {
  edm::Handle<edm::SortedCollection<HBHERecHit, edm::StrictWeakOrdering<HBHERecHit>>> hcalRecHits;
  iEvent.getByToken(HBHERecHit_Label, hcalRecHits);
  double MatchedEnergy = 0;
  const CaloGeometry* caloGeom = &iSetup.getData(geometryToken);
  if (!hcalRecHits.isValid())
    return 1000;
  const HBHERecHitCollection* hbhe = hcalRecHits.product();
  TrajectoryStateClosestToPoint t0 = track.impactPointTSCP();
  for (HBHERecHitCollection::const_iterator hbherechit = hbhe->begin(); hbherechit != hbhe->end(); hbherechit++) {
    HcalDetId id(hbherechit->detid());
    std::shared_ptr<const CaloCellGeometry> hbhe_cell = caloGeom->getGeometry(hbherechit->id());
    const GlobalPoint hitPos = hbhe_cell->getPosition();
    //Check if hit and track trajectory are in the same endcap
    if ((hitPos.eta() * t0.momentum().eta()) < 0) {
      continue;
    }
    TrajectoryStateClosestToPoint traj = track.trajectoryStateClosestToPoint(hitPos);
    math::XYZVector idPositionRoot(hitPos.x(), hitPos.y(), hitPos.z());
    math::XYZVector trajRoot(traj.position().x(), traj.position().y(), traj.position().z());
    if (ROOT::Math::VectorUtil::DeltaR(idPositionRoot, trajRoot) < coneSize) {
      MatchedEnergy += hbherechit->energy();
    }
  }

  return MatchedEnergy;
}

int HCAL::GetTransientProjectedCellsNeighbors(const HcalTopology* theHBHETopology,
                                              const CaloSubdetectorGeometry* HEGeom,
                                              HcalDetId* TrackAlignedCells,
                                              reco::TransientTrack muTrack,
                                              bool etaPlus) {
  double start, step, end;
  start = 320;
  step = 5;
  end = 530;
  int j = 0;
  HcalDetId lastClosestCell;
  for (int i = start; i < end; i += step) {
    double testPointPerp = fabs(i * tan(muTrack.track().theta()));
    double testPointX = testPointPerp * cos(muTrack.track().phi());
    double testPointY = testPointPerp * sin(muTrack.track().phi());
    double testPointZ;
    if (muTrack.track().eta() > 0) {
      testPointZ = i;
    } else {
      testPointZ = -i;
    }

    GlobalPoint testGlobalPoint = GlobalPoint(testPointX, testPointY, testPointZ);
    TrajectoryStateClosestToPoint traj = muTrack.trajectoryStateClosestToPoint(testGlobalPoint);
    HcalDetId testClosestCell = (HcalDetId)HEGeom->getClosestCell(traj.position());
    if (testClosestCell != lastClosestCell) {
      lastClosestCell = testClosestCell;
      HcalDetId IEtaNeighbor[2];
      if (etaPlus) {
        theHBHETopology->incIEta(testClosestCell, IEtaNeighbor);
      } else {
        theHBHETopology->decIEta(testClosestCell, IEtaNeighbor);
      }

      TrackAlignedCells[j] = IEtaNeighbor[0];
      j++;
      if (theHBHETopology->validHcal(IEtaNeighbor[1])) {
        TrackAlignedCells[j] = IEtaNeighbor[1];
        j++;
      }
    }
  }
  return j;
}

int HCAL::GetTransientProjectedCellsNeighborsPhi(const HcalTopology* theHBHETopology,
                                                 const CaloSubdetectorGeometry* HEGeom,
                                                 HcalDetId* TrackAlignedCells,
                                                 reco::TransientTrack muTrack,
                                                 bool phiPlus) {
  double start, step, end;
  start = 320;
  step = 5;
  end = 530;
  int j = 0;
  HcalDetId lastClosestCell;
  for (int i = start; i < end; i += step) {
    double testPointPerp = fabs(i * tan(muTrack.track().theta()));
    double testPointX = testPointPerp * cos(muTrack.track().phi());
    double testPointY = testPointPerp * sin(muTrack.track().phi());
    double testPointZ;
    if (muTrack.track().eta() > 0) {
      testPointZ = i;
    } else {
      testPointZ = -i;
    }

    GlobalPoint testGlobalPoint = GlobalPoint(testPointX, testPointY, testPointZ);
    TrajectoryStateClosestToPoint traj = muTrack.trajectoryStateClosestToPoint(testGlobalPoint);
    HcalDetId testClosestCell = (HcalDetId)HEGeom->getClosestCell(traj.position());
    if (testClosestCell != lastClosestCell) {
      lastClosestCell = testClosestCell;
      HcalDetId IPhiNeighbor;
      if (phiPlus) {
        theHBHETopology->incIPhi(testClosestCell, IPhiNeighbor);
      } else {
        theHBHETopology->decIPhi(testClosestCell, IPhiNeighbor);
      }

      TrackAlignedCells[j] = IPhiNeighbor;
      j++;
    }
  }
  return j;
}

int HCAL::GetTransientProjectedCells(const CaloSubdetectorGeometry* HEGeom,
                                     HcalDetId* TrackAlignedCells,
                                     reco::TransientTrack muTrack) {
  double start, step, end;
  start = 320;
  step = 5;
  end = 530;
  int j = 0;
  HcalDetId lastClosestCell;
  for (int i = start; i < end; i += step) {
    double testPointPerp = fabs(i * tan(muTrack.track().theta()));
    double testPointX = testPointPerp * cos(muTrack.track().phi());
    double testPointY = testPointPerp * sin(muTrack.track().phi());
    double testPointZ;
    if (muTrack.track().eta() > 0) {
      testPointZ = i;
    } else {
      testPointZ = -i;
    }

    GlobalPoint testGlobalPoint = GlobalPoint(testPointX, testPointY, testPointZ);
    TrajectoryStateClosestToPoint traj = muTrack.trajectoryStateClosestToPoint(testGlobalPoint);
    HcalDetId testClosestCell = (HcalDetId)HEGeom->getClosestCell(traj.position());

    if (testClosestCell != lastClosestCell) {
      lastClosestCell = testClosestCell;
      //Check that the track is geometrically within the eta/phi of the cell
      GlobalPoint cellGlobalPoint = HEGeom->getGeometry(testClosestCell)->getPosition();
      double etaSize = HEGeom->getGeometry(testClosestCell)->etaSpan();
      double phiSize = HEGeom->getGeometry(testClosestCell)->phiSpan();
      double cellDeta = traj.position().eta() - cellGlobalPoint.eta();
      double cellDphi = traj.position().phi() - cellGlobalPoint.phi();
      double etaEdge = etaSize / 2. - fabs(cellDeta);
      double phiEdge = phiSize / 2. - fabs(cellDphi);
      if (etaEdge < 0 || phiEdge < 0)
        continue;

      TrackAlignedCells[j] = testClosestCell;
      j = j + 1;
    }
  }
  return j;
}
