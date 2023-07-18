#include "DarkPhoton/MuAnalyzer/interface/ECAL.h"
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

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DarkPhoton/MuAnalyzer/interface/Histograms.h"
#include "Math/VectorUtil.h"
#include "TH1F.h"

#include <algorithm>
#include <iostream>

ECAL::ECAL() {
  badCellPositions.push_back(std::make_pair(-2, -2.5));
  badCellPositions.push_back(std::make_pair(1.77, -2.08));
  badCellPositions.push_back(std::make_pair(-1.55, 0.63));
  badCellPositions.push_back(std::make_pair(1.55, -1.51));
  badCellPositions.push_back(std::make_pair(1.72, 0.76));
}

double ECAL::GetIsolation(const edm::Event& iEvent,
                          const edm::EventSetup& iSetup,
                          edm::EDGetTokenT<EcalRecHitCollection> reducedEndcapRecHitCollection_Label,
                          edm::EDGetTokenT<EcalRecHitCollection> reducedBarrelRecHitCollection_Label,
                          const reco::TransientTrack track,
                          double dR) {
  TrajectoryStateClosestToPoint t0 = track.impactPointTSCP();
  for (uint8_t i = 0; i < badCellPositions.size(); i++) {
    if (deltaR(badCellPositions[i].first, badCellPositions[i].second, t0.momentum().eta(), t0.momentum().phi()) < 0.1) {
      return -1.;
    }
  }
  edm::Handle<EcalRecHitCollection> rechitsEE;
  iEvent.getByToken(reducedEndcapRecHitCollection_Label, rechitsEE);
  edm::Handle<EcalRecHitCollection> rechitsEB;
  iEvent.getByToken(reducedBarrelRecHitCollection_Label, rechitsEB);
  edm::ESHandle<CaloGeometry> TheCALOGeometry;
  iSetup.get<CaloGeometryRecord>().get(TheCALOGeometry);

  const CaloGeometry* caloGeom = TheCALOGeometry.product();
  double eDR = 0;
  for (EcalRecHitCollection::const_iterator hit = rechitsEE->begin(); hit != rechitsEE->end(); hit++) {
    const DetId id = (*hit).detid();
    const GlobalPoint hitPos = caloGeom->getSubdetectorGeometry(id)->getGeometry(id)->getPosition();
    //Check if hit and track trajectory are in the same endcap
    if ((hitPos.eta() * t0.momentum().eta()) < 0) {
      continue;
    }
    TrajectoryStateClosestToPoint traj = track.trajectoryStateClosestToPoint(hitPos);
    math::XYZVector idPositionRoot(hitPos.x(), hitPos.y(), hitPos.z());
    math::XYZVector trajRoot(traj.position().x(), traj.position().y(), traj.position().z());
    if (ROOT::Math::VectorUtil::DeltaR(idPositionRoot, trajRoot) < dR && (*hit).energy() > 0.3) {
      eDR += (*hit).energy();
    }
  }
  for (EcalRecHitCollection::const_iterator hit = rechitsEB->begin(); hit != rechitsEB->end(); hit++) {
    const DetId id = (*hit).detid();
    const GlobalPoint hitPos = caloGeom->getSubdetectorGeometry(id)->getGeometry(id)->getPosition();
    //Check if hit and track trajectory are in the same endcap
    if ((hitPos.eta() * t0.momentum().eta()) < 0) {
      continue;
    }
    TrajectoryStateClosestToPoint traj = track.trajectoryStateClosestToPoint(hitPos);
    math::XYZVector idPositionRoot(hitPos.x(), hitPos.y(), hitPos.z());
    math::XYZVector trajRoot(traj.position().x(), traj.position().y(), traj.position().z());
    if (ROOT::Math::VectorUtil::DeltaR(idPositionRoot, trajRoot) < dR && (*hit).energy() > 0.3) {
      eDR += (*hit).energy();
    }
  }
  return eDR;
}
