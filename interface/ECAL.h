#ifndef ECAL_H
#define ECAL_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

class ECAL {
public:
  ECAL();
  double GetIsolation(const edm::Event&,
                      const edm::EventSetup&,
                      edm::EDGetTokenT<EcalRecHitCollection>,
                      edm::EDGetTokenT<EcalRecHitCollection>,
                      const reco::TransientTrack,
                      double dR = 0.4);
  std::vector<std::pair<double, double>> badCellPositions;
};

#endif
