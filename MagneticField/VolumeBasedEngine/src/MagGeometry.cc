/*
 *  See header file for a description of this class.
 *
 *  \author N. Amapane - INFN Torino
 */

#include "MagneticField/VolumeBasedEngine/interface/MagGeometry.h"
#include "MagneticField/VolumeGeometry/interface/MagVolume.h"
#include "MagneticField/VolumeGeometry/interface/MagVolume6Faces.h"
#include "MagneticField/Layers/interface/MagBLayer.h"
#include "MagneticField/Layers/interface/MagESector.h"

#include "Utilities/BinningTools/interface/PeriodicBinFinderInPhi.h"

#include "FWCore/Utilities/interface/isFinite.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "MagneticField/Layers/interface/MagVerbosity.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>

using namespace std;
using namespace edm;

MagGeometry::MagGeometry(int geomVersion,
                         const std::vector<MagBLayer*>& tbl,
                         const std::vector<MagESector*>& tes,
                         const std::vector<MagVolume6Faces*>& tbv,
                         const std::vector<MagVolume6Faces*>& tev)
    : MagGeometry(geomVersion,
                  reinterpret_cast<std::vector<MagBLayer const*> const&>(tbl),
                  reinterpret_cast<std::vector<MagESector const*> const&>(tes),
                  reinterpret_cast<std::vector<MagVolume6Faces const*> const&>(tbv),
                  reinterpret_cast<std::vector<MagVolume6Faces const*> const&>(tev)) {}

MagGeometry::MagGeometry(int geomVersion,
                         const std::vector<MagBLayer const*>& tbl,
                         const std::vector<MagESector const*>& tes,
                         const std::vector<MagVolume6Faces const*>& tbv,
                         const std::vector<MagVolume6Faces const*>& tev)
    : theBLayers(tbl),
      theESectors(tes),
      theBVolumes(tbv),
      theEVolumes(tev),
      cacheLastVolume(true),
      geometryVersion(geomVersion) {
  vector<double> rBorders;
  for (const auto& layer : theBLayers) {
    //FIXME assume layers are already sorted in minR
    rBorders.push_back(layer->minR());
  }

  theBarrelBinFinder = new MagBinFinders::GeneralBinFinderInR<double>(rBorders);

  //FIXME assume sectors are already sorted in phi
  //FIXME: PeriodicBinFinderInPhi gets *center* of first bin
  int nEBins = theESectors.size();
  if (nEBins > 0) {
    float firstPhi = theESectors.front()->minPhi() + Geom::pi() / nEBins;
    theEndcapBinFinder = new PeriodicBinFinderInPhi<float>(firstPhi, nEBins);
  }

  if (!theEndcapBinFinder || theESectors.empty())
    edm::LogError("MagGeometry") << "Endcap empty";

  // Compute barrel dimensions based on geometry version
  switch (geomVersion >= 120812 ? 0 : (geomVersion >= 90812 ? 1 : 2)) {
    case 0:
      R1 = 172.400f;
      R2 = 308.735f;
      Z0 = 350.000f;
      Z1 = 633.290f;
      Z2 = 662.010f;
      break;
    case 1:
      R1 = 172.400f;
      R2 = 308.755f;
      Z0 = 350.000f;
      Z1 = 633.890f;
      Z2 = 662.010f;
      break;
    case 2:
      R1 = 172.400f;
      R2 = 308.755f;
      Z0 = 350.000f;
      Z1 = 633.290f;
      Z2 = 661.010f;
      break;
  }
}

MagGeometry::~MagGeometry() {
  delete theBarrelBinFinder;
  delete theEndcapBinFinder;

  for (const auto& layer : theBLayers)
    delete layer;

  for (const auto& sector : theESectors)
    delete sector;
}

// Return field vector at the specified global point
GlobalVector MagGeometry::fieldInTesla(const GlobalPoint& gp) const {
  if (const auto v = findVolume(gp))
    return v->fieldInTesla(gp);

  // Fall-back case: no volume found
  if (edm::isNotFinite(gp.mag())) {
    LogWarning("InvalidInput") << "Input value invalid (not a number): " << gp << endl;
  } else {
    LogWarning("MagneticField") << "MagGeometry::fieldInTesla: failed to find volume for " << gp << endl;
  }
  return GlobalVector();
}

// Linear search implementation (just for testing)
MagVolume const* MagGeometry::findVolume1(const GlobalPoint& gp, double tolerance) const {
  MagVolume6Faces const* found = nullptr;

  int errCnt = 0;

  float R = gp.perp();
  float Z = fabs(gp.z());

  if (inBarrel(R, Z)) {  // Barrel
    for (const auto& v : theBVolumes) {
      if (!v) {  //FIXME: remove this check
        cout << endl << "***ERROR: MagGeometry::findVolume: MagVolume for barrel not set" << endl;
        ++errCnt;
        if (errCnt < 3)
          continue;
        else
          break;
      }
      if (v->inside(gp, tolerance)) {
        found = v;
        break;
      }
    }
  } else {  // Endcaps
    for (const auto& v : theEVolumes) {
      if (!v) {  //FIXME: remove this check
        cout << endl << "***ERROR: MagGeometry::findVolume: MagVolume for endcap not set" << endl;
        ++errCnt;
        if (errCnt < 3)
          continue;
        else
          break;
      }
      if (v->inside(gp, tolerance)) {
        found = v;
        break;
      }
    }
  }

  return found;
}

// Use hierarchical structure for fast lookup.
MagVolume const* MagGeometry::findVolume(const GlobalPoint& gp, double tolerance) const {
  static thread_local MagVolume const* lastVolume = nullptr;

  MagVolume const* result = lastVolume;

  if (result && result->inside(gp))
    return result;

  float R = gp.perp();
  float Z = fabs(gp.z());

  if (inBarrel(R, Z)) {  // Barrel
    // Search up to 3 layers inwards. This may happen for very thin layers.
    int bin_hi = theBarrelBinFinder->binIndex(R);
    int bin_lo = max(0, bin_hi - 3);

    for (int bin = bin_hi; bin >= bin_lo; --bin)
      if ((result = theBLayers[bin]->findVolume(gp, tolerance)))
        break;
  } else {  // Endcaps
    Geom::Phi<float> phi = gp.phi();
    int bin = theEndcapBinFinder->binIndex(phi);
    result = theESectors[bin]->findVolume(gp, tolerance);
  }

  if (!result && tolerance < 0.0001) {
    // If search fails, retry with a 300 micron tolerance.
    // This is a hack for thin gaps on air-iron boundaries,
    // which will not be present anymore once surfaces are matched.
    result = findVolume(gp, 0.03);
  }

  if (cacheLastVolume)
    lastVolume = result;

  return result;
}

bool MagGeometry::inBarrel(float R, float Z) const {
  // FIXME: Get these dimensions from the builder.
  return ((Z < Z0) || (Z < Z1 && R > R1) || (Z < Z2 && R > R2));
}
