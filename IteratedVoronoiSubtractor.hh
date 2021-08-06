// Piranha Package
//
//  Questions/comments? pkomiske@mit.edu
//
//  Copyright (c) 2019-2021
//  Patrick T. Komiske III
//
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#ifndef PIRANHA_ITERATEDVORONOISUBTRACTOR_HH
#define PIRANHA_ITERATEDVORONOISUBTRACTOR_HH

// Piranha contrib
#include "DynamicVoronoiCylinder.hh"
#include "DynamicVoronoiDisk.hh"
#include "IteratedVoronoiSubtractorBase.hh"

BEGIN_PIRANHA_NAMESPACE

#ifdef DECLARE_PIRANHA_TEMPLATES
  PIRANHA_TEMPLATE_CLASS(IteratedVoronoiSubtractorBase<DynamicVoronoiCylinder>)
  PIRANHA_TEMPLATE_CLASS(IteratedVoronoiSubtractorBase<DynamicVoronoiDisk>)
#endif

///////////////////////////////////////////////////////////////////////////////
//
// IteratedVoronoiSubtractorCylinder
//   - Grooms a vertically periodic region in the plane
//   - Suitable for subtracting from the entire rapidity-azimuth plane
//   - Minimum and maximum rapidities must be provided to make region finite
//
///////////////////////////////////////////////////////////////////////////////

class IteratedVoronoiSubtractorCylinder : public IteratedVoronoiSubtractorBase<DynamicVoronoiCylinder> {
public:

  // set a maximum absolute rapidity
  IteratedVoronoiSubtractorCylinder(double abs_rap_max,
                                    SubtractionType subtype = SubtractionType::Area,
                                    double R = 1.0,
                                    bool allow_repeats = false,
                                    bool store_history = false) :
    IteratedVoronoiSubtractorCylinder(-abs_rap_max, 0, abs_rap_max, TWOPI, subtype, R, allow_repeats, store_history)
  {}

  // separate-z constructor with fully general region
  IteratedVoronoiSubtractorCylinder(double rap_min, double phi_min,
                                    double rap_max, double phi_max,
                                    SubtractionType subtype = SubtractionType::Area,
                                    double R = 1.0,
                                    bool allow_repeats = false,
                                    bool store_history = false) :
    IteratedVoronoiSubtractorBase(subtype, R, allow_repeats, store_history, rap_min, phi_min, rap_max, phi_max)
  {}

  // version taking in a jet assumed to have constituents
  PIRANHA_PSEUDOJET_CONTAINER operator()(const PseudoJet & jet) {
    individual_jet_ = true;
    if (have_bge())
      jet_rho_ = background_estimator()->rho(jet);

    std::vector<PseudoJet> consts(jet.constituents());
    preprocess(consts);
    return result(consts);
  }

  PIRANHA_PSEUDOJET_CONTAINER operator()(const std::vector<PseudoJet> & pjs) {
    individual_jet_ = false;
    preprocess(pjs);
    return result(pjs);
  }

}; // IteratedVoronoiSubtractorCylinder

///////////////////////////////////////////////////////////////////////////////
//
// IteratedVoronoiSubtractorDisk
//   - Grooms a circualr region in the plane
//   - Suitable for subtracting a conical jet, such as most anti-kT jets
//
///////////////////////////////////////////////////////////////////////////////

class IteratedVoronoiSubtractorDisk : public IteratedVoronoiSubtractorBase<DynamicVoronoiDisk> {
public:

  IteratedVoronoiSubtractorDisk(double R,
                                SubtractionType subtype = SubtractionType::Area,
                                bool allow_repeats = false,
                                bool store_history = false,
                                bool ensure_consistent_phi = true) :
    IteratedVoronoiSubtractorBase(subtype, R, allow_repeats, store_history),
    ensure_consistent_phi_(ensure_consistent_phi)
  {}

  // version taking in a jet assumed to have constituents
  PIRANHA_PSEUDOJET_CONTAINER operator()(const PseudoJet & jet) {
    individual_jet_ = true;
    if (have_bge())
      jet_rho_ = background_estimator()->rho(jet);

    vor_.set_center(jet.rap(), jet.phi());
    std::vector<PseudoJet> consts(jet.constituents());
    preprocess(consts);
    return result(consts);
  }

  // assume center is the origin
  PIRANHA_PSEUDOJET_CONTAINER operator()(const std::vector<PseudoJet> & pjs) {
    return operator()(pjs, fastjet::PtYPhiM(0, 0, 0));
  }

  // take in a non-zero center
  PIRANHA_PSEUDOJET_CONTAINER operator()(const std::vector<PseudoJet> & pjs, const PseudoJet & center) { 
    individual_jet_ = false;
    vor_.set_center(center.rap(), center.phi());
    preprocess(pjs);
    return result(pjs);
  }

private:

  // ensure that the phi fixing is done properly
  Point process_point(const PseudoJet & pj) const {
    double phi(pj.phi());
    return Point(pj.rap(), ensure_consistent_phi_ ? phi_fix(phi, vor().y0()) : phi);
  }

  bool ensure_consistent_phi_;

}; // IteratedVoronoiSubtractorDisk

END_PIRANHA_NAMESPACE

#endif // PIRANHA_ITERATEDVORONOISUBTRACTOR_HH
