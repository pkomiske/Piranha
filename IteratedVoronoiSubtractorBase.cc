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

// make templates visible
#define PIRANHA_TEMPLATE_VISIBILITY

#include "DynamicVoronoiCylinder.hh"
#include "DynamicVoronoiDisk.hh"
#include "IteratedVoronoiSubtractorBase.hh"

BEGIN_PIRANHA_NAMESPACE

///////////////////////////////////////////////////////////////////////////////
// IteratedVoronoiSubtractorBase - implementation
///////////////////////////////////////////////////////////////////////////////

template<class DynamicVoronoi>
std::string IteratedVoronoiSubtractorBase<DynamicVoronoi>::description() const {

  std::ostringstream oss;
  oss << std::setprecision(16);
  oss << "IteratedVoronoiSubtractor - ";

  if (subtraction_type() == SubtractionType::Area)
    oss << "subtract by areas\n";
  else if (subtraction_type() == SubtractionType::EMD)
    oss << "subtract by EMDs\n";
  else if (subtraction_type() == SubtractionType::AreaTrackEMD)
    oss << "subtract by areas, track EMDs\n";

  oss << "  zpt - " << zpt() << '\n'
      << "  zemd - " << zemd() << '\n'
      << "  allow_repeats - " << (allow_repeats() ? "true" : "false") << '\n'
      << "  have_bge - " << (have_bge() ? "true" : "false") << '\n'
      << "  have_rescaling - " << (have_rescaling() ? "true" : "false") << '\n'
      << "  rho_subtraction_mode - " << (rho_subtraction_mode() == RhoSubtractionMode::Additive ? "Additive" : "Fractional") << '\n'
      << "  individual_jets_use_same_rho - " << (jet_constituents_use_same_rho() ? "true" : "false") << '\n'
      << '\n'
      << vor().description();

  if (store_history()) {
    oss << "\nSubtraction history\n";
    for (const SubtractionHistory & hist : history()) {
      if (subtraction_type() == SubtractionType::Area)
        oss << "  zpt - " << hist.zpt() << '\n';
      else if (subtraction_type() == SubtractionType::EMD)
        oss << "  zpt - " << hist.zpt() << ", zemd - " << hist.zemd() << '\n';
      else if (subtraction_type() == SubtractionType::AreaTrackEMD)
        oss << "  zpt - " << hist.zpt() << ", zemd - " << hist.zemd() << ", emdtot - " << hist.emdtot() << '\n';
    }
  }

  return oss.str();
}

template<class DynamicVoronoi>
void IteratedVoronoiSubtractorBase<DynamicVoronoi>::preprocess(const std::vector<PseudoJet> & pjs) {

  // initialize internal vector of points
  points_.resize(pjs.size());
  validmask_.resize(pjs.size());
  pts_.resize(pjs.size());

  // store points if we're going to allow repeated applications
  if (allow_repeats())
    particles_ = pjs;

  // determine if we're using individual rhos per particle
  bool local_rho(have_rescaling() && !(individual_jet_ && jet_constituents_use_same_rho_));
  if (local_rho)
    rhos_.resize(pjs.size());

  // clear stored histories
  if (store_history())
    subtraction_histories_.clear();

  // set amount of subtraction so far to 0
  rhotot_ = emdtot_ = pttot_ = 0;
  nvalid_ = nremoved_ = 0;

  for (unsigned i = 0; i < pjs.size(); i++) {

    // test if point is valid, only keep nonzero pt particles
    if (vor().valid_point(points_[nvalid_] = process_point(pjs[i]))
        && (pts_[nvalid_] = pjs[i].pt()) > 0) {

      validmask_[i] = true;
      if (local_rho)
        rhos_[nvalid_] = background_estimator()->rho(pjs[i]);
      pttot_ += pts_[nvalid_++];
    }
    else validmask_[i] = false;
  }

  // resize points in case some were removed from being valid
  points_.resize(nvalid_);
  pts_.resize(nvalid_);
  if (local_rho)
    rhos_.resize(nvalid_);

  // pass points to voronoi
  vor_(points_);
}

template<class DynamicVoronoi>
std::vector<PseudoJet>
IteratedVoronoiSubtractorBase<DynamicVoronoi>::result(const std::vector<PseudoJet> & pjs) {

  assert(pts_.size() == nvalid());

  // store subtraction history
  if (store_history_)
    subtraction_histories_.emplace_back(zpt_, zemd_, nvalid() - nremoved(), nremoved());

  // determing if we have bge set
  if (have_bge()) {
    if (subtraction_type() == SubtractionType::EMD)
      THROW_PIRANHA_ERROR("subtraction type cannot be EMD if BackgroundEstimator is set");

    // we use jet_rho_ in this case
    if (individual_jet_ && jet_constituents_use_same_rho_)
      subtract_using_areas_global_rho(rhotot_ + jet_rho_);
    else if (have_rescaling()) {
      subtract_using_areas_local_rho();
    }
    else
      subtract_using_areas_global_rho(rhotot_ + background_estimator()->rho());
  }

  // no bge set, use zs
  else {
    double rhomax(rhotot_ + zpt_ * pttot_ / total_area());
    if (subtraction_type() == SubtractionType::EMD) {
      double emdmax(emdtot_ + zemd_ * pttot_);
      subtract_using_emds_global_rho(rhomax, emdmax);
    }
    else subtract_using_areas_global_rho(rhomax);
  }

  // store emd value
  if (store_history_) {
    subtraction_histories_.back().emdtot(emd());
    subtraction_histories_.back().nremoved(nremoved() - subtraction_histories_.back().nremoved());
  }

  // get and return current particles
  std::vector<PseudoJet> newPseudoJets(points_.size());
  unsigned k(0);
  for (unsigned i = 0, n = 0; i < pjs.size(); i++) {

    // check that particle was initially valid
    if (validmask_[i]) {

      // check that particle is remaining
      if (pts_[n] > 0) {
        newPseudoJets[k] = pjs[i];
        newPseudoJets[k++].reset_momentum_PtYPhiM(pts_[n], points_[n].x(), points_[n].y(), pjs[i].m());
      }
      n++;
    }
  }

  // resize newPseudoJets to the number of particles we actually have
  newPseudoJets.resize(k);

  return newPseudoJets;
}

template<class DynamicVoronoi>
void IteratedVoronoiSubtractorBase<DynamicVoronoi>::subtract_using_areas_global_rho(double rhomax) {

  //std::cout << "Subtracting using areas and global rho\n";

  // iterate until we have removed enough pt or have removed all the particles
  bool subtract(rhotot_ < rhomax);
  while (subtract && nremoved() < nvalid()) {

    // record areas if asked to store history
    if (store_history_)
      subtraction_histories_.back().areas().push_back(vor().areas());

    // iterate over pts to determine how much we remove
    double rho(std::numeric_limits<double>::infinity());
    int iremoved(VERTEX_DEFAULT_IND);
    for (unsigned i = 0; i < nvalid(); i++) {

      // particles with zero pt are considered removed
      if (pts_[i] > 0) {
        double newrho(pts_[i]/vor().area(i));
        if (newrho < rho) {
          rho = newrho;
          iremoved = i;
        }
      }
    }

    // ensure we don't go past rhomax
    if (rhotot_ + rho >= rhomax) {
      rho = rhomax - rhotot_;
      rhotot_ = rhomax;
      iremoved = VERTEX_DEFAULT_IND;
      subtract = false;
    }
    else rhotot_ += rho;

    // compute emd of this operation
    if (subtraction_type() == SubtractionType::AreaTrackEMD) {
      double emddensitytot(0);
      for (unsigned i = 0; i < nvalid(); i++) 
        if (pts_[i] > 0) 
          emddensitytot += vor().emd_density(i);
      emdtot_ += rho * emddensitytot;
    }

    // remove rho*area from each particle's pt
    bool removed_any(do_global_subtraction(rho, iremoved));
    (void) removed_any;
    assert(removed_any || !subtract);
  }
}

template<class DynamicVoronoi>
void IteratedVoronoiSubtractorBase<DynamicVoronoi>::subtract_using_emds_global_rho(double rhomax, double emdmax) {

  //std::cout << "Subtracting using emds and global rho\n";

  // iterate until we have reached that target emd value or removed all the particles
  bool subtract(emdtot_ < emdmax && rhotot_ < rhomax);
  while (subtract && nremoved() < nvalid()) {

    // record areas if asked to store history
    if (store_history_)
      subtraction_histories_.back().areas().push_back(vor().areas());

    // find minimum pt and total the emddensity
    double emddensitytot(0), rho(std::numeric_limits<double>::infinity());
    int iremoved(VERTEX_DEFAULT_IND);
    for (unsigned i = 0; i < nvalid(); i++) {

      // particles with zero pt are considered removed
      if (pts_[i] > 0) {
        double newrho(pts_[i]/vor().area(i));
        if (newrho < rho) {
          rho = newrho;
          iremoved = i;
        }
        emddensitytot += vor().emd_density(i);
      }
    }

    // check that the calculated rho doesn't exceed rhomax
    if (rho >= rhomax - rhotot_) {

      // update rho
      rho = rhomax - rhotot_;

      // check if we exceed emdmax with the new rho
      double newemd(emdtot_ + rho * emddensitytot);
      if (newemd >= emdmax) {
        rho = (emdmax - emdtot_)/emddensitytot;
        emdtot_ = emdmax;
        rhotot_ += rho;
      }

      // emd not exceeded, update just as if rhomax exceeded
      else {
        rhotot_ = rhomax;
        emdtot_ = newemd;
      }

      // we stop subtraction no matter what at this point
      iremoved = VERTEX_DEFAULT_IND;
      subtract = false;
    }

    // check that this rho doesn't exceed emdmax
    else {
      
      // emd would be exceeded, change rho and exactly update emdtot_ and stop subtraction
      double newemd(emdtot_ + rho * emddensitytot);
      if (newemd >= emdmax) {
        rho = (emdmax - emdtot_)/emddensitytot;
        emdtot_ = emdmax;
        iremoved = VERTEX_DEFAULT_IND;
        subtract = false;
      }

      // no need to change rho
      else emdtot_ = newemd;

      // update rhotot_ with rho
      rhotot_ += rho;
    }

    // remove rho*area from each particle's pt
    bool removed_any(do_global_subtraction(rho, iremoved));
    (void) removed_any;
    assert(removed_any || !subtract);
  }
}

template<class DynamicVoronoi>
bool IteratedVoronoiSubtractorBase<DynamicVoronoi>::do_global_subtraction(double rho, int iremoved) {

  // set pt of particle iremoved exactly to zero
  std::forward_list<unsigned> removals;
  if (iremoved != VERTEX_DEFAULT_IND) {
    pts_[iremoved] = 0;
    removals.push_front(iremoved);
  }

  // remove this amount from each particle's pt, keeping track of which particles to remove
  for (unsigned i = 0; i < nvalid(); i++) {
    if (pts_[i] > 0) {
      double newpt(pts_[i] - rho * vor().area(i));
      if (newpt > 0) pts_[i] = newpt;

      // zero out pt and add it to list to be removed
      else {
        pts_[i] = 0;
        removals.push_front(i);
      }
    }

  #ifdef PIRANHA_DEBUG
    else if (pts_[i] < 0)
      std::cerr << "Negative pT detected at vertex i = " << i << std::endl;
  #endif
  }

  // remove points
  for (unsigned p : removals) {
    RemovalResult rr(vor_.remove(p));
    (void) rr;
    assert(rr == RemovalResult::Success || rr == RemovalResult::Coincidence);
    nremoved_++;
  }

  return !removals.empty();
}

template<class DynamicVoronoi>
void IteratedVoronoiSubtractorBase<DynamicVoronoi>::subtract_using_areas_local_rho() {

  // iterate until rhos are all exhausted or we remove all particles
  while (nremoved() < nvalid()) {

    // record areas if asked to store history
    if (store_history_)
      subtraction_histories_.back().areas().push_back(vor().areas());

    // iterate over pts to determine how much we remove
    double rho_fraction(1);
    int iremoved(VERTEX_DEFAULT_IND);
    bool particle_removal(false);
    std::forward_list<unsigned> removals;

    // handle different rho subtraction modes
    if (rho_subtraction_mode() == RhoSubtractionMode::Additive) {

      // determine rho to subtract additively
      double rho(std::numeric_limits<double>::infinity());
      for (unsigned i = 0; i < nvalid(); i++) {

        // ignore particles that are not under consideration
        // either they've been removed (pt <= 0) or have no more rho to subtract (rho <= 0)
        if (pts_[i] <= 0 || rhos_[i] <= 0)
          continue;

        // calculate the maximum rho this particle can take before it is removed
        double max_rho_i(pts_[i]/vor().area(i));

        // check if we may potentially remove this particle
        if (max_rho_i <= rho && max_rho_i <= rhos_[i]) {
          particle_removal = true;
          iremoved = i;
          rho = max_rho_i;
        }

        // check if we hit the rho limit for this particle
        else if (rhos_[i] < rho) {
          particle_removal = false;
          iremoved = i;
          rho = rhos_[i];
        }
      }

      // check for the end of the subtraction
      if (iremoved == VERTEX_DEFAULT_IND) break;

      // update rhos and pts
      for (unsigned i = 0; i < nvalid(); i++) {

        // particle is still active if pt > 0 and rho > 0
        if (pts_[i] > 0 && rhos_[i] > 0) {
          double newpt(pts_[i] - rho * vor().area(i));
          if (newpt > 0)
            pts_[i] = newpt;

          // zero out pt and add to list for removal
          else {
            pts_[i] = 0;
            removals.push_front(i);
          }

          // update rho to reflect what we subtracted
          rhos_[i] -= rho;
        }

      #ifdef PIRANHA_DEBUG
        if (pts_[i] < 0)
          THROW_PIRANHA_ERROR("Negative pT detected at vertex i = " + std::to_string(i));
        if (rhos_[i] < 0)
          THROW_PIRANHA_ERROR("Negative rho detected at vertex i = " + std::to_string(i));
      #endif
      }

      // track EMD
      if (track_emds()) {
        double emddensitytot(0);
        for (unsigned i = 0; i < nvalid(); i++)
          emddensitytot += vor().emd_density(i);
        emdtot_ += rho * emddensitytot;
      }
    }

    // rho decreases according to the same fraction for each particle
    else if (rho_subtraction_mode() == RhoSubtractionMode::Fractional) {

      // determine rho to remove fractionally
      for (unsigned i = 0; i < nvalid(); i++) {

        // ignore particles that are not under consideration
        // either they've been removed (pt <= 0) or have no more rho to subtract (rho <= 0)
        if (pts_[i] <= 0 || rhos_[i] <= 0)
          continue;

        // calculate the maximum fraction rho this particle can take before it is removed
        double rho_fraction_i(pts_[i]/(rhos_[i]*vor().area(i)));

        // check if we may potentially remove this particle
        if (rho_fraction_i <= rho_fraction) {
          particle_removal = true;
          iremoved = i;
          rho_fraction = rho_fraction_i;
        }
      }

      // update rhos and pts
      for (unsigned i = 0; i < nvalid(); i++) {

        // particle is still active if pt > 0 and rho > 0
        if (pts_[i] > 0 && rhos_[i] > 0) {

          double rhosub(rho_fraction * rhos_[i]), newpt(pts_[i] - rhosub * vor().area(i));
          if (newpt > 0)
            pts_[i] = newpt;

          // zero out pt and add to list for removal
          else {
            pts_[i] = 0;
            removals.push_front(i);
          }

          // update rho to reflect what we subtracted
          rhos_[i] -= rhosub;
        }

      #ifdef PIRANHA_DEBUG
        if (pts_[i] < 0)
          THROW_PIRANHA_ERROR("Negative pT detected at vertex i = " + std::to_string(i));
        if (rhos_[i] < 0)
          THROW_PIRANHA_ERROR("Negative rho detected at vertex i = " + std::to_string(i));
      #endif
      }

      // track EMD
      if (track_emds()) {
        for (unsigned i = 0; i < nvalid(); i++)
          emdtot_ += rho_fraction * rhos_[i] * vor().emd_density(i);
      }
    }

    // check for invalid rho mode
    else
      THROW_PIRANHA_ERROR("Unrecognized RhoSubtractionMode");

    // handle iremoved specially to ensure that we subtract it completely
    if (particle_removal && pts_[iremoved] > 0) {
    #ifdef PIRANHA_DEBUG
      // check that we didn't miss by too much
      if (pts_[iremoved] > piranha_epsilon)
        THROW_PIRANHA_ERROR("Numerical error with pT subtraction, " + std::to_string(pts_[iremoved]));
    #endif
      pts_[iremoved] = 0;
      removals.push_front(iremoved);
    }

    // What happens if two points are to be removed that are neighbors?
    // This is a set of measure zero except in contrived situations.

    // update rhos and remove particles
    for (unsigned removal : removals) {

      // get neighbors of point to be removed
      std::vector<int> nbs(vor().neighbors(removal));

      // get areas prior to removal
      std::vector<double> nb_areas(nbs.size());
      for (unsigned i = 0; i < nbs.size(); i++) {
        //std::cout << "neighbor " << nbs[i] << std::endl;
        if (vor().vertex_is_primary_and_active(nbs[i]))
          nb_areas[i] = vor().area(nbs[i]);
        else nbs[i] = -1;
      }

      // remove particle
      RemovalResult rr(vor_.remove(removal));
      (void) rr;
      assert(rr == RemovalResult::Success || rr == RemovalResult::Coincidence);
      nremoved_++;
      
      // update rhos after removal
      double rho_removed(rhos_[removal]);
      for (unsigned i = 0; i < nbs.size(); i++) {
        int nb(nbs[i]);
        if (nb == -1) continue;

        // new rho comes from old rho and amount from removed particle
        rhos_[nb] = (nb_areas[i] * rhos_[nb] + (vor().area(nb) - nb_areas[i]) * rho_removed) / vor().area(nb);
      
        if (rhos_[nb] < piranha_epsilon) {
        #ifdef PIRANHA_DEBUG
          std::ostringstream oss;
          oss << "Rounding rho, " << rhos_[nb] << ", detected at vertex " << nb;
        #endif
          rhos_[nb] = 0;
        }
      }
      rhos_[removal] = 0;
    }

    // check end condition for fractional rho mode
    if (rho_subtraction_mode() == RhoSubtractionMode::Fractional && rho_fraction == 1)
      break;
  }
}

// ensure that explicit templates are included in compilation unit
PIRANHA_TEMPLATE_CLASS(IteratedVoronoiSubtractorBase<DynamicVoronoiCylinder>)
PIRANHA_TEMPLATE_CLASS(IteratedVoronoiSubtractorBase<DynamicVoronoiDisk>)

END_PIRANHA_NAMESPACE