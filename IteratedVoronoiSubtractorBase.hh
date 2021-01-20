// Piranha Package
//
//  Questions/comments? pkomiske@mit.edu
//
//  Copyright (c) 2019-2021
//  Patrick T. Komiske III, Eric M. Metodiev,
//  Samuel Alipour-fard, Jesse Thaler
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

#ifndef PIRANHA_ITERATEDVORONOISUBTRACTORBASE_HH
#define PIRANHA_ITERATEDVORONOISUBTRACTORBASE_HH

// C++ standard library
#include <algorithm>
#include <cassert>
#include <forward_list>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/BackgroundEstimatorBase.hh"

#include "PiranhaUtils.hh"

#define VERTEX_DEFAULT_IND -1

BEGIN_PIRANHA_NAMESPACE

///////////////////////////////////////////////////////////////////////////////
// SubtractionHistory - Stores important data about a subtraction application
///////////////////////////////////////////////////////////////////////////////

// structure to hold subtraction history
class SubtractionHistory {

  // for each point removal, the areas of each region
  std::vector<std::vector<double>> areas_;

  // the zpt and zemd of this subtraction step
  double zpt_, zemd_;

  // the total emd of the subtraction up to this point
  double emdtot_;

  // how many points we valid at the beginning of this subtraction step
  unsigned nvalid_;

  // how many points were removed in the current subtraction step
  unsigned nremoved_;

public:

  SubtractionHistory() {}
  SubtractionHistory(double zp, double ze, unsigned nv, unsigned nr) :
    zpt_(zp), zemd_(ze), nvalid_(nv), nremoved_(nr)
  {}

  // getters
  std::vector<std::vector<double>> & areas() { return areas_; }
  const std::vector<std::vector<double>> & areas() const { return areas_; }
  double zpt() const { return zpt_; }
  double zemd() const { return zemd_; }
  double emdtot() const { return emdtot_; }
  unsigned nvalid() const { return nvalid_; }
  unsigned nremoved() const { return nremoved_; }

  // setters
  void emdtot(double e) { emdtot_ = e; }
  void nremoved(unsigned nr) { nremoved_ = nr; }

}; // SubtractionHistory

///////////////////////////////////////////////////////////////////////////////
// IteratedVoronoiSubtractorBase
///////////////////////////////////////////////////////////////////////////////

template<class DynamicVoronoi>
class IteratedVoronoiSubtractorBase {
public:

  typedef typename DynamicVoronoi::Point Point;
  typedef typename DynamicVoronoi::Triangulation Triangulation;

  // this exists for SWIG
  IteratedVoronoiSubtractorBase() : bge_(nullptr) {}

  // main constructor
  template<typename... Args>
  IteratedVoronoiSubtractorBase(SubtractionType subtype, double R,
                                bool allow_repeats, bool store_history,
                                Args && ... args) :
    vor_(subtype, R, std::forward<Args>(args)...),
    zpt_(0), zemd_(0), subtype_(subtype), 
    allow_repeats_(allow_repeats), store_history_(store_history),
    total_area_(vor_.total_area()),
    emdtot_(-1),
    nvalid_(0), nremoved_(0)
  {
    set_background_estimator();
    set_jet_constituents_use_same_rho(true);
  }

  // destructor
  virtual ~IteratedVoronoiSubtractorBase() {}

  // a description of this object
  std::string description() const;

  // sets the background estimator, generally for use in pileup removal
  void set_background_estimator(BackgroundEstimatorBase * bge = nullptr) {
    bge_ = bge;
  }

  // sets whether constituents of jets use local rhos or not
  void set_jet_constituents_use_same_rho(bool b) {
    jet_constituents_use_same_rho_ = b;
  }

  // set the z values explicitly, generally for use in grooming
  // for subtype = Area, zpt = z
  // for subtype = EMD, zpt = 1, zemd = z
  void set_z(double z) {
    set_zemd(z);
    set_zpt(subtype_ == SubtractionType::EMD ? 1.0 : z);
  }
  void set_zpt(double zpt) {
    if (zpt > 1 || zpt < 0)
      throw std::invalid_argument("zpt cannot be negative or greater than one");
    zpt_ = zpt;
  }
  void set_zemd(double zemd) {
    if (zemd_ < 0)
      throw std::invalid_argument("zemd cannot be negative");
    zemd_ = zemd;
  }

  // setter functions for constructor params
  void set_subtraction_type(SubtractionType subtype) {
    subtype_ = subtype;
    vor_.set_subtraction_type(subtype);
  }
  void set_R(double R) { vor_.set_R(R); }
  void set_allow_repeats(bool allow) { allow_repeats_ = allow; }
  void set_store_history(bool store) { store_history_ = store; }

  // access stored dynamic voronoi object
  const DynamicVoronoi & vor() const { return vor_; }

  // use current settings to groom again
  std::vector<PseudoJet> reapply() {

    // check that we allow repeats, so that particles_ has been populated
    if (!allow_repeats_)
      throw std::invalid_argument("Not setup for repeated subtraction; pass true to allow_repeats.");

    return result(particles_);
  }

  // access emd value of the last subtraction
  double emd() const {
    return emdtot_;
  }

  // number of valid/removed points
  unsigned nvalid() const { return nvalid_; }
  unsigned nremoved() const { return nremoved_; }

  // access areas throughout the subtraction
  //const std::vector<SubtractionHistory> & history() const {
  const std::vector<SubtractionHistory> & history() const {
    if (!store_history_)
      THROW_PIRANHA_ERROR("subtraction history not stored; use store_history = true");
    return subtraction_histories_;
  }

  unsigned hist_size() const { return subtraction_histories_.size(); }
  unsigned hist_i_areas_size(unsigned i) const { return subtraction_histories_[i].areas().size(); }

protected:

  bool have_bge() const { return bge_ != nullptr; }
  bool have_rescaling() const { return have_bge() && bge_->rescaling_class() != nullptr; }

  // produces a point from 2D coordinates from a PseudoJet
  virtual Point process_point(const PseudoJet & pj) const {
    return Point(pj.rap(), pj.phi());
  }

  // any preprocessing that might need to be done to the PseudoJets, such as phi fixing
  // also converts pseudojets into points needed for voronoi and stores which were valid
  // stores points in voronoi object and sets pts and pttot
  void preprocess(const std::vector<PseudoJet> & pjs);

  // implement subtraction on a set of PseudoJets
  std::vector<PseudoJet> result(const std::vector<PseudoJet> & pjs);

  // instantiation of a dynamic voronoi class
  DynamicVoronoi vor_;

  // pointer to rho estimator
  BackgroundEstimatorBase * bge_;

  // track how this subtractor was called
  bool individual_jet_;

  double jet_rho_;

private:

  // rhomax is the maximum density we can remove in this step (it includes rhotot_)
  void subtract_using_areas_global_rho(double rhomax);

  // emdmax is the total amount of emd that we'd like to remove (it includes emdtot_)
  void subtract_using_emds_global_rho(double rhomax, double emdmax);

  // actually remove this much rho from every particle
  bool do_subtraction(double rho, int iremoved);

  // subtract using rhos set for each particle
  void subtract_using_areas_local_rho();

  // energy fractions available for subtraction by area or emd
  double zpt_, zemd_;

  // which type of subtraction we're doing
  SubtractionType subtype_;

  // whether or not we allow repeated applications of the Subtractor on the same set of particles
  bool allow_repeats_, store_history_, jet_constituents_use_same_rho_;

  // 2D points for each PseudoJet
  std::vector<Point> points_;

  // tracks whether a point was considered "valid" according to the DynamicVoronoi region
  std::vector<char> validmask_;

  // the pts of each PseudoJet, to be adjusted (possibly to zero) by the subtraction
  std::vector<double> pts_, rhos_;

  // this is only populated if repeats are allowed (saves redoing triangulation)
  std::vector<PseudoJet> particles_;

  // some quantities related to the current event
  double total_area_, pttot_, rhotot_, emdtot_;

  // number of PseudoJets initially valid
  unsigned nvalid_;

  // number of PseudoJets completely removed by the subtraction
  unsigned nremoved_;

  // vector of subtraction histories, in case those are being stored
  std::vector<SubtractionHistory> subtraction_histories_;

}; // IteratedVoronoiSubtractorBase

///////////////////////////////////////////////////////////////////////////////
// IteratedVoronoiSubtractorBase - implementation
///////////////////////////////////////////////////////////////////////////////

template<class T>
inline std::string
IteratedVoronoiSubtractorBase<T>::description() const {

  std::ostringstream oss;
  oss << std::setprecision(16);
  oss << "IteratedVoronoiSubtractor - ";

  if (subtype_ == SubtractionType::Area)
    oss << "subtract by areas\n";
  else if (subtype_ == SubtractionType::EMD)
    oss << "subtract by EMDs\n";
  else if (subtype_ == SubtractionType::AreaTrackEMD)
    oss << "subtract by areas, track EMDs\n";

  oss << "  zpt - " << zpt_ << '\n'
      << "  zemd - " << zemd_ << '\n'
      << "  allow_repeats - " << (allow_repeats_ ? "true" : "false") << '\n'
      << "  have_bge - " << (have_bge() ? "true" : "false") << '\n'
      << "  have_rescaling - " << (have_rescaling() ? "true" : "false") << '\n'
      << "  individual_jets_use_same_rho - " << (jet_constituents_use_same_rho_ ? "true" : "false") << '\n'
      << '\n'
      << vor_.description();

  if (store_history_) {
    oss << "\nSubtraction history\n";
    for (const SubtractionHistory & hist : subtraction_histories_) {
      if (subtype_ == SubtractionType::Area)
        oss << "  zpt - " << hist.zpt() << '\n';
      else if (subtype_ == SubtractionType::EMD)
        oss << "  zpt - " << hist.zpt() << ", zemd - " << hist.zemd() << '\n';
      else if (subtype_ == SubtractionType::AreaTrackEMD)
        oss << "  zpt - " << hist.zpt() << ", zemd - " << hist.zemd() << ", emdtot - " << hist.emdtot() << '\n';
    }
  }

  return oss.str();
}

template<class T>
inline void
IteratedVoronoiSubtractorBase<T>::preprocess(const std::vector<PseudoJet> & pjs) {

  // create vector of points to return
  points_.resize(pjs.size());
  validmask_.resize(pjs.size());
  pts_.resize(pjs.size());

  // store points if we're going to allow repeated applications
  if (allow_repeats_)
    particles_ = pjs;

  // determine if we're using individual rhos per particle
  bool local_rho(have_rescaling() && !(individual_jet_ && jet_constituents_use_same_rho_));
  if (local_rho)
    rhos_.resize(pjs.size());

  // clear stored histories
  if (store_history_)
    subtraction_histories_.clear();

  // set amount of subtraction so far to 0
  rhotot_ = emdtot_ = pttot_ = 0;
  nvalid_ = nremoved_ = 0;

  for (unsigned i = 0; i < pjs.size(); i++) {

    // test if point is valid, only keep nonzero pt particles
    if (vor_.valid_point(points_[nvalid_] = process_point(pjs[i]))
        && (pts_[nvalid_] = pjs[i].pt()) > 0) {

      validmask_[i] = true;
      if (local_rho)
        rhos_[nvalid_] = bge_->rho(pjs[i]);
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

template<class T>
inline std::vector<PseudoJet>
IteratedVoronoiSubtractorBase<T>::result(const std::vector<PseudoJet> & pjs) {

  assert(pts_.size() == nvalid());

  // store subtraction history
  if (store_history_)
    subtraction_histories_.emplace_back(zpt_, zemd_, nvalid() - nremoved(), nremoved());

  // determing if we have bge set
  if (have_bge()) {
    if (subtype_ == SubtractionType::EMD)
      THROW_PIRANHA_ERROR("subtraction type cannot be EMD if BackgroundEstimator is set");

    // we use jet_rho_ in this case
    if (individual_jet_ && jet_constituents_use_same_rho_)
      subtract_using_areas_global_rho(rhotot_ + jet_rho_);
    else if (have_rescaling()) {
      subtract_using_areas_local_rho();
    }
    else
      subtract_using_areas_global_rho(rhotot_ + bge_->rho());
  }

  // no bge set, use zs
  else {
    double rhomax(rhotot_ + zpt_ * pttot_ / total_area_);
    if (subtype_ == SubtractionType::EMD) {
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

template<class T>
inline void
IteratedVoronoiSubtractorBase<T>::subtract_using_areas_global_rho(double rhomax) {

  //std::cout << "Subtracting using areas and global rho\n";

  // iterate until we have removed enough pt or have removed all the particles
  bool subtract(rhotot_ < rhomax);
  while (subtract && nremoved() < nvalid()) {

    // record areas if asked to store history
    if (store_history_)
      subtraction_histories_.back().areas().push_back(vor_.areas());

    // iterate over pts to determine how much we remove
    double rho(std::numeric_limits<double>::infinity());
    int iremoved(VERTEX_DEFAULT_IND);
    for (unsigned i = 0; i < nvalid(); i++) {

      // particles with zero pt are considered removed
      if (pts_[i] > 0) {
        double newrho(pts_[i]/vor_.area(i));
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
    if (subtype_ == SubtractionType::AreaTrackEMD) {
      double emddensitytot(0);
      for (unsigned i = 0; i < nvalid(); i++) 
        if (pts_[i] > 0) 
          emddensitytot += vor_.emd_density(i);
      emdtot_ += rho * emddensitytot;
    }

    // remove rho*area from each particle's pt
    bool removed_any(do_subtraction(rho, iremoved));
    (void) removed_any;
    assert(removed_any || !subtract);
  }
}

template<class T>
inline void
IteratedVoronoiSubtractorBase<T>::subtract_using_emds_global_rho(double rhomax, double emdmax) {

  //std::cout << "Subtracting using emds and global rho\n";

  // iterate until we have reached that target emd value or removed all the particles
  bool subtract(emdtot_ < emdmax && rhotot_ < rhomax);
  while (subtract && nremoved() < nvalid()) {

    // record areas if asked to store history
    if (store_history_)
      subtraction_histories_.back().areas().push_back(vor_.areas());

    // find minimum pt and total the emddensity
    double emddensitytot(0), rho(std::numeric_limits<double>::infinity());
    int iremoved(VERTEX_DEFAULT_IND);
    for (unsigned i = 0; i < nvalid(); i++) {

      // particles with zero pt are considered removed
      if (pts_[i] > 0) {
        double newrho(pts_[i]/vor_.area(i));
        if (newrho < rho) {
          rho = newrho;
          iremoved = i;
        }
        emddensitytot += vor_.emd_density(i);
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
    bool removed_any(do_subtraction(rho, iremoved));
    (void) removed_any;
    assert(removed_any || !subtract);
  }
}

template<class T>
inline bool
IteratedVoronoiSubtractorBase<T>::do_subtraction(double rho, int iremoved) {

  // set pt of particle iremoved exactly to zero
  std::forward_list<unsigned> removals;
  if (iremoved != VERTEX_DEFAULT_IND) {
    pts_[iremoved] = 0;
    removals.push_front(iremoved);
  }

  // remove this amount from each particle's pt, keeping track of which particles to remove
  for (unsigned i = 0; i < nvalid(); i++) {
    if (pts_[i] > 0) {
      double newpt(pts_[i] - rho * vor_.area(i));
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

template<class T>
inline void
IteratedVoronoiSubtractorBase<T>::subtract_using_areas_local_rho() {

  //std::cout << "Subtracting using areas and local rho" << std::endl;

  // iterate until rhos are all exhausted or we remove all particles
  while (nremoved() < nvalid()) {

    // record areas if asked to store history
    if (store_history_)
      subtraction_histories_.back().areas().push_back(vor_.areas());

    // iterate over pts to determine how much we remove
    double rho(std::numeric_limits<double>::infinity());
    int iremoved(VERTEX_DEFAULT_IND);
    bool particle_removal(false);
    for (unsigned i = 0; i < nvalid(); i++) {

      // ignore particles that are not under consideration
      // either they've been removed (pt <= 0) or have no more rho to subtract (rho <= 0)
      if (pts_[i] <= 0 || rhos_[i] <= 0)
        continue;

      // calculate the maximum rho this particle can take before it is removed
      double particle_maxrho(pts_[i]/vor_.area(i));

      // check if we may potentially remove this particle
      if (particle_maxrho <= rho && particle_maxrho <= rhos_[i]) {
        particle_removal = true;
        iremoved = i;
        rho = particle_maxrho;
      }

      // check if we hit the rho limit for this particle
      else if (rhos_[i] < rho) {
        particle_removal = false;
        iremoved = i;
        rho = rhos_[i];
      }
    }

    // check for the end of the subtraction
    if (iremoved == VERTEX_DEFAULT_IND)
      break;

    //std::cout << "selected particle " << iremoved << " with rho " << rho << std::endl;

    // iterate through particles to check which ones are to be removed
    std::forward_list<unsigned> removals;
    double emddensitytot(0);
    for (unsigned i = 0; i < nvalid(); i++) {

      // particle is still active if pt > 0 and rho > 0
      if (pts_[i] > 0 && rhos_[i] > 0) {
        double newpt(pts_[i] - rho * vor_.area(i));
        if (newpt > 0)
          pts_[i] = newpt;

        // zero out pt and add to list for removal
        else {
          //std::cout << "particle " << i << " is being removed" << std::endl;
          pts_[i] = 0;
          removals.push_front(i);
        }

        // update rho to reflect what we subtracted
        rhos_[i] -= rho;

        // track emd of this operation
        if (subtype_ == SubtractionType::AreaTrackEMD)
          emddensitytot += vor_.emd_density(i);
      }

    #ifdef PIRANHA_DEBUG
      if (pts_[i] < 0)
        THROW_PIRANHA_ERROR("Negative pT detected at vertex i = " + std::to_string(i));
      if (rhos_[i] < 0)
        THROW_PIRANHA_ERROR("Negative rho detected at vertex i = " + std::to_string(i));
    #endif
    }

    // track EMD
    if (subtype_ == SubtractionType::AreaTrackEMD)
      emdtot_ += rho * emddensitytot;

    // handle iremoved specially to ensure that we subtract it completely
    if (particle_removal && pts_[iremoved] > 0) {
    #ifdef PIRANHA_DEBUG
      // check that we didn't miss by too much
      if (pts_[iremoved] > 10*std::numeric_limits<double>::epsilon())
        THROW_PIRANHA_ERROR("Numerical error with pT subtraction, " + std::to_string(pts_[iremoved]));
    #endif
      //std::cout << "particle " << iremoved << " is being removed" << std::endl;
      pts_[iremoved] = 0;
      removals.push_front(iremoved);
    }

    // update rhos and remove particles
    for (unsigned removal : removals) {

      //std::cout << "just before getting neighbors of " << removal << std::endl;

      // get neighbors of point to be removed
      std::vector<int> nbs(vor_.neighbors(removal));

      //std::cout << "just after getting neighbors of " << removal << std::endl;

      // get areas prior to removal
      std::vector<double> nb_areas(nbs.size());
      for (unsigned i = 0; i < nbs.size(); i++) {
        //std::cout << "neighbor " << nbs[i] << std::endl;
        if (vor_.vertex_is_primary_and_active(nbs[i]))
          nb_areas[i] = vor_.area(nbs[i]);
        else nbs[i] = -1;
      }

      //std::cout << "just before removing particle" << std::endl;

      // remove particle
      RemovalResult rr(vor_.remove(removal));
      (void) rr;
      assert(rr == RemovalResult::Success || rr == RemovalResult::Coincidence);
      nremoved_++;

      //std::cout << "just after removing particle" << std::endl;
      
      // update rhos after removal
      double rho_removed(rhos_[removal]);
      for (unsigned i = 0; i < nbs.size(); i++) {
        int nb(nbs[i]);
        if (nb == -1) continue;

        // new rho comes from old rho and amount from removed particle
        rhos_[nb] = (nb_areas[i] * rhos_[nb] + (vor_.area(nb) - nb_areas[i]) * rho_removed) / vor_.area(nb);
      
        if (rhos_[nb] < 10*std::numeric_limits<double>::epsilon()) {
        #ifdef PIRANHA_DEBUG
          //std::cout << "nb_areas[i] " << nb_areas[i] << " rhos[nb] " << rhos_[nb] << std::endl;
          std::ostringstream oss;
          oss << "Negative rho, " << rhos_[nb] << ", detected at vertex " << nb;
        #endif
          rhos_[nb] = 0;
        }
      }
    }

    //std::cout << "finished updating rhos" << std::endl;
  }
}

END_PIRANHA_NAMESPACE

#endif // PIRANHA_ITERATEDVORONOISUBTRACTORBASE_HH
