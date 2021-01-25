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

// a small value to use for debugging
const double piranha_epsilon = 10*std::numeric_limits<double>::epsilon();

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
  IteratedVoronoiSubtractorBase() {
    set_background_estimator(nullptr);
    set_rho_subtraction_mode(RhoSubtractionMode::Additive);
  }

  // main constructor
  template<typename... Args>
  IteratedVoronoiSubtractorBase(SubtractionType subtype, double R,
                                bool allow_repeats, bool store_history,
                                Args && ... args) :
    vor_(subtype, R, std::forward<Args>(args)...),
    zpt_(0), zemd_(0),
    emdtot_(-1),
    nvalid_(0), nremoved_(0)
  {
    // set options
    set_allow_repeats(allow_repeats);
    set_store_history(store_history);

    // set defaults
    set_background_estimator(nullptr);
    set_jet_constituents_use_same_rho(true);
  }

  // destructor
  virtual ~IteratedVoronoiSubtractorBase() {}

  // a description of this object
  std::string description() const;

  // set the z values explicitly, generally for use in grooming
  // for subtype = Area, zpt = z
  // for subtype = EMD, zpt = 1, zemd = z
  void set_z(double z) {
    set_zemd(z);
    set_zpt(subtraction_type() == SubtractionType::EMD ? 1.0 : z);
  }
  void set_zpt(double z) {
    if (z > 1 || z < 0)
      throw std::invalid_argument("zpt cannot be negative or greater than one");
    zpt_ = z;
  }
  void set_zemd(double z) {
    if (z < 0)
      throw std::invalid_argument("zemd cannot be negative");
    zemd_ = z;
  }

  // sets the background estimator, generally for use in pileup removal
  void set_background_estimator(BackgroundEstimatorBase * bge) { bge_ = bge; }
  void set_rho_subtraction_mode(RhoSubtractionMode rho_mode) { rho_mode_ = rho_mode; }

  // sets whether constituents of jets use local rhos or not
  void set_jet_constituents_use_same_rho(bool b) { jet_constituents_use_same_rho_ = b; }

  // setter functions for constructor params
  void set_subtraction_type(SubtractionType subtype) { vor_.set_subtraction_type(subtype); }
  void set_R(double R) { vor_.set_R(R); }
  void set_allow_repeats(bool allow) { allow_repeats_ = allow; }
  void set_store_history(bool store) { store_history_ = store; }

  // getter functions for options stored in this object
  double zpt() const { return zpt_; }
  double zemd() const { return zemd_; }
  BackgroundEstimatorBase * background_estimator() { return bge_; }
  const BackgroundEstimatorBase * background_estimator() const { return bge_; }
  bool jet_constituents_use_same_rho() const { return jet_constituents_use_same_rho_; }
  RhoSubtractionMode rho_subtraction_mode() const { return rho_mode_; }
  bool track_emds() const { return vor().track_emds(); }
  double total_area() const { return vor().total_area(); }

  // getter functions for constructor params
  SubtractionType subtraction_type() const { return vor().subtraction_type(); }
  double R() const { return vor().R(); }
  bool allow_repeats() const { return allow_repeats_; }
  bool store_history() const { return store_history_; }
  
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

  bool have_bge() const { return background_estimator() != nullptr; }
  bool have_rescaling() const {
    return have_bge() && background_estimator()->rescaling_class() != nullptr;
  }

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
  bool do_global_subtraction(double rho, int iremoved);

  // subtract using rhos set for each particle
  void subtract_using_areas_local_rho();

  // energy fractions available for subtraction by area or emd
  double zpt_, zemd_;

  // whether or not we allow repeated applications of the Subtractor on the same set of particles
  bool allow_repeats_, store_history_, jet_constituents_use_same_rho_;

  // control how rho is subtracted
  RhoSubtractionMode rho_mode_;

  // 2D points for each PseudoJet
  std::vector<Point> points_;

  // tracks whether a point was considered "valid" according to the DynamicVoronoi region
  std::vector<char> validmask_;

  // the pts of each PseudoJet, to be adjusted (possibly to zero) by the subtraction
  std::vector<double> pts_, rhos_;

  // this is only populated if repeats are allowed (saves redoing triangulation)
  std::vector<PseudoJet> particles_;

  // some quantities related to the current event
  double pttot_, rhotot_, emdtot_;

  // number of PseudoJets initially valid
  unsigned nvalid_;

  // number of PseudoJets completely removed by the subtraction
  unsigned nremoved_;

  // vector of subtraction histories, in case those are being stored
  std::vector<SubtractionHistory> subtraction_histories_;

}; // IteratedVoronoiSubtractorBase

END_PIRANHA_NAMESPACE

#endif // PIRANHA_ITERATEDVORONOISUBTRACTORBASE_HH
