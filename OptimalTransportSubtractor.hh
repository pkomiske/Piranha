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

#ifndef PIRANHA_OPTIMALTRANSPORTSUBTRACTOR_HH
#define PIRANHA_OPTIMALTRANSPORTSUBTRACTOR_HH

// FastJet
#include "fastjet/PseudoJet.hh"

// EventGeometry contrib
#include "fastjet/contrib/EventGeometry.hh"

// Piranha contrib
#include "PiranhaUtils.hh"

BEGIN_PIRANHA_NAMESPACE

////////////////////////////////////////////////////////////////////////////////
// GhostGridBase
////////////////////////////////////////////////////////////////////////////////

class GhostGridBase {
public:

  GhostGridBase(double drap, double dphi) : 
    drap_(drap), dphi_(dphi)
  {
    if (drap_ < 0 || dphi_ < 0)
      throw std::invalid_argument("drap and dphi must be positive");
  }
  GhostGridBase(unsigned nrap, unsigned nphi) :
    nrap_(nrap), nphi_(nphi)
  {
    if (nrap_ == 0 || nphi_ == 0)
      throw std::invalid_argument("nrap and nphi must be positive");
  }

  virtual ~GhostGridBase() {}

  virtual std::string description() const {
    std::ostringstream oss;
    oss << "  (nrap, drap) - (" << nrap_ << ", " << drap_ << ")\n"
        << "  (nphi, dphi) - (" << nphi_ << ", " << dphi_ << ")\n";
    return oss.str();
  }

  // the number of ghosts there will be
  std::size_t nghosts() const { return points_.size(); }

  // create ghosts such that their total pT is a given amount
  const std::vector<PseudoJet> & ghosts_with_total_pt(double total_ghost_pt, double rap_offset = 0, double phi_offset = 0) {
    return ghosts_with_pt(total_ghost_pt/nghosts(), rap_offset, phi_offset);
  }

  // create ghosts, each with a given pt
  const std::vector<PseudoJet> & ghosts_with_pt(double ghost_pt, double rap_offset = 0, double phi_offset = 0) {
    ghosts_.clear();
    for (const std::pair<double, double> & p : points_)
      ghosts_.push_back(PtYPhiM(ghost_pt, p.first + rap_offset, p.second + phi_offset));
    return ghosts_;
  }

protected:

  double drap_, dphi_, rap_start_, phi_start_;
  unsigned nrap_, nphi_;

  std::vector<std::pair<double, double>> points_;
  std::vector<PseudoJet> ghosts_;

  virtual bool keep_point(double rap, double phi) const { return true; }

  void construct_points() {
    points_.reserve(std::size_t(nrap_)*std::size_t(nphi_));
    for (unsigned i = 0; i < nrap_; i++) {
      double rap(rap_start_ + i*drap_);
      for (unsigned j = 0; j < nphi_; j++) {
        double phi(phi_start_ + j*dphi_);
        if (keep_point(rap, phi))
          points_.emplace_back(rap, phi);
      }
    }
  }

}; // GhostGridBase

////////////////////////////////////////////////////////////////////////////////
// GhostGridRectangle
////////////////////////////////////////////////////////////////////////////////

class GhostGridRectangle : public GhostGridBase {
public:

  // partial constructors, assume full phi range and symmetric rapidity coverage
  GhostGridRectangle(double abs_rap_max, double d) :
    GhostGridRectangle(-abs_rap_max, 0, abs_rap_max, TWOPI, d, d)
  {}
  GhostGridRectangle(double abs_rap_max, unsigned n) :
    GhostGridRectangle(-abs_rap_max, 0, abs_rap_max, TWOPI, n, n)
  {}

  // partial constructors, assume full phi range and specified rapidity coverage
  GhostGridRectangle(double rap_min, double rap_max, double drap, double dphi) :
    GhostGridRectangle(rap_min, 0, rap_max, TWOPI, drap, dphi)
  {}
  GhostGridRectangle(double rap_min, double rap_max, unsigned nrap, unsigned nphi) :
    GhostGridRectangle(rap_min, 0, rap_max, TWOPI, nrap, nphi)
  {}

  // full constructors
  GhostGridRectangle(double rap_min, double phi_min, double rap_max, double phi_max,
                               double drap, double dphi) :
    GhostGridBase(drap, dphi)
  {
    nrap_ = unsigned((rap_max - rap_min)/drap) + 1;
    nphi_ = unsigned((phi_max - phi_min)/dphi) + 1;
    setup(rap_min, phi_min, rap_max, phi_max);
  }
  GhostGridRectangle(double rap_min, double phi_min, double rap_max, double phi_max,
                               unsigned nrap, unsigned nphi) :
    GhostGridBase(nrap, nphi)
  {
    drap_ = std::min((rap_max - rap_min)/(nrap - 1), double(1000));
    dphi_ = std::min((phi_max - phi_min)/(nphi - 1), double(1000));
    setup(rap_min, phi_min, rap_max, phi_max);
  }

  std::string description() const {
    std::ostringstream oss;
    oss << "GhostGridRectangle\n"
        << "  (rap_min, phi_min) - (" << rap_min_ << ", " << phi_min_ << ")\n"
        << "  (rap_max, phi_max) - (" << rap_max_ << ", " << phi_max_ << ")\n"
        << GhostGridBase::description();
    return oss.str();
  }

private:

  double rap_min_, phi_min_, rap_max_, phi_max_;

  void setup(double rap_min, double phi_min, double rap_max, double phi_max) {
    if (rap_min >= rap_max)
      throw std::invalid_argument("rap_min must be less than rap_max");
    if (phi_min >= phi_max)
      throw std::invalid_argument("phi_min must be less than phi_max");

    rap_min_ = rap_min;
    phi_min_ = phi_min;
    rap_max_ = rap_max;
    phi_max_ = phi_max;
    rap_start_ = (rap_max_ + rap_min_)/2 - (nrap_ - 1)*drap_/2;
    phi_start_ = (phi_max_ + phi_min_)/2 - (nphi_ - 1)*dphi_/2;

    construct_points();
  }

}; // GhostGridRectangle

////////////////////////////////////////////////////////////////////////////////
// GhostGridDisk
////////////////////////////////////////////////////////////////////////////////

// Note: the disk will always be centered at (0, 0). To get a translated disk,
//       use the extra arguments to the `ghosts_...` functions of base class.
class GhostGridDisk : public GhostGridBase {
public:

  // partial constructors, symmetric treatment of rap and phi
  GhostGridDisk(double R, double d) : GhostGridDisk(R, d, d) {}
  GhostGridDisk(double R, unsigned n) : GhostGridDisk(R, n, n) {}

  // full constructors, specify center and different treatment of rap and phi
  GhostGridDisk(double R, double drap, double dphi) :
    GhostGridBase(drap, dphi)
  {
    nrap_ = unsigned(2*R/drap) + 1;
    nphi_ = unsigned(2*R/dphi) + 1;
    setup(R);
  }
  GhostGridDisk(double R, unsigned nrap, unsigned nphi) :
    GhostGridBase(nrap, nphi)
  {
    drap_ = std::min(2*R/(nrap - 1), double(1000));
    dphi_ = std::min(2*R/(nphi - 1), double(1000));
    setup(R);
  }

  std::string description() const {
    std::ostringstream oss;
    oss << "GhostGridDisk\n"
        << "  R - " << R_ << '\n'
        << GhostGridBase::description();
    return oss.str();
  }

private:

  double R_, R2_;

  bool keep_point(double rap, double phi) const {
    return (rap*rap + phi*phi <= R2_);
  }

  void setup(double R) {
    if (R <= 0)
      throw std::invalid_argument("R must be positive");

    R_ = R;
    R2_ = R*R;
    rap_start_ = -drap_*(nrap_ - 1)/2;
    phi_start_ = -dphi_*(nphi_ - 1)/2;

    construct_points();
  }

}; // GhostGridDisk

////////////////////////////////////////////////////////////////////////////////
// OptimalTransportSubtractor
////////////////////////////////////////////////////////////////////////////////

template<class EMD>
class OptimalTransportSubtractor {
public:

  // constructor from an EMD object and an GhostGrid object
  OptimalTransportSubtractor(double z, const EMD & emd, const GhostGridDisk & grid) :
    OptimalTransportSubtractor(z, emd, grid, 0)
  {}

  OptimalTransportSubtractor(double z, const EMD & emd, const GhostGridRectangle & grid) :
    OptimalTransportSubtractor(z, emd, grid, 0)
  {}

  // return a description of the object
  std::string description() const {
    std::ostringstream oss;
    oss << "OptimalTransportSubtractor\n"
        << "  z - " << z_ << '\n'
        << '\n'
        << grid_ptr_->description() << '\n'
        << emd_obj_.description();
    return oss.str();
  }

  // access underlying EMD object
  const EMD & emd_obj() const { return emd_obj_; }
  double total_subtracted() const { return total_subtracted_; }
  
  // operate on a single PseudoJet
  std::vector<PseudoJet> operator()(const PseudoJet & jet, double min_weight_to_keep = 1e-14) {

    if (!jet.has_constituents())
      throw std::runtime_error("jet must have constituents in order to subtract");

    // get constituents
    std::vector<PseudoJet> jet_consts(jet.constituents());

    // get ghosts
    double total_ghost_weight(z_*get_total_weight(jet_consts));
    std::vector<PseudoJet> ghosts(grid_ptr_->ghosts_with_total_pt(total_ghost_weight, jet.rap(), jet.phi()));

    // do the subtracting
    return subtract(jet_consts, ghosts, min_weight_to_keep);
  }

  // operator on a vector of PseudoJets
  std::vector<PseudoJet> operator()(const std::vector<PseudoJet> & pjs,
                                    const PseudoJet & offset = PtYPhiM(0, 0, 0),
                                    double min_weight_to_keep = 1e-14) {

    // get ghosts
    double total_ghost_weight(z_*get_total_weight(pjs));
    std::vector<PseudoJet> ghosts(grid_ptr_->ghosts_with_total_pt(total_ghost_weight, offset.rap(), offset.phi()));

    // do the subtracting
    return subtract(pjs, ghosts, min_weight_to_keep);
  }

private:

  template<class GhostGrid>
  OptimalTransportSubtractor(double z, const EMD & emd, const GhostGrid & grid, int ignored) :
    emd_obj_(emd),
    grid_ptr_(new GhostGrid(grid)),
    z_(z), total_subtracted_(0)
  {
    if (z < 0 || z > 1)
      throw std::invalid_argument("z must be in [0,1]");
    if (emd_obj_.norm())
      throw std::invalid_argument("EMD object should have norm = false");
  }

  double get_total_weight(const std::vector<PseudoJet> & pjs) const {
    double total_weight(0);
    for (const PseudoJet & pj : pjs)
      total_weight += EMD::ParticleWeight::weight(pj);
    return total_weight;
  }

  std::vector<PseudoJet> subtract(const std::vector<PseudoJet> & pjs,
                               const std::vector<PseudoJet> & ghosts,
                               double min_weight_to_keep) {

    // run emd computation
    emd_obj_(pjs, ghosts);

    // verify that extra particle went to the ghosts
    if (emd_obj_.extra() == emd::ExtraParticle::Zero)
      throw std::runtime_error("event should not have gotten an extra particle");

    // get flow vector
    std::vector<double> flows(emd_obj_.flows());

    // subtract pt from each PseudoJet
    std::vector<PseudoJet> subtracted_pjs;
    subtracted_pjs.reserve(pjs.size());
    total_subtracted_ = 0;
    for (unsigned i = 0; i < pjs.size(); i++) {

      // tally pt to subtract
      double wsub_i(0);
      for (unsigned j = 0, in1 = i*emd_obj_.n1(); j < ghosts.size(); j++)
        wsub_i += flows[in1 + j];
      total_subtracted_ += wsub_i;

      double new_w(EMD::ParticleWeight::weight(pjs[i]) - wsub_i);
      if (new_w > min_weight_to_keep) {
        subtracted_pjs.push_back(pjs[i]);
        EMD::ParticleWeight::set_weight(subtracted_pjs.back(), new_w);
      }
    }

    return subtracted_pjs;
  }

  EMD emd_obj_;
  std::shared_ptr<GhostGridBase> grid_ptr_;
  double z_, total_subtracted_;

}; // OptimalTransportSubtractor

END_PIRANHA_NAMESPACE

#endif // PIRANHA_OPTIMALTRANSPORTSUBTRACTOR_HH