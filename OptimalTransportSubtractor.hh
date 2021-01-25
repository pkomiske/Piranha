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

// this avoids unnecessary compilation of EMD templates
#ifndef SWIG_PREPROCESSOR
DECLARE_EMD_TEMPLATES
#endif

BEGIN_PIRANHA_NAMESPACE

////////////////////////////////////////////////////////////////////////////////
// GhostGridBase
////////////////////////////////////////////////////////////////////////////////

class GhostGridBase {
public:

  GhostGridBase(double drap, double dphi) : 
    drap_(drap), dphi_(dphi), constructed_points_(false)
  {
    if (drap_ < 0 || dphi_ < 0)
      throw std::invalid_argument("drap and dphi must be positive");
  }
  GhostGridBase(unsigned nrap, unsigned nphi) :
    nrap_(nrap), nphi_(nphi), constructed_points_(false)
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

  // access parameters
  double drap() const { return drap_; }
  double dphi() const { return dphi_; }
  unsigned nrap() const { return nrap_; }
  unsigned nphi() const { return nphi_; }

  // access points/ghosts
  const std::vector<std::pair<double, double>> & points() {
    if (!constructed_points_)
      construct_points();

    return points_;
  }
  const std::vector<PseudoJet> & ghosts() const { return ghosts_; }

  // the number of ghosts there will be
  std::size_t nghosts() { return points().size(); }

  // create ghosts such that their total weight is a given amount
  template<class ParticleWeight>
  const std::vector<PseudoJet> & ghosts_with_total_weight(double total_ghost_weight, double rap_offset = 0, double phi_offset = 0) {
    return ghosts_with_individual_weight<ParticleWeight>(total_ghost_weight/nghosts(), rap_offset, phi_offset);
  }

  // create ghosts, each with a given pt
  template<class ParticleWeight>
  const std::vector<PseudoJet> & ghosts_with_individual_weight(double ghost_weight, double rap_offset = 0, double phi_offset = 0) {
    ghosts_.clear();
    for (const std::pair<double, double> & p : points()) {
      ghosts_.push_back(PtYPhiM(0, p.first + rap_offset, p.second + phi_offset));
      ParticleWeight::set_weight(ghosts_.back(), ghost_weight);
    }
    return ghosts();
  }

protected:

  double drap_, dphi_, rap_start_, phi_start_;
  unsigned nrap_, nphi_;
  bool constructed_points_;

  std::vector<std::pair<double, double>> points_;
  std::vector<PseudoJet> ghosts_;

  virtual bool keep_point(double rap, double phi) const { return true; }
  void construct_points();

}; // GhostGridBase

////////////////////////////////////////////////////////////////////////////////
// GhostGridRectangle
////////////////////////////////////////////////////////////////////////////////

class GhostGridRectangle : public GhostGridBase {

  // store the corners of the rectangle
  std::pair<double, double> lower_left_, upper_right_;

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

  std::string description() const;

  // access region as pair of points
  std::pair<std::pair<double,double>,std::pair<double,double>> region() const {
    return std::make_pair(lower_left_, upper_right_);
  }

private:

  void setup(double rap_min, double phi_min, double rap_max, double phi_max);

}; // GhostGridRectangle

////////////////////////////////////////////////////////////////////////////////
// GhostGridDisk
////////////////////////////////////////////////////////////////////////////////

// Note: the disk will always be centered at (0, 0). To get a translated disk,
//       use the extra arguments to the `ghosts_...` functions of base class.
class GhostGridDisk : public GhostGridBase {

  // store radius and radius squared
  double R_, R2_;

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

  std::string description() const;
  double R() const { return R_; }

private:

  bool keep_point(double rap, double phi) const { return (rap*rap + phi*phi <= R2_); }
  void setup(double R);

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
  std::string description() const;

  // access underlying objects
  const EMD & emd_obj() const { return emd_obj_; }
  const 
  double z() const { return z_; }
  double total_subtracted() const { return total_subtracted_; }
  
  // operate on a single PseudoJet
  std::vector<PseudoJet> operator()(const PseudoJet & jet, double min_weight_to_keep = 1e-14) {

    // get constituents
    if (!jet.has_constituents())
      throw std::runtime_error("jet must have constituents in order to subtract");
    std::vector<PseudoJet> jet_consts(jet.constituents());

    // do the subtracting
    return subtract(jet_consts, construct_ghosts(jet_consts, jet.rap(), jet.phi()), min_weight_to_keep);
  }

  // operator on a vector of PseudoJets
  std::vector<PseudoJet> operator()(const std::vector<PseudoJet> & pjs,
                                    const PseudoJet & offset = PtYPhiM(0, 0, 0),
                                    double min_weight_to_keep = 1e-14) {

    // do the subtracting
    return subtract(pjs, construct_ghosts(pjs, offset.rap(), offset.phi()), min_weight_to_keep);
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
    if (emd_obj().norm())
      throw std::invalid_argument("EMD object should have norm = false");
  }

  std::vector<PseudoJet> construct_ghosts(const std::vector<PseudoJet> & pjs, double rap_off, double phi_off) const;
  std::vector<PseudoJet> subtract(const std::vector<PseudoJet> & pjs,
                                  const std::vector<PseudoJet> & ghosts,
                                  double min_weight_to_keep);

  EMD emd_obj_;
  std::shared_ptr<GhostGridBase> grid_ptr_;
  double z_, total_subtracted_;

}; // OptimalTransportSubtractor

END_PIRANHA_NAMESPACE

#endif // PIRANHA_OPTIMALTRANSPORTSUBTRACTOR_HH
