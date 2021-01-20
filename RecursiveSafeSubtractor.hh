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

#ifndef PIRANHA_RECURSIVESAFESUBTRACTOR_HH
#define PIRANHA_RECURSIVESAFESUBTRACTOR_HH

// C++ standard library
#include <forward_list>
#include <sstream>
#include <string>
#include <vector>

// FastJet
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/Recluster.hh"

// Piranha contrib
#include "PiranhaUtils.hh"

#define INVALID_PTTOT -1

BEGIN_PIRANHA_NAMESPACE

class RecursiveSafeSubtractor {
public:

  // constructor setting default values of z and f and (re)clustering algorithm
  RecursiveSafeSubtractor(double z, double f = 1, JetAlgorithm jet_alg=cambridge_algorithm) :
    jet_alg_(jet_alg),
    recluster_(jet_alg_),
    default_z_(z), default_f_(f), pttot_(INVALID_PTTOT)
  {
    if (z < 0 || z > 1) throw std::invalid_argument("z must be in the range [0,1]");
    if (f < 0 || f > 1) throw std::invalid_argument("f must be in the range [0,1]");
  }

  // describes the object
  std::string description() const;

  // internally set jet and groom according to default parameters
  std::forward_list<PseudoJet> operator()(const PseudoJet & jet);

  // internally set jet (given as a vector of PseudoJets) and groom according to default parameters
  std::forward_list<PseudoJet> operator()(const std::vector<PseudoJet> & pjs);

  // groom the internally set jet according to the parameters z and f
  // this can be used to avoid multiple reclusterings/initializings of the same jet
  std::forward_list<PseudoJet> apply(double z, double f);

private:

  // initialize and run for the internally stored jet
  std::forward_list<PseudoJet> init_and_run();

  // recursive function that tallies the pt of the particles of this jet
  double initialize_pt_sums(const PseudoJet & jet);

  // recursive function that applies the safe drop by subtracting ptsub from the jet
  std::forward_list<PseudoJet> do_RecursiveSafeSubtractor(const PseudoJet & jet, double ptsub) const;

  // the jet algorithm to recluster with
  JetAlgorithm jet_alg_;

  // object that facilitates reclustering
  Recluster recluster_;

  // the internally stored reclustered jet
  PseudoJet reclustered_jet_;

  // cached ptsum values for each point in clustering tree
  std::vector<double> ptsums_;

  // stored parameters
  double default_z_, default_f_, pttot_, f_;
};

END_PIRANHA_NAMESPACE

#endif // PIRANHA_RECURSIVESAFESUBTRACTOR_HH