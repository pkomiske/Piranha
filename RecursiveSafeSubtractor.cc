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

#include <stdexcept>

#include "fastjet/ClusterSequence.hh"

#include "RecursiveSafeSubtractor.hh"

BEGIN_PIRANHA_NAMESPACE

// describes the object
std::string RecursiveSafeSubtractor::description() const {
  std::ostringstream oss;
  std::string jet_alg_descr;
  try { jet_alg_descr = JetDefinition::algorithm_description(jet_alg_); }
  catch (Error & e) {
    throw std::runtime_error(e.message());
  }

  oss << "RecursiveSafeSubtractor\n"
      << "  z - " << default_z_ << '\n'
      << "  f - " << default_f_ << '\n'
      << "  jet_alg - " << JetDefinition::algorithm_description(jet_alg_) << '\n';
  return oss.str();
}

// sets RecursiveSafeSubtractor to operate on this jet
std::forward_list<PseudoJet> RecursiveSafeSubtractor::operator()(const PseudoJet & jet) {

  // recluster the jet
  if (!jet.has_valid_cs())
    throw std::runtime_error("jet must have an associated ClusterSequence");
  reclustered_jet_ = recluster_(jet);

  return init_and_run();
}

// sets RecursiveSafeSubtractor to operate on a vector of PseudoJets
std::forward_list<PseudoJet> RecursiveSafeSubtractor::operator()(const std::vector<PseudoJet> & pjs) {
  
  ClusterSequence * cs = new ClusterSequence(pjs, JetDefinition(jet_alg_, JetDefinition::max_allowable_R));
  reclustered_jet_ = cs->inclusive_jets()[0];
  cs->delete_self_when_unused();

  return init_and_run();
}

// initialize and run for the internally stored jet
std::forward_list<PseudoJet> RecursiveSafeSubtractor::init_and_run() {

  // initialize ptsums
  ptsums_.resize(2*reclustered_jet_.validated_cs()->n_particles());
  pttot_ = initialize_pt_sums(reclustered_jet_);

  // groom according to default parameters
  f_ = default_f_;
  return do_RecursiveSafeSubtractor(reclustered_jet_, default_z_ * pttot_);
}

// grooms the internally set jet according to z and f parameters
std::forward_list<PseudoJet> RecursiveSafeSubtractor::apply(double z, double f) {
  if (z < 0 || z > 1) throw std::invalid_argument("z must be in the range [0,1]");
  if (f < 0 || f > 1) throw std::invalid_argument("f must be in the range [0,1]");
  if (pttot_ == INVALID_PTTOT)
    throw std::runtime_error("cannot use `apply` without first calling RecursiveSafeSubtractor directly on a jet or vector of PseudoJets");

  f_ = f;
  return do_RecursiveSafeSubtractor(reclustered_jet_, z * pttot_);
}

// recursive method, returns the sum of pts of the pieces of this pseudojet
double RecursiveSafeSubtractor::initialize_pt_sums(const PseudoJet & pj) {
  PseudoJet piece1, piece2;
  if (!pj.has_parents(piece1, piece2))
    return ptsums_[pj.cluster_hist_index()] = pj.pt();
  return ptsums_[pj.cluster_hist_index()] = initialize_pt_sums(piece1) + initialize_pt_sums(piece2);
}

// recursive method, applies RecursiveSafeSubtractor algorithm
std::forward_list<PseudoJet> RecursiveSafeSubtractor::do_RecursiveSafeSubtractor(const PseudoJet & pj, double ptsub) const {

  // check for boundary case where ptsub is zero
  if (ptsub <= 0) {
    if (pj.has_constituents()) {
      std::vector<PseudoJet> consts(pj.constituents());
      return std::forward_list<PseudoJet>(consts.begin(), consts.end());
    }
    else return std::forward_list<PseudoJet>{pj};
  }

  // recurse to parents of this PseudoJet, if none then subtract the pt
  PseudoJet piece1, piece2;
  if (!pj.has_parents(piece1, piece2))
    return std::forward_list<PseudoJet>{PtYPhiM(pj.pt() - ptsub, pj.rap(), pj.phi(), pj.m())};;

  // lookup ptsums of each parent that were cached in initialization
  double ptsum1(ptsums_[piece1.cluster_hist_index()]), ptsum2(ptsums_[piece2.cluster_hist_index()]);

  // sort out which branch is harder
  bool piece2harder(ptsum2 > ptsum1);
  double maxpt(ptsum1), minpt(ptsum2);
  if (piece2harder) {
    maxpt = ptsum2;
    minpt = ptsum1;
  }
  const PseudoJet & hard(piece2harder ? piece2 : piece1), & soft(piece2harder ? piece1 : piece2);
  double f(minpt == maxpt ? 0.5 : f_), ptsubsoft(f*ptsub), ptsubhard((1-f)*ptsub);

  // soft prong dies entirely
  if (minpt <= ptsubsoft)
    return do_RecursiveSafeSubtractor(hard, ptsub - minpt);

  // hard prong dies entirely (can happen if f < 0.5)$
  if (maxpt <= ptsubhard)
    return do_RecursiveSafeSubtractor(soft, ptsub - maxpt);

  // run RecursiveSafeSubtractor on the hard and soft branches, splice, and return result
  std::forward_list<PseudoJet> pieces(do_RecursiveSafeSubtractor(soft, ptsubsoft));
  pieces.splice_after(pieces.before_begin(), do_RecursiveSafeSubtractor(hard, ptsubhard));
  return pieces;
}

END_PIRANHA_NAMESPACE
