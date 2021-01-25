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
template<class ParticleWeight>
std::string RecursiveSafeSubtractor<ParticleWeight>::description() const {
  std::ostringstream oss;
  std::string jet_alg_descr;
  try { jet_alg_descr = JetDefinition::algorithm_description(jet_alg_); }
  catch (Error & e) {
    throw std::runtime_error(e.message());
  }

  oss << "RecursiveSafeSubtractor\n"
      << "  z - " << default_z_ << '\n'
      << "  f - " << default_f_ << '\n'
      << "  jet_alg - " << JetDefinition::algorithm_description(jet_alg_) << '\n'
      << "  particle_weight - " << ParticleWeight::name() << '\n';
  return oss.str();
}

// internally set jet and groom according to default parameters
template<class ParticleWeight>
std::forward_list<PseudoJet>
RecursiveSafeSubtractor<ParticleWeight>::operator()(const PseudoJet & jet) {

  // recluster the jet
  if (!jet.has_valid_cs())
    throw std::runtime_error("jet must have an associated ClusterSequence");
  reclustered_jet_ = recluster_(jet);

  return init_and_run();
}

// internally set jet (given as a vector of PseudoJets) and groom according to default parameters
template<class ParticleWeight>
std::forward_list<PseudoJet>
RecursiveSafeSubtractor<ParticleWeight>::operator()(const std::vector<PseudoJet> & pjs) {
  ClusterSequence * cs = new ClusterSequence(pjs, JetDefinition(jet_alg_, JetDefinition::max_allowable_R));
  reclustered_jet_ = cs->inclusive_jets()[0];
  cs->delete_self_when_unused();

  return init_and_run();
}

// groom the internally set jet according to the parameters z and f
// this can be used to avoid multiple reclusterings/initializings of the same jet
template<class ParticleWeight>
std::forward_list<PseudoJet> RecursiveSafeSubtractor<ParticleWeight>::apply(double z, double f) {
  validate_params(z, f);
  if (total_weight_ == INVALID_TOTAL_WEIGHT)
    throw std::runtime_error("cannot use `apply` without first calling RecursiveSafeSubtractor directly on a jet or vector of PseudoJets");

  f_ = f;
  return do_recursive_subtraction(reclustered_jet_, z * total_weight_);
}

// initialize and run for the internally stored jet
template<class ParticleWeight>
std::forward_list<PseudoJet> RecursiveSafeSubtractor<ParticleWeight>::init_and_run() {

  // initialize weight_sums
  weight_sums_.resize(2*reclustered_jet_.validated_cs()->n_particles());
  total_weight_ = initialize_weight_sums(reclustered_jet_);

  // groom according to default parameters
  f_ = default_f_;
  return do_recursive_subtraction(reclustered_jet_, default_z_ * total_weight_);
}

// recursive function that tallies the pt of the particles of this jet
template<class ParticleWeight>
double RecursiveSafeSubtractor<ParticleWeight>::initialize_weight_sums(const PseudoJet & pj) {
  PseudoJet piece1, piece2;
  if (!pj.has_parents(piece1, piece2))
    return weight_sums_[pj.cluster_hist_index()] = ParticleWeight::weight(pj);
  return weight_sums_[pj.cluster_hist_index()] = initialize_weight_sums(piece1) + initialize_weight_sums(piece2);
}

// recursive function that applies the safe drop by subtracting weight_sub from the jet
template<class ParticleWeight>
std::forward_list<PseudoJet>
RecursiveSafeSubtractor<ParticleWeight>::do_recursive_subtraction(const PseudoJet & pj, double weight_sub) const {

  // check for boundary case where weight_sub is zero
  if (weight_sub <= 0) {
    if (pj.has_constituents()) {
      std::vector<PseudoJet> consts(pj.constituents());
      return std::forward_list<PseudoJet>(consts.begin(), consts.end());
    }
    else return std::forward_list<PseudoJet>{pj};
  }

  // recurse to parents of this PseudoJet, if none then subtract the weight
  PseudoJet piece1, piece2;
  if (!pj.has_parents(piece1, piece2)) {
    std::forward_list<PseudoJet> pjlist{pj};
    ParticleWeight::set_weight(pjlist.front(), ParticleWeight::weight(pj) - weight_sub);
    return pjlist;
  }

  // lookup weight_sums of each parent that were cached in initialization
  double weightsum1(weight_sums_[piece1.cluster_hist_index()]), weightsum2(weight_sums_[piece2.cluster_hist_index()]);

  // sort out which branch is harder
  bool piece2harder(weightsum2 > weightsum1);
  double maxweight(weightsum1), minweight(weightsum2);
  if (piece2harder) {
    maxweight = weightsum2;
    minweight = weightsum1;
  }
  const PseudoJet & hard(piece2harder ? piece2 : piece1), & soft(piece2harder ? piece1 : piece2);
  double f(minweight == maxweight ? 0.5 : f_), weight_subsoft(f*weight_sub), weight_subhard((1-f)*weight_sub);

  // soft prong dies entirely
  if (minweight <= weight_subsoft)
    return do_recursive_subtraction(hard, weight_sub - minweight);

  // hard prong dies entirely (can happen if f < 0.5)
  if (maxweight <= weight_subhard)
    return do_recursive_subtraction(soft, weight_sub - maxweight);

  // run RecursiveSafeSubtractor on the hard and soft branches, splice, and return result
  std::forward_list<PseudoJet> pieces(do_recursive_subtraction(soft, weight_subsoft));
  pieces.splice_after(pieces.before_begin(), do_recursive_subtraction(hard, weight_subhard));
  return pieces;
}

// specify explicit templates
template class RecursiveSafeSubtractor<emd::TransverseMomentum>;
template class RecursiveSafeSubtractor<emd::TransverseEnergy>;
template class RecursiveSafeSubtractor<emd::Energy>;
template class RecursiveSafeSubtractor<emd::Momentum>;

END_PIRANHA_NAMESPACE
