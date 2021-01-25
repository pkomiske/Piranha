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

#include "OptimalTransportSubtractor.hh"

BEGIN_PIRANHA_NAMESPACE

void GhostGridBase::construct_points() {
  points_.reserve(std::size_t(nrap_)*std::size_t(nphi_));
  for (unsigned i = 0; i < nrap_; i++) {
    double rap(rap_start_ + i*drap_);
    for (unsigned j = 0; j < nphi_; j++) {
      double phi(phi_start_ + j*dphi_);
      if (keep_point(rap, phi))
        points_.emplace_back(rap, phi);
    }
  }
  constructed_points_ = true;
}

std::string GhostGridRectangle::description() const {
  std::ostringstream oss;
  oss << "GhostGridRectangle\n"
      << "  (rap_min, phi_min) - (" << lower_left_.first << ", " << lower_left_.second << ")\n"
      << "  (rap_max, phi_max) - (" << upper_right_.first << ", " << upper_right_.second << ")\n"
      << GhostGridBase::description();
  return oss.str();
}

void GhostGridRectangle::setup(double rap_min, double phi_min, double rap_max, double phi_max) {
  if (rap_min >= rap_max)
    throw std::invalid_argument("rap_min must be less than rap_max");
  if (phi_min >= phi_max)
    throw std::invalid_argument("phi_min must be less than phi_max");

  lower_left_ = std::make_pair(rap_min, phi_min);
  upper_right_ = std::make_pair(rap_max, phi_max);
  rap_start_ = (rap_max + rap_min)/2 - (nrap_ - 1)*drap_/2;
  phi_start_ = (phi_max + phi_min)/2 - (nphi_ - 1)*dphi_/2;
}

std::string GhostGridDisk::description() const {
  std::ostringstream oss;
  oss << "GhostGridDisk\n"
      << "  R - " << R_ << '\n'
      << GhostGridBase::description();
  return oss.str();
}

void GhostGridDisk::setup(double R) {
  if (R <= 0)
    throw std::invalid_argument("R must be positive");

  R_ = R;
  R2_ = R_*R_;
  rap_start_ = -drap_*(nrap_ - 1)/2;
  phi_start_ = -dphi_*(nphi_ - 1)/2;
}

// return a description of the object
template<class EMD>
std::string OptimalTransportSubtractor<EMD>::description() const {
  std::ostringstream oss;
  oss << "OptimalTransportSubtractor\n"
      << "  z - " << z_ << '\n'
      << '\n'
      << grid_ptr_->description() << '\n'
      << emd_obj().description();
  return oss.str();
}

template<class EMD>
std::vector<PseudoJet> OptimalTransportSubtractor<EMD>::construct_ghosts(const std::vector<PseudoJet> & pjs, double rap_off, double phi_off) const {
  double total_weight(0);
  for (const PseudoJet & pj : pjs)
    total_weight += EMD::ParticleWeight::weight(pj);

  return grid_ptr_->ghosts_with_total_weight<typename EMD::ParticleWeight>(z()*total_weight, rap_off, phi_off);
}

template<class EMD>
std::vector<PseudoJet>
OptimalTransportSubtractor<EMD>::subtract(const std::vector<PseudoJet> & pjs,
                                          const std::vector<PseudoJet> & ghosts,
                                          double min_weight_to_keep) {

  // run emd computation
  emd_obj_(pjs, ghosts);

  // verify that extra particle went to the ghosts
  if (emd_obj().extra() == emd::ExtraParticle::Zero)
    throw std::runtime_error("event should not have gotten an extra particle");

  // subtract pt from each PseudoJet
  std::vector<PseudoJet> subtracted_pjs;
  subtracted_pjs.reserve(pjs.size());
  total_subtracted_ = 0;
  for (unsigned i = 0; i < pjs.size(); i++) {

    // tally pt to subtract
    double wsub_i(0);
    for (unsigned j = 0, in1 = i*emd_obj().n1(); j < ghosts.size(); j++)
      wsub_i += emd_obj().flow(in1 + j);
    total_subtracted_ += wsub_i;

    double new_w(EMD::ParticleWeight::weight(pjs[i]) - wsub_i);
    if (new_w > min_weight_to_keep) {
      subtracted_pjs.push_back(pjs[i]);
      EMD::ParticleWeight::set_weight(subtracted_pjs.back(), new_w);
    }
  }

  return subtracted_pjs;
}

// explicit template instantiations
template class OptimalTransportSubtractor<emd::EMD<emd::TransverseMomentum, emd::DeltaR>>;
template class OptimalTransportSubtractor<emd::EMD<emd::TransverseMomentum, emd::HadronicDot>>;
template class OptimalTransportSubtractor<emd::EMD<emd::TransverseMomentum, emd::HadronicDotMassive>>;
template class OptimalTransportSubtractor<emd::EMD<emd::TransverseEnergy, emd::DeltaR>>;
template class OptimalTransportSubtractor<emd::EMD<emd::TransverseEnergy, emd::HadronicDot>>;
template class OptimalTransportSubtractor<emd::EMD<emd::TransverseEnergy, emd::HadronicDotMassive>>;
template class OptimalTransportSubtractor<emd::EMD<emd::Momentum, emd::EEDot>>;
template class OptimalTransportSubtractor<emd::EMD<emd::Momentum, emd::EEDotMassive>>;
template class OptimalTransportSubtractor<emd::EMD<emd::Momentum, emd::EEArcLength>>;
template class OptimalTransportSubtractor<emd::EMD<emd::Momentum, emd::EEArcLengthMassive>>;
template class OptimalTransportSubtractor<emd::EMD<emd::Energy, emd::EEDot>>;
template class OptimalTransportSubtractor<emd::EMD<emd::Energy, emd::EEDotMassive>>;
template class OptimalTransportSubtractor<emd::EMD<emd::Energy, emd::EEArcLength>>;
template class OptimalTransportSubtractor<emd::EMD<emd::Energy, emd::EEArcLengthMassive>>;


END_PIRANHA_NAMESPACE