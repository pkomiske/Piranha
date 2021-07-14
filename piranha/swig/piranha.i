// -*- C++ -*-
%define PIRANHA_DOCSTRING
"# Piranha FastJet Contrib

Python interface to the Piranha FastJet contrib package.
"
%enddef

%module(docstring=PIRANHA_DOCSTRING) piranha

#define PIRANHANAMESPACE fastjet::contrib::piranha

// C++ standard library wrappers
%include <exception.i>
%include <std_list.i>
%include <std_pair.i>
%include <std_string.i>
%include <std_vector.i>

// this makes SWIG aware of the types contained in the main fastjet library
// but does not generate new wrappers for them here
#ifdef FASTJET_PREFIX
  %import FASTJET_PREFIX/share/fastjet/pyinterface/fastjet.i
#else
  //%import EventGeometry/PyFJCore/pyfjcore/swig/pyfjcore.i
#endif

// make library aware of EventGeometry fjcontrib wrappers
%import eventgeometry/swig/eventgeometry.i

// converts fastjet::Error into a FastJetError Python exception
FASTJET_ERRORS_AS_PYTHON_EXCEPTIONS(Piranha)

// include FastJetError in python module
%pythoncode {
from fastjet import FastJetError
}

// turn off exception handling for now, since fastjet::Error is not thrown here
//%exception;

// include headers in source file
%{
#ifndef SWIG
#define SWIG
#endif

//#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
//#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "Piranha.hh"

using namespace fastjet::contrib;
%}

// pair templates
%template(pairDouble) std::pair<double, double>;
%template(pairPairDouble) std::pair<std::pair<double, double>, std::pair<double, double>>;

// vector templates
%template(vectorVectorDouble) std::vector<std::vector<double>>;
%template(vectorSubtractionHistory) std::vector<PIRANHANAMESPACE::SubtractionHistory>;
%template(vectorVectorPseudoJet) std::vector<std::vector<fastjet::PseudoJet>>;
%template(vectorPairDouble) std::vector<std::pair<double, double>>;
%template(vectorInt) std::vector<int>;

// ignored classes/methods
%ignore PIRANHANAMESPACE::PiranhaError;
%ignore PIRANHANAMESPACE::RecursiveSafeSubtractor::operator()(const PseudoJet & jet);
%ignore PIRANHANAMESPACE::RecursiveSafeSubtractor::operator()(const std::vector<PseudoJet> & pjs);
%ignore PIRANHANAMESPACE::RecursiveSafeSubtractor::apply(double z, double f);

// wrap core utils
%include "PiranhaUtils.hh"

// wrap DynamicVoronoi code
%include "DynamicVoronoiBase.hh"
%include "DynamicVoronoiCylinder.hh"
%include "DynamicVoronoiDisk.hh"

// generate grooming wrappers
%include "IteratedVoronoiSubtractorBase.hh"
%include "OptimalTransportSubtractor.hh"
%include "RecursiveSafeSubtractor.hh"

%define PIRANHA_OPTIMAL_TRANSPORT_SUBTRACTOR_TEMPLATE(Weight, Distance)
  %template(OptimalTransportSubtractor##Weight##Distance)
        OptimalTransportSubtractor<eventgeometry::EMD<double, eventgeometry::Weight, eventgeometry::Distance>>;
%enddef

// explicit templates
BEGIN_PIRANHA_NAMESPACE

  %extend IteratedVoronoiSubtractorBase { ADD_STR_FROM_DESCRIPTION() }
  %extend GhostGridDisk { ADD_STR_FROM_DESCRIPTION() }
  %extend GhostGridRectangle { ADD_STR_FROM_DESCRIPTION() }
  %extend OptimalTransportSubtractor { ADD_STR_FROM_DESCRIPTION() }

  %template(IteratedVoronoiSubtractorDiskBase) IteratedVoronoiSubtractorBase<DynamicVoronoiDisk>;
  %template(IteratedVoronoiSubtractorCylinderBase) IteratedVoronoiSubtractorBase<DynamicVoronoiCylinder>;

  PIRANHA_OPTIMAL_TRANSPORT_SUBTRACTOR_TEMPLATE(TransverseMomentum, DeltaR)
  PIRANHA_OPTIMAL_TRANSPORT_SUBTRACTOR_TEMPLATE(TransverseMomentum, HadronicDot)
  PIRANHA_OPTIMAL_TRANSPORT_SUBTRACTOR_TEMPLATE(TransverseMomentum, HadronicDotMassive)
  PIRANHA_OPTIMAL_TRANSPORT_SUBTRACTOR_TEMPLATE(TransverseEnergy, DeltaR)
  PIRANHA_OPTIMAL_TRANSPORT_SUBTRACTOR_TEMPLATE(TransverseEnergy, HadronicDot)
  PIRANHA_OPTIMAL_TRANSPORT_SUBTRACTOR_TEMPLATE(TransverseEnergy, HadronicDotMassive)
  PIRANHA_OPTIMAL_TRANSPORT_SUBTRACTOR_TEMPLATE(Momentum, EEDot)
  PIRANHA_OPTIMAL_TRANSPORT_SUBTRACTOR_TEMPLATE(Momentum, EEDotMassive)
  PIRANHA_OPTIMAL_TRANSPORT_SUBTRACTOR_TEMPLATE(Momentum, EEArcLength)
  PIRANHA_OPTIMAL_TRANSPORT_SUBTRACTOR_TEMPLATE(Momentum, EEArcLengthMassive)
  PIRANHA_OPTIMAL_TRANSPORT_SUBTRACTOR_TEMPLATE(Energy, EEDot)
  PIRANHA_OPTIMAL_TRANSPORT_SUBTRACTOR_TEMPLATE(Energy, EEDotMassive)
  PIRANHA_OPTIMAL_TRANSPORT_SUBTRACTOR_TEMPLATE(Energy, EEArcLength)
  PIRANHA_OPTIMAL_TRANSPORT_SUBTRACTOR_TEMPLATE(Energy, EEArcLengthMassive)

/*%template(OptimalTransportSubtractorTransverseMomentumDeltaR) OptimalTransportSubtractor<eventgeometry::EMD<eventgeometry::TransverseMomentum, eventgeometry::DeltaR>>;
%template(OptimalTransportSubtractorTransverseMomentumHadronicDot) OptimalTransportSubtractor<eventgeometry::EMD<eventgeometry::TransverseMomentum, eventgeometry::HadronicDot>>;
%template(OptimalTransportSubtractorTransverseMomentumHadronicDotMassive) OptimalTransportSubtractor<eventgeometry::EMD<eventgeometry::TransverseMomentum, eventgeometry::HadronicDotMassive>>;
%template(OptimalTransportSubtractorTransverseEnergyDeltaR) OptimalTransportSubtractor<eventgeometry::EMD<eventgeometry::TransverseEnergy, eventgeometry::DeltaR>>;
%template(OptimalTransportSubtractorTransverseEnergyHadronicDot) OptimalTransportSubtractor<eventgeometry::EMD<eventgeometry::TransverseEnergy, eventgeometry::HadronicDot>>;
%template(OptimalTransportSubtractorTransverseEnergyHadronicDotMassive) OptimalTransportSubtractor<eventgeometry::EMD<eventgeometry::TransverseEnergy, eventgeometry::HadronicDotMassive>>;
%template(OptimalTransportSubtractorMomentumEEDot) OptimalTransportSubtractor<eventgeometry::EMD<eventgeometry::Momentum, eventgeometry::EEDot>>;
%template(OptimalTransportSubtractorMomentumEEDotMassive) OptimalTransportSubtractor<eventgeometry::EMD<eventgeometry::Momentum, eventgeometry::EEDotMassive>>;
%template(OptimalTransportSubtractorMomentumEEArcLength) OptimalTransportSubtractor<eventgeometry::EMD<eventgeometry::Momentum, eventgeometry::EEArcLength>>;
%template(OptimalTransportSubtractorMomentumEEArcLengthMassive) OptimalTransportSubtractor<eventgeometry::EMD<eventgeometry::Momentum, eventgeometry::EEArcLengthMassive>>;
%template(OptimalTransportSubtractorEnergyEEDot) OptimalTransportSubtractor<eventgeometry::EMD<eventgeometry::Energy, eventgeometry::EEDot>>;
%template(OptimalTransportSubtractorEnergyEEDotMassive) OptimalTransportSubtractor<eventgeometry::EMD<eventgeometry::Energy, eventgeometry::EEDotMassive>>;
%template(OptimalTransportSubtractorEnergyEEArcLength) OptimalTransportSubtractor<eventgeometry::EMD<eventgeometry::Energy, eventgeometry::EEArcLength>>;
%template(OptimalTransportSubtractorEnergyEEArcLengthMassive) OptimalTransportSubtractor<eventgeometry::EMD<eventgeometry::Energy, eventgeometry::EEArcLengthMassive>>;*/

END_PIRANHA_NAMESPACE

// include this after explicit IVSBase templates
%include "IteratedVoronoiSubtractor.hh"

BEGIN_PIRANHA_NAMESPACE

  %extend RecursiveSafeSubtractor {
    ADD_STR_FROM_DESCRIPTION()

    std::vector<PseudoJet> operator()(const fastjet::PseudoJet & jet, int _ = 0) {
      std::forward_list<fastjet::PseudoJet> result((*$self)(jet));
      return std::vector<fastjet::PseudoJet>(result.begin(), result.end());
    }
    std::vector<PseudoJet> operator()(const std::vector<fastjet::PseudoJet> & pjs, int _ = 0) {
      std::forward_list<fastjet::PseudoJet> result((*$self)(pjs));
      return std::vector<fastjet::PseudoJet>(result.begin(), result.end());
    }
    std::vector<PseudoJet> apply(double z, double f, int _ = 0) {
      std::forward_list<fastjet::PseudoJet> result($self->apply(z, f));
      return std::vector<fastjet::PseudoJet>(result.begin(), result.end()); 
    }
  }

  %template(RecursiveSafeSubtractorTransverseMomentum) RecursiveSafeSubtractor<eventgeometry::TransverseMomentum<double>>;
  %template(RecursiveSafeSubtractorTransverseEnergy) RecursiveSafeSubtractor<eventgeometry::TransverseEnergy<double>>;
  %template(RecursiveSafeSubtractorEnergy) RecursiveSafeSubtractor<eventgeometry::Energy<double>>;
  %template(RecursiveSafeSubtractorMomentum) RecursiveSafeSubtractor<eventgeometry::Momentum<double>>;

END_PIRANHA_NAMESPACE

// add convenience functions for accessing templated OptimalTransportSubtractor classes
%pythoncode %{

def OptimalTransportSubtractor(*args, weight='TransverseMomentum', pairwise_distance='DeltaR', **kwargs):

    if weight == 'TransverseMomentum':
        if pairwise_distance == 'DeltaR':
            return OptimalTransportSubtractorTransverseMomentumDeltaR(*args, **kwargs)
        elif pairwise_distance == 'HadronicDot':
            return OptimalTransportSubtractorTransverseMomentumHadronicDot(*args, **kwargs)
        elif pairwise_distance == 'HadronicDotMassive':
            return OptimalTransportSubtractorTransverseMomentumHadronicDotMassive(*args, **kwargs)
        else:
            raise TypeError('pairwise distance `{}` not recognized'.format(pairwise_distance))

    elif weight == 'TransverseEnergy':
        if pairwise_distance == 'DeltaR':
            return OptimalTransportSubtractorTransverseEnergyDeltaR(*args, **kwargs)
        elif pairwise_distance == 'HadronicDot':
            return OptimalTransportSubtractorTransverseEnergyHadronicDot(*args, **kwargs)
        elif pairwise_distance == 'HadronicDotMassive':
            return OptimalTransportSubtractorTransverseEnergyHadronicDotMassive(*args, **kwargs)
        else:
            raise TypeError('pairwise distance `{}` not recognized'.format(pairwise_distance))

    elif weight == 'Energy':
        if pairwise_distance == 'EEDot':
            return OptimalTransportSubtractorEnergyEEDot(*args, **kwargs)
        elif pairwise_distance == 'EEDotMassive':
            return OptimalTransportSubtractorEnergyEEDotMassive(*args, **kwargs)
        elif pairwise_distance == 'EEArcLength':
            return OptimalTransportSubtractorEnergyEEArcLength(*args, **kwargs)
        elif pairwise_distance == 'EEArcLengthMassive':
            return OptimalTransportSubtractorEnergyEEArcLengthMassive(*args, **kwargs)
        else:
            raise TypeError('pairwise distance `{}` not recognized'.format(pairwise_distance))

    elif weight == 'Momentum':
        if pairwise_distance == 'EEDot':
            return OptimalTransportSubtractorMomentumEEDot(*args, **kwargs)
        elif pairwise_distance == 'EEDotMassive':
            return OptimalTransportSubtractorMomentumEEDotMassive(*args, **kwargs)
        elif pairwise_distance == 'EEArcLength':
            return OptimalTransportSubtractorMomentumEEArcLength(*args, **kwargs)
        elif pairwise_distance == 'EEArcLengthMassive':
            return OptimalTransportSubtractorMomentumEEArcLengthMassive(*args, **kwargs)
        else:
            raise TypeError('pairwise distance `{}` not recognized'.format(pairwise_distance))

    else:
        raise TypeError('weight `{}` not recognized'.format(weight))

def RecursiveSafeSubtractor(*args, weight='TransverseMomentum', **kwargs):

    if weight == 'TransverseMomentum':
        return RecursiveSafeSubtractorTransverseMomentum(*args, **kwargs)
    elif weight == 'TransverseEnergy':
        return RecursiveSafeSubtractorTransverseEnergy(*args, **kwargs)
    elif weight == 'Energy':
        return RecursiveSafeSubtractorEnergy(*args, **kwargs)
    elif weight == 'Momentum':
        return RecursiveSafeSubtractorMomentum(*args, **kwargs)
    else:
        raise TypeError('weight `{}` not recognized'.format(weight))

%}