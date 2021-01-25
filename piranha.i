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
%import fastjet/pyinterface/fastjet.i

// make library aware of EventGeometry fjcontrib wrappers
%import fastjet/contrib/eventgeometry.i

// converts fastjet::Error into a FastJetError Python exception
FASTJET_ERRORS_AS_PYTHON_EXCEPTIONS(Piranha)

// include FastJetError in python module
%pythoncode {
from fastjet import FastJetError

__version__ = '1.0.0a0'
}

// turn off exception handling for now, since fastjet::Error is not thrown here
%exception;

// include headers in source file
%{
#ifndef SWIG
#define SWIG
#endif

#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "Piranha.hh"

// macros for exception handling
#define CATCH_STD_EXCEPTION catch (std::exception & e) { SWIG_exception(SWIG_SystemError, e.what()); } \
  catch (fastjet::Error & e) { PyErr_SetString(FastJetError_, e.message().c_str()); SWIG_fail; }
#define CATCH_STD_INVALID_ARGUMENT catch (std::invalid_argument & e) { SWIG_exception(SWIG_ValueError, e.what()); }
#define CATCH_STD_RUNTIME_ERROR catch (std::runtime_error & e) { SWIG_exception(SWIG_RuntimeError, e.what()); }
#define CATCH_STD_LOGIC_ERROR catch (std::logic_error & e) { SWIG_exception(SWIG_RuntimeError, e.what()); }
#define CATCH_STD_OUT_OF_RANGE catch (std::out_of_range & e) { SWIG_exception(SWIG_IndexError, e.what()); }

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

// basic exception handling for any function
%exception {
  try { $action }
  CATCH_STD_RUNTIME_ERROR
  CATCH_STD_EXCEPTION
}

// IVG exception handling
%exception PIRANHANAMESPACE::IteratedVoronoiSubtractorBase::IteratedVoronoiSubtractorBase {
  try { $action }
  CATCH_STD_INVALID_ARGUMENT
  CATCH_STD_RUNTIME_ERROR
  CATCH_STD_EXCEPTION
}
%exception PIRANHANAMESPACE::IteratedVoronoiSubtractorBase::reapply {
  try { $action }
  CATCH_STD_INVALID_ARGUMENT
  CATCH_STD_RUNTIME_ERROR
  CATCH_STD_EXCEPTION
}

// RecursiveSafeSubtractor exception handling
%exception PIRANHANAMESPACE::RecursiveSafeSubtractor::RecursiveSafeSubtractor {
  try { $action }
  CATCH_STD_INVALID_ARGUMENT
  CATCH_STD_EXCEPTION
}
%exception PIRANHANAMESPACE::RecursiveSafeSubtractor::apply {
  try { $action }
  CATCH_STD_INVALID_ARGUMENT
  CATCH_STD_RUNTIME_ERROR
  CATCH_STD_EXCEPTION
}

// OptimalTransport exception handling
%exception PIRANHANAMESPACE::GhostGridDisk::GhostGridDisk {
  try { $action }
  CATCH_STD_INVALID_ARGUMENT
  CATCH_STD_EXCEPTION
}
%exception PIRANHANAMESPACE::GhostGridRectangle::GhostGridRectangle {
  try { $action }
  CATCH_STD_INVALID_ARGUMENT
  CATCH_STD_EXCEPTION
}
%exception PIRANHANAMESPACE::OptimalTransportSubtractor::OptimalTransportSubtractor {
  try { $action }
  CATCH_STD_INVALID_ARGUMENT
  CATCH_STD_EXCEPTION
}

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

// make printable
%extend PIRANHANAMESPACE::IteratedVoronoiSubtractorBase { ADD_STR_FROM_DESCRIPTION }
%extend PIRANHANAMESPACE::GhostGridDisk { ADD_STR_FROM_DESCRIPTION }
%extend PIRANHANAMESPACE::GhostGridRectangle { ADD_STR_FROM_DESCRIPTION }
%extend PIRANHANAMESPACE::OptimalTransportSubtractor { ADD_STR_FROM_DESCRIPTION }

// explicit templates
BEGIN_PIRANHA_NAMESPACE
%template(IteratedVoronoiSubtractorDiskBase) IteratedVoronoiSubtractorBase<DynamicVoronoiDisk>;
%template(IteratedVoronoiSubtractorCylinderBase) IteratedVoronoiSubtractorBase<DynamicVoronoiCylinder>;
%template(OptimalTransportSubtractorTransverseMomentumDeltaR) OptimalTransportSubtractor<emd::EMD<emd::TransverseMomentum, emd::DeltaR>>;
%template(OptimalTransportSubtractorTransverseMomentumHadronicDot) OptimalTransportSubtractor<emd::EMD<emd::TransverseMomentum, emd::HadronicDot>>;
%template(OptimalTransportSubtractorTransverseMomentumHadronicDotMassive) OptimalTransportSubtractor<emd::EMD<emd::TransverseMomentum, emd::HadronicDotMassive>>;
%template(OptimalTransportSubtractorTransverseEnergyDeltaR) OptimalTransportSubtractor<emd::EMD<emd::TransverseEnergy, emd::DeltaR>>;
%template(OptimalTransportSubtractorTransverseEnergyHadronicDot) OptimalTransportSubtractor<emd::EMD<emd::TransverseEnergy, emd::HadronicDot>>;
%template(OptimalTransportSubtractorTransverseEnergyHadronicDotMassive) OptimalTransportSubtractor<emd::EMD<emd::TransverseEnergy, emd::HadronicDotMassive>>;
%template(OptimalTransportSubtractorMomentumEEDot) OptimalTransportSubtractor<emd::EMD<emd::Momentum, emd::EEDot>>;
%template(OptimalTransportSubtractorMomentumEEDotMassive) OptimalTransportSubtractor<emd::EMD<emd::Momentum, emd::EEDotMassive>>;
%template(OptimalTransportSubtractorMomentumEEArcLength) OptimalTransportSubtractor<emd::EMD<emd::Momentum, emd::EEArcLength>>;
%template(OptimalTransportSubtractorMomentumEEArcLengthMassive) OptimalTransportSubtractor<emd::EMD<emd::Momentum, emd::EEArcLengthMassive>>;
%template(OptimalTransportSubtractorEnergyEEDot) OptimalTransportSubtractor<emd::EMD<emd::Energy, emd::EEDot>>;
%template(OptimalTransportSubtractorEnergyEEDotMassive) OptimalTransportSubtractor<emd::EMD<emd::Energy, emd::EEDotMassive>>;
%template(OptimalTransportSubtractorEnergyEEArcLength) OptimalTransportSubtractor<emd::EMD<emd::Energy, emd::EEArcLength>>;
%template(OptimalTransportSubtractorEnergyEEArcLengthMassive) OptimalTransportSubtractor<emd::EMD<emd::Energy, emd::EEArcLengthMassive>>;
END_PIRANHA_NAMESPACE

// include this after explicit IVSBase templates
%include "IteratedVoronoiSubtractor.hh"

%extend PIRANHANAMESPACE::RecursiveSafeSubtractor {
  ADD_STR_FROM_DESCRIPTION

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

BEGIN_PIRANHA_NAMESPACE
%template(RecursiveSafeSubtractorTransverseMomentum) RecursiveSafeSubtractor<emd::TransverseMomentum>;
%template(RecursiveSafeSubtractorTransverseEnergy) RecursiveSafeSubtractor<emd::TransverseEnergy>;
%template(RecursiveSafeSubtractorEnergy) RecursiveSafeSubtractor<emd::Energy>;
%template(RecursiveSafeSubtractorMomentum) RecursiveSafeSubtractor<emd::Momentum>;
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