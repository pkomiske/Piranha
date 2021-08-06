// -*- C++ -*-
%define PIRANHA_DOCSTRING
"# Piranha FastJet Contrib

Python interface to the Piranha FastJet contrib package.
"
%enddef

%module(docstring=PIRANHA_DOCSTRING) piranha

#define PIRANHA_NAMESPACE fastjet::contrib::piranha

// C++ standard library wrappers
%include <exception.i>
%include <std_list.i>
%include <std_pair.i>
%include <std_string.i>
%include <std_vector.i>

// ensure pyfjcore usage is carried along
#ifdef PIRANHA_USE_PYFJCORE
# define EVENTGEOMETRY_USE_PYFJCORE
#endif

%pythonbegin %{
import eventgeometry
%}

// import eventgeometry (handles importing pyfjcore or fastjet)
%import eventgeometry/swig/eventgeometry.i

// converts fastjet::Error into a FastJetError Python exception
FASTJET_ERRORS_AS_PYTHON_EXCEPTIONS(piranha)

// include FastJetError in python module
%pythoncode %{
  from eventgeometry import FastJetError
%}

// define as macro for use in contrib files
%define PIRANHA_ERRORS_AS_PYTHON_EXCEPTIONS(module)
%{
// Python class for representing errors from FastJet
static PyObject * PiranhaError_;
%}

// this gets placed in the SWIG_init function
%init %{
  // setup error class
  char * msg_piranha = (char *) calloc(strlen(`module`)+15, sizeof(char));
  strcpy(msg_piranha, `module`);
  strcat(msg_piranha, ".PiranhaError");
  PiranhaError_ = PyErr_NewException(msg_piranha, NULL, NULL);
  Py_INCREF(PiranhaError_);
  if (PyModule_AddObject(m, "PiranhaError", PiranhaError_) < 0) {
    Py_DECREF(m);
    Py_DECREF(PiranhaError_);
  }
%}
%enddef

PIRANHA_ERRORS_AS_PYTHON_EXCEPTIONS(piranha)

// include headers in source file
%{
#ifndef SWIG
# define SWIG
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
%template(vectorSubtractionHistory) std::vector<PIRANHA_NAMESPACE::SubtractionHistory>;
%template(vectorPairDouble) std::vector<std::pair<double, double>>;
%template(vectorInt) std::vector<int>;

// ignored classes/methods
%ignore PIRANHA_NAMESPACE::PiranhaError;

// exception handling that catches PiranhaError
%exception {
  try { $action }
  catch (piranha::PiranhaError & e) {
    PyErr_SetString(PiranhaError_, e.what());
    SWIG_fail;
  }
  SWIG_CATCH_STDEXCEPT
  catch (...) {
    SWIG_exception_fail(SWIG_UnknownError, "unknown exception");
  }
}

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

// explicit templates
BEGIN_PIRANHA_NAMESPACE

#ifdef PIRANHA_USE_PYFJCORE

  %extend IteratedVoronoiSubtractorCylinder {
    PIRANHA_PSEUDOJET_CONTAINER operator()(const PseudoJetContainer & pjc) {
      return $self->operator()(pjc.as_vector());
    }
    PIRANHA_PSEUDOJET_CONTAINER operator()(const PseudoJetContainer & pjc) {
      return $self->operator()(pjc.as_vector());
    }
  }

  %extend IteratedVoronoiSubtractorDisk {
    PIRANHA_PSEUDOJET_CONTAINER operator()(const PseudoJetContainer & pjc) {
      return $self->operator()(pjc.as_vector());
    }
    PIRANHA_PSEUDOJET_CONTAINER operator()(const PseudoJetContainer & pjc, const PseudoJet & center) {
      return $self->operator()(pjc, center);
    }
  }

  %extend OptimalTransportSubtractor {
    PIRANHA_PSEUDOJET_CONTAINER operator()(const PseudoJetContainer & pjc,
                                           const PseudoJet & offset = PtYPhiM(0, 0, 0),
                                           double min_weight_to_keep = 1e-14) {
      return $self->operator()(pjc, offset, min_weight_to_keep);
    }
  }

#endif // PIRANHA_USE_PYFJCORE

  %extend GhostGridDisk { ADD_REPR_FROM_DESCRIPTION }
  %extend GhostGridRectangle { ADD_REPR_FROM_DESCRIPTION }
  %extend IteratedVoronoiSubtractorCylinder { ADD_REPR_FROM_DESCRIPTION }
  %extend IteratedVoronoiSubtractorDisk { ADD_REPR_FROM_DESCRIPTION }
  %extend OptimalTransportSubtractor { ADD_REPR_FROM_DESCRIPTION }

  %template(IteratedVoronoiSubtractorDiskBase) IteratedVoronoiSubtractorBase<DynamicVoronoiDisk>;
  %template(IteratedVoronoiSubtractorCylinderBase) IteratedVoronoiSubtractorBase<DynamicVoronoiCylinder>;

%define PIRANHA_OPTIMAL_TRANSPORT_SUBTRACTOR_TEMPLATE(Weight, Distance)
  %template(OptimalTransportSubtractor##Weight##Distance)
        OptimalTransportSubtractor<eventgeometry::EMD<double, eventgeometry::Weight, eventgeometry::Distance>>;
%enddef

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

END_PIRANHA_NAMESPACE

// include this after explicit IVSBase templates
%include "IteratedVoronoiSubtractor.hh"

BEGIN_PIRANHA_NAMESPACE

  %extend RecursiveSafeSubtractor {
    ADD_REPR_FROM_DESCRIPTION

    // only accept PseudoJetContainer if compiling with pyfjcore
    #ifdef PIRANHA_USE_PYFJCORE
      PIRANHA_PSEUDOJET_CONTAINER operator()(const PseudoJetContainer & pjc) {
        std::forward_list<fastjet::PseudoJet> result($self->operator()(pjc.as_vector()));
        return std::vector<fastjet::PseudoJet>(result.begin(), result.end());
      }
    #endif

    // accept vector of PseudoJets
    PIRANHA_PSEUDOJET_CONTAINER operator()(const std::vector<PseudoJet> & pjs) {
      std::forward_list<fastjet::PseudoJet> result($self->operator()(pjs));
      return std::vector<fastjet::PseudoJet>(result.begin(), result.end());
    }

    // accept single PseudoJet
    PIRANHA_PSEUDOJET_CONTAINER operator()(const fastjet::PseudoJet & jet) {
      std::forward_list<fastjet::PseudoJet> result($self->operator()(jet));
      return std::vector<fastjet::PseudoJet>(result.begin(), result.end());
    }

    PIRANHA_PSEUDOJET_CONTAINER apply(double z, double f) {
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