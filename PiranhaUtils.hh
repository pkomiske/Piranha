// Piranha Package
//
//  Questions/comments? pkomiske@mit.edu
//
//  Copyright (c) 2019-2021
//  Patrick T. Komiske III
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

#ifndef PIRANHA_PIRANHAUTILS_HH
#define PIRANHA_PIRANHAUTILS_HH

// C++ standard library
#include <stdexcept>
#include <string>

#ifndef VERBOSE
# define VERBOSE 1
#endif

#define DEBUG

// defined constants
#define INFINITE_VERTEX_ID -1
#define DEFAULT_VORONOI_QUANTITY -1.0
#define REMOVED_COINCIDENCE -1

// namespace for this package
#define BEGIN_PIRANHA_NAMESPACE namespace fastjet { namespace contrib { namespace piranha {
#define END_PIRANHA_NAMESPACE } } }

BEGIN_PIRANHA_NAMESPACE

const double PI = 3.14159265358979323846;
const double TWOPI = 2*PI;

#ifndef PIRANHA_USE_FJCORE
  // this function is included in pyfjcore
  inline double phi_fix(double phi, double ref_phi) {
    double diff(phi - ref_phi);
    if (diff > PI) phi -= TWOPI;
    else if (diff < -PI) phi += TWOPI;
    return phi; 
  }
#endif

// enum to hold type of quantity to be monitored during grooming
enum class SubtractionType { Area, AreaTrackEMD, EMD };
enum class RhoSubtractionMode { Additive, Fractional };

// internal enum to hold result of removing a point
enum class RemovalResult { Success, Coincidence };

// custom error class
class PiranhaError : public std::runtime_error {
public:

  PiranhaError(const std::string & message) : std::runtime_error(message) {}
  PiranhaError(const std::string & message, const std::string & file, int line) : 
    std::runtime_error(file + " - line " + std::to_string(line) + ": " + message) {}
};

#define THROW_PIRANHA_ERROR(MSG) throw PiranhaError(MSG, __FILE__, __LINE__)

END_PIRANHA_NAMESPACE

#endif // PIRANHA_PIRANHAUTILS_HH
