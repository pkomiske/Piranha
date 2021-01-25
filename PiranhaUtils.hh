#ifndef EVENTGEOMETRY_EVENTGEOMETRYUTILS_HH
#define EVENTGEOMETRY_EVENTGEOMETRYUTILS_HH

// C++ standard library
#include <stdexcept>
#include <string>

#include "fastjet/internal/base.hh"

#ifndef VERBOSE
#define VERBOSE 1
#endif

#define DEBUG

// defined constants
#define INFINITE_VERTEX_ID -1
#define DEFAULT_VORONOI_QUANTITY -1.0
#define REMOVED_COINCIDENCE -1

// namespace for this package
#define BEGIN_PIRANHA_NAMESPACE FASTJET_BEGIN_NAMESPACE namespace contrib { namespace piranha {
#define END_PIRANHA_NAMESPACE } } FASTJET_END_NAMESPACE

BEGIN_PIRANHA_NAMESPACE

const double PI = 3.14159265358979323846;
const double TWOPI = 2*PI;
inline double phi_fix(double phi, double ref_phi) {
  double diff(phi - ref_phi);
  if (diff > PI) phi -= TWOPI;
  else if (diff < -PI) phi += TWOPI;
  return phi; 
}

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

#endif // EVENTGEOMETRY_EVENTGEOMETRYUTILS_HH
