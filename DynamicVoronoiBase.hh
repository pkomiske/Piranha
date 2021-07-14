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

#ifndef PIRANHA_DYNAMICVORONOIBASE_HH
#define PIRANHA_DYNAMICVORONOIBASE_HH

// C++ standard library
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

// CGAL Delaunay headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_id_2.h>

// CGAL sorting headers
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_2.h>
#include <CGAL/property_map.h>

#include "PiranhaUtils.hh"

BEGIN_PIRANHA_NAMESPACE

///////////////////////////////////////////////////////////////////////////////
// Helper functions
///////////////////////////////////////////////////////////////////////////////

// function to compute squared euclidean distance between two points
template<class Point>
inline double dist2(const Point & p0, const Point & p1) {
  double dx(p0.x() - p1.x()), dy(p0.y() - p1.y());
  return dx*dx + dy*dy;
}

///////////////////////////////////////////////////////////////////////////////
// DynamicVoronoiBase - Interfaces with CGAL and useful data structures
///////////////////////////////////////////////////////////////////////////////

// base class for dynamic voronoi instantiations
class DynamicVoronoiBase {
public:

  // CGAL kernel typedef
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

  // CGAL triangulation typedefs
  typedef CGAL::Triangulation_vertex_base_with_id_2<K> Tvb2;
  typedef CGAL::Triangulation_face_base_2<K> Tfb2;
  typedef CGAL::Triangulation_data_structure_2<Tvb2,Tfb2> Tds2;
  typedef CGAL::Delaunay_triangulation_2<K,Tds2> Triangulation;

  //typedef typename K::Point_2 Point;
  typedef Triangulation::Face_handle Face_handle;
  typedef Triangulation::Face_circulator Face_circulator;
  typedef Triangulation::Vertex_handle Vertex_handle;
  typedef Triangulation::Vertex_circulator Vertex_circulator;
  typedef Triangulation::Point Point;

  typedef CGAL::First_of_pair_property_map<std::pair<Point,int>> Property_map;
  typedef CGAL::Spatial_sort_traits_adapter_2<K,Property_map> Search_traits_2;
  typedef std::vector<std::pair<Point,int>> IndexedPointsVec;

  // allow IVSubtractorBase to access class internals
  //template<class T>
  //friend class IteratedVoronoiSubtractorBase;

protected:

  // CGAL objects
  Triangulation triangulation_;
  Search_traits_2 searchTraits_;

  // for holding how many points there are
  unsigned nInitPoints_, nPrimaryDelaunayVerts_;

  // data storage
  IndexedPointsVec indexedPoints_;
  std::vector<Vertex_handle> delaunayVerts_;
  std::vector<double> voronoiAreas_, voronoiEMDDensities_;
  double total_area_;

  // map from finite faces to their circumcenters, for caching voronoi vertices
  std::unordered_map<Face_handle,std::pair<Point,bool>> voronoiVerts_;

  // vector to store pairs of coinciding vertices
  std::vector<int> coincidences_;

  // the R parameter in the EMD definition
  double R_;

  // whether or not we will be tracking EMDs along with areas
  SubtractionType subtype_;
  bool track_emds_, track_intersection_vertices_;

  // intersection points
  Point p0_, p1_;

  // validity of the intersection points
  bool p0_valid_, p1_valid_;

public:

  // exists for swig
  DynamicVoronoiBase() {}

  // explicit constructor
  DynamicVoronoiBase(SubtractionType subtype, double R, bool track_intersection_vertices) :
    track_intersection_vertices_(track_intersection_vertices)
  {
    set_subtraction_type(subtype);
    set_R(R);
  }

  // destructor
  virtual ~DynamicVoronoiBase() {}

  // return a decription of this object
  std::string description() const;

  // setter functions
  void set_subtraction_type(SubtractionType subtype) {
    subtype_ = subtype;
    track_emds_ = (subtype_ == SubtractionType::EMD) || 
                  (subtype_ == SubtractionType::AreaTrackEMD);
  }
  virtual void set_R(double R) {
    if (R <= 0)
      throw std::invalid_argument("R must be positive");
    R_ = R;
  }

  // getter functions
  SubtractionType subtraction_type() const { return subtype_; }
  bool track_emds() const { return track_emds_; }
  double R() const { return R_; }

  // processes and stores points
  void operator()(const std::vector<std::pair<double, double>> & coords) {
    std::vector<Point> points(coords.size());
    unsigned i(0);
    for (const auto & xy : coords)
      points[i++] = Point(xy.first, xy.second);
    operator()(points);
  }
  void operator()(const std::vector<Point> & points) {

    // set the number of points
    nInitPoints_ = points.size();

    // reset some data structures
    reset();

    // construct the indexed points
    construct_indexed_points(points);

    // construct triangulation
    construct_triangulation();
  }

  // function to determine if point is valid
  virtual bool valid_point(const Point & p) const = 0;
  bool valid_point(double x, double y) const { return valid_point(Point(x, y)); }

  // get number of unique Delaunay vertices that were created
  unsigned number_of_primary_delaunay_vertices() const { return nPrimaryDelaunayVerts_; }
  
  // area and EMD access functions
  double total_area() const { return total_area_; }
  const std::vector<double> & areas() const { return voronoiAreas_; }
  double area(unsigned i) const { return voronoiAreas_[i]; }

  // access coincidences
  const std::vector<int> & coincidences() const { return coincidences_; }

  // EMD density is EMD/rho. It has units of [length]^2.
  double emd_density(unsigned i) const { return voronoiEMDDensities_[i]; }

  // get neighbors of a particular vertex
  std::vector<int> neighbors(unsigned i) const;

  // - a vertex is primary if its id is in 0,...,nInitPoints_-1
  // - a vertex is active if it owns its region, that is, it is not part of
  //   a coincidence chain or if it is, its the vertex that has the vertex
  //   handle associated to it
  bool vertex_is_primary_and_active(int i) const {
    return unsigned(i) < nInitPoints_ && i >= 0 && delaunayVerts_[i]->id() == i;
  }

protected:

  // methods to return the name and parameters of the derived class
  virtual std::string name() const = 0;
  virtual std::string parameters() const = 0;

  // functions that construct the triangulation and associated data structures
  virtual void construct_indexed_points(const std::vector<Point> &) = 0;

  // returns true if the point is in the acceptance region
  virtual bool in_acceptance(const Point & p) const = 0;

  // do intersection of a segment with the acceptance region of this object
  virtual void intersect_acceptance(const std::pair<Point,bool> &, 
                                    const std::pair<Point,bool> &) = 0;

  // processes region and updates voronoi verts and areas/emds
  void process_region(const Vertex_handle &);

  // handle any post-processing in process region, do nothing by default
  virtual void finish_process_region(const std::vector<std::pair<Point,bool>> & sivs, 
                                     Vertex_handle vh) {}

  // the EMD density of a triangular region
  double emd_density_triangle(const Point &, const Point &, const Point &) const;

  // clear incident faces to vertex from the voronoiVerts mapping
  void clear_voronoi_vert(Vertex_handle vh) {
    Face_circulator fc(triangulation_.incident_faces(vh)), done(fc);
    if (fc != nullptr) do voronoiVerts_.erase(fc); while (++fc != done);
  }

private:

  // resets some data structures
  void reset();

  // builds the triangulation from the already set indexedPoints
  void construct_triangulation();

  // assumes that exactly one point is valid and returns it
  // throws an error if this is not the case
  const Point & retrieve_single_intersection() const {
    if (p0_valid_ && !p1_valid_) return p0_;
    if (!p0_valid_ && p1_valid_) return p1_;

    int n(int(p0_valid_) + int(p1_valid_));
    throw PiranhaError("Expected one intersection but found " + std::to_string(n));
    return p0_;
  }

}; // DynamicVoronoiBase

END_PIRANHA_NAMESPACE

#endif // PIRANHA_DYNAMICVORONOIBASE_HH
