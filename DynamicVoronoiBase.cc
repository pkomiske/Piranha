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

#include <CGAL/Polygon_2.h>

#include "DynamicVoronoiBase.hh"

BEGIN_PIRANHA_NAMESPACE

using Point = DynamicVoronoiBase::Point;
using K = DynamicVoronoiBase::K;

// return a decription of this object
std::string DynamicVoronoiBase::description() const {
  std::ostringstream oss;
  oss << name() << '\n'
      << "  - Tracking region areas\n";

  if (track_emds())
  oss << "  - Tracking region EMDs\n";
  if (track_intersection_vertices_)
  oss << "  - Tracking intersection vertices\n";

  oss << '\n'
      << parameters();

  return oss.str();
}

void DynamicVoronoiBase::reset() {

  // clear triangulation and set infinite vertex id
  triangulation_.clear();
  triangulation_.infinite_vertex()->id() = INFINITE_VERTEX_ID;

  // clear unordered_map of Voronoi vertices
  // use the fact that points in 2D Voronoi diagrams have on average 6 neighbors
  voronoiVerts_.clear();
  voronoiVerts_.reserve(6*nInitPoints_);

  // initialize number of points and area/emd vectors
  coincidences_.resize(nInitPoints_);
  voronoiAreas_.resize(nInitPoints_);
  if (track_emds_) voronoiEMDDensities_.resize(nInitPoints_);
}

void DynamicVoronoiBase::construct_triangulation() {

  // ensure that delaunayVerts_ can hold a vertex for each indexedPoint
  delaunayVerts_.resize(indexedPoints_.size());

  // spatially sort points
  CGAL::spatial_sort(indexedPoints_.begin(), indexedPoints_.end(), searchTraits_);
  
  // get vertex handles
  // nPrimaryDelaunayVerts_ will for this loop be the number of total Delaunay verts
  // this is corrected at the end by subtracting nNonPrimaryVerts
  Face_handle f;
  unsigned nNonPrimaryVerts(0);
  nPrimaryDelaunayVerts_ = 0;
  for (const auto & p : indexedPoints_) {

    // insert vertex, update face, and store vertex handle
    Vertex_handle vh(triangulation_.insert(p.first, f));
    f = vh->face();
    delaunayVerts_[p.second] = vh;

    // determine if point was primary, coincident (size of triangulation did not increase)
    bool primary(unsigned(p.second) < nInitPoints_),
         coincident(triangulation_.number_of_vertices() == nPrimaryDelaunayVerts_);

    // update vertex id if not coincident
    if (!coincident) {
      vh->id() = p.second;
      nPrimaryDelaunayVerts_++;
      if (primary) coincidences_[p.second] = p.second;
      else nNonPrimaryVerts++;
    }

    // handle coincidence
    else {

      // only track coincidence chains among primary vertices
      // uses the same scheme as DnnPlace.cc from fastjet
      // set the area and emd to be zero for the coincident point
      if (primary) {
        coincidences_[p.second] = coincidences_[vh->id()];
        coincidences_[vh->id()] = p.second;
        voronoiAreas_[p.second] = 0;
        if (track_emds_) voronoiEMDDensities_[p.second] = 0;
      }
    }
  }

  // correct nPrimaryDelaunayVerts_ for the non-primary vertices we inserted
  nPrimaryDelaunayVerts_ -= nNonPrimaryVerts;

  // initialize regions (use same order as before to take advantage of locality)
  // avoid recomputation for coincident vertices
  for (const auto & p : indexedPoints_) {

    // check that vertex is has an associated region
    if (vertex_is_primary_and_active(p.second)) 
      process_region(delaunayVerts_[p.second]);
  }

#ifdef PIRANHA_DEBUG
  // check that areas sum properly
  double area(0);
  for (double a : voronoiAreas_) area += a;
  if (fabs(area - total_area()) > 1e-8 && area > 1e-8) {
    std::ostringstream oss;
    oss << std::setprecision(10) << "Total area differs from expected by " << area - total_area();
    std::cerr << __FILE__ << " - line " << __LINE__ << ": " << oss.str() << std::endl;
  }
#endif
}

// get the neighbors of a particular vertex
std::vector<int> DynamicVoronoiBase::neighbors(unsigned i) const {
  
  // check for valid removal
#ifdef PIRANHA_DEBUG
  if (!vertex_is_primary_and_active(i))
    THROW_PIRANHA_ERROR("Point to be removed " + std::to_string(i) + "is not primary and/or active");
#endif

  // get vertex handle of this point
  Vertex_handle vh(delaunayVerts_[i]);
  (void) vh; // supress unused variable warning
#ifdef PIRANHA_DEBUG
  if (vh == nullptr)
    THROW_PIRANHA_ERROR("Invalid vertex handle");
#endif

  // acquire neighbor indices
  std::vector<int> nbs;
  nbs.reserve(16);
  Vertex_circulator vc(triangulation_.incident_vertices(delaunayVerts_[i])), done(vc);
  if (vc != nullptr)
    do nbs.push_back(vc->id());
    while (++vc != done);

  // shouldn't ever get here
  else THROW_PIRANHA_ERROR("Vertex unexpectedly had no incident vertices");

  return nbs;
}

// the following formula implements the EMD density contribution from a triangle
// assumed to have uniform piranhas and flowing to one of the vertices
// let the side lengths be a, b, c with c opposite the vertex where the piranhas go
// EMD_triangle_density = A(c(a+b)((a-b)^2 + c^2) + 16A^2 atanh(c/(a+b)))/(6c^3)
// where A = sqrt((a+b+c)(-a+b+c)(a-b+c)(a+b-c))/4 by Heron's formula
// let p0 be the vertex that the piranhas flow to
double DynamicVoronoiBase::emd_density_triangle(const Point & p0, const Point & p1, const Point & p2) const {
  double A(K::Triangle_2(p0, p1, p2).area()), 
         a(sqrt(dist2(p0, p1))), b(sqrt(dist2(p0, p2))), aplusb(a + b), aminusb(a - b),
         c2(dist2(p1, p2)), c(sqrt(c2));

  // some numerical checks
  if (c2 <= std::numeric_limits<double>::epsilon() || c >= aplusb) return 0;

  // implement formula
  return A*(c*aplusb*(aminusb*aminusb + c2) + 16*A*A*atanh(c/aplusb))/(6*c2*c*R_);
}

// processes region and updates voronoi verts and areas/emds
void DynamicVoronoiBase::process_region(const Vertex_handle & vh) {

  // get vertex id
  int vhid(vh->id());

  // should never get a non-primary vertex here
  // a non-primary vertex is one that is either coincident or has an id outside of 0,...,nInitPoints_-1
#ifdef PIRANHA_DEBUG
  if (!vertex_is_primary_and_active(vhid))
    THROW_PIRANHA_ERROR("processing region for non-primary or active vertex");
#endif

  // iterate over faces and get vertices of voronoi region
  std::vector<std::pair<Point,bool>> regionVerts;
  regionVerts.reserve(10);
  Face_circulator fc(triangulation_.incident_faces(vh)), done(fc);
  bool allValid(true);

  if (fc != nullptr) do {

    // should never have an infinite face here
  #ifdef PIRANHA_DEBUG
    if (triangulation_.is_infinite(fc))
      THROW_PIRANHA_ERROR("Unexpected infinite face encountered");
  #endif

    // find/compute voronoi vertex and store it in regionVerts
    auto v(voronoiVerts_.find(fc));
    if (v == voronoiVerts_.end()) {
      Point p(triangulation_.dual(fc));
      regionVerts.push_back(std::make_pair(p, in_acceptance(p)));
      voronoiVerts_[fc] = regionVerts.back();
      allValid &= regionVerts.back().second;
    }
    else {
      regionVerts.push_back(v->second);
      allValid &= v->second.second;
    }
  } while (++fc != done);

  // shouldn't ever get here
  else THROW_PIRANHA_ERROR("Vertex unexpectedly had no incident faces");

  // if all points are valid, just create a polygon and get the area that way
  // this will be a common case so we handle it first
  CGAL::Polygon_2<K> polygon;
  const Point & point(vh->point());
  double emd(0);
  if (allValid) {

    // compute area
    for (const auto & rv : regionVerts) polygon.push_back(rv.first);
    voronoiAreas_[vhid] = polygon.area();

    // compute emds
    if (track_emds_) {
      unsigned i(0);
      for (unsigned s = regionVerts.size() - 1; i < s; i++)
        emd += emd_density_triangle(point, regionVerts[i].first, regionVerts[i+1].first);
      voronoiEMDDensities_[vhid] = emd + emd_density_triangle(point, regionVerts[i].first, regionVerts[0].first);
    }
    return;
  }

  // we have some vertices outside of the acceptance region and maybe some inside (not guaranteed)
  // traverse the vertices of the region and construct polygon and circular segment vector
  // segIntVerts: holds points of intersection of edges with acceptance bounday, the bool values
  // are whether or not they are the first vertex of the segment
  std::vector<std::pair<Point,bool>> segIntVerts;
  for (unsigned i = 0, s = regionVerts.size(); i < s; i++) {

    // get references to the two vertices under consideration
    const std::pair<Point,bool> & vi(regionVerts[i]), & vj(regionVerts[(i + 1) % s]);

    // if both points are inside of the acceptance region, we don't have to check for intersection
    if (vi.second && vj.second) {

      // the policy will be to include only the first vertex of an internal pair
      // so that they can be safely chained together
      polygon.push_back(vi.first);
      if (track_emds_)
        emd += emd_density_triangle(point, vi.first, vj.first);
      continue;
    }

    // get intersection of this segment with the acceptance region
    intersect_acceptance(vi, vj);

    // entering edge
    if (!vi.second && vj.second) {
      const Point & segIntersect(retrieve_single_intersection());

      // the boundary intersection is false since it is the second vertex of that segment
      // include only the boundary intersection in the polygon since the internal 
      // vertex will be handled when it is first
      polygon.push_back(segIntersect);
      if (track_emds_)
        emd += emd_density_triangle(point, segIntersect, vj.first);
      if (track_intersection_vertices_)
        segIntVerts.push_back(std::make_pair(segIntersect, false));
    }

    // leaving edge
    else if (vi.second && !vj.second) {
      const Point & segIntersect(retrieve_single_intersection());

      // add both the first (internal) vertex and the boundary intersection (since we won't see it again)
      polygon.push_back(vi.first);
      polygon.push_back(segIntersect);
      if (track_emds_)
        emd += emd_density_triangle(point, vi.first, segIntersect);
      if (track_intersection_vertices_)
        segIntVerts.push_back(std::make_pair(segIntersect, true));
    }

    // edge passes completely through the acceptance region
    else if (p0_valid_ || p1_valid_) {

      // the points have been oriented properly by the intersector
      polygon.push_back(p0_);
      polygon.push_back(p1_);
      if (track_emds_)
        emd += emd_density_triangle(point, p0_, p1_);
      if (track_intersection_vertices_) {
        segIntVerts.push_back(std::make_pair(p0_, false));
        segIntVerts.push_back(std::make_pair(p1_, true));
      }
    }
  }

  // store area of polygon and emd
  voronoiAreas_[vhid] = polygon.area();
  if (track_emds_)
    voronoiEMDDensities_[vhid] = emd;

  // delegate to derived class for what to do now
  if (track_intersection_vertices_)
    finish_process_region(segIntVerts, vh);
}

END_PIRANHA_NAMESPACE