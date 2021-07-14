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

#include <iomanip>
#include <unordered_set>

#include <CGAL/intersections.h>

#include "DynamicVoronoiCylinder.hh"

BEGIN_PIRANHA_NAMESPACE

// Cylinder will be periodic in the y direction
DynamicVoronoiCylinder::DynamicVoronoiCylinder(SubtractionType subtype, double R,
                                               double xmin, double ymin,
                                               double xmax, double ymax) :
  DynamicVoronoiBase(subtype, R, false),

  // initialize coordinates
  xmin_(xmin), ymin_(ymin), xmax_(xmax), ymax_(ymax),
  xsize_(xmax_ - xmin_),
  ysize_(ymax_ - ymin_),

  // initialize transformations
  translateDown_(CGAL::TRANSLATION, K::Vector_2(0, -ysize_)),
  translateUp_(CGAL::TRANSLATION, K::Vector_2(0, ysize_))
{
  if (xsize_ <= 0 || ysize_ <= 0)
    throw std::invalid_argument("rectangle should have non-zero area");

  // initialize rectangles
  primary_region_ = K::Iso_rectangle_2(xmin_, ymin_, xmax_, ymax_);
  acceptance_ = K::Iso_rectangle_2(xmin_, ymin_ - 2*ysize_, xmax_, ymax_ + 2*ysize_);
  total_area_ = primary_region_.area();

  // there is some wiggle room in these choices (by design)
  double maxDist(std::sqrt(xsize_*xsize_ + ysize_*ysize_)), halfy(ymin_ + ysize_/2);
  leftEndCap_ = Point(xmin_ - maxDist, halfy);
  rightEndCap_ = Point(xmax_ + maxDist, halfy);
}

void DynamicVoronoiCylinder::construct_indexed_points(const std::vector<Point> & points) {

  // additional quantities
  two_nInitPoints_ = 2*nInitPoints_;
  three_nInitPoints_ = 3*nInitPoints_;

  // resize vectors
  indexedPoints_.resize(three_nInitPoints_ + 2);

  // give each point and index and also create the points in the upper and lower panels
  // 0,...,nInitPoints_-1 will be the given points
  // nInitPoints_,...,2*nInitPoints_-1 will be the points in the lower panel
  // 2*nInitPoints_,...,3*nInitPoints_-1 will be the points in the upper panel
  for (unsigned i = 0; i < nInitPoints_; i++) {
    unsigned downind(nInitPoints_ + i), upind(two_nInitPoints_ + i);
    indexedPoints_[i] = std::make_pair(points[i], i);
    indexedPoints_[downind] = std::make_pair(translateDown_(points[i]), downind);
    indexedPoints_[upind] = std::make_pair(translateUp_(points[i]), upind);
  }

  // add two endcap points
  indexedPoints_[three_nInitPoints_] = std::make_pair(leftEndCap_, three_nInitPoints_);
  indexedPoints_[three_nInitPoints_ + 1] = std::make_pair(rightEndCap_, three_nInitPoints_ + 1);
}

std::vector<int> DynamicVoronoiCylinder::neighbors(unsigned i) const {

  // get raw neighbors
  std::vector<int> nbs(DynamicVoronoiBase::neighbors(i));

  // adjust translated points to primary region
  for (int & nb : nbs) {
    if (0 <= nb) {
      if (nb < int(nInitPoints_)) continue;
      else if (nb < int(two_nInitPoints_)) nb -= nInitPoints_;
      else if (nb < int(three_nInitPoints_)) nb -= two_nInitPoints_;
    }
  }

  return nbs;
}

void DynamicVoronoiCylinder::intersect_acceptance(const std::pair<Point,bool> & start, 
                                                  const std::pair<Point,bool> & end) {

  // reset validities
  p0_valid_ = p1_valid_ = false;

  // get segment for these points
  K::Segment_2 seg(start.first, end.first);
  CGAL::Object intersect(CGAL::intersection(seg, acceptance_));

  // check for empty object, meaning no intersection
  if (intersect.empty())
    return;

  // cast to segment
  const K::Segment_2 * segInt(CGAL::object_cast<K::Segment_2>(&intersect));
  if (segInt) {
    p0_ = segInt->source();
    p1_ = segInt->target();

    // validity of the endpoints is determined by whether 
    // start and end were in or out of acceptance
    p0_valid_ = !start.second;
    p1_valid_ = !end.second;
  }
  else {

    // cast to point
    const Point * pointInt(CGAL::object_cast<Point>(&intersect));
    if (pointInt) {
      p0_ = *pointInt;
      p0_valid_ = true;
    }
    else THROW_PIRANHA_ERROR("bad casts to segment and point");
  }
}

RemovalResult DynamicVoronoiCylinder::remove(int i) {

  // check for valid removal
#ifdef PIRANHA_DEBUG
  if (!vertex_is_primary_and_active(i))
    THROW_PIRANHA_ERROR("Point to be removed " + std::to_string(i) + 
                        "is not primary and/or active");
#endif

  // get vertex handles of this point and the translated points
  int id(nInitPoints_ + i), iu(2*nInitPoints_ + i);
  Vertex_handle vh(delaunayVerts_[i]),
                vhd(delaunayVerts_[id]),
                vhu(delaunayVerts_[iu]);

  // check that the translated points own their own vertices
  // (only relevant if there was a coincidence previously)
#ifdef PIRANHA_DEBUG
  if (vhd->id() != id || vhu->id() != iu)
    THROW_PIRANHA_ERROR("Translated points do not own their own vertices");
#endif

  // get next coincidence value
  int next_i(coincidences_[i]);

  // check that we haven't already removed the point
#ifdef PIRANHA_DEBUG
  if (vh == nullptr || next_i == REMOVED_COINCIDENCE)
    THROW_PIRANHA_ERROR("Previously removed point " + std::to_string(i));
#endif

  // check if this vertex is part of a coincidence chain
  if (next_i != i) {

    // coincidences_[i] holds the "previous" coincidence, get the "next" coincident vertex
    while (coincidences_[next_i] != i)
      next_i = coincidences_[next_i];

    // - check if this is the index that the vertex handle is 
    //   associated to currently, update to next_i
    // - update the quantity vectors for next_i to reflect those 
    //   for the associated index
    if (vh->id() == i) {
      voronoiAreas_[next_i] = voronoiAreas_[i];
      if (track_emds_)
        voronoiEMDDensities_[next_i] = voronoiEMDDensities_[i];
      vh->id() = next_i;
    }

    // set area and emd density for the removed particle
    voronoiAreas_[i] = 0;
    if (track_emds_)
      voronoiEMDDensities_[i] = 0;

    // update chain of coincidences
    coincidences_[next_i] = coincidences_[i];

    // mark coincidence as removed
    coincidences_[i] = REMOVED_COINCIDENCE;

    // update translated vertices to have new index
    vhd->id() = nInitPoints_ + next_i;
    vhu->id() = 2*nInitPoints_ + next_i;

    return RemovalResult::Coincidence;
  }

  // set area and emd density for the removed particle
  voronoiAreas_[i] = 0;
  if (track_emds_)
    voronoiEMDDensities_[i] = 0;

  // clear incident faces from our storing of the voronoiVerts
  // store neighbors of this point for later updating their regions
  // use a set to avoid updating the same region twice
  std::unordered_set<Vertex_handle> primary_nbs;
  for (Vertex_handle v : {vh, vhd, vhu}) {
    clear_voronoi_vert(v);

    // circulate over the neighbors and add them if they're primary
    Vertex_circulator vc(triangulation_.incident_vertices(v)), done(vc);
    if (vc != nullptr) {
      do {
        if (unsigned(vc->id()) < nInitPoints_ && vc->id() != i)
          primary_nbs.insert(vc);
      } while (++vc != done);
    }

    // shouldn't ever get here
    else THROW_PIRANHA_ERROR("Vertex unexpectedly had no incident vertices");
  }

  // actually remove vertices
  triangulation_.remove(vh);
  triangulation_.remove(vhd);
  triangulation_.remove(vhu);

  // invalidate vertex handle for this vertex
  delaunayVerts_[i] = delaunayVerts_[id] = delaunayVerts_[iu] = nullptr;

  // update each primary neighbor's region
  for (Vertex_handle nb : primary_nbs)
    process_region(nb);

#ifdef PIRANHA_DEBUG
  // check that areas sum properly
  double area(0);
  for (double a : voronoiAreas_) area += a;
  if (fabs(area - total_area()) > 1e-8 && area > 1e-8) {
    std::ostringstream oss;
    oss << std::setprecision(10) << "Total area differs from expected by " 
        << area - total_area();
    std::cerr << __FILE__ << " - line " << __LINE__ << ": " 
              << oss.str() << std::endl;
  }
#endif

  return RemovalResult::Success;
}

END_PIRANHA_NAMESPACE