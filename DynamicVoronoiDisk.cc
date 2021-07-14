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

#include <cassert>

// boost headers for elliptic integrals
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

#include "DynamicVoronoiDisk.hh"

BEGIN_PIRANHA_NAMESPACE

void DynamicVoronoiDisk::construct_indexed_points(const std::vector<Point> & points) {

  // resize vectors
  indexedPoints_.resize(nInitPoints_ + 4);

  // create corners, using 3R should guarantee that none of the circle is assigned 
  // to the corners assuming that there is at least one point in the circle
  // insert corners at the end of indexedPoints_
  indexedPoints_[nInitPoints_]   = std::make_pair(Point(x0() - threeR_, y0() - threeR_), nInitPoints_);
  indexedPoints_[nInitPoints_+1] = std::make_pair(Point(x0() + threeR_, y0() - threeR_), nInitPoints_+1);
  indexedPoints_[nInitPoints_+2] = std::make_pair(Point(x0() + threeR_, y0() + threeR_), nInitPoints_+2);
  indexedPoints_[nInitPoints_+3] = std::make_pair(Point(x0() - threeR_, y0() + threeR_), nInitPoints_+3);

  // index points
  for (unsigned i = 0; i < nInitPoints_; i++)
    indexedPoints_[i] = std::make_pair(points[i], i);
}

void DynamicVoronoiDisk::intersect_acceptance(const std::pair<Point,bool> & start,
                                              const std::pair<Point,bool> & end) {

  // reset validities
  p0_valid_ = p1_valid_ = false;

  const Point & s(start.first), & e(end.first);
  K::Line_2 line(s, e);
  double a(line.a()), b(line.b()), c(line.c()), a2(a*a), b2(b*b);

  // compute argument of square root and check that it's positive
  double ax0(a*x0()), by0(b*y0()), 
         tmp0(ax0 + by0 + c), tmp1(a2 + b2), 
         tmp2(tmp1*R2_ - tmp0*tmp0);

  // line doesn't intersect at all
  if (tmp2 < 0.) return;

  // compute intermediate values
  double sqrtval(std::sqrt(tmp2)), tmp3(b2*x0() - a*(c + by0)), tmp4(std::fabs(b)*sqrtval),
         tmp5(-b*(c + ax0) + a2*y0()), tmp6(std::signbit(b) ? a*sqrtval : -a*sqrtval),

         // compute intersection points
         xplus((tmp3 + tmp4)/tmp1), yplus((tmp5 + tmp6)/tmp1),
         xminus((tmp3 - tmp4)/tmp1), yminus((tmp5 - tmp6)/tmp1),

         // compute xrange and yrange
         xmin(std::min(s.x(), e.x())), xmax(std::max(s.x(), e.x())),
         ymin(std::min(s.y(), e.y())), ymax(std::max(s.y(), e.y()));

  // results to return
  p0_ = Point(xplus, yplus);
  p1_ = Point(xminus, yminus);

  // check that xplus is in xrange and yplus is in yrange
  if (xmin <= xplus && xplus <= xmax && ymin <= yplus && yplus <= ymax) 
    p0_valid_ = true;

  // check that xminus is in xrange and yminus is in yrange
  if (xmin <= xminus && xminus <= xmax && ymin <= yminus && yminus <= ymax) 
    p1_valid_ = true;

  // orient points if both were valid intersections
  if (p0_valid_ && p1_valid_) {

    // if they're not aligned, swap p0 and p1
    if ((e.x() - s.x())*(xminus - xplus) + (e.y() - s.y())*(yminus - yplus) < 0) {
      Point p2(p0_);
      p0_ = p1_;
      p1_ = p2;
    }
  }
}

// finish processing a region by adding contributions from circular segments
void DynamicVoronoiDisk::finish_process_region(const std::vector<std::pair<Point,bool>> & segIntVerts, Vertex_handle vh) {

  unsigned s(segIntVerts.size()), vhid(vh->id());
  const Point & point(vh->point());

  // if we had zero intersections of any kind, there will be no points
  if (segIntVerts.size() == 0) {
    voronoiAreas_[vhid] = total_area_;
    if (track_emds_)
      voronoiEMDDensities_[vhid] = emd_density_entire_circle(point);
    return;
  }

  // get the area of each of the circular segments
  double a(0);
  for (unsigned i = (segIntVerts[0].second ? 0 : 1); i < s; i += 2)
    a += circular_segment_area(segIntVerts[i].first, segIntVerts[(i + 1) % s].first);
  voronoiAreas_[vhid] += a; 

  // get emd contribution from each of the circular segments
  double emd(0);
  if (track_emds_) {
    for (unsigned i = (segIntVerts[0].second ? 0 : 1); i < s; i += 2) {
      const Point & p1(segIntVerts[i].first), & p2(segIntVerts[(i + 1) % s].first);
      if (CGAL::orientation(point, p1, p2) == CGAL::LEFT_TURN)
        emd += emd_density_circular_wedge(point, p1, p2);
      else emd += emd_density_entire_circle(point) - emd_density_circular_wedge(point, p2, p1);
    }
    voronoiEMDDensities_[vhid] += emd;
  }
}

RemovalResult DynamicVoronoiDisk::remove(int i) {

  // check for valid removal
#ifdef PIRANHA_DEBUG
  if (!vertex_is_primary_and_active(i))
    THROW_PIRANHA_ERROR("Point to be removed " + std::to_string(i) + "is not primary and/or active");
#endif

  // get vertex handle of this point
  Vertex_handle vh(delaunayVerts_[i]);

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

    // check if this is the index that the vertex handle is associated to currently, update to next_i
    // update the quantity vectors for next_i to reflect those for the associated index
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

    return RemovalResult::Coincidence;
  }

  // set area and emd density for the removed particle
  voronoiAreas_[i] = 0;
  if (track_emds_)
    voronoiEMDDensities_[i] = 0;

  // get incident faces and remove them from our storing of the voronoiVerts
  clear_voronoi_vert(vh);

  // store neighbors of this point for later updating their regions
  // circulate over incident vertices and record their id
  std::vector<Vertex_handle> primary_nbs;
  Vertex_circulator vc(triangulation_.incident_vertices(vh)), done(vc);
  if (vc != nullptr) {
    do {
      if (unsigned(vc->id()) < nInitPoints_)
        primary_nbs.push_back(vc);
    } while (++vc != done);
  }

  // shouldn't ever get here
  else THROW_PIRANHA_ERROR("Vertex unexpectedly had no incident vertices");

  // actually remove vertex
  triangulation_.remove(vh);

  // invalidate vertex handle for this vertex
  delaunayVerts_[i] = nullptr;

  // update each primary neighbor's region
  for (Vertex_handle nb : primary_nbs)
    process_region(nb);

#ifdef PIRANHA_DEBUG
  // check that areas sum properly
  double area(0);
  for (double a : voronoiAreas_) area += a;
  if (fabs(area - total_area()) > 1e-8 && area > 1e-8) {
    std::ostringstream oss;
    oss << std::setprecision(10) << "Total area differs from piR^2 by " << area - total_area();
    std::cerr << __FILE__ << " - line " << __LINE__ << ": " << oss.str() << std::endl;
  }
#endif

  return RemovalResult::Success;
}

// p0 and p1 should be points on the circle
double DynamicVoronoiDisk::circular_segment_area(const Point & p0, const Point & p1) const {
  double aover2R(std::min(sqrt(dist2(p0, p1))/twoR_, 1.0)),
         segment_area(R2_ * (asin(aover2R) - aover2R * sqrt(1 - aover2R*aover2R)));

  // determine if we need adjust for flipped region
  bool flip(CGAL::orientation(p0, p1, center_) == CGAL::RIGHT_TURN);
#ifdef PIRANHA_DEBUG
  assert(K::Line_2(p0, p1).has_on_negative_side(center_) == flip);
#endif
  if (flip) segment_area = total_area_ - segment_area;
  return std::max(segment_area, 0.0);
}

// NOTE: This function explicitly assumes that the points p0, p1, p2 are oriented counter-clockwise!
double DynamicVoronoiDisk::emd_density_circular_wedge(const Point & p0, const Point & p1, const Point & p2) const {

  double dp0x(p0.x() - x0()), dp0y(p0.y() - y0()),
         dp1x(p1.x() - x0()), dp1y(p1.y() - y0()),
         dp2x(p2.x() - x0()), dp2y(p2.y() - y0()),
         d2(dp0x*dp0x + dp0y*dp0y),
         d(sqrt(d2));

  // ensure the proper signs for th1 and sinth1
  bool p1below(dp0x*dp1y < dp1x*dp0y), p2below(dp0x*dp2y < dp2x*dp0y);
  if (d > std::numeric_limits<double>::epsilon()) {
    double dR(d*R_), x(std::min(d/R_, 1.0)), x2(std::min(d2/R2_, 1.0)),
           costh1(std::max(std::min((dp0x*dp1x + dp0y*dp1y)/dR, 1.0), -1.0)), // dot product
           costh2(std::max(std::min((dp0x*dp2x + dp0y*dp2y)/dR, 1.0), -1.0)), // dot product
           cos2th1(costh1*costh1), cos2th2(costh2*costh2),
           sin2th1(1 - cos2th1), sin2th2(1 - cos2th2),
           sinth1(sqrt(sin2th1)), sinth2(sqrt(sin2th2)),
           th1(acos(costh1)), th2(acos(costh2));

    // if p2 is "below" the c-p0 line, need to update th2 and sinth2
    // always put th2 in [0,2pi]
    if (p2below) {
      sinth2 *= -1;
      th2 = TWOPI - th2;
    }

    // if p1 is "below" the c-p0 line, need to update th1 and sinth1
    if (p1below) {
      sinth1 *= -1;

      // th1 may need to be in [-pi,pi] or [0,2pi] to avoid a branch cut passing through the integration region
      // if p2 is below also, then branch cut is at 0, otherwise it's at pi
      th1 = (p2below ? TWOPI - th1 : -th1);
    }

  #ifdef PIRANHA_DEBUG
    if (th2 < th1) {
      std::cerr << __FILE__ << " - line " << __LINE__ << ": "
                << "Incorrect angles encountered\n"
                << "center: " << center_ << '\n'
                << "p0: " << p0 << '\n'
                << "p1: " << p1 << '\n'
                << "p2: " << p2 << '\n'
                << "th1 th2: " << th1 << ' ' << th2 << std::endl;
      //THROW_PIRANHA_ERROR("Incorrect angles");
    }
  #endif

    // implement formula
    return R2over9_*((x2 + 7)*(boost::math::ellint_2(x, th2) - boost::math::ellint_2(x, th1))
                     -4*(1 - x2)*(boost::math::ellint_1(x, th2) - boost::math::ellint_1(x, th1))
                     -9*x*(sinth2 - sinth1)
                     +4*x2*(sinth2*costh2*sqrt(1 - x2*sin2th2) - sinth1*costh1*sqrt(1 - x2*sin2th1))
                     +x*x2*(sin2th2*sinth2 - sin2th1*sinth1 - 3*(cos2th2*sinth2 - cos2th1*sinth1))
                    );
  }
  else {
    double th1(atan2(dp1y, dp1x)), th2(atan2(dp2y, dp2x));

    // check that there is no cut in between th1 and th2
    // this occurrs when p0vec x p1vec > 0 and p0vec x p2vec < 0
    if (!p1below && p2below) th2 += TWOPI;

  #ifdef PIRANHA_DEBUG
    if (th2 < th1)
      std::cerr << __FILE__ << " - line " << __LINE__ << ": Incorrect angles encountered" << std::endl;
  #endif

    // formula in the limit d -> 0 is very simple
    return 3*R2over9_*(th2 - th1);
  }
}

double DynamicVoronoiDisk::emd_density_entire_circle(const Point & p0) const {
  double d2(dist2(p0, center_)), x2(std::min(d2/R2_, 1.0)), x(sqrt(x2));
  return 4*R2over9_*((7 + x2)*boost::math::ellint_2(x) + 4*(x2 - 1)*boost::math::ellint_1(x));
}

END_PIRANHA_NAMESPACE