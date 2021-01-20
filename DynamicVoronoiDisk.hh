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

#ifndef PIRANHA_DYNAMICVORONOIDISK_HH
#define PIRANHA_DYNAMICVORONOIDISK_HH

#include <iostream>

#include "DynamicVoronoiBase.hh"

BEGIN_PIRANHA_NAMESPACE

///////////////////////////////////////////////////////////////////////////////
// DynamicVoronoiDisk - Points lie inside of a circle in the plane
///////////////////////////////////////////////////////////////////////////////

class DynamicVoronoiDisk : public DynamicVoronoiBase {
private:

  // variables held by this class
  double R2_, R2over9_, piR2_, twoR_, threeR_, x0_, y0_;
  Point center_;

public:

  // exists for SWIG
  DynamicVoronoiDisk() {}

  // constructor
  DynamicVoronoiDisk(SubtractionType subtype, double R) :
    DynamicVoronoiBase(subtype, R, true),
    x0_(0), y0_(0), 
    center_(x0_, y0_)
  {
    set_R(R);
  }

  // overloads setting R
  void set_R(double R) {
    DynamicVoronoiBase::set_R(R);
    R2_ = R_*R_;
    R2over9_ = R2_/9;
    piR2_ = PI*R2_;
    twoR_ = 2*R_;
    threeR_ = 3*R_;
  }

  // function for setting the center of the circle
  void set_center(double x0, double y0) {
    if (x0 != x0_ || y0 != y0_) {
      x0_ = x0;  
      y0_ = y0;
      center_ = Point(x0_, y0_);
    }
  }

  // functions for accessing the center
  double x0() const { return x0_; }
  double y0() const { return y0_; }

  // gets the total area of the circular region
  double total_area() const { return piR2_; }

  // function deciding if point is valid
  bool valid_point(const Point & p) const { return dist2(p, center_) <= R2_; }

  // removes point and updates quantities appropriately
  RemovalResult remove(int);

// these functions overload the base class ones
protected:

  // methods to return the name and parameters of this class
  std::string name() const { return "DynamicVoronoiDisk"; }
  std::string parameters() const {
    std::ostringstream oss;
    oss << std::setprecision(16)
        << "Circular primary region\n"
        << "  Center - (" << x0() << ", " << y0() << ")\n"
        << "  Radius - " << R_ << '\n'
        << "  Area - " << total_area() << '\n';

    return oss.str();
  }

  // setup triangulation with these points
  void construct_indexed_points(const std::vector<Point> &);

  // returns true if the point is in the acceptance region, which here is the same as the primary region
  bool in_acceptance(const Point & p) const { return valid_point(p); }

  // compute the intersection of the circle with the segment defined by these points
  void intersect_acceptance(const std::pair<Point,bool> &, const std::pair<Point,bool> &);

  // finish processing a region by adding contributions from circular segments
  void finish_process_region(const std::vector<std::pair<Point,bool>> & segIntVerts, Vertex_handle vh);

private:

  // p0 and p1 should be points on the circle
  double circular_segment_area(const Point &, const Point &) const;

  // p0 is as above, p1 and p2 must lie on the circle, p1 is the source and p2 is the target
  // let d be the distance between p0 and the center of the circle
  // let the angle p0-center-p1 = theta1 and p0-center-p2 = theta2
  double emd_density_circular_wedge(const Point & p0, const Point & p1, const Point & p2) const;

  // the emd density for the entire circle
  // should give the same value as for th2 = 2pi, th1 = 0 in the above function
  double emd_density_entire_circle(const Point & p0) const;

}; // DynamicVoronoiDisk

END_PIRANHA_NAMESPACE

#endif // PIRANHA_DYNAMICVORONOIDISK_HH
