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

#ifndef PIRANHA_DYNAMICVORONOICYLINDER_HH
#define PIRANHA_DYNAMICVORONOICYLINDER_HH

#include <CGAL/Aff_transformation_2.h>

#include "DynamicVoronoiBase.hh"

BEGIN_PIRANHA_NAMESPACE

///////////////////////////////////////////////////////////////////////////////
// DynamicVoronoiCylinder - Points lie on a cylinder (period assumed in y direction)
///////////////////////////////////////////////////////////////////////////////

class DynamicVoronoiCylinder : public DynamicVoronoiBase {
private:

  typedef CGAL::Aff_transformation_2<K> Transformation;

  // coordinate of the region
  double xmin_, ymin_, xmax_, ymax_, xsize_, ysize_;

  // additional numbers of points
  unsigned two_nInitPoints_, three_nInitPoints_;

  // rectangles representing the primary region and the acceptance region
  K::Iso_rectangle_2 primary_region_, acceptance_;

  // transformations to shift points down and up
  Transformation translateDown_, translateUp_;
  
  // points added to the ends to ensure the relevant part of the triangulation is finite
  Point leftEndCap_, rightEndCap_;

public:

  // exists for swig
  DynamicVoronoiCylinder() {}

  // Cylinder will be periodic in the y direction
  DynamicVoronoiCylinder(SubtractionType subtype, double R,
                         double xmin, double ymin,
                         double xmax, double ymax);

  // gets the total area of the rectangle
  double total_area() const { return primary_region_.area(); }

  // here, point is valid if it is in the primary region
  bool valid_point(const Point & p) const { return primary_region_.has_on_bounded_side(p); }

  // get neighbors (only in the primary region) of given point
  std::vector<int> neighbors(unsigned i) const;

  // removes point and updates quantities appropriately
  RemovalResult remove(int i);

// these functions overload the base class ones
protected:

  // methods to return the name and parameters of this class
  std::string name() const { return "DynamicVoronoiCylinder"; }
  std::string parameters() const {
    std::ostringstream oss;
    oss << std::setprecision(16)
        << "Rectangular primary region - periodic in vertical direction\n"
        << "  Lower-left - (" << xmin_ << ", " << ymin_ << ")\n"
        << "  Upper-right - (" << xmax_ << ", " << ymax_ << ")\n"
        << "  Xsize - " << xsize_ << '\n'
        << "  Ysize - " << ysize_ << '\n'
        << "  Area - " << total_area() << '\n';

    return oss.str();
  }

  // setup triangulation with these points
  void construct_indexed_points(const std::vector<Point> &);

  // returns true if the point is in the acceptance region, which here means it has a valid x coordinate
  bool in_acceptance(const Point & p) const {
    return xmin_ <= p.x() && p.x() <= xmax_;
  }

  // compute the intersection of the acceptance rectangle with the segment defined by these points
  void intersect_acceptance(const std::pair<Point,bool> &, const std::pair<Point,bool> &);

}; // DynamicVoronoiCylinder

END_PIRANHA_NAMESPACE

#endif // PIRANHA_DYNAMICVORONOICYLINDER_HH