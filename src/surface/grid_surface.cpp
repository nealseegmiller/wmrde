#include <wmrde/surface/grid_surface.h>

//to read in data from csv file
#include <iostream>
#include <fstream>

namespace wmrde
{

GridSurface::GridSurface(
    const GridVector& X,
    const GridVector& Y,
    const std::vector<Real>& Z) //elevations
{
  assert(Z.size() == X.numel*Y.numel);

  std::array<GridVector,2> XY = {X, Y};
  Vec3 invalid_value = Vec3::Constant(RealNan());
  G_ = GriddedInterpolant<Vec3,2>(XY, invalid_value);

  //set the interpolant values,
  //interp elevation and gradient of elevation for fast calculation of surface normals
  std::vector<Vec3> V(Z.size());
  for (int i = 0; i < X.numel; i++)
  {
    for (int j = 0; j < Y.numel; j++)
    {
      int idx = G_.sub2ind(i,j);
      Vec3 v;
      v[0] = Z[idx]; //elevation

      //dz/dx
      if (i+1 < X.numel) { v[1] = (Z[G_.sub2ind(i+1,j)] - v[0])/X.spacing; } //use forward difference if possible
      else if (i-1 >= 0) { v[1] = (v[0] - Z[G_.sub2ind(i-1,j)])/X.spacing; }
      //dz/dy
      if (j+1 < Y.numel) { v[2] = (Z[G_.sub2ind(i,j+1)] - v[0])/Y.spacing; } //use forward difference if possible
      else if (j-1 >= 0) { v[2] = (v[0] - Z[G_.sub2ind(i,j-1)])/Y.spacing; }

      G_[idx] = v;
    }
  }
}

//assumes comma delimited txt file of the format:
//X.llim, X.spacing, X.numel, Y.llim, Y.spacing, Y.numel, Z[0], Z[1], ..., Z[n]
GridSurface::GridSurface(const std::string& filename)
{
  // declare file stream: http://www.cplusplus.com/reference/iostream/ifstream/
  std::ifstream file ( filename.c_str() );
  std::string value;
  if ( file.is_open() )
  {
    //read in the X, Y grid vectors
    std::array<GridVector,2> XY;
    for (int i=0; i<2; i++)
    {
      //read a string until next comma: http://www.cplusplus.com/reference/string/getline/
      std::getline (file, value, ','); XY[i].llim = (Real) atof(value.c_str());
      std::getline (file, value, ','); XY[i].spacing = (Real) atof(value.c_str());
      std::getline (file, value, ','); XY[i].numel = atoi(value.c_str());
    }
    GridVector& X = XY[0];
    GridVector& Y = XY[1];

    //read the elevation data
    std::vector<Real> Z(X.numel*Y.numel);
    size_t idx = 0;
    while ( file.good() ) {
      assert( idx < Z.size() ); //check if too many elevation values in file
      std::getline ( file, value, ',' );
      Z[idx] = (Real) atof(value.c_str());
      idx++;
    }
    file.close();

    assert(idx == Z.size()); //check if all elevation values have been set

    *this = GridSurface(X,Y,Z);
  }
}

bool GridSurface::getHeight(const Vec3& pt, Real& height) const
{
  std::array<Real,2> xq = {pt[0], pt[1]};
  Vec3 v = G_.interpolate(xq);
  height = v[0];
  return std::isnan(height);
}

bool GridSurface::getDistance(const Vec3& pt, Real& distance, Vec3& normal) const
{
  std::array<Real,2> xq = {pt[0], pt[1]};
  Vec3 v = G_.interpolate(xq);
  Real height = v[0];
  normal << -v[1], -v[2], 1.0; //dz/dx, dz/dy, 1.0
  normal.normalize();

  distance = (pt[2] - height)*normal[2];

  return std::isnan(distance);
}

bool GridSurface::getNormal(const Vec3& pt, Vec3& normal) const
{
  std::array<Real,2> xq = {pt[0], pt[1]};
  Vec3 v = G_.interpolate(xq);
  normal << -v[1], -v[2], 1.0; //dz/dx, dz/dy, 1.0
  normal.normalize();

  return std::isnan(Real(normal[0]));
}

} //namespace
