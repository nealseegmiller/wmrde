#include <wmrde/surface/grid_surface.h>

#include <wmrde/rosout.h>

//to read in data from csv file
#include <iostream>
#include <fstream>

namespace wmrde
{

Vec3 GridSurface::getPoint(const int i, const int j) const
{
  Vec3 pt;
  pt[0] = G_.getX(0).at(i);
  pt[1] = G_.getX(1).at(j);
  pt[2] = G_[G_.sub2ind(i,j)][0];
  return pt;
}

GridSurface::GridSurface() : initialized_(false) {}

GridSurface::GridSurface(
    const GridVector& X,
    const GridVector& Y,
    const std::vector<Real>& Z) //elevations
  : initialized_(false)
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
      //use central difference if possible
      int ip = std::min(i+1, X.numel-1);
      int im = std::max(0, i-1);
      int jp = std::min(j+1, X.numel-1);
      int jm = std::max(0, j-1);

      Vec3 v;
      v[0] = Z[G_.sub2ind(i,j)]; //z
      v[1] = (Z[G_.sub2ind(ip,j)] - Z[G_.sub2ind(im,j)]) / (2.0 * X.spacing); //dz/dx
      v[2] = (Z[G_.sub2ind(i,jp)] - Z[G_.sub2ind(i,jm)]) / (2.0 * Y.spacing); //dz/dy

      G_[G_.sub2ind(i,j)] = v;
    }
  }

  initialized_ = true;
}

//assumes comma delimited txt file of the format:
//X.llim, X.spacing, X.numel, Y.llim, Y.spacing, Y.numel, Z[0], Z[1], ..., Z[n]
GridSurface::GridSurface(const std::string& filename)
  : initialized_(false)
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
      ROS_DEBUG("%s: llim = %f, spacing = %f, numel = %d", i == 0 ? "X" : "Y", XY[i].llim, XY[i].spacing, XY[i].numel); //DEBUGGING
    }
    GridVector& X = XY[0];
    GridVector& Y = XY[1];

    //read the elevation data
    std::vector<Real> Z(X.numel*Y.numel);
    size_t idx = 0;
    while ( file.good() )
    {
      assert( idx < Z.size() ); //check if too many elevation values in file
      std::getline ( file, value, ',' );
      Z[idx] = (Real) atof(value.c_str());
//      ROS_DEBUG("Z[%zu] = %f", idx, Z[idx]); //DEBUGGING
      idx++;
    }
    file.close();

    assert(idx == Z.size()); //check if all elevation values have been set

    *this = GridSurface(X,Y,Z);
    initialized_ = true;
  }
  else
  {
    ROS_ERROR("Failed to open file '%s'", filename.c_str());
  }
}

bool GridSurface::getHeight(const Vec3& pt, Real& height) const
{
  std::array<Real,2> xq = {pt[0], pt[1]};
  Vec3 v = G_.interpolate(xq);
  height = v[0];
  return !std::isnan(height);
}

bool GridSurface::getDistance(const Vec3& pt, Real& distance, Vec3& normal) const
{
  std::array<Real,2> xq = {pt[0], pt[1]};
  Vec3 v = G_.interpolate(xq);
  Real height = v[0];
  normal << -v[1], -v[2], 1.0; //dz/dx, dz/dy, 1.0
  normal.normalize();

  distance = (pt[2] - height)*normal[2];

  return !std::isnan(distance);
}

bool GridSurface::getNormal(const Vec3& pt, Vec3& normal) const
{
  std::array<Real,2> xq = {pt[0], pt[1]};
  Vec3 v = G_.interpolate(xq);
  normal << -v[1], -v[2], 1.0; //dz/dx, dz/dy, 1.0
  normal.normalize();

  return !std::isnan(Real(normal[0]));
}

} //namespace
