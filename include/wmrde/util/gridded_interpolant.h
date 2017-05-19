#ifndef _WMRDE_GRIDDED_INTERPOLANT_H_
#define _WMRDE_GRIDDED_INTERPOLANT_H_


#include <array>
#include <iostream>
#include <assert.h>

#include <wmrde/common.h>
#include <wmrde/util/xs_Float.h> //cast float to int

namespace wmrde
{

struct GridVectorIndex
{
  int floor;
  Real fraction;

  GridVectorIndex() :
    floor(0),
    fraction(0.0)
  {}
  GridVectorIndex(int _floor, Real _fraction) :
    floor(_floor),
    fraction(_fraction)
  {}
};

//regularly spaced grid vector
struct GridVector
{
/*!
 * \struct GridVector
 * A struct for compactly representing a regularly-spaced vector of points
 */
  Real llim; //lower limit
  Real spacing;
  int numel; //number of elements

  GridVector() :
    llim(0.0),
    spacing(0.0),
    numel(0)
  {}

  /*!
   * Constructs the grid vector
   * \param llim The lower limit on point location
   * \param spacing The spacing between points
   * \param numel The number of points
   */
  GridVector(Real _llim, Real _spacing, Real _numel) :
    llim(_llim),
    spacing(_spacing),
    numel(_numel)
  {}

  Real at(const int idx) const { return llim + spacing*idx; }
  GridVectorIndex getIndex(const Real query) const
  {
    GridVectorIndex idx;
//    idx.floor = std::floor((query-llim)/spacing); //implicit cast to int
    idx.floor = xs_FloorToInt((query-llim)/spacing); //faster?
    idx.fraction = (query - (llim + idx.floor*spacing)) / spacing; // range of [0,1)
    return idx;
  }
};

//ValueType can be float, double or a multichannel object as long as
//arithmetic operators in interp1d are defined
template <typename ValueType, int N>
class GriddedInterpolant
{
/*!
 * \class GriddedInterpolant
 * A class for performing linear interpolation on a regularly-spaced, N-dimensional grid
 */
 public:

  GriddedInterpolant() {}

  /*!
   * Constructs the gridded interpolant
   * \param X The N-dimensional grid, specified by an array of N vectors of regularly spaced points.
   * \param V The values at every grid point. Must have at least as many elements as points in grid.
   */
  GriddedInterpolant(
      const std::array<GridVector,N>& X,
      const std::vector<ValueType>& V,
      const ValueType invalid_value)
  {
    //copy to member variables
    X_ = X;
    V_ = V;
    invalid_value_ = invalid_value;

    initGridSize();
    assert(V.size() == grid_size_[N-1]); //check the size of value array
  }

  GriddedInterpolant(
      const std::array<GridVector,N>& X,
      const ValueType invalid_value)
  {
    //copy to member variables
    X_ = X;
    invalid_value_ = invalid_value;

    initGridSize();

    //initialize the value array
    //values may be set later using the [] operator
    V_.clear();
    V_.resize(grid_size_[N-1]);
  }

  ~GriddedInterpolant() {}

  void initGridSize()
  {
    //compute grid size for each dimension
    grid_size_[0] = X_[0].numel;
    for (int dim = 1; dim < N; dim++) //loop over dimension
    {
      grid_size_[dim] = grid_size_[dim-1]*X_[dim].numel;
    }
  }

  /*!
   * Convert subscripts (indices in each grid vector) to index in value vector V
   * Overloaded for 2, 3, and N dimensions
   */
  inline int sub2ind(int i, int j) const { return i + j*grid_size_[0]; }
  inline int sub2ind(int i, int j, int k) const { return i + j*grid_size_[0] + k*grid_size_[1]; }
  inline int sub2ind(const int sub[N]) const
  {
    int idx = sub[0];
    for (int i = 1; i < N; i++) { idx += sub[i]*grid_size_[i-1]; }
    return idx;
  }

  //access functions
  const GridVector& getX(const int dim) const { return X_[dim]; }

  /*!
   * Get a mutable reference to value in V
   * Compute index using sub2ind(). Will segfault if idx is out of bounds.
   */
  inline ValueType& operator[](int idx) { return V_[idx]; }
  inline ValueType operator[](int idx) const { return V_[idx]; } //const overload

  /*!
   * Interpolate on the grid at a query point
   * \param xq The N-dimensional query point
   * \return The interpolated value. NaN if out of bounds.
   */
  ValueType interpolate(
      const std::array<Real,N>& xq) const //query point
  {
    std::array<GridVectorIndex,N> inds;
    for (int i = 0; i < N; i++) //loop over dimensions
    {
      inds[i] = X_[i].getIndex(xq[i]);

//      printf("inds[%d].floor = %d, .fraction = %f\n", i, inds[i].floor, inds[i].fraction); //DEBUGGING

      //TODO, check this
      //necessary to interpolate at upper limit
      if (inds[i].floor == X_[i].numel-1 &&
          inds[i].fraction == 0.0)
      {
        inds[i].floor -= 1;
        inds[i].fraction = 1.0;
      }

      //check if out of bounds
      if (inds[i].floor < 0 ||
          inds[i].floor+1 > X_[i].numel-1)
      {
        return invalid_value_;
      }
    }
    return recursiveinterp(inds, V_.data(), N-1);
  }

 private:
  std::array<GridVector,N> X_; //the grid points
  std::vector<ValueType> V_; //value at every grid point
  std::array<int,N> grid_size_; //grid size for each dimension
  ValueType invalid_value_; //TODO, set to NaN if ValueType is a float

  inline ValueType interp1d(const Real fraction, ValueType v0, ValueType v1) const
  {
    return v0 + fraction*(v1 - v0); //arithmetic operators must be defined for ValueType
  }

  //recursive subfunction
  ValueType recursiveinterp(
      const std::array<GridVectorIndex,N>& inds, //indices of query point for each GridVector
      const ValueType* V_ptr, // a pointer into block of V data of size grid_size_[dim]
      const int dim) const //dimension to interpolate on (decrements with each recursion)
  {
    const GridVectorIndex& idx = inds[dim];

#if 1
    //TODO, faster?
    if (dim == 0)
    {
      return interp1d(idx.fraction, V_ptr[idx.floor], V_ptr[idx.floor+1]);
    }
    else
    {
      int offset = grid_size_[dim-1];
      return interp1d(idx.fraction,
          recursiveinterp(inds, V_ptr + idx.floor*offset, dim-1),
          recursiveinterp(inds, V_ptr + (idx.floor+1)*offset, dim-1));
    }
#else
    //DEBUGGING
    ValueType v0, v1;
    if (dim == 0)
    {
      v0 = V_ptr[idx.floor];
      v1 = V_ptr[idx.floor+1];
    }
    else
    {
      int offset = grid_size_[dim-1];
      v0 = recursiveinterp(inds, V_ptr + idx.floor*offset, dim-1);
      v1 = recursiveinterp(inds, V_ptr + (idx.floor+1)*offset, dim-1);
    }
    ValueType v = interp1d(idx.fraction, v0, v1);
    std::cout << std::string(2*(N-dim),' ') << "v0 = " << v0 << ", v1 = " << v1 << ", fraction = " << idx.fraction << ", v = " << v << std::endl;
    return v;
#endif

  }
};

} //namespace

#endif
