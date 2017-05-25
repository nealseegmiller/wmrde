//Surface.h
//abstract Surface class
//all pure (=0) virtual functions must be defined by derived classes

#ifndef _WMRDE_SURFACE_H_
#define _WMRDE_SURFACE_H_

#include <vector>
#include <memory> //for unique_ptr
#include <wmrde/algebra/linalg3.h> //for Vec3

namespace wmrde
{


class Surface
{
/*!
 * \class Surface
 * An abstract class for representing 2.5 dimensional surfaces,
 * i.e. for which every x,y location has exactly one height.
 * Use this class for collision checking between wheels and ground in wmrde simulations.
 */
 public:

  virtual ~Surface() {};

  /*!
   * Calculate the surface height at an x,y location
   * \param pt The x,y location at which to calculate surface height.
   * \param height The surface height at the location.
   * \return true if valid
   */
  virtual bool getHeight(const Vec3& pt, Real& height) const = 0;

  /*!
   * Calculate the surface normal at an x,y location
   * \param pt The x,y location at which to calculate surface normal.
   * \param normal The surface normal vector.
   * \return true if valid
   */
  virtual bool getNormal(const Vec3& pt, Vec3& normal) const = 0;

  /*!
   * Calculate the distance of a point from the surface
   * \param pt The point to which distance from the surface is calculated.
   * \param distance The distance of the point from the surface, measured along the surface normal.
   *           positive indicates the point is above the surface, negative indicates below.
   * \param normal The surface normal vector
   * \return true if valid
   */
  virtual bool getDistance(const Vec3& pt, Real& distance, Vec3& normal) const = 0;
};

class Surfaces
{
/*!
 * \class Surfaces
 * A class to manage multiple Surface objects.
 */
 public:

  /*!
   * Copy object of derived Surface type into Surfaces
   */
  template<typename SurfaceType>
  void addSurface(const SurfaceType& surf)
  {
    data_.emplace_back(new SurfaceType(surf)); //copy constructor
  }
  int numSurfaces() const { return data_.size(); }
  void clear() { data_.clear(); }

  /*!
   * Get Surface pointer by index. will Segmentation fault if index is invalid!
   */
  const Surface* getSurfacePtr(const size_t idx) const
  {
    return data_[idx].get();
  }

  /*!
   * Calls getHeight on all Surface objects and returns the max value
   * \return The index of the Surface object that corresponds to height value. -1 if invalid
   */
  int getHeight(const Vec3& pt, Real& height) const ;

  /*!
   * Calls getDistance on all Surface objects and returns the min value,
   * as well as the corresponding surface normal for that min distance value
   * \return The index of Surface object that corresponds to the distance value. -1 if invalid
   */
  int getDistance(const Vec3& pt, Real& distance, Vec3& normal) const ;

 private:
  std::vector<std::unique_ptr<Surface> > data_; //! A vector of Surface objects to iterate over in function calls
  //use unique_ptr to avoid object slicing
};

} //namespace

#endif //_WMRDE_SURFACE_H_



