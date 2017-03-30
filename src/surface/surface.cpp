#include <limits>
#include <wmrde/surface/surface.h>

namespace wmrde
{

int Surfaces::getHeight(const Vec3& pt, Real& height) const
{
	height = -std::numeric_limits<Real>::infinity();
	int surf_idx = -1; //init to invalid
	for (size_t i=0; i < data_.size(); i++)
	{
		Real this_height;
		if (data_[i]->getHeight(pt, this_height) && //valid
		    this_height > height) //tallest
		{
		  height = this_height;
		  surf_idx = i;
		}
	}
	return surf_idx;
}


int Surfaces::getDistance(const Vec3& pt, Real& distance, Vec3& normal) const
{
  distance = std::numeric_limits<Real>::infinity();
  int surf_idx = -1;
	for (size_t i=0; i < data_.size(); i++)
	{
	  Real this_distance;
	  Vec3 this_normal;
	  if (data_[i]->getDistance(pt, this_distance, this_normal) &&
	      this_distance < distance)
	  {
	    distance = this_distance;
	    normal = this_normal;
	    surf_idx = i;
	  }
	}
	return surf_idx;
}

} //namespace





