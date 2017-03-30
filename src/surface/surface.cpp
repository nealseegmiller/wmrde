#include <limits>
#include <wmrde/surface/surface.h>

inline void copy3(const Real src[3], Real dst[3])
{
  dst[0] = src[0];
  dst[1] = src[1];
  dst[2] = src[2];
}

namespace wmrde
{

int Surfaces::getHeight(const Real pt[3], Real& height) const
{
	height = -std::numeric_limits<Real>::infinity();
	int surf_idx = -1; //init to invalid
	for (size_t i=0; i < data_.size(); i++)
	{
		Real height_;
		if (data_[i]->getHeight(pt, height_) && //valid
		    height_ > height) //tallest
		{
		  height = height_;
		  surf_idx = i;
		}
	}
	return surf_idx;
}


int Surfaces::getDistance(const Real pt[3], Real& distance, Real normal[3]) const
{
  distance = std::numeric_limits<Real>::infinity();
  int surf_idx = -1;
	for (size_t i=0; i < data_.size(); i++)
	{
	  Real distance_;
	  Real normal_[3];
	  if (data_[i]->getDistance(pt, distance_, normal_) &&
	      distance_ < distance)
	  {
	    distance = distance_;
	    copy3(normal_, normal);
	    surf_idx = i;
	  }
	}
	return surf_idx;
}

} //namespace





