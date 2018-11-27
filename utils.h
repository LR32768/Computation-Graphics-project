/* Define handful utils data structures and methods here
 * Author: Rui Lu
 * Last edit: 2018/11/26, 10:28
 */
#ifndef UTILS_H_
#define UTILS_H_
#include "maths.h"
#include <algorithm>
#include <cmath>
#include <cstdio>

//Hit point structure
///record the hitting point 
struct Hitpoint{


};


//Box: bounding box
struct Box{
	Vec bounds[2];
	//bounds[0] is the min axis for box and bounds[1] is the max

	void update(const Vec &v)
	{
		if (v.x < bounds[0].x) bounds[0].x = v.x;
		if (v.y < bounds[0].y) bounds[0].y = v.y;
		if (v.z < bounds[0].z) bounds[0].z = v.z;
		if (v.x > bounds[1].x) bounds[1].x = v.x;
		if (v.y > bounds[1].y) bounds[1].y = v.y;
		if (v.z > bounds[1].z) bounds[1].z = v.z;
	}

	//fast bounding box intersection, return whether intersect or not
	double inter(const Ray &r) const
	{
		double 
		tmin  = (bounds[r.sign[0]].x - r.o.x)   * r.inv_d.x,
		tmax  = (bounds[1-r.sign[0]].x - r.o.x) * r.inv_d.x,
		tymin = (bounds[r.sign[1]].y - r.o.y)   * r.inv_d.y,
		tymax = (bounds[1-r.sign[1]].y - r.o.y) * r.inv_d.y; 

		if ((tmin > tymax) || (tymin > tmax)) return false; 
		if (tymin > tmin) tmin = tymin; 
		if (tymax < tmax) tmax = tymax; 

		double
		tzmin = (bounds[r.sign[2]].z - r.o.z) * r.inv_d.z,
		tzmax = (bounds[1-r.sign[2]].z - r.o.z) * r.inv_d.z; 

		if ((tmin > tzmax) || (tzmin > tmax)) return false; 
		if (tzmin > tmin) tmin = tzmin; 
		if (tzmax < tmax) tmax = tzmax; 

		return true; 
	}
};

#endif
