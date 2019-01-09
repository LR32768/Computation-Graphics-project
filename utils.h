/* Define handful utils data structures and methods here
 * Author: Rui Lu
 * Last edit: 2018/11/26, 10:28
 */
#ifndef UTILS_H_
#define UTILS_H_
#include "maths.h"
#include "objects.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

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

static std::vector<Object*> tri_obj;

// Loading a triangle mesh model from objpath and assign it in the space with  
// position bias of b, strething parameter k
void load_obj(char* objpath, Vec b = Vec(), Vec k = Vec(1,1,1), int threshold = 131072)
///	object path, bias, stretch, number of meshes that want to read-in
{
	printf("Start loading triangle mesh model...\n");
	FILE* fobj = fopen(objpath, "r");
	char buf[255], s[10];
	fgets(buf, 255, fobj);
	fgets(buf, 255, fobj);
	Vec vpoints[67000];
	int vtail = 0, ftail = 0;

	while (~fscanf(fobj, "%s", s)) {
		if (strcmp(s, "v") == 0) {
			double x, y, z;
			fscanf(fobj, "%lf%lf%lf", &x, &y, &z);
			//printf("reading point %lf,%lf,%lf\n", x, y, z);
			vpoints[vtail++] = Vec(x, y, z).mult(k) + b;
		}
		else if (strcmp(s, "f") == 0) {
			int na, nb, nc, nd;
			fscanf(fobj, "%d/%*d/%*d %d/%*d/%*d %d/%*d/%*d %d/%*d/%*d", &na, &nb, &nc, &nd);
			--na; --nb; --nc; --nd;
			//printf("reading quadra no %d,%d,%d,%d\n", na, nb, nc, nd);
			//change this line's reflection type into REFR after testing and debugging
			Object* tri1 = new TriangleObj(vpoints[na], vpoints[nb], vpoints[nc], Vec(), Vec(), DIFF); 
			Object* tri2 = new TriangleObj(vpoints[nc], vpoints[nd], vpoints[na], Vec(), Vec(), DIFF);
			ftail++;
			tri_obj.push_back(tri1);
			tri_obj.push_back(tri2);
		}
		else if (strcmp(s, "vt") == 0) {
			fscanf(fobj, "%*f%*f");
		}
		else if (strcmp(s, "vn") == 0) {
			fscanf(fobj, "%*f%*f%*f");
		}
		if (ftail >= threshold)
			break;
	}
	printf("Loading finished! Totally %d faces.\n", ftail);
}

#endif
