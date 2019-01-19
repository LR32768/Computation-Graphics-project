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

//A link list for the hitpoint
struct hpList {
	Hitpoint* hp;
	hpList* next=NULL;
	hpList(Hitpoint* _hp = NULL, hpList* _next = NULL) : hp(_hp), next(_next) {}
};

//add one hitpoint to the list
hpList* ListAdd(Hitpoint* _hp, hpList* oldlist) {
	hpList* p = new hpList(_hp, oldlist);
	return p;
}

// A hash table data structure for finding hitting points recorded somewhere in the space
class SpaceHash {
public:
	std::vector<std::vector<Hitpoint*> > HitList;
	double stride = .2;
	int tablesize;
	SpaceHash () {}
	void initial (int size = 10000000) {
		HitList.clear();
		HitList.resize(size);
		tablesize = size;
	}
	unsigned int hash_int(const int a, const int b, const int c) { 
	//simple hash function for int
		unsigned int pr1 = 73850693, pr2 = 83492791, pr3 = 19349663;
		unsigned key = ((pr1 * a) ^ (pr2 * b) ^ (pr3 * c)) % tablesize;
		return key;
	}
	unsigned int hash_pos(const Vec pos) { 
	//simple hash function	
		int dx = floor(pos.x / stride), dy = floor(pos.y / stride), dz = floor(pos.z / stride);
		return hash_int(dx, dy, dz);
	}
	std::vector<Hitpoint*>& lookup(const Vec pos) {
		int key = hash_pos(pos);
		return HitList[key];
	}
	void insert(Hitpoint* newp)
	{
		//printf("get hit point at %lf,%lf,%lf\n", newp->pos.x, newp->pos.y, newp->pos.z);
		//printf("x range %d to %d\n", int((newp->pos.x - newp->radius) / stride), int((newp->pos.x + newp->radius) / stride));
		//printf("y range %d to %d\n", int((newp->pos.y - newp->radius) / stride), int((newp->pos.y + newp->radius) / stride));
		//printf("z range %d to %d\n", int((newp->pos.z - newp->radius) / stride), int((newp->pos.z + newp->radius) / stride));
		for(int _dx = int((newp->pos.x - newp->radius) / stride);  _dx <= int((newp->pos.x + newp->radius) / stride); _dx ++)
			for(int _dy = int((newp->pos.y - newp->radius) / stride);  _dy <= int((newp->pos.y + newp->radius) / stride); _dy ++)
				for(int _dz = int((newp->pos.z - newp->radius) / stride);  _dz <= int((newp->pos.z + newp->radius) / stride); _dz ++)
					{
						unsigned int key = hash_int(_dx, _dy, _dz);
						//printf("The key is %d\n", key);
						HitList[key].push_back(newp);
					}
		//printf("insert finished!\n");
	}
};

//Box: bounding box
struct Box{
	Vec bounds[2];
	//bounds[0] is the min axis for box and bounds[1] is the max
	Box() : bounds({Vec(INF, INF, INF), Vec(-INF, -INF, -INF)}) {}

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
			Object* tri1 = new TriangleObj(vpoints[na], vpoints[nb], vpoints[nc], Vec(), Vec(), REFR); 
			Object* tri2 = new TriangleObj(vpoints[nc], vpoints[nd], vpoints[na], Vec(), Vec(), REFR);
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


//Generating random photon's starting position and direction
void genp(Ray* pr, Vec* f, double ratio) {
	*f = Vec(2500, 2500, 2500) * M_PI * 7.0;
	double p = 2. * M_PI * rand01(), t = 2. * acos(sqrt(1.-rand01()));
	double st = sin(t);
	pr->d = Vec(cos(p)*st,cos(t),sin(p)*st);
	Vec dev = Vec(20 * rand01() - 10, 0, 20 * rand01() - 10);
	pr->o = Vec(50, 60, 85) + dev; //to be adjusted into random starting posistion
}

void newgenp(Ray *pr, Vec* f, double ratio) {
	*f = Vec(2500, 2500, 2500) * M_PI * 5.0;
	double theta = M_PI * .2 * rand01(), phi = M_PI * 2. * rand01();
	pr->d = Vec(sin(phi) * sin(theta), cos(theta), cos(phi) * sin(theta)).norm();
	double theta2 = M_PI * 2. * rand01(), rho = rand01() * 5.;
	pr->o = Vec(70, 80, 40) + Vec(sin(theta2) * rho, 0, cos(theta2) * rho);
}

#endif
