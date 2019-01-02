/* Define different objects and its properties
 * Author: Rui Lu 
 * Last Edit: 2019/1/1 18:33
*/
#ifndef OBJECTS_H_
#define OBJECTS_H_

#include "maths.h"
#include <cstring>

enum Refl_t { DIFF, SPEC, REFR }; 


class Texture {
public:
	Vec color, emission;
	Refl_t refl;
	double brdf;
	std::string filename;
	int w, h, c;
	//unsigned char *buf;
	Texture(std::string tname, double b, Vec col, Vec e, Refl_t r): filename(tname), brdf(b), color(col), emission(e), refl(r) {
	}

	/*
	Texture(const Texture&t): brdf(t.brdf), filename(t.filename), color(t.color), emission(t.emission), refl(t.refl) {
		if (t.filename != "")
			buf = stbi_load(filename.c_str(), &w, &h, &c, 0);
		else
			buf = NULL;
	}
	std::pair<Refl_t, P3> getcol(ld a, ld b) {
		if (buf == NULL)
			return std::make_pair(refl, color);
		int pw = (int(a * w) % w + w) % w, ph = (int(b * h) % h + h) % h;
		int idx = ph * w * c + pw * c;
		int x = buf[idx + 0], y = buf[idx + 1], z = buf[idx + 2];
		// printf("find point %d %d %lf %lf\n", ph, pw,a,b);
		if (x == 233 && y == 233 && z == 233) {
			return std::make_pair(SPEC, P3(1, 1, 1)*.999);
		}
		return std::make_pair(refl, P3(x, y, z) / 255.);
	}*/
};

class obj {
public:
	Texture text;
	Object(Texture t): texture(t) {}
	Object(Refl_t refl, Vec color, Vec emission, double brdf, std::string tname):
		texture(tname, brdf, color, emission, refl) {}
	virtual double intersect(const Ray &r) { printf("Not defined intersect yet!\n"); }
	virtual Vec get_norm(const Vec &pos) { printf("Not defined norm vector yet!\n"); }
	//To be added for Bezier object...
};


class SphereObj: public Object {
public:
	Vec o; 
	double rad;
	SphereObj(Vec o_, double r_, Texture t):
		o(o_), r(r_), Object(t) {}
	SphereObj(Vec o_, double r_, Refl_t refl, double brdf = 1.5, Vec color = Vec(), Vec emission = Vec(), std::string tname = ""):
		o(o_), r(r_), Object(refl, color, emission, brdf, tname) {}
	virtual double intersect(const Ray &r) const 
	{ // returns distance, 0 if nohit 
		Vec op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0 
		double t, b=op.dot(r.d), det=b*b-op.dot(op)+rad*rad; 
		if (det<0) return 0; else det=sqrt(det); 
		return (t = b-det) > EPS ? t : ((t = b+det) > EPS ? t : 0); 
	}
	virtual Vec get_norm(Vec & pos) {
		double d = (pos - o).dot(pos - o) - rad;
		return d < EPS ? (p - o).norm() : Vec(); // return 
	}
};

class PlaneObj: public Object{
	Vec p, n;  //special point and normal vector
	PlaneObj (Vec p_, Vec n_, Texture t):
		p(p_), n(n_), Object(t) {}
	PlaneObj
};
#endif
