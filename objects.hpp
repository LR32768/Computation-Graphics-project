/* Define different objects and its properties
 * Author: Rui Lu 
 * Last Edit: 2019/1/2 19:35
 */
#ifndef OBJECTS_H_
#define OBJECTS_H_

#include "maths.h"
#include <cstring>
#include <string>
#include <vector>
#include <cassert>

enum Refl_t { DIFF, SPEC, REFR }; 


class Texture {
	public:
		Vec color, emission;
		Refl_t refl;
		double brdf;
		std::string filename;
		int w, h;
		std::vector<std::vector<Vec> > pattern;
		//unsigned char *buf;
		Texture(std::string tname, double b, Vec col, Vec e, Refl_t r): filename(tname), brdf(b), color(col), emission(e), refl(r) 
		{ }
		Texture(const char* tname, double b = 1.5): filename(tname), brdf(b) {
			printf("starting loading textures...\n");
			FILE* f = fopen(tname, "r");
			char tmp[100];
			fgets(tmp, 100, f);
			fgets(tmp, 100, f);
			fscanf(f,"%d%d%*d\n", &w, &h);
			printf("we have %d*%d", w, h);
			pattern.resize(h);
			for (int i = 0; i < h; i ++) {
				pattern[i].resize(w);
				unsigned char ta, tb, tc;
				for (int j = 0; j < w; j++)	{
					fread(&ta, sizeof(char), 1, f);
					fread(&tb, sizeof(char), 1, f);
					fread(&tc, sizeof(char), 1, f);
					pattern[i][j] = Vec(double(ta) / 255.0, double(tb) / 255.0, double(tc) / 255.0);
				}
			}
			printf("Texture%d*%d loaded finish!\n", w, h);
		}


};

class Object {
	public:
		Texture text;
		Object(Texture t): text(t) {}
		Object(Refl_t refl, Vec color, Vec emission, double brdf, std::string tname):
			text(tname, brdf, color, emission, refl) {}
		virtual double intersect(const Ray &r) { printf("Not defined intersect yet!\n"); }
		virtual Vec get_norm(const Vec &pos) const { printf("Not defined norm vector yet!\n"); }
		virtual Vec get_color(const Vec &pos) const { printf("Not defined get_color yet!\n"); }
		//To be added for Bezier object...
};


class SphereObj: public Object {
	public:
		Vec o; 
		double rad;
		SphereObj(Vec o_, double r_, Texture t):
			o(o_), rad(r_), Object(t) {}
		SphereObj(double r_, Vec o_, Vec emission = Vec(), Vec color = Vec(),  Refl_t refl = DIFF, double brdf = 1.5, std::string tname = ""):
			o(o_), rad(r_), Object(refl, color, emission, brdf, tname) {}
		virtual double intersect(const Ray &r) 
		{ // returns distance, 0 if nohit 
			Vec op = o-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0 
			double t, b = op.dot(r.d), det=b*b - op.dot(op) + rad*rad; 
			if (det<0) return 0; else det=sqrt(det); 
			return (t = b-det) > EPS ? t : ((t = b+det) > EPS ? t : 0); 
		}
		virtual Vec get_norm(const Vec& pos) const
		{
			double d = fabs(vnorm(pos - o) - rad);
			return d < EPS ? (pos - o).norm() : Vec(); // return 
		}
		virtual Vec get_color(const Vec& pos) const
		{
			if (text.filename == "") 
				return text.color;
			else 
				return Vec();
		}
};

class TriangleObj: public Object {
	public:
		Vec a, b, c, n;  // The position for three end points, anti-clockwise
		Vec ta, tb, tc; // The corresponding point on the texture map
		TriangleObj (Vec a_, Vec b_, Vec c_, Vec ta_, Vec tb_, Vec tc_, Texture t): 
			a(a_), b(b_), c(c_), ta(ta_), tb(tb_), tc(tc_), Object(t)
		{ Vec _n = (a-b) % (c-b); n = _n.norm(); }
		TriangleObj (Vec a_, Vec b_, Vec c_,  Vec emission = Vec(), Vec color = Vec(),  Refl_t refl = DIFF, 
				double brdf = 1.5, std::string tname = ""): a(a_), b(b_), c(c_), Object(refl, color, emission, brdf, tname)
		{ Vec _n = (a-b) % (c-b); n = _n.norm(); }
		virtual double intersect(const Ray &r)
		{
			Vec E1 = b - a;
			Vec E2 = c - a;
			Vec P = r.d % E2;
			double det = E1.dot(P);
			Vec T;
			if( det >0 ) { T = r.o - a; } else { T = a - r.o; det = -det; }
			if( det < EPS )	return false;
			double u = T.dot(P);
			if( u < 0.0f || u > det )	return false;
			Vec Q = T % E1;
			double v = r.d.dot(Q);
			if( v < 0.0f || u + v > det )	return false;
			return E2.dot(Q) / det;
		}
		virtual Vec get_norm(const Vec & pos) const { return n; }
		virtual Vec get_color(const Vec& pos) const
		{
			if (text.filename == "") 
				return text.color;
			Vec pa = pos - a, pb = pos - b, pc = pos - c;
			double rx = vnorm(pa % pb) / (vnorm(pa % pb) + vnorm(pc % pa));
			Vec ter = b * rx + c * (1 - rx);
			double ry = vnorm(pa) / vnorm(ter - a);
			Vec _ter = tb * rx + tc * (1 - rx);
			Vec ans = _ter * ry + ta * (1 - ry);
			assert(ans.x >= 0 && ans.x < text.h && ans.y >= 0 && ans.y < text.w);
			return text.pattern[int(ans.x)][int(ans.y)];
		}
};
#endif
