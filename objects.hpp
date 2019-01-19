/* Define different objects and its properties
 * Author: Rui Lu 
 * Last Edit: 2019/1/15 19:35
 */
#ifndef OBJECTS_H_
#define OBJECTS_H_

#include "maths.h"
#include <cstring>
#include <string>
#include <vector>
#include <cassert>

enum Refl_t { DIFF, SPEC, REFR, GLOSS }; 

//Hit point structure
///record the hitting point 
struct Hitpoint {
	Vec pos, n, col, flux;
	double radius = .1;
	int num, vx, vy;
	Hitpoint(Vec _pos = Vec(), Vec _n = Vec(), Vec _col=Vec(), int _vx=0, int _vy=0):
		pos(_pos), n(_n), col(_col), vx(_vx), vy(_vy), num(0), flux(Vec()) {}
};

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
			refl = DIFF;
			printf("starting loading textures...\n");
			FILE* f = fopen(tname, "r");
			char tmp[100];
			fgets(tmp, 100, f);
			fgets(tmp, 100, f);
			fscanf(f,"%d%d%*d\n", &w, &h);
			printf("we have %d*%d image as texture.\n", w, h);
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
			printf("Texture %d*%d loaded finish!\n", w, h);
		}


};

class Object {
	public:
		Texture text;
		Object(Texture t): text(t) {}
		Object(Refl_t refl, Vec color, Vec emission, double brdf, std::string tname):
			text(tname, brdf, color, emission, refl) {}
		virtual double intersect(const Ray &r) { printf("Not defined intersect yet!\n"); }
		virtual Vec get_norm(Vec &pos) const { printf("Not defined norm vector yet!\n"); }
		virtual Vec get_color(Vec &pos) const { printf("Not defined get_color yet!\n"); }
		//To be added for Bezier object...
};


class SphereObj: public Object {
	public:
		Vec o; 
		double rad;
		SphereObj(double r_, Vec o_, Texture t):
			o(o_), rad(r_), Object(t) {}
		SphereObj(double r_, Vec o_, Vec emission = Vec(), Vec color = Vec(),  Refl_t refl = DIFF, double brdf = 1.5, std::string tname = ""):
			o(o_), rad(r_), Object(refl, color, emission, brdf, tname) {}
		virtual double intersect(const Ray &r) {
			// returns distance, 0 if nohit 
			Vec op = o-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0 
			double t, b = op.dot(r.d), det=b*b - op.dot(op) + rad*rad; 
			if (det<0) return 0; else det=sqrt(det); 
			return (t = b-det) > EPS ? t : ((t = b+det) > EPS ? t : -1); 
		}
		virtual Vec get_norm(Vec& pos) const {
			double d = fabs(vnorm(pos - o) - rad);
			return d < EPS ? (pos - o).norm() : Vec(); // return 
		}
		virtual Vec get_color(Vec& pos) const {
			if (text.filename == "") 
				return text.color;
			else {
				Vec op = pos - o;
				double phi = acos(op.y / rad), theta = atan2(op.z, op.x);
				//printf("%lf and %lf\n", phi, theta);
				theta = theta < 0 ? theta + M_PI : theta;
				int dx = int(text.h * phi / M_PI), dy = int(text.w * theta / 2 / M_PI);
				//printf("%d,%d/%d\n", dx, dy, text.w);
				dx = dx < 0 ? 0 : dx >= text.h ? text.h-1 : dx;
				dy = dy < 0 ? 0 : dy >= text.w ? text.w-1 : dy;
				return text.pattern[dx][dy];
			}
		}
};

class TriangleObj: public Object {
	public:
		Vec a, b, c, n;  // The position for three end points, anti-clockwise
		Vec ta, tb, tc; // The corresponding point on the texture map
		TriangleObj (Vec a_, Vec b_, Vec c_, Vec ta_, Vec tb_, Vec tc_, Texture t): 
			a(a_), b(b_), c(c_), ta(ta_), tb(tb_), tc(tc_), Object(t)
	{ Vec _n = (c-b) % (a-b); n = _n.norm(); }
		TriangleObj (Vec a_, Vec b_, Vec c_,  Vec emission = Vec(), Vec color = Vec(),  Refl_t refl = DIFF, 
				double brdf = 1.5, std::string tname = ""): a(a_), b(b_), c(c_), Object(refl, color, emission, brdf, tname)
		{ Vec _n = (a-b) % (c-b); n = _n.norm(); }
		virtual double intersect(const Ray &r) {
			Vec E1 = b - a, E2 = c - a, P = r.d % E2, T;
			double det = E1.dot(P);
			if( det > 0 ) { T = r.o - a; } else { T = a - r.o; det = -det; }
			if( det < EPS )	return -1;
			double u = T.dot(P);
			if( u < 0.0f || u > det )	return -1;
			Vec Q = T % E1;
			double v = r.d.dot(Q);
			if( v < 0.0f || u + v > det )	return -1;
			double ans = E2.dot(Q) / det;
			return ans > EPS ? ans : -1;
		}
		virtual Vec get_norm(Vec & pos) const { return n; }
		virtual Vec get_color(Vec& pos) const {
			if (text.filename == "") 
				return text.color;
			Vec pa = a-pos, pb = b-pos, pc = c-pos;
			double rx = vnorm(pa % (b-a)) / (vnorm(pa % (b-a)) + vnorm((c-a) % pa));
			Vec ter = c * rx + b * (1.0 - rx);
			double ry = vnorm(pa) / vnorm(ter - a);
			Vec _ter = tc * rx + tb * (1 - rx);
			Vec ans = _ter * ry + ta * (1 - ry);
			assert(ans.x >= 0 && ans.x < text.h && ans.y >= 0 && ans.y < text.w);
			return text.pattern[int(ans.x)][int(ans.y)];
		}
		Vec get_interpolation(Vec& pos) const {
			if (text.filename == "") 
				return text.color;
			Vec pa = a-pos, pb = b-pos, pc = c-pos;
			double rx = vnorm(pa % (b-a)) / (vnorm(pa % (b-a)) + vnorm((c-a) % pa));
			Vec ter = c * rx + b * (1.0 - rx);
			double ry = vnorm(pa) / vnorm(ter - a);
			Vec _ter = tc * rx + tb * (1 - rx);
			Vec ans = _ter * ry + ta * (1 - ry);
			assert(ans.x >= 0 && ans.x < text.h && ans.y >= 0 && ans.y < text.w);
			return ans;
		}
		/*Vec map(const Triangle &dst, const Vec &p) {
		  double t1, t2;
		  Vec px, PX, P;
		  t1 = fabs((e1 % (p - v)).len() / (e2 % (p - v)).len());
		  t1 = t1 / (1. + t1);
		  px = v + e1 * (1 - t1) + e2 * t1;
		  PX = dst.v + dst.e1 * (1 - t1) + dst.e2 * t1;
		  t2 = (p - v).len() / (px - v).len();
		  P = dst.v * (1 - t2) + PX * t2;
		  return P;
		  }*/
};

class BCurve2D
{
	public:
		double *dx, *dy, max, height, max2, r, num, tt1=0, tt2=1;
		int n;
		struct D {
			double t0, t1, width, y0, y1, width2;
		} data[20];

		BCurve2D(double* px, double* py, int n_, int num_, double r_): num(num_), n(n_), r(r_) {
			dx = new double[n];
			dy = new double[n];
			assert(fabs(py[0]) <= 1e-6);
			--n;
			// preproces
			for(int i = 0; i <= n; ++i) {
				dx[i] = px[0];
				dy[i] = py[0];
				for (int j = 0; j <= n - i; ++j) {
					px[j] = px[j + 1] - px[j];
					py[j] = py[j + 1] - py[j];
				}
			}
			double n_down = 1, fac = 1, nxt = n;
			for (int i = 0; i <= n; ++i, --nxt) {
				fac = fac * (i == 0 ? 1 : i);
				dx[i] = dx[i] * n_down / fac;
				dy[i] = dy[i] * n_down / fac;
				n_down *= nxt;
			}
			max = 0;
			double interval = 1. / (num - 1), c = 0;
			for (int cnt = 0; cnt <= num; c += interval, ++cnt) {
				data[cnt].width = 0;
				data[cnt].t0 = c - r > 0. ? c - r : 0.;
				data[cnt].t1 = c + r > 1. ? c + r : 0.;
				data[cnt].y0 = getpos(data[cnt].t0).y;
				data[cnt].y1 = getpos(data[cnt].t1).y;
				for (double t = data[cnt].t0; t <= data[cnt].t1; t += 1e-5) {
					Vec pos = getpos(t);
					if (data[cnt].width < pos.x)
						data[cnt].width = pos.x;
				}
				if (max < data[cnt].width)
					max = data[cnt].width;
				data[cnt].width += 1e-6;
				data[cnt].width2 = pow(data[cnt].width, 2);
			}
			max += 1e-6;
			max2 = max * max;
			height = getpos(1).y;
		}

		//compute postition given parameter t
		Vec getpos(double t) const {
			double ans_x = 0, ans_y = 0, t_pow = 1;
			for (int i = 0; i <= n; ++i) {
				ans_x += dx[i] * t_pow;
				ans_y += dy[i] * t_pow;
				t_pow *= t;
			}
			return Vec(ans_x, ans_y);
		}

		//compute tagent direction
		Vec getdir(double t) const {
			double ans_x = 0, ans_y = 0, t_pow = 1;
			for(int i = 1; i <= n; ++i) {
				ans_x += dx[i] * i * t_pow;
				ans_y += dy[i] * i * t_pow;
				t_pow *= t;
			}
			return Vec(ans_x, ans_y);
		}
		//compute second derivative
		Vec getdir2(double t) const {
			double ans_x = 0, ans_y = 0, t_pow = 1;
			for(int i = 2; i <= n; ++i) {
				ans_x += dx[i] * i * (i - 1) * t_pow;
				ans_y += dy[i] * i * (i - 1) * t_pow;
				t_pow *= t;
			}
			return Vec(ans_x, ans_y);
		}
};

class BezierObj: public Object {
	public:
		BCurve2D curve;
		Vec pos; //The center point at bottom

		BezierObj(Vec pos_, BCurve2D c_, Texture t): pos(pos_), curve(c_), Object(t) {}
		BezierObj(Vec pos_, BCurve2D c_, Vec emission = Vec(), Vec color = Vec(), Refl_t refl=DIFF, 
				double brdf = 1.5, std::string tname = ""):
			pos(pos_), curve(c_), Object(refl, color, emission, brdf, tname) {}

		double solve_t(double yc) const { 
			// solve y(t)=yc by tangent iteration
			double t = .5, ft, dft;
			for (int i = 10; i--; ) {
				if (t < 0) t = 0;
				else if (t > 1) t = 1;
				ft = curve.getpos(t).y - yc, dft = curve.getdir(t).y;
				if (fabs(ft) < 1e-6)
					return t;
				t -= ft / dft;
			}
			return -1;
		}

		double get_sphere_intersect(Ray ray, Vec o, double r) {
			Vec ro = o - ray.o;
			double b = ray.d.dot(ro);
			double d = b * b - ro.dot(ro) + r * r;
			if (d < 0) return -1;
			else d = sqrt(d);
			double t = b - d > 1e-6 ? b - d :(b + d > 1e-6? b + d : -1);
			return t < 0 ? -1 : t;
		}

		virtual double intersect(const Ray &ray) {
			double final_dis = (1<<20);
			if (fabs(ray.d.y) < 5e-4) {
				double dis_to_axis = vnorm(Vec(pos.x, ray.o.y, pos.z) - ray.o);
				double hit = (ray.o + ray.d * dis_to_axis).y;
				if (hit < pos.y + 1e-6 || hit > pos.y + curve.height - 1e-6)
					return -1;
				// solve function pos.y+y(t)=ray.o.y to get x(t)
				double t = solve_t(hit - pos.y);
				if (t < 0 || t > 1) return -1;

				Vec loc = curve.getpos(t);
				double ft = pos.y + loc.y - hit;
				if (fabs(ft) > 1e-6)	return -1;

				final_dis = get_sphere_intersect(ray, Vec(pos.x, pos.y + loc.y, pos.z), loc.x);
				if (final_dis < 0) return -1;

				Vec inter_p = ray.o + ray.d * final_dis;
				if (fabs((inter_p - Vec(pos.x, inter_p.y, pos.z)).Norm2() - loc.x*loc.x) > 1e-1) return -1;

				hit = inter_p.y;
				if (hit < pos.y + 1e-6 || hit > pos.y + curve.height - 1e-6) return -1;
				t = solve_t(hit - pos.y);
				loc = curve.getpos(t);
				ft = pos.y + loc.y - hit;
				if (fabs(ft) > 1e-6) return -1;
				final_dis = get_sphere_intersect(ray, Vec(pos.x, pos.y + loc.y, pos.z), loc.x);
				if (final_dis < 0) return -1;
				inter_p = ray.o + ray.d * final_dis;
				if (fabs((inter_p - Vec(pos.x, hit, pos.z)).Norm2() - (loc.x)*(loc.x)) > 1e-2) return -1;
				return final_dis;
			}
			double a = 0, b = 0, c = 0, t1, t2;
			// (xo-x'+xd/yd*(y-yo))^2 -> (t1+t2*y)^2
			t1 = ray.o.x - pos.x - ray.d.x / ray.d.y * ray.o.y; t2 = ray.d.x / ray.d.y;
			a += t2 * t2; b += 2 * t1 * t2; c += t1 * t1;
			// (zo-z'+zd/yd*(y-yo))^2 -> (t1+t2*y)^2
			t1 = ray.o.z - pos.z - ray.d.z / ray.d.y * ray.o.y; t2 = ray.d.z / ray.d.y;
			a += t2 * t2; b += 2 * t1 * t2; c += t1*t1;
			// ay^2+by+c -> a'(y-b')^2+c'
			c = c - b * b / 4 / a;
			b = -b / 2 / a - pos.y;
			// printf("%lf %lf %lf\n",a,b,c);
			if (0 <= b && b <= curve.height && c > curve.max2
					|| (b < 0 || b > curve.height) && std::min(b*b, (curve.height - b)*(curve.height-b)) * a + c > curve.max2) //no intersect
				return -1;
			for(int ind = 0; ind <= curve.num; ++ind) {
				double t0 = curve.data[ind].t0, t1 = curve.data[ind].t1;
				check(t0, t1, (t0 + t1 + t0) / 3, ray, a, b, c, final_dis);
				check(t0, t1, (t1 + t0 + t1) / 3, ray, a, b, c, final_dis);
			}
			check(curve.tt1-1e-2, curve.tt1+1e-2, curve.tt1, ray, a, b, c, final_dis);
			check(curve.tt2-1e-2, curve.tt2+1e-2, curve.tt2, ray, a, b, c, final_dis);
			return (final_dis < (1<<20) / 2) ? final_dis : -1;
		}

		bool check(double low, double upp, double init, Ray ray, double a, double b, double c, double& final_dis) {
			double t = newton(init, a, b, c, low, upp);
			if (t <= 0 || t >= 1)
				return false;
			Vec loc = curve.getpos(t);
			double x = loc.x, y = loc.y;
			double ft = x - sqrt(a * (y-b) * (y-b) + c);
			if (fabs(ft) > 1e-6)
				return false;
			// calc t for ray
			double dis = (pos.y + y - ray.o.y) / ray.d.y;
			if (dis < 1e-6)
				return false;
			Vec inter_p = ray.o + ray.d * dis;
			if (fabs((Vec(pos.x, pos.y + y, pos.z) - inter_p).Norm2() - x * x) > 1e-6)
				return false;
			if (dis < final_dis) {
				final_dis = dis;
				return true;
			}
			return false;
		}

		double newton(double t, double a, double b, double c, double low = 1e-6, double upp = 1-1e-6) {
			double ft, dft, x, y, dx, dy, sq;
			Vec loc, dir;

			for (int i = 10; i--; ) {
				if (t < 0) t = low;
				if (t > 1) t = upp;
				loc = curve.getpos(t), dir = curve.getdir(t);
				x = loc.x, dx = dir.x;
				y = loc.y, dy = dir.y;
				// printf("%lf %lf %lf\n",t,x,y);
				sq = sqrt(a * (y - b) * (y - b) + c);
				ft = x - sq;
				dft = dx - a * (y - b) * dy / sq;
				if (fabs(ft) < 1e-6)
					return t;
				t -= ft / dft;
			}
			return -1;
		}



		Vec Bezier_axis_trans(Vec& inter_p) const {
			double t = solve_t(inter_p.y - pos.y);
			double u = atan2(inter_p.z - pos.z, inter_p.x - pos.x); // between -pi ~ pi
			if (u < 0)
				u += 2 * M_PI;
			return Vec(u, t);
		}

		virtual Vec get_norm(Vec& pos) const {
			Vec tmp = Bezier_axis_trans(pos);
			Vec dir = curve.getdir(tmp.y);
			Vec d_surface = Vec(cos(tmp.x)*dir.x, dir.y, sin(tmp.x)*dir.x);
			Vec d_circ = Vec(-sin(tmp.x), 0, cos(tmp.x));
			return (d_circ % d_surface).norm();
		}

		virtual Vec get_color(Vec& pos) const {
			if (text.filename == "") 
				return text.color;
			Vec ut = Bezier_axis_trans(pos);
			int dy = int((ut.x / (2 * M_PI)) * text.w), dx = int(ut.y * text.h);
			//printf("%lf,%lf\n",ut.x,ut.y);
			dx = dx < 0 ? 0 : dx >= text.h ? text.h : dx; 
			dy = dy < 0 ? 0 : dy >= text.w ? text.w : dy;
			return text.pattern[dx][dy];
		}
};


#endif

