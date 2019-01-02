/* Math tools and basic functions here  
 * Author: Rui Lu
 * Last Edit: 2018/11/26 13:11
 */
#ifndef MATHS_H_
#define MATHS_H_
#include <cmath>
#define INF 1e20
#define EPS 1e-4

struct Vec {
	double x, y, z;
	Vec(double x_ = 0, double y_ = 0, double z_ = 0) {x = x_; y = y_; z = z_; }
	Vec operator+ (const Vec& b) const { return Vec(x+b.x, y+b.y, z+b.z); }
	Vec operator- (const Vec& b) const { return Vec(x-b.x, y-b.y, z-b.z); }
	Vec operator* (const double b) const { return Vec(x*b, y*b, z*b); }
	Vec	mult(const Vec& b) const { return Vec(x*b.x, y*b.y, z*b.z); }
	Vec& norm() { return *this = *this * (1 / sqrt(x*x + y*y + z*z)); }
	//dot product
	double dot(const Vec& b) const {return x*b.x + y*b.y + z*b.z; } 
	//cross product
	Vec operator% (Vec& b) { return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x); }
};

struct Ray {
	Vec o, d; //origin and direction
	Vec inv_d; //inverse direction for fast intersection
	int sign[3]; //sign for inverser direction
	Ray(Vec o_, Vec d_) : o(o_), d(d_) 
	{
		inv_d = Vec(1.0 / d_.x, 1.0 / d_.y, 1.0 / d_.z);
		sign[0] = (inv_d.x < 0);
		sign[1] = (inv_d.y < 0);
		sign[2] = (inv_d.z < 0);
	} 
};


static unsigned long rx = 123456789, ry = 362436069, rz = 521288629;
unsigned long frand(void) 
	///A fast random generator based on LFSR with period 2^96-1
{	
	unsigned long t;
	rx ^= rx << 16;
	rx ^= rx >> 5;
	rx ^= rx << 1;
	t = rx;
	rx = ry;
	ry = rz;
	rz = t ^ rx ^ ry;
	return rz;
}

double rand01(void)
{
	unsigned long MAX = (1 << 20);
	return double (frand() % MAX) / double(MAX);
}

//clamp function, clamp every double into [0, 1]
inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; } 

//round every colour channel lightness into 0-255 value
inline int toInt(double x){ return int(pow(clamp(x),1/2.2)*255+.5); }
 
#endif
