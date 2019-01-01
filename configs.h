#ifndef CONFIGS_H_
#define CONFIGS_H_
#include "maths.h"
enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance() 
struct Sphere { 
	double rad;       // radius 
	Vec p, e, c;      // position, emission, color 
	Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive) 
	Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_): 
		rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {} 
	double intersect(const Ray &r) const { // returns distance, 0 if nohit 
		Vec op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0 
		double t, eps=1e-4, b=op.dot(r.d), det=b*b-op.dot(op)+rad*rad; 
		if (det<0) return 0; else det=sqrt(det); 
		return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0); 
	} 
}; 

Sphere spheres[] = {//Scene: radius, position, emission, color, material 
	Sphere(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.25,.75,.25),DIFF),//Left 
	Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF),//Rght 
	Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFF),//Back 
	Sphere(1e5, Vec(50,40.8,-1e5+170), Vec(),Vec(.75,.75,.75),DIFF),//Frnt 
	Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.1,.1),DIFF),//Botm 
	Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top 
	Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(1,1,1)*.999, SPEC),//Mirr 
	Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.999, REFR),//Glas 
	Sphere(16.5,Vec(73,43,30),         Vec(),Vec(.75,.25,.25), DIFF),//Diff
	Sphere(600, Vec(50,681.6-.27,81.6),Vec(12,12,12),  Vec(), DIFF) //Lite 
}; 


//define camera position and diretion
Ray cam(Vec(40,62,285.6), Vec(0.08,-0.12612,-1).norm()); // cam pos, dir
#endif
