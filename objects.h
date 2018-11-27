/* Define different kind of objects here
 * Author: Rui Lu
 * Last Edit: 2018/11/19 12:11
 */
#ifndef OBJECTS_H_
#define OBJECTS_H_
#include "maths.h"
#include <cmath>
#include <cstdlib>
#include <algorithm>

// Define different kind of reflection
// DIFFuse, SPECular, REFlection, LIGHt, GROSsy
enum Refl_t { DIFF, SPEC, REFR, LIGH, GROS };

struct Sphere { 
	double rad;       // radius 
	Vec p, e, c;      // position, emission, color 
	Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive) 
	Sphere (double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_) : 
		rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {} 
	double intersect (const Ray &r) const 
	{ 
		// returns distance, 0 if nohit 
		Vec op = p - r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0 
		double t, b = op.dot(r.d), det = b*b - op.dot(op) + rad*rad; 
		if (det<0) 
			return 0; 
		else det = sqrt(det); 
		return (t = b-det) > EPS ? t : ((t = b+det) > EPS ? t : 0); 
	} 
}; 

#endif
