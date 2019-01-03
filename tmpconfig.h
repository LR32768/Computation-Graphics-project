#ifndef TMPCONFIGS_H_
#define TMPCONFIGS_H_
#include "objects.hpp"
#include "maths.h"

Object* spheres[] = {//Scene: radius, position, emission, color, material 
	new SphereObj(1e5, Vec( 1e5+1,40.8,81.6), Vec(), Vec(.25,.75,.25), DIFF),//Left 
	new SphereObj(1e5, Vec(-1e5+99,40.8,81.6),Vec(), Vec(.25,.25,.75), DIFF),//Rght 
	new SphereObj(1e5, Vec(50,40.8, 1e5),     Vec(), Vec(.75,.75,.75), DIFF),//Back 
	new SphereObj(1e5, Vec(50,40.8,-1e5+170), Vec(), Vec(.75,.75,.75), DIFF),//Frnt 
	new SphereObj(1e5, Vec(50, 1e5, 81.6),    Vec(), Vec(.75,.1,.1),   DIFF),//Botm 
	new SphereObj(1e5, Vec(50,-1e5+81.6,81.6),Vec(), Vec(.75,.75,.75), DIFF),//Top 
	new SphereObj(14,  Vec(27,14,37),         Vec(), Vec(1,1,1)*.999,  SPEC),//Mirr Ball 
	new SphereObj(14,  Vec(73,14,68),         Vec(), Vec(1,1,1)*.999,  REFR),//Glas Ball
	new SphereObj(14,  Vec(73,43,20),         Vec(), Vec(.75,.25,.25), DIFF),//Diff Ball
	new TriangleObj(Vec(20,30,14), Vec(50,30,14), Vec(20,55,14),  Vec(100, 100, 0), Vec(100, 600, 0), Vec(500, 350, 0), Texture("rainbow.ppm")),
	new SphereObj(600, Vec(50,681.6-.27,81.6),Vec(12,12,12),    Vec(), DIFF) //Lite 
}; 


//define camera position and diretion
Ray cam(Vec(40,62,285.6), Vec(0.08,-0.12612,-1).norm()); // cam pos, dir
#endif

