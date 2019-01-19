#ifndef TMPCONFIGS_H_
#define TMPCONFIGS_H_
#include "objects.hpp"
#include "maths.h"

const double bz_x = 0.8, bz_y = 1.2;
double control_x[] = {4.2/bz_x,15.5/bz_x,1.3/bz_x,10./bz_x};
double control_y[] = {0./bz_y,10/bz_y,30./bz_y,40./bz_y};
double slide = 3950-2500;
BCurve2D bezier(control_x, control_y, 4, 20, .365);

Object* scenes[] = {//Scene: radius, position, emission, color, material  
	new TriangleObj(Vec(0,-1e3,-1e3)  ,   Vec(0, 2e3, -1e3),	  Vec(0, 2e3, 3e3),     Vec(), Vec(.75,.75,.75), DIFF), //new Left1
	new TriangleObj(Vec(0, 2e3, 3e3) ,    Vec(0, -1e3, 3e3),	  Vec(0,-1e3,-1e3),     Vec(), Vec(.75,.75,.75), DIFF), //new Left2
	new TriangleObj(Vec(100), Vec(100, 82, 150), Vec(100, 82),   Vec(380), Vec(0,600), Vec(), Texture("windows.ppm")), //new Right1
	new TriangleObj(Vec(100, 82, 150), Vec(100), Vec(100,0,150), Vec(0,600), Vec(380), Vec(380,600), Texture("windows.ppm")), //new Right2
	new TriangleObj(Vec(), Vec(100), Vec(100, 82), Vec(82*25,slide), Vec(82*25,100*25+slide), Vec(0,100*25+slide), Texture("background.ppm")), //new Back1
	new TriangleObj(Vec(100, 82) ,  Vec(0, 100), Vec(),    Vec(0,100*25+slide), Vec(0,slide), Vec(82*25,slide), Texture("background.ppm")), //new Back2
	new TriangleObj(Vec(-1e3, -1e3, 150), Vec(1.5e3, 2e3, 150),   Vec(1.5e3, -1e3, 150),Vec(), Vec(.75,.75,.75), DIFF), //new Frnt1
	new TriangleObj(Vec(1.5e3, 2e3, 150), Vec(-1e3, -1e3, 150),   Vec(-1e3, 2e3, 150),  Vec(), Vec(.75,.75,.75), DIFF), //new Frnt2
	new TriangleObj(Vec(), Vec(0, 0, 170),  Vec(100, 0, 170), Vec(), Vec(850,0), Vec(850,500), Texture("woodfloor.ppm")), //new Bottom1
	new TriangleObj(Vec(), Vec(100, 0, 170),Vec(100, 0, 0),   Vec(),  Vec(850,500), Vec(0,500), Texture("woodfloor.ppm")), //new Bottom2
	new TriangleObj(Vec(-1e3,81.6,-1e3),  Vec(1.5e3, 81.6, 3e3),  Vec(-1e3, 81.6, 3e3), Vec(), Vec(.75,.75,.75), DIFF), //new Top1
	new TriangleObj(Vec(-1e3,81.6,-1e3),  Vec(1.5e3, 81.6, -1e3), Vec(1.5e3, 81.6, 3e3),Vec(), Vec(.75,.75,.75), DIFF), //new Top2
	new SphereObj(8,  Vec(27,8,37),     Vec(), Vec(1,1,1)*.999,  SPEC),//Mirr Ball 
	new SphereObj(7,  Vec(73,7,73),   Vec(), Vec(1,.1,1)*.999,  REFR),//Glas Ball
	new SphereObj(6,  Vec(44,6,59),     Vec(), Vec(.75,.25,.25), REFR),//Diff Ball
	new SphereObj(6,  Vec(58,64,80),    Vec(), Vec(.25,.75,.25), REFR),//Moon
	new BezierObj(Vec(77, 0, 25), bezier, Texture("vase.ppm")),
	//new TriangleObj(Vec(20,30,14), Vec(50,30,14), Vec(20,55,14),  Vec(100, 400, 0), Vec(550, 400, 0), Vec(100, 550, 0), Texture("rainbow.ppm")),
	//new SphereObj(600, Vec(50,681.6-.27,81.6),Vec(12,12,12),  Vec(), DIFF) //Lite 
}; 


//define camera position and diretion
Ray cam(Vec(20,44,235.6), Vec(0.269,-0.07,-1).norm()); // cam pos, dir
#endif


