#include "maths.h"
#include "utils.h"
#include "objects.h"
#include "configs.h"
#include <cstdio>
#include <iostream>
#include <cmath>
#include <cstring>

inline bool intersect(const Ray &r, double &t, int &id){ 
	double n = sizeof(spheres) / sizeof(Sphere), d; 
	for(int i = int(n); i--; ) 
		if((d = spheres[i].intersect(r)) && d<t) { t=d; id=i; } 
	return t < INF; 
}

Vec radiance(const Ray &r, int depth){ 
	double t;                               // distance to intersection 
	int id=0;                               // id of intersected object 
	if (!intersect(r, t, id)) return Vec(); // if miss, return black 
	const Sphere &obj = spheres[id];        // the hit object 
	Vec x=r.o+r.d*t, n=(x-obj.p).norm(), nl=n.dot(r.d)<0?n:n*-1, f=obj.c; 
	double p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl 
	if (++depth>5) if (rand01()<p) f=f*(1/p); else return obj.e; //R.R. 
	if (obj.refl == DIFF){                  // Ideal DIFFUSE reflection 
		double r1=2*M_PI*rand01(), r2=rand01(), r2s=sqrt(r2); 
		Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1):Vec(1))%w).norm(), v=w%u; 
		Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm(); 
		return obj.e + f.mult(radiance(Ray(x,d),depth)); 
	} else 
	if (obj.refl == SPEC)            // Ideal SPECULAR reflection 
		return obj.e + f.mult(radiance(Ray(x,r.d-n*2*n.dot(r.d)),depth)); 
	Ray reflRay(x, r.d-n*2*n.dot(r.d));     // Ideal dielectric REFRACTION 
	bool into = n.dot(nl)>0;                // Ray from outside going in? 
	double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=r.d.dot(nl), cos2t; 
	if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    // Total internal reflection 
		return obj.e + f.mult(radiance(reflRay,depth)); 
	Vec tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm(); 
	double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n)); 
	double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P); 
	Vec ans = obj.e + f.mult(depth>2 ? (rand01()<P ?   // Russian roulette 
				radiance(reflRay,depth)*RP:radiance(Ray(x,tdir),depth)*TP) : 
			radiance(reflRay,depth)*Re+radiance(Ray(x,tdir),depth)*Tr); 
	//debugging
	//printf("%f,%f,%f\n", ans.x, ans.y, ans.z);
	return ans;
}

void print_to_ppm(Vec* c, int w, int h)
{
	FILE *f = fopen("output_image.ppm","w");
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255); 
	for (int i=0; i<w*h; i++) 
	{
		fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
		//printf("%d,%d,%d\n",toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
	}
}

int main(int argc, char *argv[])
{
	int w=1024, h=768, samps = argc==2 ? atoi(argv[1])/4 : 1; // # samples 
	Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).norm()); // cam pos, dir 
	Vec cx=Vec(w*.5135/h), cy=(cx%cam.d).norm()*.5135, r, *c=new Vec[w*h]; 
//#pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP 
	for (int y=0; y<h; y++){                       // Loop over image rows 
		fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1)); 
		for (unsigned short x=0; x<w; x++)   // Loop cols 
			for (int sy = 0, i=(h-y-1)*w+x; sy<2; sy++)     // 2x2 subpixel rows 
				for (int sx = 0; sx < 2; sx++, r=Vec()){        // 2x2 subpixel cols 
					for (int s=0; s < samps; s++){ 
						double r1=2*rand01(), dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1); 
						double r2=2*rand01(), dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2); 
						Vec d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) + 
							cy*( ( (sy+.5 + dy)/2 + y)/h - .5) + cam.d; 
						//printf("%f\n", radiance(Ray(cam.o+d*140,d.norm()),0).x*(1./samps));
						r = r + radiance(Ray(cam.o+d*140,d.norm()),0)*(1./samps); 
					} // Camera rays are pushed ^^^^^ forward to start in interior 
					c[i] = c[i] + Vec(clamp(r.x),clamp(r.y),clamp(r.z))*.25; 
					//std::cout << r.x << " " << r.y << " " << r.z << std::endl;
				} 
	}
	print_to_ppm(c, w, h);
	return 0;
}
