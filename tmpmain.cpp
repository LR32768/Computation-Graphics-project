#include <cstdlib>
#include <string>
#include <cstdio> 
#include <vector>
#include "maths.h"
#include "tmpconfig.h"
#include "objects.hpp"
#include "utils.h"

extern std::vector<Object*> tri_obj; //triangle meshes
unsigned short visx;
int visy, hptail = 0;
hpList* hlist;
SpaceHash Htable;

inline bool intersect(const Ray &r, double &t, int &id) { 
	double d, inf = t = INF; 
	int n = int(sizeof(scenes) / sizeof(SphereObj*));
	for(int i = n; i--; ) {
		d = scenes[i]->intersect(r);
		if(d > 0 && d < t)
		{ t = d; id = i; }
	}
	for(int i = 0; i < tri_obj.size(); i++) {
		d = tri_obj[i]->intersect(r);
		if(d > 0 && d < t)
		{ t = d; id = n + i; }
	}
	return t < INF; 
}

void radiance(const Ray &r, int depth, int type, Vec flx) {
	//Core code for algorithm, 
	//type = 0 means eye ray, otherwise photon ray

	double t;                               // distance to intersection 
	int id = 0, num = int(sizeof(scenes) / sizeof(SphereObj*));         // id of intersected object 
	if (!intersect(r, t, id) || depth >= 8) return; // if miss, return black
	//printf("ray from %d,%d hit at object %d\n", visx, visy, id);
	depth++;

	const Object *obj = id < num ? scenes[id] : tri_obj[id - num];      // the hit object 
	Vec x = r.o + r.d*t, n = obj->get_norm(x), nl = n.dot(r.d) < 0 ? n : n*-1, f = obj->get_color(x); 
	double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z;

	if (obj->text.refl == DIFF)   // Ideal DIFFUSE reflection 
		if(!type) { // eye ray, then store the Hitpoint position
			Hitpoint* tmph = new Hitpoint(x, n, f.mult(flx), visx, visy);
			hlist = ListAdd(tmph, hlist);
			Htable.insert(hlist->hp);
			return;
		}
		else { //photon ray, update the neighborhood hitpoint's brightness and pass on
			std::vector<Hitpoint*>& candi = Htable.lookup(x);	
			for(int tt = 0; tt < candi.size(); tt ++)
				if(candi[tt]->n.dot(n) > 1e-3 && vnorm(candi[tt]->pos - x) < candi[tt]->radius) {
					double g = (candi[tt]->num * GAMMA + GAMMA) / (candi[tt]->num * GAMMA + 1.0);
					candi[tt]->num++;
					candi[tt]->radius *= sqrt(g);
					//update radius, lr=.8
					candi[tt]->flux = (candi[tt]->flux + f.mult(flx)) * g;
				}
			double r1 = 2*M_PI*rand01(), r2 = rand01(), r2s = sqrt(r2); 
			Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0,1) : Vec(1))%w).norm(), v = w % u; 
			Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm(); 
			if(rand01() < p)
				radiance(Ray(x, d), depth, type, flx.mult(f) * (1.0/p) + obj->text.emission); 
		} 
	else 
		if (obj->text.refl == SPEC) // Ideal SPECULAR reflection 
			radiance(Ray(x, r.d - n*2*n.dot(r.d)), depth, type, f.mult(flx) + obj->text.emission); 
		else { 
			Ray lr(x, r.d - n*2.0*n.dot(r.d)); bool into = (n.dot(nl) > 0.0);
			double nc = 1.0, nt = 1.5, nnt = into ? nc/nt:nt/nc, ddn = r.d.dot(nl), cos2t;
			if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0) 
				return radiance(lr, depth, type, flx);
			Vec td = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
			double a = nt-nc, b = nt+nc, R0 = a*a/(b*b), c = 1-(into?-ddn:td.dot(n));
			double Re=R0+(1-R0)*c*c*c*c*c,P=Re;Ray rr(x,td); Vec fa=f.mult(flx);
			if(!type) {
				radiance(lr, depth, type, fa*Re);
				radiance(rr, depth, type, fa*(1.0-Re));
			}
			else {
				(rand01()<P) ? 
					radiance(lr, depth, type, flx) : radiance(rr, depth, type, flx);
			}
		}
}

int main(int argc, char *argv[]) {

	int w = 2048, h = 1536, samps = argc==2 ? atoi(argv[1]) / 1000 : 1000000; // # samples 
	load_obj("Ocean.obj", Vec(50, 40, 85), Vec(60, 60, 60), 4);
	//TODO if you want to load triangle meshes, enable last row

	//Loading meshes
	//tri_obj = load_obj("Ocean.obj", Vec(50, 40, 85), Ve:c(10, 10, 10), 16);
	printf("%d,\n", tri_obj.size());
	printf("%s.\n", tri_obj[0]->text.filename.c_str());
	//for (int i = 0; i < 16; i++)
	//	printf("the i'th coordinate as %lf,%lf,%lf\n", tri_obj[i].a.x, tri_obj[i].a.y, tri_obj[i].a.z);

	Vec cx = Vec(w*.5135 / h), cy = (cx%cam.d).norm()*.5135, r, *c = new Vec[w*h]; 

	for(int iter = 0; iter < 10; iter++) {
		Htable.initial(); //Initialize Htable
		hpList* p = hlist;
		while (p != NULL) { hpList* tmp = p->next; delete p; p = tmp; }
		hlist = NULL;
		printf("\n%d round pass...\n", iter);
		//Construct the eye-ray hit point
		for (visy = 0; visy < h; ++visy) {
			fprintf(stderr, "\rHitPointPass %5.2f%%", 100.0 * visy / (h-1));
			for (visx = 0; visx < w; ++visx) {
				double r1=2*rand01(), dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1); 
				double r2=2*rand01(), dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2); 
				Vec d = cx*(((1 + dx) + visx)/w - .5) + 
					cy*(-((1 + dy) + visy)/h + .5) + cam.d; 
				radiance(Ray(cam.o + d * 140, d.norm()), 0, 0,  Vec(1, 1, 1));
			}
		}
		fprintf(stderr, "\n");
		#pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP 
		//Now wrting SPPM
		for (int i = 0; i < samps/10; i++)	{
			double ratio = 100.0 * (i + 1) / samps;
			fprintf(stderr, "\rPhotonPass %5.2f%%", ratio);
			int m = 1000 * i; Ray r; Vec f;
			for (int j = 0; j < 1000; j++) {
				newgenp(&r, &f, i / samps); 
				radiance(r, 0, 1, f);
			}
		}

		double shat = .05; //visual depth constant

		//printf("safe here...\n");
		hpList* lst = hlist;
		while (lst->next != NULL) {
			//printf("(%d,%d) flux (%lf,%lf,%lf) color (%lf,%lf,%lf)\n", hp.vx, hp.vy, hp.flux.x, hp.flux.y, 
			//		hp.flux.z, hp.col.x, hp.col.y, hp.col.z);
			Hitpoint hp_ = *(lst->hp);
			int index = hp_.vy * w + hp_.vx;
			c[index] = c[index] + hp_.flux.mult(hp_.col) * (1.0 / (hp_.radius*hp_.radius*M_PI * samps * 1000.0));
			lst = lst->next;
		}
	}
	FILE *f = fopen("tmp_image.ppm", "w");         // Write image to PPM file. 
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255); 
	for (int i = 0; i < w * h; i++) 
		fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z)); 
} 
