#ifndef PARTICLE_H
#define PARTICLE_H

#include"VectorMath.h"
class Particle
{
public:
	Particle()=default;
	Particle(Particle *particle);
	Particle& operator=(const Particle *particle);

	void update()
	{
		vecScaleAdd(x, x, _clock.timeScale, v);
		vecScaleAdd(v, v, _clock.timeScale, f);
	}
public:
	static double meanDiameterScale;
	//properties of particle
	int id;
	int type;
	int density;
	double3 x;   //position
	double3 v;   //velocity
	double3 f;   //force
	double3 xImage;

	//double m;            //mass
	double r;            //radius
	double e;            //energy
	double p;            //pressure
	double z;            //number of contact

	double DiameterScale;      

	~Particle()=default;
};
double Particle::meanDiameterScale=0.0;
extern Particle *particle;
Particle *particle=NULL;
#endif
