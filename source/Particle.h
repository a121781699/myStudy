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
	    vecScaleAdd(x,x,0.001,v);
	    vecScaleAdd(v,v,0.001,f);
	}
public:
	static double meanDiameterScale;
	//properties of particle
	double id;
	double type;
	double density;
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
