#ifndef PARTICLE_H
#define PARTICLE_H

#include"VectorMath.h"

#define calcFi(f, prefactor, x, sigma)     \
{                                          \
    vecScale();                            \
    vecSub();                              \
    vecScale();                            \
}
void calcForcei(double3& f, double const prefactor, double3 x, double const sigma)
{
	double3 vecUnit = { 1.0,1.0,1.0 };
	vecScale(x, 1.0 / sigma, x);
	vecSub(f, vecUnit, x);
	vecScale(f, prefactor, f);
}


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
	static void calcForce(void(*func)(double3&, double const, double3, double const), Particle *particle, int nAtoms);
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

void Particle::calcForce(void(*func)(double3&,double const,double3,double const), Particle *particle, int nAtoms)
{
	for (int i = 0; i < nAtoms; ++i)
	{
		double3 fsum = { 0,0,0 };
		for (int j = 0; j < nAtoms; ++j)
		{
			double3 f, x;
			vecSub(x, particle[j].x, particle[i].x);
			double sigma = (particle[i].DiameterScale + particle[j].DiameterScale)*Particle::meanDiameterScale / 2.0;
			func(f, -2 * var.epsilon / sigma, x, sigma);
			vecAdd(fsum, fsum, f);
			//vecSub(particle[i].f, particle[i].f, f);
		}
		particle[i].f = fsum;
	}
}
#endif
