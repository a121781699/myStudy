#ifndef PARTICLE_H
#define PARTICLE_H

#include"VectorMath.h"

#define calcFi(f, prefactor, x, vec)     \
{                                          \
    vecScale();                            \
    vecSub();                              \
    vecScale();                            \
}
void calcForcei(double3& f, double const prefactor, double const x, double3 const vec)
{
	double3 vec_unit;
	vecUnit(vec_unit, vec);
	vecScale(f, prefactor*(1 - x), vec_unit);
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
	static void calcForce(void(*func)(double3&, double const, double const, double3 const), Particle *particle, int nAtoms);
	static void eulerInteger(Particle *particle, int nAtoms);
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

void Particle::calcForce(void(*func)(double3&,double const,double const,double3 const), Particle *particle, int nAtoms)
{
	double xNorm;
	double3 f, x;
	for (int i = 0; i < nAtoms; ++i)
	{
		double3 fsum = { 0,0,0 };
		for (int j = 0; j < nAtoms; ++j)
		{
			if (j == i) continue;
			vecSub(x, particle[j].x, particle[i].x);
			vecNorm(xNorm, x);
			double sigma = (particle[i].DiameterScale + particle[j].DiameterScale)*Particle::meanDiameterScale / 2.0;
			if (sigma <= xNorm) continue;
			double factor = -2 * var.epsilon / sigma;
			func(f, factor, 1 - xNorm / sigma, x);
			vecAdd(fsum, fsum, f);
			//vecSub(particle[i].f, particle[i].f, f);
		}
		particle[i].f = fsum;
	}
}

void Particle::eulerInteger(Particle *particle, int nAtoms)
{
	for (int i = 0; i < nAtoms; ++i)
		particle[i].update();
}
#endif
