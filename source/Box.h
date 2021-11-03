#ifndef BOX_H
#define BOX_H

#include"VectorMath.h"
#include"Particle.h"

#include"VectorMath.h"
#include"Particle.h"

class Box
{
public:
	Box() = default;
	void initialize();
	void setBasicInfo();
	void inflate_volume(double target_volFrac);
	void InitShear(double delg);
	void XzShearStrain();
	void read_positions(const char* file, int ignore_line);
	void record_positions(const char* file);

	double getE();
	~Box() = default;
public:
	int nAtoms;
	//box parameter
	Hvoigt6 Length;
	double boxVolume;
	double epsilon;

	double3 boxLo;
	double3 boxHi;
	Hvoigt6 invLength;

	double E;
};
extern Box *box;
Box *box = NULL;

double Box::getE()
{
	if (IsCalcuE) return E;
	E = 0;
	for (int i = 0; i < nAtoms; ++i)
	{
		double e = 0;
		for (int j = 0; j < nAtoms; ++j)
		{
			if (j == i) continue;
			double sigma = (particle[i].DiameterScale + particle[j].DiameterScale)*Particle::meanDiameterScale / 2.0;
			double3 x;
			double x_norm, factor;
			vecSub(x, particle[i].x, particle[j].x);
			vecNorm(x_norm, x);
			if (x_norm >= sigma) continue;
			factor = var.epsilon / 2.0;
			e += factor * pow((1.0 - x_norm / sigma), 2);
		}
		E += e / 2.0;
	}
	IsCalcuE = true;
	//E = 1.0;
	return E;
}
#endif
