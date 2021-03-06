#ifndef BOX_H
#define BOX_H

#include"VectorMath.h"
#include"Particle.h"

class Box
{
public:
	Box()=default;
	void initialize();
	void setBasicInfo();
	void inflate_volume(double target_volFrac);
	void InitShear(double delg);
	void XzShearStrain();
	void read_positions(const char* file, int ignore_line);
	void record_positions(const char* file);

	~Box()=default;
public:
	int nAtoms;
	//box parameter
	double3 boxLo;
	double3 boxHi;

	Hvoigt6 Length;  
	double boxVolume;
	double epsilon;
	Hvoigt6 invLength;  
};
extern Box *box;
Box *box = NULL;

void defVar()
{
	_clock.timeUnit = 1.0;
	_clock.timeScale = 0.005*_clock.timeUnit;

	var.epsilon = 1.0;
}
#endif
