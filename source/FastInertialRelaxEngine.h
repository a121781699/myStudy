#ifndef FASTINERTIALRELAXENGINE
#define FASTINERTIALRELAXENGINE
#include"../../source/VectorMath.h"
#include"../../source/Particle.h"
#include"../../source/Box.h"
#include"EulerIntegral.h"

class FIRE :public EulerIntegral
{
public:
	FIRE(int argc, char *argv[]) :
		EulerIntegral(argc, argv)
	{
		readCmdInfo(&p, &b, argv[1]);
		dt_max = 100 * dt_step;
		alpha = alpha0;
		steps = 0;
	}

	void energyMinimize();
	void calcFandV();

	~FIRE();

protected:
	const double alpha0 = 0.1;
	const double f_inc = 1.1;
	const double f_dec = 0.5;
	const double f_alpha = 0.99;
	const int N_min = 5;

	double dt_max;
	double alpha;
	int steps;
private:
	Particle *p;
	Box *b;
};

void FIRE::energyMinimize()
{


}
void FIRE::calcFandV()
{
	calcForce(calcForcei);
	double Pmomentum = 0.0;
	double V_abs = 0.0;
	for (int i = 0; i < b->nAtoms; i++)
	{
		double fDotv, v2;
		vecDot(fDotv, p[i].v, p[i].f);
		vecDot(v2, p[i].v, p[i].v);
		Pmomentum += fDotv;
		V_abs += v2;
	}

	if (Pmomentum > 0)
	{
		steps++;
		if (steps > N_min)
		{
			dt_step = ((dt_step*f_inc) > dt_max ? dt_max : (dt_step*f_inc));
			alpha *= f_alpha;
		}

		double3 f_scale;
		V_abs = sqrt(V_abs);
		for (int i = 0; i < b->nAtoms; i++)
		{
			vecScale(p[i].v, 1 - alpha, p[i].v);
			vecScale(f_scale, alpha*V_abs, p[i].f);
			vecAdd(p[i].v, f_scale, p[i].v);
		}
	}
	else
	{
		dt_step *= f_dec;
		for (int i = 0; i < b->nAtoms; i++)
			p[i].v = { 0,0,0 };
		alpha = alpha0;
		steps = 0;
	}
}

FIRE::~FIRE()
{
	delete[]p;
	delete b;
}
#endif // !FASTINERTIALRELAXENGINE

