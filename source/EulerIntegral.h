#ifndef EULERINTE
#define EULERINTE
#include"../../source/VectorMath.h"
#include"../../source/Particle.h"
#include"../../source/Box.h"

class EulerIntegral
{
public:
	EulerIntegral(int argc, char *argv[])
	{
		readCmdInfo(&p, &b, argv[1]);
		dt_step = _clock.timeScale;
	}

	void readCmdInfo(Particle **particle, Box **box, char *filename);
	void melocularMove()
	{
		calcForce(calcForcei);
		for (int i = 0; i < b->nAtoms; i++)
		{
			vecScaleAdd(p[i].x, p[i].x, dt_step, p[i].v);
			vecScaleAdd(p[i].v, p[i].v, dt_step, p[i].f);
		}
	}
	void adjustImg()
	{
		for (int i = 0; i < b->nAtoms; ++i)
		{
			if (b->boxLo.x > p[i].x.x) p[i].xImage.x = 1;
			else if (b->boxHi.x < p[i].x.x) p[i].xImage.x = -1;
			else p[i].xImage.x = 0;
			if (b->boxLo.y > p[i].x.y) p[i].xImage.y = 1;
			else if (b->boxHi.y > p[i].x.y) p[i].xImage.y = -1;
			else p[i].xImage.y = 0;
			if (b->boxLo.z > p[i].x.z) p[i].xImage.z = 1;
			else if (b->boxHi.z > p[i].x.z) p[i].xImage.z = -1;
			else p[i].xImage.z = 0;
		}

		imgFlag = false;
		for (int i = 0; i < b->nAtoms; ++i)
		{
			if (p[0].xImage.x != 0) imgFlag = true;
			if (p[0].xImage.y != 0) imgFlag = true;
			if (p[0].xImage.z != 0) imgFlag = true;
		}

		if (imgFlag)
		{
			for (int i = 0; i < b->nAtoms; ++i)
			{
				p[i].x.x = p[i].x.x + (b->boxHi.x - b->boxLo.x)*p[i].xImage.x;
				p[i].x.y = p[i].x.y + (b->boxHi.y - b->boxLo.y)*p[i].xImage.y;
				p[i].x.z = p[i].x.z + (b->boxHi.z - b->boxLo.z)*p[i].xImage.z;
			}
		}
	}
	void calcForce(void(*func)(double3&, double const, double const, double3 const));
	void printInfo();
	void fprintInfo(char const *filename);
	virtual ~EulerIntegral();
protected:
	static bool imgFlag;
	double dt_step;
private:
	Particle *p;
	Box *b;
};
bool EulerIntegral::imgFlag = false;
void EulerIntegral::readCmdInfo(Particle **particle, Box **box, char *filename)
{
	FILE *fp = fopen(filename, "r");
	if (!fp)
	{
		fprintf(stderr, "Can't open this file!\n");
		return;
	}

	char buff[256];
	for (int i = 0; i < 2; ++i)
		fgets(buff, 256, fp);

	*box = new Box;
	fscanf(fp, "%d", &(*box)->nAtoms); fgets(buff, 256, fp);
	*particle = new Particle[(*box)->nAtoms];

	for (int i = 0; i < 2; ++i)
		fgets(buff, 256, fp);

	//size of box
	fscanf(fp, "%lf", &(*box)->boxLo.x); fscanf(fp, "%lf", &(*box)->boxHi.x); fgets(buff, 256, fp);
	fscanf(fp, "%lf", &(*box)->boxLo.y); fscanf(fp, "%lf", &(*box)->boxHi.y); fgets(buff, 256, fp);
	fscanf(fp, "%lf", &(*box)->boxLo.z); fscanf(fp, "%lf", &(*box)->boxHi.z); fgets(buff, 256, fp);

	for (int i = 0; i < 4; ++i)
		fgets(buff, 256, fp);

	//particles
	for (int i = 0; i < (*box)->nAtoms; ++i)
	{
		fscanf(fp, "%d", &(*particle)[i].id);   fgetc(fp);
		fscanf(fp, "%d", &(*particle)[i].type); fgetc(fp);
		fscanf(fp, "%lf", &(*particle)[i].DiameterScale); fgetc(fp);
		fscanf(fp, "%d", &(*particle)[i].density); fgetc(fp);
		fscanf(fp, "%lf", &(*particle)[i].x.x);  fgetc(fp);
		fscanf(fp, "%lf", &(*particle)[i].x.y);  fgetc(fp);
		fscanf(fp, "%lf", &(*particle)[i].x.z);  fgetc(fp);
		fscanf(fp, "%lf", &(*particle)[i].xImage.x); fgetc(fp);
		fscanf(fp, "%lf", &(*particle)[i].xImage.y); fgetc(fp);
		fscanf(fp, "%lf", &(*particle)[i].xImage.z); fgetc(fp);
	}

	for (int i = 0; i < (*box)->nAtoms; ++i)
		(*particle)[i].v = { 0,0,0 };
	for (int i = 0; i < (*box)->nAtoms; ++i)
		(*particle)[i].f = { 0,0,0 };

	Particle::meanDiameterScale = 0.0;
	for (int i = 0; i < (*box)->nAtoms; ++i)
		Particle::meanDiameterScale +=
		((*particle)[i].DiameterScale / (*box)->nAtoms);

	for (int i = 0; i < (*box)->nAtoms; ++i)
		(*particle)[i].DiameterScale =
		(*particle)[i].DiameterScale / Particle::meanDiameterScale;

	defVar();
}

void EulerIntegral::calcForce(void(*func)(double3&, double const, double const, double3 const))
{
	double xNorm;
	double3 f, x;
	for (int i = 0; i < b->nAtoms; ++i)
	{
		double3 fsum = { 0,0,0 };
		for (int j = 0; j < b->nAtoms; ++j)
		{
			if (j == i) continue;
			vecSub(x, p[j].x, p[i].x);
			vecNorm(xNorm, x);
			double sigma = (p[i].DiameterScale + p[j].DiameterScale)*Particle::meanDiameterScale / 2.0;
			if (sigma <= xNorm) continue;
			double factor = -2 * var.epsilon / sigma;
			func(f, factor, 1 - xNorm / sigma, x);
			vecAdd(fsum, fsum, f);
			//vecSub(particle[i].f, particle[i].f, f);
		}
		p[i].f = fsum;
	}
}

void EulerIntegral::printInfo() 
{
	fprintf(stdout, "nAtoms: %d\n"
		"xLo: %lf yLo: %lf zLo: %lf\n",
		b->nAtoms,
		b->boxLo.x, b->boxLo.y, b->boxLo.z);
	for (int i = 0; i < b->nAtoms; ++i)
		fprintf(stdout, "%lf %lf %lf\n",
			p[i].x.x, p[i].x.y, p[i].x.z);

	//Particle::calcForce(calcForcei, particle, box->nAtoms);
	//Particle::eulerInteger(particle, box->nAtoms);

	fprintf(stdout, "nAtoms: %d\n", b->nAtoms);
	for (int i = 0; i < b->nAtoms; ++i)
		fprintf(stdout, "%lf %lf %lf\n",
			p[i].f.x, p[i].f.y, p[i].f.z);
}

void EulerIntegral::fprintInfo(char const *filename)
{
	FILE *fp = fopen(filename, "w");

	fprintf(fp, "nAtoms: %d\n"
		"xLo: %lf yLo: %lf zLo: %lf\n",
		b->nAtoms,
		b->boxLo.x, b->boxLo.y, b->boxLo.z);
	for (int i = 0; i < b->nAtoms; ++i)
		fprintf(fp, "%lf %lf %lf\n",
			p[i].x.x, p[i].x.y, p[i].x.z);

	//Particle::calcForce(calcForcei, particle, box->nAtoms);
	//Particle::eulerInteger(particle, box->nAtoms);

	fprintf(fp, "nAtoms: %d\n", b->nAtoms);
	for (int i = 0; i < b->nAtoms; ++i)
		fprintf(fp, "%lf %lf %lf\n",
			p[i].f.x, p[i].f.y, p[i].f.z);
}

EulerIntegral::~EulerIntegral()
{
	delete[]p;
	delete b;
}
#endif // !EULERINTE

