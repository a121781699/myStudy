#include"../Project1/Project1/EulerIntegral.h"
#include"../Project1/Project1/FastInertialRelaxEngine.h"
void foo()
{
    fprintf(stdout,"nAtoms: %d\n",box->nAtoms);
    for(int i=0;i<box->nAtoms;++i) 
	fprintf(stdout,"%lf %lf %lf\n",
		particle[i].x.x,particle[i].x.y,particle[i].x.z);

	//Particle::calcForce(calcForcei, particle, box->nAtoms);
	//Particle::eulerInteger(particle, box->nAtoms);

	fprintf(stdout, "nAtoms: %d\n", box->nAtoms);
	for (int i = 0; i < box->nAtoms; ++i)
		fprintf(stdout, "%lf %lf %lf\n",
			particle[i].f.x, particle[i].f.y, particle[i].f.z);
}

int main(int argc, char *argv[])
{
	FIRE assembly(argc, argv);
    //readCmdInfo(&box,&particle,argv[1]);

	//assembly.calcForce(calcForcei);
	assembly.melocularMove();
	assembly.printInfo();
    //foo();

    //safe_exit();
    return 0;
}
