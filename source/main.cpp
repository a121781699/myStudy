#include<time.h>

#include"VectorMath.h"
#include"Particle.h"
#include"Box.h"
#include"SimFunc.h"

void foo()
{
    fprintf(stdout,"nAtoms: %d\n",box->nAtoms);
    for(int i=0;i<box->nAtoms;++i) 
	fprintf(stdout,"%lf %lf %lf\n",
		particle[i].x.x,particle[i].x.y,particle[i].x.z);
}

int main(int argc, char *argv[])
{
    readCmdInfo(&box,&particle,argv[1]);

    foo();

    safe_exit();
    return 0;
}