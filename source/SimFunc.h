#ifndef SIMFUNC_H
#define SIMFUNC_H

#include"VectorMath.h"
#include"Particle.h"
#include"Box.h"

void readCmdInfo(Box **box, Particle **particle, char *filename)
{
     FILE *fp=fopen(filename,"r");
     if(!fp)
     {
         fprintf(stderr,"Can't open this file!\n");
	 return;
     }

     char buff[256];
     for(int i=0;i<2;++i)
         fgets(buff,256,fp);

     *box = new Box;
     fscanf(fp, "%d", &(*box)->nAtoms); fgets(buff,256,fp);
     *particle = new Particle[(*box)->nAtoms];

	 for (int i = 0; i < 2; ++i)
		 fgets(buff, 256, fp);

	 //size of box
	 fscanf(fp, "lf", &(*box)->boxLo.x); fscanf(fp, "lf", &(*box)->boxHi.x); fgets(buff, 256, fp);
	 fscanf(fp, "lf", &(*box)->boxLo.y); fscanf(fp, "lf", &(*box)->boxHi.y); fgets(buff, 256, fp);
	 fscanf(fp, "lf", &(*box)->boxLo.z); fscanf(fp, "lf", &(*box)->boxHi.z); fgets(buff, 256, fp);

	 for (int i = 0; i < 4; ++i)
		 fgets(buff, 256, fp);

	 //particles
     for(int i=0;i<(*box)->nAtoms;++i)
     {
	     fscanf(fp,"%d",&(*particle)[i].id);   fgetc(fp);
         fscanf(fp,"%d",&(*particle)[i].type); fgetc(fp);
         fscanf(fp,"%lf",&(*particle)[i].DiameterScale); fgetc(fp);
         fscanf(fp,"%d",&(*particle)[i].density); fgetc(fp);
         fscanf(fp,"%lf",&(*particle)[i].x.x);  fgetc(fp);
         fscanf(fp,"%lf",&(*particle)[i].x.y);  fgetc(fp);
         fscanf(fp,"%lf",&(*particle)[i].x.z);  fgetc(fp);
         fscanf(fp,"%lf",&(*particle)[i].xImage.x); fgetc(fp);
         fscanf(fp,"%lf",&(*particle)[i].xImage.y); fgetc(fp);
         fscanf(fp,"%lf",&(*particle)[i].xImage.z); fgetc(fp);
     }
      
	 for (int i = 0; i < (*box)->nAtoms; ++i)
		 (*particle)[i].v = { 0,0,0 };
	 for (int i = 0; i < (*box)->nAtoms; ++i)
		 (*particle)[i].f = { 0,0,0 };

     Particle::meanDiameterScale=0.0;
     for(int i=0;i<(*box)->nAtoms;++i)
	  Particle::meanDiameterScale+=
		  ((*particle)[i].DiameterScale/(*box)->nAtoms);

     for(int i=0;i<(*box)->nAtoms;++i)
		 (*particle)[i].DiameterScale=
		  (*particle)[i].DiameterScale/Particle::meanDiameterScale;

	 defVar();
}

void safe_exit()
{
    delete box;
    delete []particle;
}

#endif
