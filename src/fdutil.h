/* This file is automatically generated. DO NOT EDIT! */

#ifndef _fdutil_h
#define _fdutil_h


typedef struct fdm2 *fdm2d;


typedef struct fdm3 *fdm3d;


typedef struct lcoef2 *lint2d;


typedef struct lcoef3 *lint3d;


typedef struct abc2 *abcone2d;


typedef struct abc3 *abcone3d;


typedef struct spon *sponge;


typedef struct ofg *ofg2d;


struct fdm2{
    int nb;
    int   nz,nzpad;
    int   nx,nxpad;
    float oz,ozpad;
    float ox,oxpad;
    float dz;
    float dx;
    bool verb;
    bool top;
    bool free;
    int ompchunk;
};


struct fdm3{
    int nb;
    int   nz,nzpad;
    int   nx,nxpad;
    int   ny,nypad;
    float oz,ozpad;
    float ox,oxpad;
    float oy,oypad;
    float dz;
    float dx;
    float dy;
    bool verb;
    bool free;
    int ompchunk;
};


struct lcoef2{
    int n;
    float *w00;
    float *w01;
    float *w10;
    float *w11;
    int *jz;
    int *jx;
};


struct lcoef3{
    int n;
    float *w000;
    float *w001;
    float *w010;
    float *w011;
    float *w100;
    float *w101;
    float *w110;
    float *w111;
    int *jz;
    int *jx;
    int *jy;
};


struct abc2{
    bool free;
    float *bzl;
    float *bzh;
    float *bxl;
    float *bxh;
};


struct abc3{
    bool free;
    float**bzl;
    float**bzh;
    float**bxl;
    float**bxh;
    float**byl;
    float**byh;
};


struct spon{
    float *w;
};


struct ofg{
    float **tt;
};


/*------------------------------------------------------------*/
fdm2d fdutil_init(bool verb_,
		  bool free_,
		  sf_axis az_,
		  sf_axis ax_,
		  int     nb_,
		  int ompchunk_);
/*< init fdm utilities >*/


/*------------------------------------------------------------*/
fdm3d fdutil3d_init(bool verb_,
		    bool free_,
		    sf_axis az_,
		    sf_axis ax_,
		    sf_axis ay_,
		    int     nb_,
		    int ompchunk_);
/*< init fdm utilities >*/


/*------------------------------------------------------------*/
ofg2d offgrid_init(fdm2d fdm);
/*< init off-grid interpolation >*/


/*------------------------------------------------------------*/
void offgridfor(float **ti,
		ofg2d  ofg,
		fdm2d  fdm);
/*< forward off-grid interpolation (in place) >*/


/*------------------------------------------------------------*/
void offgridadj(float **ti,
		ofg2d  ofg,
		fdm2d  fdm);
/*< adjoint off-grid interpolation (in place) >*/


/*------------------------------------------------------------*/
void expand(float** a,
	    float** b,
	    fdm2d fdm);
/*< expand domain >*/


/*------------------------------------------------------------*/
void expand3d(float ***a,
	      float ***b,
	      fdm3d  fdm);
/*< expand domain >*/


/*------------------------------------------------------------*/
void cut2d(float**  a,
	   float**  b,
	   fdm2d  fdm,
	   sf_axis c1,
	   sf_axis c2);
/*< cut a rectangular wavefield subset >*/


/*------------------------------------------------------------*/
void cut3d(float*** a,
	   float*** b,
	   fdm3d  fdm,
	   sf_axis c1,
	   sf_axis c2,
	   sf_axis c3);
/*< cut a rectangular wavefield subset >*/


/*------------------------------------------------------------*/
void bfill(float** b,
	   fdm2d fdm);
/*< fill boundaries >*/


/*------------------------------------------------------------*/
lint2d lint2d_make(int    na,
		   pt2d*  aa,
		   fdm2d fdm);
/*< init 2D linear interpolation >*/


/*------------------------------------------------------------*/
lint3d lint3d_make(int    na,
		   pt3d*  aa,
		   fdm3d fdm);
/*< init 3D linear interpolation >*/


/*------------------------------------------------------------*/
void lint2d_hold(float**uu,
		 float *ww,
		 lint2d ca);
/*< hold fixed value in field >*/


/*------------------------------------------------------------*/
void lint2d_inject(float**uu,
		   float *ww,
		   lint2d ca);
/*< inject into wavefield >*/


/*------------------------------------------------------------*/
void lint3d_inject(float***uu,
		   float  *ww,
		   lint3d  ca);
/*< inject into wavefield >*/


/*------------------------------------------------------------*/
void lint2d_inject1(float**uu,
		    float  ww,
		    lint2d ca);
/*< inject into wavefield >*/


/*------------------------------------------------------------*/
void lint3d_inject1(float***uu,
		    float   ww,
		    lint3d  ca);
/*< inject into wavefield >*/


/*------------------------------------------------------------*/
void lint2d_extract(float**uu,
		    float* dd,
		    lint2d ca);
/*< extract from wavefield >*/


void lint3d_extract(float***uu,
		    float  *dd,
		    lint3d  ca);
/*< extract from wavefield >*/


/*------------------------------------------------------------*/
void fdbell_init(int n);
/*< init bell taper >*/


/*------------------------------------------------------------*/
void fdbell3d_init(int n);
/*< init bell taper >*/


/*------------------------------------------------------------*/
void lint2d_bell(float**uu,
		 float *ww,
		 lint2d ca);
/*< apply bell taper >*/


/*------------------------------------------------------------*/
void lint3d_bell(float***uu,
		 float  *ww,
		 lint3d  ca);
/*< apply bell taper >*/


/*------------------------------------------------------------*/
abcone2d abcone2d_make(int     nop,
		       float    dt,
		       float**  vv,
		       bool   free,
		       fdm2d   fdm);
/*< init 2D ABC >*/


/*------------------------------------------------------------*/
abcone3d abcone3d_make(int     nop,
		       float    dt,
		       float ***vv,
		       bool   free,
		       fdm3d   fdm);
/*< init 3D ABC >*/


/*------------------------------------------------------------*/
void abcone2d_apply(float**   uo,
		    float**   um,
		    int      nop,
		    abcone2d abc,
		    fdm2d    fdm);
/*< apply 2D ABC >*/


/*------------------------------------------------------------*/
void abcone3d_apply(float  ***uo,
		    float  ***um,
		    int      nop,
		    abcone3d abc,
		    fdm3d    fdm);
/*< apply 3D ABC >*/


/*------------------------------------------------------------*/
sponge sponge_make(int nb);
/*< init boundary sponge >*/


/*------------------------------------------------------------*/
void sponge2d_apply(float**   uu,
		    sponge   spo,
		    fdm2d    fdm);
/*< apply boundary sponge >*/


void sponge2d_apply_test(float**   uu,
		    sponge   spo,
		    fdm2d    fdm);
/*< apply boundary sponge >*/


/*------------------------------------------------------------*/
void sponge3d_apply(float  ***uu,
		    sponge   spo,
		    fdm3d    fdm);
/*< apply boundary sponge >*/

#endif
