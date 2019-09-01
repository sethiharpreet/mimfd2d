#include <rsf.h>
#include <math.h>

#include "fdutil.h"
#include "mfd2d_utils.h"

#define NOP 4 /* derivative operator half-size */

int main(int argc, char* argv[]){

  bool verb ;  /* verbose  */
  bool ssou;   /* stress source */
  bool snap;   /* snapshot */
  bool fsrf;   /* free surface */
  bool dabc;   /* Absorbing boundary */

  int  jsnap,ntsnap,jdata; /* wavefield parameters*/

  /* I/O files */
  sf_file Fwav  = NULL; /* wavelet   */
  sf_file Fwfl  = NULL; /* wavefield */
  sf_file Fsou  = NULL; /* sources   */
  sf_file Frec  = NULL; /* receivers */

  sf_file Fvel = NULL; /* velocity 1st layer */
  sf_file Fden = NULL; /* density  1st layer */

  sf_file Fccc  = NULL; /* stiffness 2nd layer */
  sf_file Fden1  = NULL; /* density 2nd layer */

  sf_file Fdat  = NULL; /* data      */


  /* cube axes */
  sf_axis at,ax,az,az1;
  sf_axis as,ar,ac;

  int     nt,nz,nx,ns,nr,nc,nb; /* Layer 1*/
  int     it,iz,ix;
  float   dt,dz,dx,dtx,dtz;
  int     nzpad,nxpad;
  int     nz1,nx1;  /* Layer 2 */
  int    nz1pad;

  /* FDM structure */
  fdm2d    fdm=NULL, fdm1=NULL;
  abcone2d abcp=NULL,abcs=NULL;
  sponge   spo=NULL;


  /* I/O arrays */
  float***ww=NULL;           /* wavelet   */
  pt2d   *ss=NULL;           /* sources   */
  pt2d   *rr=NULL;           /* receivers */
  float **dd=NULL;           /* data      */

  float  **tt  = NULL;       /* Layer 1 */
  float  *  *tt1 = NULL;       /* Layer 2 */


  /* Layer 1 */

  /* acoustic */

  float ** vp = NULL;
  float ** ro = NULL;

  /* Layer 2 */

  /* orthorombic footprint - 4 coefficients */
  /* c11 c13
     .   c33
             c55 */

  float **c11 = NULL;
  float **c33 = NULL;
  float **c55 = NULL;
  float **c13 = NULL;

  float **ro1 = NULL;


  /* Temporary variables */
  float **t11, **t12, **t21, **t22;

  /* Layer 1*/

  float **vx_ff, **vz_ff;
  float **vx_vv, **vz_vv;

  float **p_vf, **p_fv;

  /* Layer 2*/
  float **vx1_ff, **vz1_ff;
  float **vx1_vv, **vz1_vv;

  float **sxx1_vf, **sxz1_vf, **szz1_vf;
  float **sxx1_fv, **sxz1_fv, **szz1_fv;

  // combined grids
  float **wfl1, **wfl2;

  /* linear interpolation weights/indices */
  lint2d cs,cr;

  /* Gaussian bell */
  int nbell;

  /* wavefield cut params */
  sf_axis   acz=NULL, acx=NULL;
  int       nqz,nqx;
  float     oqz,oqx;
  float     dqz,dqx;
  float     **uc=NULL;

  /*------------------------------------------------------------*/
  /* init RSF */
  sf_init(argc,argv);

  /* execution flags */
  if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
  if(! sf_getbool("snap",&snap)) snap=false; /* wavefield snapshots flag */
  if(! sf_getbool("free",&fsrf)) fsrf=false; /* free surface flag */
  if(! sf_getbool("ssou",&ssou)) ssou=false; /* stress source */
  if(! sf_getbool("dabc",&dabc)) dabc=false; /* absorbing BC */

  /* I/O files */
  Fwav = sf_input ("in" ); /* wavelet   */

  Fvel = sf_input("vel"); /* velocity 1st layer */
  Fden = sf_input("den"); /* density 1st layer */

  Fccc = sf_input ("ccc"); /* stiffness 2nd layer */
  Fden1 = sf_input ("den1"); /* density  2nd layer */

  Fsou = sf_input ("sou"); /* sources   */
  Frec = sf_input ("rec"); /* receivers */
  Fdat = sf_output("out"); /* data      */
  if(snap){
      Fwfl = sf_output("wfl"); /* wavefield */
  }

  /* axes */
  at = sf_iaxa(Fwav,3); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* time */
  ax = sf_iaxa(Fccc,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* space */
  az = sf_iaxa(Fccc,1); sf_setlabel(az1,"z"); if(verb) sf_raxa(az); /* depth: Layer 2 */

  az1 = sf_iaxa(Fvel,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az1); /* depth: Layer 1 */

  as = sf_iaxa(Fsou,2); sf_setlabel(as,"s"); if(verb) sf_raxa(as); /* sources */
  ar = sf_iaxa(Frec,2); sf_setlabel(ar,"r"); if(verb) sf_raxa(ar); /* receivers */

  nt = sf_n(at); dt = sf_d(at);
  nz = sf_n(az); dz = sf_d(az);
  nx = sf_n(ax); dx = sf_d(ax);

  nz1 = sf_n(az1);  /* Layer 1 dx1=dx */

  ns = sf_n(as);
  nr = sf_n(ar);

  /*------------------------------------------------------------*/
  /* other execution parameters */
  if(! sf_getint("nbell",&nbell)) nbell=5;  /* bell size */
  if(verb) sf_warning("nbell=%d",nbell);
  if(! sf_getint("jdata",&jdata)) jdata=1;
  if(snap) {  /* save wavefield every *jsnap* time steps */
    if(! sf_getint("jsnap",&jsnap)) jsnap=nt;
  }
  /*------------------------------------------------------------*/

  /* 2D vector components */
  nc=2;
  ac=sf_maxa(nc,0,1);
  int ncs=3; /* 3-component source */
  /*------------------------------------------------------------*/
  /* setup output data header */
  sf_oaxa(Fdat,ar,1);
  sf_oaxa(Fdat,ac,2);

  sf_setn(at,nt/jdata);
  sf_setd(at,dt*jdata);
  sf_oaxa(Fdat,at,3);

  /*------------------------------------------------------------*/
  /* expand domain for FD operators and ABC */
  if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP;

  fdm  = fdutil_init(verb,fsrf,true,az,ax,nb,1); /* Layer 1 --> Top */
  fdm1 = fdutil_init(verb,false,false,az1,ax,nb,1); /* Layer 2 --> Bottom */

  fdbell_init(nbell);
  /* Layer 1*/
  nzpad=fdm->nzpad;
  nxpad=fdm->nxpad;

  /* Layer 2*/
  nz1pad= fdm1->nzpad;

  /*------------------------------------------------------------*/

  /* setup output wavefield header */
  if(snap) {
    if(!sf_getint  ("nqz",&nqz)) nqz=sf_n(az)+sf_n(az1); // nz+nz1
    if(!sf_getint  ("nqx",&nqx)) nqx=sf_n(ax);

    if(!sf_getfloat("oqz",&oqz)) oqz=sf_o(az);
    if(!sf_getfloat("oqx",&oqx)) oqx=sf_o(ax);

    dqz=sf_d(az); // dz=dz1
    dqx=sf_d(ax);

    acz = sf_maxa(nqz,oqz,dqz); sf_raxa(acz);
    acx = sf_maxa(nqx,oqx,dqx); sf_raxa(acx);
    uc=sf_floatalloc2(sf_n(acz),sf_n(acx));

    ntsnap=0;
    for(it=0; it<nt; it++) {
        if(it%jsnap==0) ntsnap++;
    }

    sf_setn(at,  ntsnap);
    sf_setd(at,dt*jsnap);
    if(verb) sf_raxa(at);

      /* wavefields */
    sf_oaxa(Fwfl,acz,1);
    sf_oaxa(Fwfl,acx,2);
    sf_oaxa(Fwfl,ac, 3);
    sf_oaxa(Fwfl,at, 4);
  }

  /* source array */
  ww=sf_floatalloc3(ns,ncs,nt);
  sf_floatread(ww[0][0],nt*ncs*ns,Fwav);

  /* data array */
  dd=sf_floatalloc2(nr,nc);

  /*------------------------------------------------------------*/
  /* setup source/receiver coordinates */
  ss = (pt2d*) sf_alloc(ns,sizeof(*ss));
  rr = (pt2d*) sf_alloc(nr,sizeof(*rr));

  pt2dread1(Fsou,ss,ns,2); /* read (x,z) coordinates */
  pt2dread1(Frec,rr,nr,2); /* read (x,z) coordinates */

  cs = lint2d_make(ns,ss,fdm);
  cr = lint2d_make(nr,rr,fdm);

  /*------------------------------------------------------------*/
  /* setup FD coefficients */
  dtz = dt/dz;
  dtx = dt/dx;

  /*------------------------------------------------------------*/

  /* Temporary variable for input paramter reading */
  tt  = sf_floatalloc2(nz,nx);  /* Layer 1 */
  tt1 = sf_floatalloc2(nz1,nx); /* Layer 2 */


  /* Layer 1*/
  ro = sf_floatalloc2(nzpad+2,nxpad+2);
  vp = sf_floatalloc2(nzpad+2,nxpad+2);

  /* input density */
  sf_floatread(tt[0],nz*nx,Fden); expand2d_top(tt,ro,fdm);

  /* input density */
  sf_floatread(tt[0],nz*nx,Fvel); expand2d_top(tt,vp,fdm);


  free(*tt); free(tt);

  /* Layer 2*/
  ro1   = sf_floatalloc2(nz1pad+2,nxpad+2);
  c11   = sf_floatalloc2(nz1pad+2,nxpad+2);
  c33   = sf_floatalloc2(nz1pad+2,nxpad+2);
  c55   = sf_floatalloc2(nz1pad+2,nxpad+2);
  c13   = sf_floatalloc2(nz1pad+2,nxpad+2);

  /* Temporary variables */
  t11 = sf_floatalloc2(nz1pad+1,nxpad+1);
  t12 = sf_floatalloc2(nz1pad+1,nxpad+2);
  t21 = sf_floatalloc2(nz1pad+2,nxpad+1);
  t22 = sf_floatalloc2(nz1pad+2,nxpad+2);

  /* input density */
  sf_floatread(tt1[0],nz*nx,Fden1);    expand2d_bottom(tt1,ro1,fdm1);

  /* input stiffness */
  sf_floatread(tt1[0],nz*nx,Fccc);    expand2d_bottom(tt1,c11,fdm1);
  sf_floatread(tt1[0],nz*nx,Fccc);    expand2d_bottom(tt1,c33,fdm1);
  sf_floatread(tt1[0],nz*nx,Fccc);    expand2d_bottom(tt1,c55,fdm1);
  sf_floatread(tt1[0],nz*nx,Fccc);    expand2d_bottom(tt1,c13,fdm1);

  free(*tt1); free(tt1);

  /* Field Variables */

  /* Layer 1 */

  /* Velocity */
  /* [f,f] grid */
  vx_ff  = sf_floatalloc2(nzpad+2,nxpad+2);
  vz_ff  = sf_floatalloc2(nzpad+2,nxpad+2);

  /* [v,v] grid */
  vx_vv  = sf_floatalloc2(nzpad+1,nxpad+1);
  vz_vv  = sf_floatalloc2(nzpad+1,nxpad+1);

  /* Pressure */
  /* [v,f] grid */
  p_vf   = sf_floatalloc2(nzpad+1,nxpad+2);

  /* [f,v] */
  p_fv   = sf_floatalloc2(nzpad+2,nxpad+1);

  /*------------------------------------------------*/

  /* Layer 2 */

  /* Velocity */
  /* [f,f] grid */
  vx1_ff   = sf_floatalloc2(nz1pad+2,nxpad+2);
  vz1_ff   = sf_floatalloc2(nz1pad+2,nxpad+2);

  /* [v,v] grid */
  vx1_vv   = sf_floatalloc2(nz1pad+1,nxpad+1);
  vz1_vv   = sf_floatalloc2(nz1pad+1,nxpad+1);


  /* Stress */
  /* [v,f] grid */
  sxx1_vf  = sf_floatalloc2(nz1pad+1,nxpad+2);
  sxz1_vf  = sf_floatalloc2(nz1pad+1,nxpad+2);
  szz1_vf  = sf_floatalloc2(nz1pad+1,nxpad+2);

  /* [f,v] grid */

  sxx1_fv  = sf_floatalloc2(nz1pad+2,nxpad+1);
  sxz1_fv  = sf_floatalloc2(nz1pad+2,nxpad+1);
  szz1_fv  = sf_floatalloc2(nz1pad+2,nxpad+1);

  /* combined grid */

  wfl1    = sf_floatalloc2(nz1pad+nzpad,nxpad);
  wfl2    = sf_floatalloc2(nz1pad+nzpad,nxpad);

  /*------------------------------------------------------------*/
  /*
   *  MAIN LOOP
   */
  /*------------------------------------------------------------*/
  if(verb) fprintf(stderr,"\n");

  for (it=0; it<nt; it++) {

    if(verb) fprintf(stderr,"\b\b\b\b\b%d",it);

    /* Layer 1 */
    update_pressure_fv(p_fv, vx_ff, vz1_vv, t21, ro1, vp, dtx, dtz, nzpad, nxpad);
    update_pressure_vf(p_vf, vz1_ff, vx1_vv, t12, ro1, vp, dtx, dtz, nzpad, nxpad);

    update_ac_velocity_ff(vz_ff, vx_ff, p_vf, p_fv, t22, ro, dtx, dtz, nzpad, nxpad);
    update_ac_velocity_vv(vz_vv, vx_vv, p_vf, p_fv, t11, ro, dtx, dtz, nzpad, nxpad);

    /* Layer 2 */
    update_stress_vf(vz1_ff, vx1_ff, vz1_vv, vx1_vv, sxx1_vf, sxz1_vf, szz1_vf, ro1, c11, c13,   \
                     c33, c55, t12, dtx, dtz, nz1pad, nxpad);

    update_stress_fv(vz1_ff, vx1_ff, vz1_vv, vx1_vv, sxx1_fv, sxz1_fv, szz1_fv, ro1, c11, c13,   \
                     c33, c55, t12, dtx, dtz, nz1pad, nxpad);


    /* Coupling conditions */

    /* source injection */

  }




}
