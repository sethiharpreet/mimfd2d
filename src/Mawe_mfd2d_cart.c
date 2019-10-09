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
  bool dabc;   /* absorbing boundary */

  int  jsnap,ntsnap,jdata; /* wavefield parameters*/

  /* I/O files */
  sf_file Fwav  = NULL; /* wavelet   */
  sf_file Fwfl  = NULL; /* wavefield */
  sf_file Fsou  = NULL; /* sources   */
  sf_file Frec  = NULL; /* receivers */

  sf_file Fvel  = NULL;
  sf_file Fden  = NULL;

  sf_file Fdat  = NULL; /* data      */


  /* cube axes */
  sf_axis at,ax,az;
  sf_axis as,ar,ac;

  int     nt,nz,nx,ns,nr,nc,nb;
  int     it,iz,ix;
  float   dt,dz,dx,idz,idx,dtx,dtz;
  int     nzpad,nxpad;

  /* FDM structure */
  fdm2d    fdm=NULL;

  /* I/O arrays */
  float***ww=NULL;           /* wavelet   */
  pt2d   *ss=NULL;           /* sources   */
  pt2d   *rr=NULL;           /* receivers */
  float **dd=NULL;           /* data      */

  float  **tt  = NULL;       /* Layer 1 */

  /* acoustic */

  float ** vp = NULL;
  float ** ro = NULL;

  /* Temporary variables */
  float **t11a, **t12a, **t21a, **t22a;

  /* Layer 1*/
  float **vx_ff, **vz_ff;
  float **vx_vv, **vz_vv;

  float **p_vf, **p_fv;

  // combined wavefields
  float **wflx, **wflz;

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
  if(! sf_getbool("dabc",&dabc)) dabc=false; /* absorbing boundary */

  /* I/O files */
  Fwav = sf_input ("in"); /* wavelet   */

  Fvel = sf_input("vel"); /* velocity */
  Fden = sf_input("den"); /* density  */

  Fsou = sf_input ("sou"); /* sources   */
  Frec = sf_input ("rec"); /* receivers */
  Fdat = sf_output("out"); /* data      */
  if(snap){
      Fwfl = sf_output("wfl"); /* wavefield */
  }

  /* axes */
  at = sf_iaxa(Fwav,3); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* time */

  az = sf_iaxa(Fvel,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* depth */
  ax = sf_iaxa(Fvel,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* space */

  as = sf_iaxa(Fsou,2); sf_setlabel(as,"s"); if(verb) sf_raxa(as); /* sources */
  ar = sf_iaxa(Frec,2); sf_setlabel(ar,"r"); if(verb) sf_raxa(ar); /* receivers */

  nt = sf_n(at); dt = sf_d(at);
  nz = sf_n(az); dz = sf_d(az);
  nx = sf_n(ax); dx = sf_d(ax);

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

  fdm  = fdutil_init(verb,fsrf,az,ax,nb,1); /* Layer 1 --> Top */

  fdbell_init(nbell);
  /* Layer 1*/
  nzpad=fdm->nzpad;
  nxpad=fdm->nxpad;

  /*------------------------------------------------------------*/

  /* setup output wavefield header */
  if(snap) {
    if(!sf_getint  ("nqz",&nqz)) nqz=nzpad;
    if(!sf_getint  ("nqx",&nqx)) nqx=nxpad;

    if(!sf_getfloat("oqz",&oqz)) oqz=sf_o(az);
    if(!sf_getfloat("oqx",&oqx)) oqx=sf_o(ax);

    dqz=sf_d(az);
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
  tt  = sf_floatalloc2(nz,nx);

  /* Layer 1*/
  ro = sf_floatalloc2(nzpad+2,nxpad+2);
  vp = sf_floatalloc2(nzpad+2,nxpad+2);

  /* Temporary variables */
  t11a = sf_floatalloc2(nzpad+1,nxpad+1); memset(t11a[0],0,(nzpad+1)*(nxpad+1)*sizeof(float));
  t12a = sf_floatalloc2(nzpad+1,nxpad+2); memset(t12a[0],0,(nzpad+1)*(nxpad+2)*sizeof(float));
  t21a = sf_floatalloc2(nzpad+2,nxpad+1); memset(t21a[0],0,(nzpad+2)*(nxpad+1)*sizeof(float));
  t22a = sf_floatalloc2(nzpad+2,nxpad+2); memset(t22a[0],0,(nzpad+2)*(nxpad+2)*sizeof(float));

  /* input density */
  sf_floatread(tt[0],nz*nx,Fden); expand2d(tt,ro,fdm);

  /* input density */
  sf_floatread(tt[0],nz*nx,Fvel); expand2d(tt,vp,fdm);


  free(*tt); free(tt);


  /* Field Variables */

  /* Velocity */
  /* [f,f] grid */
  vx_ff  = sf_floatalloc2(nzpad+2,nxpad+2); memset(vx_ff[0],0,(nzpad+2)*(nxpad+2)*sizeof(float));
  vz_ff  = sf_floatalloc2(nzpad+2,nxpad+2); memset(vz_ff[0],0,(nzpad+2)*(nxpad+2)*sizeof(float));

  /* [v,v] grid */
  vx_vv  = sf_floatalloc2(nzpad+1,nxpad+1); memset(vx_vv[0],0,(nzpad+1)*(nxpad+1)*sizeof(float));
  vz_vv  = sf_floatalloc2(nzpad+1,nxpad+1); memset(vz_vv[0],0,(nzpad+1)*(nxpad+1)*sizeof(float));

  /* Pressure */
  /* [v,f] grid */
  p_vf   = sf_floatalloc2(nzpad+1,nxpad+2); memset(p_vf[0],0,(nzpad+1)*(nxpad+2)*sizeof(float));

  /* [f,v] */
  p_fv   = sf_floatalloc2(nzpad+2,nxpad+1); memset(p_fv[0],0,(nzpad+2)*(nxpad+1)*sizeof(float));

  /*------------------------------------------------*/

  /* combined grid */

  wflz    = sf_floatalloc2(nzpad,nxpad); memset(wflz[0],0,nzpad*nxpad*sizeof(float));
  wflx    = sf_floatalloc2(nzpad,nxpad); memset(wflx[0],0,nzpad*nxpad*sizeof(float));


  /*------------------------------------------------------------*/
  /*
   *  MAIN LOOP
   */
  /*------------------------------------------------------------*/
  if(verb) fprintf(stderr,"\n");

  for (it=0; it<nt; it++) {

    if(verb) fprintf(stderr,"\b\b\b\b\b%d",it);

    /* Update pressure */
    update_pressure_fv(p_fv, vx_ff, vz_vv, t21a, ro, vp, dtx, dtz, nzpad, nxpad);
    update_pressure_vf(p_vf, vz_ff, vx_vv, t12a, ro, vp, dtx, dtz, nzpad, nxpad);

    /* Source Injection */
    inject_vel_src_cart(vz_ff, vz_vv, ww, cs, it);

    /* Update velocity */
    update_ac_velocity_ff(vz_ff, vx_ff, p_vf, p_fv, t22a, ro, dtx, dtz, nzpad, nxpad);
    update_ac_velocity_vv(vz_vv, vx_vv, p_vf, p_fv, t11a, ro, dtx, dtz, nzpad, nxpad);


    //fprintf(stderr, "%d\n", nx);
    /* */
    if(snap && it%jsnap==0) {
      /* extract wavefield on fsg */
      extract_wfld_fsg(wflz, vz_ff, vz_vv, nxpad, nzpad);
      //cut2dfree(wflz,uc,fdm,acz,acx);

      sf_floatwrite(wflz[0],nzpad*nxpad,Fwfl);

      extract_wfld_fsg(wflx,vx_ff, vx_vv, nxpad, nzpad);
      //cut2dfree(wflx,uc,fdm,acz,acx);

      sf_floatwrite(wflx[0],nzpad*nxpad,Fwfl);

      }

  }

}
