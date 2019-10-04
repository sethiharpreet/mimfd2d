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
  sf_file Fccc  = NULL; /* stiffness */
  sf_file Fden  = NULL; /* density   */
  sf_file Fdat  = NULL; /* data      */


  /* cube axes */
  sf_axis at,ax,az;
  sf_axis as,ar,ac;

  int     nt,nz,nx,ns,nr,nc,nb; /* Layer 1*/
  int     it,iz,ix;
  float   dt,dz,dx,idz,idx,dtx,dtz;
  int     nzpad,nxpad;

  /* FDM structure */
  fdm2d    fdm=NULL;
  abcone2d abcp=NULL,abcs=NULL;
  sponge   spo=NULL;


  /* I/O arrays */
  float***ww=NULL;           /* wavelet   */
  pt2d   *ss=NULL;           /* sources   */
  pt2d   *rr=NULL;           /* receivers */
  float **dd=NULL;           /* data      */

  float  **tt  = NULL;

  /* orthorombic footprint - 4 coefficients */
  /* c11 c13
     .   c33
             c55 */

  float **c11 = NULL;
  float **c33 = NULL;
  float **c55 = NULL;
  float **c13 = NULL;
  float **vp=NULL,**vs=NULL;

  float **ro = NULL;

  /* Temporary variables */
  float **t11, **t12, **t21, **t22;

  float **vx_ff, **vz_ff;
  float **vx_vv, **vz_vv;

  float **sxx_vf, **sxz_vf, **szz_vf;
  float **sxx_fv, **sxz_fv, **szz_fv;

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
  if(! sf_getbool("dabc",&dabc)) dabc=false; /* absorbing BC */

  /* I/O files */
  Fwav = sf_input ("in" ); /* wavelet   */
  Fccc = sf_input ("ccc"); /* stiffness */
  Fden = sf_input ("den"); /* density   */
  Fsou = sf_input ("sou"); /* sources   */
  Frec = sf_input ("rec"); /* receivers */
  Fdat = sf_output("out"); /* data      */
  if(snap){
      Fwfl = sf_output("wfl"); /* wavefield */
  }

  /* axes */
  at = sf_iaxa(Fwav,3); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* time */
  ax = sf_iaxa(Fccc,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* space */
  az = sf_iaxa(Fccc,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* depth */

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

  fdm  = fdutil_init(verb,fsrf,az,ax,nb,1);

  fdbell_init(nbell);
  /* Layer */
  nzpad=fdm->nzpad;
  nxpad=fdm->nxpad;


  /*------------------------------------------------------------*/

  /* setup output wavefield header */
  if(snap) {
    if(!sf_getint  ("nqz",&nqz)) nqz=sf_n(az);
    if(!sf_getint  ("nqx",&nqx)) nqx=sf_n(ax);

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
  idz = 1/dz;
  idx = 1/dx;

  /*------------------------------------------------------------*/

  /* Temporary variable for input paramter reading */
  tt    = sf_floatalloc2(nz,nx);

  ro    = sf_floatalloc2(nzpad+2,nxpad+2);
  c11   = sf_floatalloc2(nzpad+2,nxpad+2);
  c33   = sf_floatalloc2(nzpad+2,nxpad+2);
  c55   = sf_floatalloc2(nzpad+2,nxpad+2);
  c13   = sf_floatalloc2(nzpad+2,nxpad+2);

  /* Temporary variables */
  t11 = sf_floatalloc2(nzpad+1,nxpad+1); memset(t11[0],0,(nzpad+1)*(nxpad+1)*sizeof(float));
  t12 = sf_floatalloc2(nzpad+1,nxpad+2); memset(t12[0],0,(nzpad+1)*(nxpad+2)*sizeof(float));
  t21 = sf_floatalloc2(nzpad+2,nxpad+1); memset(t21[0],0,(nzpad+2)*(nxpad+1)*sizeof(float));
  t22 = sf_floatalloc2(nzpad+2,nxpad+2); memset(t22[0],0,(nzpad+2)*(nxpad+2)*sizeof(float));

  /* input density */
  sf_floatread(tt[0],nz*nx,Fden);    expand2d(tt,ro,fdm);

  /* input stiffness */
  sf_floatread(tt[0],nz*nx,Fccc);    expand2d(tt,c11,fdm);
  sf_floatread(tt[0],nz*nx,Fccc);    expand2d(tt,c33,fdm);
  sf_floatread(tt[0],nz*nx,Fccc);    expand2d(tt,c55,fdm);
  sf_floatread(tt[0],nz*nx,Fccc);    expand2d(tt,c13,fdm);

  free(*tt); free(tt);


  /* Absorbing boundary conditions */
  /*------------------------------------------------------------*/
//   if(dabc) {
// /* one-way abc setup   */
//     vp = sf_floatalloc2(nzpad+2,nxpad+2);
//     vs = sf_floatalloc2(nzpad+2,nxpad+2);
//     for(ix=0; ix<nxpad+2; ix++) {
//       for(iz=0; iz<nzpad+2; iz++) {
//         vp[ix][iz] = sqrt( c11[ix][iz]/ro[ix][iz] );
//         vs[ix][iz] = sqrt( c13[ix][iz]/ro[ix][iz] );
//       }
//     }
//     abcp = abcone2d_make(NOP,dt,vp,fsrf,fdm);
//     abcs = abcone2d_make(NOP,dt,vs,fsrf,fdm);
//     free(*vp); free(vp);
//     free(*vs); free(vs);
//
//     /* sponge abc setup */
//     spo = sponge_make(fdm->nb);
//     }

  /* Field Variables */

  /* Velocity */
  /* [f,f] grid */
  vx_ff   = sf_floatalloc2(nzpad+2,nxpad+2); memset(vx_ff[0],0,(nzpad+2)*(nxpad+2)*sizeof(float));
  vz_ff   = sf_floatalloc2(nzpad+2,nxpad+2); memset(vz_ff[0],0,(nzpad+2)*(nxpad+2)*sizeof(float));

  /* [v,v] grid */
  vx_vv   = sf_floatalloc2(nzpad+1,nxpad+1); memset(vx_vv[0],0,(nzpad+1)*(nxpad+1)*sizeof(float));
  vz_vv   = sf_floatalloc2(nzpad+1,nxpad+1); memset(vz_vv[0],0,(nzpad+1)*(nxpad+1)*sizeof(float));


  /* Stress */
  /* [v,f] grid */
  sxx_vf  = sf_floatalloc2(nzpad+1,nxpad+2); memset(sxx_vf[0],0,(nzpad+1)*(nxpad+2)*sizeof(float));
  sxz_vf  = sf_floatalloc2(nzpad+1,nxpad+2); memset(sxz_vf[0],0,(nzpad+1)*(nxpad+2)*sizeof(float));
  szz_vf  = sf_floatalloc2(nzpad+1,nxpad+2); memset(szz_vf[0],0,(nzpad+1)*(nxpad+2)*sizeof(float));

  /* [f,v] grid */

  sxx_fv  = sf_floatalloc2(nzpad+2,nxpad+1); memset(sxx_fv[0],0,(nzpad+2)*(nxpad+1)*sizeof(float));
  sxz_fv  = sf_floatalloc2(nzpad+2,nxpad+1); memset(sxz_fv[0],0,(nzpad+2)*(nxpad+1)*sizeof(float));
  szz_fv  = sf_floatalloc2(nzpad+2,nxpad+1); memset(szz_fv[0],0,(nzpad+2)*(nxpad+1)*sizeof(float));
  /* [x,z] output wavefields*/
  wflx    = sf_floatalloc2(nzpad,nxpad); memset(wflx[0],0,nzpad*nxpad*sizeof(float));
  wflz    = sf_floatalloc2(nzpad,nxpad); memset(wflz[0],0,nzpad*nxpad*sizeof(float));

  /*------------------------------------------------------------*/
  /*
   *  MAIN LOOP
   */
  /*------------------------------------------------------------*/
  if(verb) fprintf(stderr,"\n");

  for (it=0; it<nt; it++) {

    if(verb) fprintf(stderr,"\b\b\b\b\b%d",it);

    update_stress_vf(vz_ff, vx_ff, vz_vv, vx_vv, sxx_vf, sxz_vf, szz_vf, ro, c11, c13,   \
                     c33, c55, t12, dtx, dtz, nzpad, nxpad);


    update_stress_fv(vz_ff, vx_ff, vz_vv, vx_vv, sxx_fv, sxz_fv, szz_fv, ro, c11, c13,   \
                     c33, c55, t21, dtx, dtz, nzpad, nxpad);

    // if(ssou) inject_ps_src_cart(szz_fv, szz_vf, ww, cs, it);
    inject_vel_src_cart(vz_ff, vz_vv, ww, cs, it);


    update_velocity_ff(vx_ff, vz_ff, sxx_fv, sxz_fv, sxz_vf, szz_vf, ro, \
                       t22, dtx, dtz, nzpad, nxpad);

    update_velocity_vv(vx_vv, vz_vv, szz_fv, sxz_fv, sxz_vf, sxx_vf, ro, \
                       t11, dtx, dtz, nzpad, nxpad);

    //
    // if(dabc) {
    //
    //   /* sponge ABC */
    //   sponge2d_apply(vx_ff,spo,fdm);
    //   sponge2d_apply(vx_vv,spo,fdm);
    //
    //   sponge2d_apply(vz_ff,spo,fdm);
    //   sponge2d_apply(vz_vv,spo,fdm);
    //     }

    if(snap && it%jsnap==0) {
      /* extract wavefield on fsg */
      extract_wfld_fsg(wflz, vz_ff, vz_vv, nxpad, nzpad);
      cut2dfree(wflz,uc,fdm,acz,acx);

      sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);

      extract_wfld_fsg(wflx,vx_ff, vx_vv, nxpad, nzpad);
      cut2dfree(wflx,uc,fdm,acz,acx);

      sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);

      }

  }

  /* deallocate arrays */
  free(**ww); free(*ww); free(ww);
  free(ss);
  free(rr);
  free(*dd);  free(dd);

  free(*ro);  free(ro);
  free(*c11); free(c11);
  free(*c33); free(c33);
  free(*c55); free(c55);
  free(*c13); free(c13);

  free(*vz_ff); free(vz_ff);
  free(*vx_ff); free(vx_ff);
  free(*vz_vv); free(vz_vv);
  free(*vx_vv); free(vx_vv);

  free(*sxx_vf); free(sxx_vf);
  free(*sxz_vf); free(sxz_vf);
  free(*szz_vf); free(szz_vf);

  free(*sxx_fv); free(sxx_fv);
  free(*sxz_fv); free(sxz_fv);
  free(*szz_fv); free(szz_fv);

  free(*wflx); free(wflx);
  free(*wflz); free(wflz);


}
