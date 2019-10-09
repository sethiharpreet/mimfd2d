#include <rsf.h>
#include <math.h>

#include "fdutil.h"
#include "mfd2d_utils.h"
#include "pml.h"

#define NOP 4 /* derivative operator half-size */

int main(int argc, char* argv[]){

  bool verb ;  /* verbose  */
  bool ssou;   /* stress source */
  bool snap;   /* snapshot */
  bool fsrf;   /* free surface */
  bool cpml;   /* Absorbing boundary */

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

  int     nt,nz,nx,ns,nr,nc,npml;
  int     it,iz,ix;
  float   dt,dz,dx,idz,idx,dtx,dtz;
  int     nzpad,nxpad;

  /* PML coeff*/
  float maxvel, fpml;

  /* FDM structure */
  fdm2d    fdm=NULL;

  /* I/O arrays */
  float***ww=NULL;           /* wavelet   */
  pt2d   *ss=NULL;           /* sources   */
  pt2d   *rr=NULL;           /* receivers */
  float **dd=NULL;           /* data      */

  float  **tt  = NULL;

  /* orthorombic footprint - 6 coefficients */
  /* c11 c13 c15
     .   c33 c35
             c55 */

  float **c11 = NULL;
  float **c33 = NULL;
  float **c55 = NULL;
  float **c13 = NULL;
  float **c15 = NULL;
  float **c35 = NULL;

  float **ro = NULL;

  /* Temporary variables */
  float **t11, **t12, **t21, **t22;

  /* C-PML variables */
  float **conv_vz_vv_vf, **conv_vx_vv_vf, **conv_vz_ff_vf, **conv_vx_ff_vf;
  float **conv_vz_ff_fv, **conv_vx_ff_fv, **conv_vz_vv_fv, **conv_vx_vv_fv;
  float **conv_sxx_fv_ff, **conv_szz_vf_ff, **conv_sxz_vf_ff, **conv_sxz_fv_ff;
  float **conv_szz_fv_vv, **conv_sxz_vf_vv, **conv_sxx_vf_vv, **conv_sxz_fv_vv;


  float *K_x, *K_x_half, *b_x, *b_x_half, *a_x, *a_x_half;
  float *K_z, *K_z_half, *b_z, *b_z_half, *a_z, *a_z_half;



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
  if(! sf_getbool("pml",&cpml))  cpml=false; /* absorbing BC */

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
  if( !sf_getint("nb",&npml) || npml<NOP) npml=30;

  fdm  = fdutil_init(verb,fsrf,az,ax,npml,1);

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

  /* PML coeff */
  fpml = 200.0;
  maxvel = 3.0;

  /*------------------------------------------------------------*/

  /* Temporary variable for input paramter reading */
  tt    = sf_floatalloc2(nz,nx);

  ro    = sf_floatalloc2(nzpad+2,nxpad+2);
  c11   = sf_floatalloc2(nzpad+2,nxpad+2);
  c33   = sf_floatalloc2(nzpad+2,nxpad+2);
  c55   = sf_floatalloc2(nzpad+2,nxpad+2);
  c13   = sf_floatalloc2(nzpad+2,nxpad+2);
  c15   = sf_floatalloc2(nzpad+2,nxpad+2);
  c35   = sf_floatalloc2(nzpad+2,nxpad+2);

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
  sf_floatread(tt[0],nz*nx,Fccc);    expand2d(tt,c15,fdm);
  sf_floatread(tt[0],nz*nx,Fccc);    expand2d(tt,c35,fdm);

  free(*tt); free(tt);


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


  /* C-PML */

  /* Coefficients */
  K_x = sf_floatalloc(2*npml);
  K_z = sf_floatalloc(2*npml);
  b_x = sf_floatalloc(2*npml);
  b_z = sf_floatalloc(2*npml);
  a_x = sf_floatalloc(2*npml);
  a_z = sf_floatalloc(2*npml);


  K_x_half = sf_floatalloc(2*npml);
  K_z_half = sf_floatalloc(2*npml);
  b_x_half = sf_floatalloc(2*npml);
  b_z_half = sf_floatalloc(2*npml);
  a_x_half = sf_floatalloc(2*npml);
  a_z_half = sf_floatalloc(2*npml);


  /* FSG convolutional variables */


  /* Stress */
  conv_sxx_fv_ff = sf_floatalloc2(nzpad+2,2*npml); memset(conv_sxx_fv_ff[0],0,(nzpad+2)*(2*npml)*sizeof(float));
  conv_szz_vf_ff = sf_floatalloc2(npml,nxpad+2);   memset(conv_szz_vf_ff[0],0,(npml)*(nxpad+2)*sizeof(float));
  conv_sxz_vf_ff = sf_floatalloc2(npml,nxpad+2);   memset(conv_sxz_vf_ff[0],0,(npml)*(nxpad+2)*sizeof(float));
  conv_sxz_fv_ff = sf_floatalloc2(nzpad+2,2*npml); memset(conv_sxz_fv_ff[0],0,(nzpad+2)*(2*npml)*sizeof(float));

  conv_szz_fv_vv = sf_floatalloc2(npml,nxpad+1);   memset(conv_szz_fv_vv[0],0,(npml)*(nxpad+1)*sizeof(float));
  conv_sxz_vf_vv = sf_floatalloc2(nzpad+1,2*npml); memset(conv_sxz_vf_vv[0],0,(nzpad+1)*(2*npml)*sizeof(float));
  conv_sxx_vf_vv = sf_floatalloc2(nzpad+1,2*npml); memset(conv_sxx_vf_vv[0],0,(nzpad+1)*(2*npml)*sizeof(float));
  conv_sxz_fv_vv = sf_floatalloc2(npml,nxpad+1);   memset(conv_szz_fv_vv[0],0,(npml)*(nxpad+1)*sizeof(float));

  /* Velocities */
  conv_vz_vv_vf = sf_floatalloc2(nzpad+1,2*npml); memset(conv_vz_vv_vf[0],0,(nzpad+1)*(2*npml)*sizeof(float));
  conv_vx_vv_vf = sf_floatalloc2(nzpad+1,2*npml); memset(conv_vx_vv_vf[0],0,(nzpad+1)*(2*npml)*sizeof(float));
  conv_vz_ff_vf = sf_floatalloc2(npml,nxpad+2);   memset(conv_vz_ff_vf[0],0,(npml)*(nxpad+2)*sizeof(float));
  conv_vx_ff_vf = sf_floatalloc2(npml,nxpad+2);   memset(conv_vx_ff_vf[0],0,(npml)*(nxpad+2)*sizeof(float));

  conv_vz_ff_fv = sf_floatalloc2(nzpad+2,2*npml); memset(conv_vz_ff_fv[0],0,(nzpad+2)*(2*npml)*sizeof(float));
  conv_vx_ff_fv = sf_floatalloc2(nzpad+2,2*npml); memset(conv_vx_ff_fv[0],0,(nzpad+2)*(2*npml)*sizeof(float));
  conv_vz_vv_fv = sf_floatalloc2(npml,nxpad+1);   memset(conv_vz_vv_fv[0],0,(npml)*(nxpad+1)*sizeof(float));
  conv_vx_vv_fv = sf_floatalloc2(npml,nxpad+1);   memset(conv_vx_vv_fv[0],0,(npml)*(nxpad+1)*sizeof(float));


  /* PML coeff initilaize */
  pml_coeff(a_x, a_x_half, b_x, b_x_half, K_x, K_x_half, a_z, a_z_half, b_z, b_z_half, K_z, K_z_half, \
            maxvel, fpml, dx, dz, dt, npml, nx, nz);
  //

  /*------------------------------------------------------------*/
  /*
   *  MAIN LOOP
   */
  /*------------------------------------------------------------*/
  if(verb) fprintf(stderr,"\n");

  for (it=0; it<nt; it++) {

    if(verb) fprintf(stderr,"\b\b\b\b\b%d",it);


    update_stress_vf_tti_PML(vz_ff, vx_ff, vz_vv, vx_vv, sxx_vf, sxz_vf, szz_vf, ro, c11, c13, c33, c55, c15, c35, \
                             t12, dtx, dtz, nzpad, nxpad, conv_vz_vv_vf, conv_vx_vv_vf, conv_vz_ff_vf, conv_vx_ff_vf, \
                             K_z, b_z, a_z, K_x_half, b_x_half, a_x_half, npml, fsrf);

    update_stress_fv_tti_PML(vz_ff, vx_ff, vz_vv, vx_vv, sxx_fv, sxz_fv, szz_fv, ro, c11, c13, c33, c55, c15, c35,\
                             t21, dtx, dtz, nzpad, nxpad, conv_vz_ff_fv, conv_vx_ff_fv, conv_vz_vv_fv, conv_vx_vv_fv, \
                             K_x, b_x, a_x, K_z_half, b_z_half, a_z_half, npml, fsrf);


    // if(ssou) inject_ps_src_cart(szz_fv, szz_vf, ww, cs, it);
    inject_vel_src_cart(vz_ff, vz_vv, ww, cs, it);


    update_velocity_ff_PML(vx_ff, vz_ff, sxx_fv, sxz_fv, sxz_vf, szz_vf, ro, t22, dtx, dtz, nzpad, nxpad, \
                           conv_sxx_fv_ff, conv_szz_vf_ff, conv_sxz_vf_ff, conv_sxz_fv_ff, K_x_half, b_x_half, a_x_half,\
                           K_z_half, b_z_half, a_z_half, npml, fsrf);


    update_velocity_vv_PML(vx_vv, vz_vv, szz_fv, sxz_fv, sxz_vf, sxx_vf, ro, t11, dtx, dtz, nzpad, nxpad, \
                           conv_szz_fv_vv, conv_sxz_vf_vv, conv_sxx_vf_vv, conv_sxz_fv_vv, K_x, b_x, a_x, \
                           K_z, b_z, a_z, npml, fsrf);



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
  free(*c15); free(c15);
  free(*c35); free(c35);

  free(*t11); free(t11);
  free(*t12); free(t12);
  free(*t21); free(t21);
  free(*t22); free(t22);

  /* PML */

  free(*conv_vz_vv_vf); free(conv_vz_vv_vf);
  free(*conv_vx_vv_vf); free(conv_vx_vv_vf);
  free(*conv_vz_ff_vf); free(conv_vz_ff_vf);
  free(*conv_vx_ff_vf); free(conv_vx_ff_vf);

  free(*conv_vz_ff_fv); free(conv_vz_ff_fv);
  free(*conv_vx_ff_fv); free(conv_vx_ff_fv);
  free(*conv_vz_vv_fv); free(conv_vz_vv_fv);
  free(*conv_vx_vv_fv); free(conv_vx_vv_fv);

  free(*conv_sxx_fv_ff); free(conv_sxx_fv_ff);
  free(*conv_sxz_fv_ff); free(conv_sxz_fv_ff);
  free(*conv_szz_vf_ff); free(conv_szz_vf_ff);
  free(*conv_sxz_vf_ff); free(conv_sxz_vf_ff);

  free(*conv_szz_fv_vv); free(conv_szz_fv_vv);
  free(*conv_sxz_vf_vv); free(conv_sxz_vf_vv);
  free(*conv_sxx_vf_vv); free(conv_sxx_vf_vv);
  free(*conv_sxz_fv_vv); free(conv_sxz_fv_vv);

  free(K_x); free(K_x_half); free(K_z); free(K_z_half);
  free(a_x); free(a_x_half); free(a_z); free(a_z_half);
  free(b_x); free(b_x_half); free(b_z); free(b_z_half);

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
