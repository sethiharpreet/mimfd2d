#include <rsf.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fdutil.h"
#include "ewefd2d_mfd_utils.h"
#include "mfd_utils_2d.h"


int main(int argc, char* argv[]){

  bool verb ;  /* verbose  */
  bool stress; /* stress source */

  int  jsnap,ntsnap,jdata; /* wavefield parameters*/

  /* I/O files */
  sf_file Fwav  = NULL; /* wavelet   */
  sf_file Fwfl  = NULL; /* wavefield */

  sf_file Fvel1 = NULL; /* velocity 1st layer*/
  sf_file Fden1 = NULL; /* density 1st layer */

  sf_file Fccc  = NULL; /* stiffness 2nd layer */
  sf_file Fden  = NULL; /* density 2nd layer */

  /* Layer 2*/
  sf_axis at,az,ax,ac; /* cube axes */

  /* Layer 1*/
  sf_axis az1; /* z-axis */

  int     nt,nz,nx,nc,nb; // Layer 2
  int     nz1             // Layer 1
  int     nzpad,nxpad;
  float   dt,dz,dx,idz,idx,k,s;

  int it,is1,is2,wfld_j; /* source location */



  /* Initialize RSF */
  sf_init(argc,argv);

  if(! sf_getbool("verb",&verb))    verb=0;  /* verbose */
  if(! sf_getbool("stress",&verb))  stress=1; /* pressure/stress or velocity source */
  if(! sf_getint("jdata",jdata)     jdata=1;


  float **tt  = NULL;
  float **tt1 = NULL;

  /* Laqyer 1 */
  /* Acoustic */
  float **vel1, **rho1;

  /* Layer 2 */
  /* orthorombic footprint - 4 coefficients */
  /* c11 c13
     .   c33
             c55 */

  float **c11 = NULL;
  float **c33 = NULL;
  float **c55 = NULL;
  float **c13 = NULL;

  float **rho;

  float *D,*Q; /* Divergence */
  float *G,*P; /* Gradience */
  float *B;    /* Boundary  */

  /* Wavelet */
  float **ww;

  /* Temporary variables */
  float **t11, **t12, **t21, **t22;


  /* Layer 1*/
  float **vx_ff1, **vz_ff1, **vx_vv1, **vz_vv1, **p_vf1, **p_fv1;


  /* Layer 2*/
  float **vx_ff, **vz_ff, **vx_vv, **vx_ff;
  float **sxx_vf, **sxz_vf, **szz_vf, **sxx_fv, **sxz_fv, **szz_fv;

  /* Wavefield */
  float ***wfld_io;

  // Boundary padding
  nzpad = nz + nb;
  nxpad = nx + 2*nb;


  /* setup I/O files */
  Fwav    = sf_input("in");    /* wavelet */

  /* Medium variables */

  /* Layer 1 */
  Fvel1   = sf_input("vel1");  /* velocity */
  Fden1   = sf_input("rho1");  /* density */

  /* Layer 2 */
  Fccc    = sf_input("ccc");  /* stiffness */
  Fden    = sf_input("rho");  /* density */

  Fwfl= sf_ouput("out"); /* Wavefield */

  /* Read/Write axes */
  at  = sf_iaxa(Fwav,1);   nt  = sf_n(at); dt = sf_d(at);
  az  = sf_iaxa(Fccc,1);   nz  = sf_n(az); dz = sf_d(az);
  ax  = sf_iaxa(Fccc,2);   nx  = sf_n(ax); dx = sf_d(ax);

  az1 = sf_iaxa(Fvel1,1);  nz1 = sf_n(az1);  /* Layer 1*/

  if(verb){
    fprintf('Time nt : %d, and  dt: %f \n',nt,dt);
    fprintf('Axis1 nz: %d and dz: %f  \n',nz,dz);
    fprintf('Axis2 nx: %d and dx: %f : \n',nx,dx);
  }

  /* Wavefield parmaeters */
  if(! sf_getint("jdata",&jdata)) jdata=1;
  if(! sf_getint("jsnap",&jsnap)) jsnap=nt;


  sf_oaxa(Fwfl,az,1);
  sf_oaxa(Fwfl,ax,2);

  sf_setn(at,  ntsnap);
  sf_setd(at,dt*jsnap);

  sf_setn(ac,2);
  sf_seto(ac,0);
  sf_setd(ac,1);
  sf_oaxa(Fwfl,ac,3);

  sf_setn(at,nt/jdata);
  sf_setd(at,dt*jdata);
  sf_oaxa(Fwfl,at,4);

  /* Source Locations */
  if(! sf_getint("is1",is1)) is1=nz1/2;
  if(! sf_getint("is2",is2)) is2=nx/2;

  /* Wavelet */
  ww  = sf_floatalloc2(nt,2); sf_floatread(ww[0],nt*2,Fw);

  tt1 = sf_floatalloc2(nz1,nx);
  tt  = sf_floatalloc2(nz,nx);

  /* Layer 1*/
  rho1 = sf_floatalloc2(nz1+2,nxpad+2);
  vel1 = sf_floatalloc2(nz1+2,nxpad+2);

  /* input density */
  sf_floatread(tt1[0],nz1*nx,Fden1);     expand(tt1,rho1,nz1,nxpad+2);

  /* input density */
  sf_floatread(tt1[0],nz1*nx,Fvel1);     expand(tt1,vel1,nz1,nxpad+2);


  free(*tt1); free(tt1);

  /* Layer 2*/
  rho   = sf_floatalloc2(nzpad+2,nxpad+2);
  c11   = sf_floatalloc2(nzpad+2,nxpad+2);
  c33   = sf_floatalloc2(nzpad+2,nxpad+2);
  c55   = sf_floatalloc2(nzpad+2,nxpad+2);
  c13   = sf_floatalloc2(nzpad+2,nxpad+2);

  /* Temporary variables */
  t11 = sf_floatalloc2(nzpad+1,nxpad+1);
  t12 = sf_floatalloc2(nzpad+1,nxpad+2);
  t21 = sf_floatalloc2(nzpad+2,nxpad+1);
  t22 = sf_floatalloc2(nzpad+2,nxpad+2);

  /* input density */
  sf_floatread(tt[0],nz*nx,Fden);     expand(tt,rho,nzpad,nxpad);

  /* input stiffness */
  sf_floatread(tt[0],nz*nx,Fccc);    expand(tt,c11,nzpad+2,nxpad+2);
  sf_floatread(tt[0],nz*nx,Fccc);    expand(tt,c33,nzpad+2,nxpad+2);
  sf_floatread(tt[0],nz*nx,Fccc);    expand(tt,c55,nzpad+2,nxpad+2);
  sf_floatread(tt[0],nz*nx,Fccc);    expand(tt,c13,nzpad+2,nxpad+2);

  free(*tt); free(tt);


  /* Layer 1*/
  float **vx_ff1, **vz_ff1, **vx_vv1, **vz_vv1;
  float **p_vf1, **p_fv1;


  /* Layer 2*/
  float **vx_ff, **vz_ff, **vx_vv, **vx_ff;
  float **sxx_vf, **sxz_vf, **szz_vf, **sxx_fv, **sxz_fv, **szz_fv;


  /* Field Variables */

  /* Layer 1 */

  /* Velocity */
  /* [f,f] grid */
  vx_ff1  = sf_floatalloc2(nzpad+2,nxpad+2);
  vz_ff1  = sf_floatalloc2(nzpad+2,nxpad+2);

  /* [v,v] grid */
  vx_vv1  = sf_floatalloc2(nzpad+1,nxpad+1);
  vz_vv1  = sf_floatalloc2(nzpad+1,nxpad+1);

  /* Pressure */
  /* [v,f] grid */
  p_vf1   = sf_floatalloc2(nzpad+1,nxpad+2);

  /* [f,v] */
  p_fv1   = sf_floatalloc2(nzpad+2,nxpad+1);


  /* Layer 2 */

  /* Velocity */
  /* [f,f] grid */
  vx_ff   = sf_floatalloc2(nzpad+2,nxpad+2);
  vz_ff   = sf_floatalloc2(nzpad+2,nxpad+2);

  /* [v,v] grid */
  vx_vv   = sf_floatalloc2(nzpad+1,nxpad+1);
  vz_vv   = sf_floatalloc2(nzpad+1,nxpad+1);


  /* Stress */
  /* [v,f] grid */
  sxx_vf  = sf_floatalloc2(nzpad+1,nxpad+2);
  sxz_vf  = sf_floatalloc2(nzpad+1,nxpad+2);
  szz_vf  = sf_floatalloc2(nzpad+1,nxpad+2);

  /* [f,v] grid */

  sxx_fv  = sf_floatalloc2(nzpad+2,nxpad+1);
  sxz_fv  = sf_floatalloc2(nzpad+2,nxpad+1);
  szz_fv  = sf_floatalloc2(nzpad+2,nxpad+1);


  /* I/O wfld */
  wfld_io=sf_floatalloc3(nz,nx,2);


  /* Operator Initialize */
  mfd_op_init(D,G,P,Q,B,k);


  for(it=0; it<nt; it++){

    if(verb && (it-1)%25==0) fprintf(stderr, "Timestep %d of %d \n", it-1,nt-1);

      /* Layer 1 */
      update_pressure_fv(p_fv1,vx_ff1,vz_vv1,t21,rho1,vel1,D,G,nzpad,nxpad,k,s);
      update_pressure_vf(p_vf1,vz_ff1,vx_vv1,t12,rho1,vel1,D,G,nzpad,nxpad,k,s);

      /* Layer 2 */


      update_stress_vf(vx_ff, vz_ff, vx_vv, vz_vv, sxx_vf, sxz_vf, szz_vf, t12, c11, c33, c55, c13, nzpad, nxpad);
      update_stress_fv(vel_ff,vel_vv,sig_fv,t21,c11,c33,c55,c13);


    }

    /* set stress condition */
    update_sm_cart(sig_fv);

    /* inject wavelet in stress components */
    if(stress) inject_stress_src_fsg_cart(sig_fv,sig_vf,vel_ff,vel_vv,ww,is1,is2,it,ccc);


    /* update velocity in field components */
    update_velocity_ff(vel_ff,sig_vf,sig_fv,rr);
    update_velocity_vv(vel_vv,sig_vf,sig_fv,rr);


    /* update velocity in field components */
    update_velocity_ff(vel_ff,sig_vf,sig_fv,rr);
    update_velocity_vv(vel_vv,sig_vf,sig_fv,rr);


    /* update vm */
    update_vm_cart(vel_ff,vel_vv,ll,mm);

    /* inject wavelet in velocity components */
    if(!stress) inject_vel_src_fsg_cart(vel_ff,vel_vv,ww,is1,is2,it);

    /* Extract data */
    if(it%wfld_j==0){
      extract_wfld_fsg_cart(vel_ff,vel_vv,wfld_io);
      sf_floatwrite(wfld_io,nz*nx*2,Fo);
      }


    }

    deallocate(rr,mm,ll,sig_fv,vel_ff,sig_vf,vel_vv,wfld_io,ww);

    return(0);

}
