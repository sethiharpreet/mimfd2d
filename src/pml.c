/* Convolutional Perfectly Matched Layer */
/* Reference: An unsplit convolutional perfectly matched layer improved
at grazing incidence for the seismic wave equation : Komatitsch and Martin 2007 */

#include <rsf.h>
#include <stdio.h>
#include <stdlib.h>
#include "fdutil.h"
#include "mfd2d_utils.h"



void pml_coeff(float *a_x, float *a_x_half, float *b_x, float *b_x_half, float *K_x, float *K_x_half, \
               float *a_z, float *a_z_half, float *b_z, float *b_z_half, float *K_z, float *K_z_half, \
               float quasi_cp_max, float fpml, float dx, float dz, float dt, int npml, int nx, int nz)
{

  const float npower = 2.0; /*  power to compute d0 profile */
  const float k_max_PML = 1.0; /* from Gedney page 8.11 */
  const float PI = 3.141592653589793;
  const float alpha_max_PML = 2.0 * PI * (fpml/2.0); /* from festa and Vilotte */

  float thickness_PML_x, thickness_PML_z, xoriginleft, xoriginright, zoriginbottom, zorigintop;
  float Rcoef , d0_x, d0_z, xval, zval, abscissa, abscissa_norm;

  float d_x, d_z, alpha_prime_x, alpha_prime_z;
  float d_x_half, d_z_half, alpha_prime_x_half, alpha_prime_z_half;


  /* thickness of the PML layer in meters */
  thickness_PML_x = npml*dx;
  thickness_PML_z = npml*dz;

  /* reflection coefficient (INRIA report section 6.1) */
  Rcoef = 0.00000001;
  //Rcoef = 0.001;

	/* local variables */
	int i, h;

  /* define profile of absorption in PML region */

  /* compute d0 from INRIA report section 6.1 */
  d0_x = - (npower + 1) * quasi_cp_max * log(Rcoef) / (2.0 * thickness_PML_x);
  d0_z = - (npower + 1) * quasi_cp_max * log(Rcoef) / (2.0 * thickness_PML_z);

  /* damping in the X direction */
  /* -------------------------- */

  /* origin of the PML layer (position of right edge minus thickness, in meters) */
  xoriginleft = thickness_PML_x;
  xoriginright = (nx-1)*dx - thickness_PML_x;



  /******** left boundary *********/

  for (i=0;i<npml;i++){

    a_x[i]      = 0.0;
    a_x_half[i] = 0.0;
    b_x[i]      = 0.0;
    b_x_half[i] = 0.0;
    K_x[i]      = 1.0;
    K_x_half[i] = 1.0;


    xval = dx * (float)i;


    /* define damping profile at the grid points */
    abscissa = xoriginleft - xval;

    if(abscissa >= 0.0){

      abscissa_norm = abscissa/thickness_PML_x;
      d_x = d0_x * pow(abscissa_norm,npower);

      /* this taken from Gedney page 8.2 */
      K_x[i] = 1.0 + (k_max_PML - 1.0) * pow(abscissa_norm,npower);
      alpha_prime_x = alpha_max_PML * (1.0 - abscissa_norm);

      }



    /* define damping profile at half the grid points */
    abscissa = xoriginleft - (xval + dx/2.0);

    if(abscissa >= 0.0){

      abscissa_norm = abscissa/thickness_PML_x;
      d_x_half = d0_x * pow(abscissa_norm,npower);

      /* this taken from Gedney page 8.2 */
      K_x_half[i] = 1.0 + (k_max_PML - 1.0) * pow(abscissa_norm,npower);
      alpha_prime_x_half = alpha_max_PML * (1.0 - abscissa_norm);

      }

     /* just in case, for -5 at the end */
     if(alpha_prime_x < 0.0){ alpha_prime_x = 0.0;}
     if(alpha_prime_x_half < 0.0) { alpha_prime_x_half = 0.0;}

     b_x[i]      = exp(- (d_x / K_x[i] + alpha_prime_x) * dt);
     b_x_half[i] = exp(- (d_x_half / K_x_half[i] + alpha_prime_x_half) * dt);

     /* avoid division by zero outside the PML */
     if(fabs(d_x) > 1.0e-6){ a_x[i] = d_x * (b_x[i] - 1.0) / (K_x[i] * (d_x + K_x[i] * alpha_prime_x));}
     if(fabs(d_x_half) > 1.0e-6){ a_x_half[i] = d_x_half * (b_x_half[i] - 1.0) / (K_x_half[i] * (d_x_half + K_x_half[i] * alpha_prime_x_half));}

      }
      /* end of left boundary */



    /******** right boundary *********/

  for (i=nx-npml;i<nx;i++){

    h=i-nx+2*npml;
    a_x[h]      = 0.0;
    a_x_half[h] = 0.0;
    b_x[h]      = 0.0;
    b_x_half[h] = 0.0;
    K_x[h]      = 1.0;
	  K_x_half[h] = 1.0;


	  xval = dx*(float)i;



    /* define damping profile at the grid points */
    abscissa = xval - xoriginright;

    if(abscissa >= 0.0){

      abscissa_norm = abscissa / thickness_PML_x;
      d_x = d0_x * pow(abscissa_norm,npower);

      /* this taken from Gedney page 8.2 */
      K_x[h] = 1.0 + (k_max_PML - 1.0) * pow(abscissa_norm,npower);
      alpha_prime_x = alpha_max_PML * (1.0 - abscissa_norm);

      }

    /* define damping profile at half the grid points */
    abscissa = xval + dx/2.0 - xoriginright;

    if(abscissa >= 0.0){

      abscissa_norm = abscissa / thickness_PML_x;
      d_x_half = d0_x * pow(abscissa_norm ,npower);

      /* this taken from Gedney page 8.2 */
      K_x_half[h] = 1.0 + (k_max_PML - 1.0) * pow(abscissa_norm,npower);
      alpha_prime_x_half = alpha_max_PML * (1.0 - abscissa_norm);

      }

      /* just in case, for -5 at the end */
      if(alpha_prime_x < 0.0){ alpha_prime_x = 0.0;}
      if(alpha_prime_x_half < 0.0) {alpha_prime_x_half = 0.0;}

      b_x[h]      = exp(- (d_x / K_x[h] + alpha_prime_x) * dt);
      b_x_half[h] = exp(- (d_x_half / K_x_half[h] + alpha_prime_x_half) * dt);

      /* avoid division by zero outside the PML */
      if(fabs(d_x) > 1.0e-6){ a_x[h] = d_x * (b_x[h] - 1.0) / (K_x[h] * (d_x + K_x[h] * alpha_prime_x));}
      if(fabs(d_x_half) > 1.0e-6){ a_x_half[h] = d_x_half * (b_x_half[h] - 1.0) / (K_x_half[h] * (d_x_half + K_x_half[h] * alpha_prime_x_half));}

       }
       /* end of right boundary */


  /* damping in the Z direction */
  /* -------------------------- */

  /* origin of the PML layer (position of right edge minus thickness, in meters) */
  zorigintop = thickness_PML_z;
  zoriginbottom = (nz-1)*dz - thickness_PML_z;

  for(i=0;i<npml;i++){

    a_z[i]      = 0.0;
    a_z_half[i] = 0.0;
    b_z[i]      = 0.0;
    b_z_half[i] = 0.0;
    K_z[i]      = 1.0;
    K_z_half[i] = 1.0;

    zval = dz*i;

  /******** top boundary *********/

    /* define damping profile at the grid points */
    abscissa = zorigintop - zval;

    if(abscissa >= 0.0){

      abscissa_norm = abscissa / thickness_PML_z;
      d_z = d0_z * pow(abscissa_norm,npower);

      /* this taken from Gedney page 8.2 */
      K_z[i] = 1.0 + (k_max_PML - 1.0) * pow(abscissa_norm,npower);
      alpha_prime_z = alpha_max_PML * (1.0 - abscissa_norm);

      }

    /* define damping profile at half the grid points */
    abscissa = zorigintop - (zval + dz/2.0);

    if(abscissa >= 0.0){

      abscissa_norm = abscissa / thickness_PML_z;
      d_z_half = d0_z * pow(abscissa_norm,npower);

      /* this taken from Gedney page 8.2 */
      K_z_half[i] = 1.0 + (k_max_PML - 1.0) * pow(abscissa_norm,npower);
      alpha_prime_z_half = alpha_max_PML * (1.0 - abscissa_norm);

      }

    b_z[i]      = exp(- (d_z / K_z[i] + alpha_prime_z) * dt);
    b_z_half[i] = exp(- (d_z_half / K_z_half[i] + alpha_prime_z_half) * dt);

    /* avoid division by zero outside the PML */
    if(fabs(d_z) > 1.0e-6){ a_z[i] = d_z * (b_z[i] - 1.0) / (K_z[i] * (d_z + K_z[i] * alpha_prime_z));}
    if(fabs(d_z_half) > 1.0e-6){ a_z_half[i] = d_z_half * (b_z_half[i] - 1.0) / (K_z_half[i] * (d_z_half + K_z_half[i] * alpha_prime_z_half));}

  } /* end of top boundary */

/******** bottom boundary *********/
  for (i=nz-npml;i<nz;i++){

    h=i-nz+2*npml;

    a_z[h]      = 0.0;
    a_z_half[h] = 0.0;
    b_z[h]      = 0.0;
    b_z_half[h] = 0.0;
    K_z[h]      = 1.0;
    K_z_half[h] = 1.0;
    zval = dz * i;

    /* define damping profile at the grid points */
    abscissa = zval - zoriginbottom;

    if(abscissa >= 0.0){

      abscissa_norm = abscissa / thickness_PML_z;
      d_z = d0_z * pow(abscissa_norm,npower);

      /* this taken from Gedney page 8.2 */
      K_z[h] = 1.0 + (k_max_PML - 1.0) * pow(abscissa_norm,npower);
      alpha_prime_z = alpha_max_PML * (1.0 - abscissa_norm);

        }

    /* define damping profile at half the grid points */
    abscissa = zval + dz/2.0 - zoriginbottom;

    if(abscissa >= 0.0){

        abscissa_norm = abscissa / thickness_PML_z;
        d_z_half = d0_z * pow(abscissa_norm,npower);

        /* this taken from Gedney page 8.2 */
        K_z_half[h] = 1.0 + (k_max_PML - 1.0) * pow(abscissa_norm,npower);
        alpha_prime_z_half = alpha_max_PML * (1.0 - abscissa_norm);

          }

     b_z[h] = exp(- (d_z / K_z[h] + alpha_prime_z) * dt);
     b_z_half[h] = exp(- (d_z_half / K_z_half[h] + alpha_prime_z_half) * dt);


     /* avoid division by zero outside the PML */
     if(fabs(d_z) > 1.0e-6){ a_z[h] = d_z * (b_z[h] - 1.0) / (K_z[h] * (d_z + K_z[h] * alpha_prime_z));}
     if(fabs(d_z_half) > 1.0e-6){ a_z_half[h] = d_z_half * (b_z_half[h] - 1.0) / (K_z_half[h] * (d_z_half + K_z_half[h] * alpha_prime_z_half));}

   } /* end of bottom boundary */

}

/* PML FSG update */

/* Acoustic */
void update_pressure_fv_PML(float **p_fv, float **vx_ff, float **vz_vv, float **t21, float **ro, float **vp, float dtx, float dtz, int nzpad, int nxpad,\
                            float **conv_vzvv_fv, float **conv_vxff_fv, float *K_x, float *b_x, float *a_x, float *K_z_half, float *b_z_half, float *a_z_half, \
                            int npml, bool fsrf, bool cpld)
/* vx[f,f] -------> p[f,v]   */
/* vz[v,v] -------> p[f,v]   */
{
  int iz,ix,ib;

  memset(t21[0],0,(nzpad+2)*(nxpad+1)*sizeof(float));

  D1_CC(t21,vz_vv,nzpad+1,nxpad+1);

  for (ix=0; ix<nxpad+1; ix++){
    for(iz=0; iz<nzpad+2; iz++){

      /* top boundary */
      if(!(fsrf) && (iz < npml)){
        conv_vzvv_fv[ix][iz] = b_z_half[iz] * conv_vzvv_fv[ix][iz] + a_z_half[iz] * t21[ix][iz];
        t21[ix][iz] = t21[ix][iz]/K_z_half[iz] + conv_vzvv_fv[ix][iz];
      }

      /* bottom boundary */
      if(!(cpld) && (iz > nzpad+2-npml)){
        ib = (iz - nzpad - 2 + 2*npml);
        conv_vzvv_fv[ix][ib] = b_z_half[ib] * conv_vzvv_fv[ix][ib] + a_z_half[ib] * t21[ix][iz];
        t21[ix][iz] = t21[ix][iz]/K_z_half[ib] + conv_vzvv_fv[ix][ib];
      }

      /* Update pressure */
      p_fv[ix][iz]-= dtz*ro[ix][iz]*vp[ix][iz]*vp[ix][iz]*t21[ix][iz];

    }
  }

  memset(t21[0],0,(nzpad+2)*(nxpad+1)*sizeof(float));

  G2_CC(t21,vx_ff,nzpad+2,nxpad+1);

  for (ix=0; ix<nxpad+1; ix++){
    for(iz=0; iz<nzpad+2; iz++){

      /* left boundary */
      if(ix < npml){
        conv_vxff_fv[ix][iz] = b_x[ix] * conv_vxff_fv[ix][iz]  + a_x[ix] * t21[ix][iz];
        t21[ix][iz] = t21[ix][iz] / K_x[ix] + conv_vxff_fv[ix][iz];
      }

      /* right boundary */
      if(ix > nxpad+1-npml){
        ib = (ix - nxpad -1 + 2*npml);
        conv_vxff_fv[ib][iz] = b_x[ib] * conv_vxff_fv[ib][iz]  + a_x[ib] * t21[ix][iz];
        t21[ix][iz] = t21[ix][iz] / K_x[ib] + conv_vxff_fv[ib][iz];
      }

      /* Update pressure */
      p_fv[ix][iz]-= dtx*ro[ix][iz]*vp[ix][iz]*vp[ix][iz]*t21[ix][iz];

    }
  }

}

void update_pressure_vf_PML(float **p_vf, float **vz_ff, float **vx_vv, float **t12, float **ro, float **vp, float dtx, float dtz, int nzpad, int nxpad,\
                            float **conv_vxvv_vf, float **conv_vzff_vf, float *K_z, float *b_z, float *a_z, float *K_x_half, float *b_x_half, float *a_x_half, \
                            int npml, bool fsrf, bool cpld)
/* vx[f,f] -------> p[f,v]   */
/* vz[v,v] -------> p[f,v]   */
{
  int iz,ix,ib;

  memset(t12[0],0,(nzpad+1)*(nxpad+2)*sizeof(float));

  G1_CC(t12,vz_ff,nzpad+1,nxpad+2);

  for (ix=0; ix<nxpad+2; ix++){
    for(iz=0; iz<nzpad+1; iz++){

      /* top boundary */
      if(!(fsrf) && (iz < npml)){
        conv_vzff_vf[ix][iz] = b_z[iz] * conv_vzff_vf[ix][iz] + a_z[iz] * t12[ix][iz];
        t12[ix][iz] = t12[ix][iz]/K_z[iz] + conv_vzff_vf[ix][iz];
      }

      /* bottom boundary */
      if(!(cpld) && (iz > nzpad+1-npml)){
        ib = (iz - nzpad - 1 + 2*npml);
        conv_vzff_vf[ix][ib] = b_z[ib] * conv_vzff_vf[ix][ib] + a_z[ib] * t12[ix][iz];
        t12[ix][iz] = t12[ix][iz]/K_z[ib] + conv_vzff_vf[ix][ib];
      }

      p_vf[ix][iz]-= dtz*ro[ix][iz]*vp[ix][iz]*vp[ix][iz]*t12[ix][iz];
    }
  }

  memset(t12[0],0,(nzpad+1)*(nxpad+2)*sizeof(float));


  D2_CC(t12,vx_vv,nzpad+1,nxpad+1);

  for (ix=0; ix<nxpad+2; ix++){
    for(iz=0; iz<nzpad+1; iz++){

      /* left boundary */
      if(ix < npml){
        conv_vxvv_vf[ix][iz] = b_x_half[ix] * conv_vxvv_vf[ix][iz]  + a_x_half[ix] * t12[ix][iz];
        t12[ix][iz] = t12[ix][iz] / K_x_half[ix] + conv_vxvv_vf[ix][iz];
      }

      /* right boundary */
      if(ix > nxpad+2-npml){
        ib = (ix - nxpad - 2 + 2*npml);
        conv_vxvv_vf[ib][iz] = b_x_half[ib] * conv_vxvv_vf[ib][iz]  + a_x_half[ib] * t12[ix][iz];
        t12[ix][iz] = t12[ix][iz] / K_x_half[ib] + conv_vxvv_vf[ib][iz];
      }

      p_vf[ix][iz]-= dtx*ro[ix][iz]*vp[ix][iz]*vp[ix][iz]*t12[ix][iz];
    }
  }

}

/* Acoustic velocity Update */

void update_ac_velocity_ff_PML(float **vz_ff, float **vx_ff, float ** p_vf, float **p_fv, float **t22, float **ro, float dtx, float dtz, int nzpad, int nxpad,\
                               float **conv_pvf_ff, float **conv_pfv_ff, float *K_x_half, float *K_z_half, float *b_x_half, float *b_z_half, float *a_x_half, float *a_z_half, \
                               int npml, bool fsrf, bool cpld)
/* p[f,v] -------> vx[f,f]   */
/* p[v,f] -------> vz[f,f]   */
{

  int iz,ix,ib;

  memset(t22[0],0,(nzpad+2)*(nxpad+2)*sizeof(float));

  D1_CC(t22,p_vf,nzpad+1,nxpad+2);

  for (ix=0; ix<nxpad+2; ix++){
    for(iz=0; iz<nzpad+2; iz++){

      /* top boundary */
      if(!(fsrf) && (iz < npml)){
        conv_pvf_ff[ix][iz] = b_z_half[iz] * conv_pvf_ff[ix][iz] + a_z_half[iz] * t22[ix][iz];
        t22[ix][iz] = t22[ix][iz]/K_z_half[iz] + conv_pvf_ff[ix][iz];
      }

      /* bottom boundary */
      if(!(cpld) && (iz > nzpad+2-npml)){
        ib = (iz - nzpad - 2 + 2*npml);
        conv_pvf_ff[ix][ib] = b_z_half[ib] * conv_pvf_ff[ix][ib] + a_z_half[ib] * t22[ix][iz];
        t22[ix][iz] = t22[ix][iz]/K_z_half[ib] + conv_pvf_ff[ix][ib];
      }

      vz_ff[ix][iz]-= dtz*t22[ix][iz]/ro[ix][iz];
    }
  }

  memset(t22[0],0,(nzpad+2)*(nxpad+2)*sizeof(float));

  D2_CC(t22,p_fv,nzpad+2,nxpad+1);

  for (ix=0; ix<nxpad+2; ix++){
    for(iz=0; iz<nzpad+2; iz++){

      /* left boundary */
      if(ix < npml){
        conv_pfv_ff[ix][iz] = b_x_half[ix] * conv_pfv_ff[ix][iz]  + a_x_half[ix] * t22[ix][iz];
        t22[ix][iz] = t22[ix][iz] / K_x_half[ix] + conv_pfv_ff[ix][iz];
      }

      /* right boundary */
      if(ix > nxpad+2-npml){
        ib = (ix - nxpad - 2 + 2*npml);
        conv_pfv_ff[ib][iz] = b_x_half[ib] * conv_pfv_ff[ib][iz]  + a_x_half[ib] * t22[ix][iz];
        t22[ix][iz] = t22[ix][iz] / K_x_half[ib] + conv_pfv_ff[ib][iz];
      }

      vx_ff[ix][iz]-= dtx*t22[ix][iz]/ro[ix][iz];
    }
  }

}



void update_ac_velocity_vv_PML(float **vz_vv, float **vx_vv, float **p_vf, float **p_fv, float **t11, float **ro, float dtx, float dtz, int nzpad, int nxpad, \
                               float **conv_pfv_vv, float **conv_pvf_vv, float *K_x, float *K_z, float *b_x, float *b_z, float *a_x, float *a_z, \
                               float *K_x_half, float *K_z_half, float *b_x_half, float *b_z_half, float *a_x_half, float *a_z_half, int npml, bool fsrf, bool cpld)
/* p[f,v] -------> vz[v,v]   */
/* p[v,f] -------> vx[v,v]   */
{

  int iz,ix,ib;

  memset(t11[0],0,(nzpad+1)*(nxpad+1)*sizeof(float));

  G1_CC(t11,p_fv,nzpad+1,nxpad+1);

  for (ix=0; ix<nxpad+1; ix++){
    for(iz=0; iz<nzpad+1; iz++){

      /* top boundary */
      if(!(fsrf) && (iz < npml)){
        conv_pfv_vv[ix][iz] = b_z[iz] * conv_pfv_vv[ix][iz] + a_z[iz] * t11[ix][iz];
        t11[ix][iz] = t11[ix][iz]/K_z[iz] + conv_pfv_vv[ix][iz];
      }

      /* bottom boundary */
      if(!(cpld) && (iz > nzpad+1-npml)){
        ib = (iz - nzpad - 1 + 2*npml);
        conv_pfv_vv[ix][ib] = b_z[ib] * conv_pfv_vv[ix][ib] + a_z[ib] * t11[ix][iz];
        t11[ix][iz] = t11[ix][iz]/K_z[ib] + conv_pfv_vv[ix][ib];
      }

      vz_vv[ix][iz]-= dtz*t11[ix][iz]/ro[ix][iz];
    }
  }

  memset(t11[0],0,(nzpad+1)*(nxpad+1)*sizeof(float));

  G2_CC(t11,p_vf,nzpad+1,nxpad+1);

  for (ix=0; ix<nxpad+1; ix++){
    for(iz=0; iz<nzpad+1; iz++){

      /* left boundary */
      if(ix < npml){
        conv_pvf_vv[ix][iz] = b_x[ix] * conv_pvf_vv[ix][iz]  + a_x[ix] * t11[ix][iz];
        t11[ix][iz] = t11[ix][iz] / K_x[ix] + conv_pvf_vv[ix][iz];
      }

      /* right boundary */
      if(ix > nxpad+1-npml){
        ib = (ix - nxpad - 1 + 2*npml);
        conv_pvf_vv[ib][iz] = b_x[ib] * conv_pvf_vv[ib][iz]  + a_x[ib] * t11[ix][iz];
        t11[ix][iz] = t11[ix][iz] / K_x[ib] + conv_pvf_vv[ib][iz];
      }

      vx_vv[ix][iz]-= dtx*t11[ix][iz]/ro[ix][iz];

    }
  }
}



/* PML FSG  TTI update */



/***************************************************************
/ Elastic (Anisotropic) orthorombic footprint - 6 coefficients
/ c11 c13   c15
/   .   c33 c35
/           c55
/  TTI
***************************************************************/


/****************************************
/ UPDATE STRESSES FOR [v,f] GRIDS - CART : TTI
*****************************************/


void update_stress_vf_tti_PML(float **vz_ff,  float **vx_ff,  float **vz_vv,  float **vx_vv, float **sxx_vf, float **sxz_vf, float **szz_vf, float **rho,   \
                              float **c11,    float **c13,    float **c33,    float **c55,   float **c15,    float **c35,    float **t12,    float dtx, float dtz,   \
                              int nzpad, int nxpad,  float **conv_vz_vv, float **conv_vx_vv, float **conv_vz_ff, float **conv_vx_ff, float *K_x, float *K_z, float *b_x, \
                              float *b_z, float *a_x, float *a_z, float *K_x_half, float *K_z_half, float *b_x_half, float *b_z_half, float *a_x_half, float *a_z_half, int npml)
/* */
/* */
{

  int iz,ix;

  memset(t12[0],0,(nzpad+1)*(nxpad+2)*sizeof(float));

  // compute DV1DE2 and update components
  // divergence operator

  D2_CC(t12,vz_vv,nzpad+1,nxpad+1);

  for(ix=0; ix<nxpad+2; ix++){
    for(iz=0; iz<nzpad+1; iz++){


      //S11 component
      // s_zz = c35*dvzdx
      szz_vf[ix][iz]+= dtx*c35[ix][iz]*t12[ix][iz];

      //S22 component
      // s_xx = c15*dvzdx
      sxx_vf[ix][iz]+= dtx*c15[ix][iz]*t12[ix][iz];

      // S12 or S21 component
      // s_zx = c55*dvzdx
      sxz_vf[ix][iz]+= dtx*c55[ix][iz]*t12[ix][iz];

      }
    }

  memset(t12[0],0,(nzpad+1)*(nxpad+2)*sizeof(float));

  // compute DV2DE2 and update components

  // Divergence
  D2_CC(t12,vx_vv,nzpad+1,nxpad+1);

  for(ix=0; ix<nxpad+2; ix++){
    for(iz=0; iz<nzpad+1; iz++){

      // S11 component
      // s_zz=c13*dvxdx
      szz_vf[ix][iz]+= dtx*c13[ix][iz]*t12[ix][iz];

      // S22 component
      // s_xx=c11*dvxdx
      sxx_vf[ix][iz]+= dtx*c11[ix][iz]*t12[ix][iz];

      //S12 or S21 components
      // s_xz = c15*dvxdx
      sxz_vf[ix][iz]+= dtx*c15[ix][iz]*t12[ix][iz];

        }
      }

  memset(t12[0],0,(nzpad+1)*(nxpad+2)*sizeof(float));


  // compute DV1DE1 and update components
  G1_CC(t12,vz_ff,nzpad+1,nxpad+2);

  for(ix=0; ix<nxpad+2; ix++){
    for(iz=0; iz<nzpad+1; iz++){

      // S11 component
      // s_zz=c33*dvzdz
      szz_vf[ix][iz]+= dtz*c33[ix][iz]*t12[ix][iz];

      // S22 component
      // s_xx=c13*dvzdz
      sxx_vf[ix][iz]+= dtz*c13[ix][iz]*t12[ix][iz];

      // S12 or S21 component
      // s_xz=c35*dvzdz
      sxz_vf[ix][iz]+= dtz*c35[ix][iz]*t12[ix][iz];
    }
  }

  // compute DV2DE1 and update components

  memset(t12[0],0,(nzpad+1)*(nxpad+2)*sizeof(float));

  G1_CC(t12,vx_ff,nzpad+1,nxpad+2);

  for(ix=0; ix<nxpad+2; ix++){
    for(iz=0; iz<nzpad+1; iz++){

      // S11 component
      // s_zz = c35*dvxdz
      szz_vf[ix][iz]+= dtz*c35[ix][iz]*t12[ix][iz];

      // S22 component
      // s_xx = c15*dvxdz
      sxx_vf[ix][iz]+= dtz*c15[ix][iz]*t12[ix][iz];

      // S12 or S21 component
      // s_zx=c55*dvxdz (dv2de1 part)
      sxz_vf[ix][iz]+= dtz*c55[ix][iz]*t12[ix][iz];

      }
    }

  }


  /****************************************
  / UPDATE STRESSES FOR [f,v] GRIDS - CART : TTI
  *****************************************/

  void update_stress_fv_tti_PML(float **vz_ff,  float **vx_ff,  float **vz_vv,  float **vx_vv,          \
                            float **sxx_fv, float **sxz_fv, float **szz_fv, float **rho,            \
                            float **c11,    float **c13,    float **c33,    float **c55,            \
                            float **c15,    float **c35,    float **t21,    float dtx,  float dtz,  \
                            int nzpad,      int nxpad)
  {
      int iz,ix;

      memset(t21[0],0,(nzpad+2)*(nxpad+1)*sizeof(float));

      //compute DV1DE2 and update components

      //Gradient
      G2_CC(t21,vz_ff,nzpad+2,nxpad+1);

      for(ix=0; ix<nxpad+1; ix++){
        for(iz=0; iz<nzpad+2; iz++){

          // S11 component
          // s_zz = c35*dvzdx
          szz_fv[ix][iz]+= dtx*c35[ix][iz]*t21[ix][iz];

          // S22 component
          // s_xx = c15*dvzdx
          sxx_fv[ix][iz]+= dtx*c15[ix][iz]*t21[ix][iz];

          // S12 or S21 component
          // s_zx=c55*dvzdx (dv1de2 part)
          sxz_fv[ix][iz]+= dtx*c55[ix][iz]*t21[ix][iz];

          }
        }


      memset(t21[0],0,(nzpad+2)*(nxpad+1)*sizeof(float));

      // compute DV2DE2 and update components

      // Gradient
      G2_CC(t21,vx_ff,nzpad+2,nxpad+1);

      for(ix=0; ix<nxpad+1; ix++){
        for(iz=0; iz<nzpad+2; iz++){

          // S11 component
          // s_zz=c13*dvxdx
          szz_fv[ix][iz]+= dtx*c13[ix][iz]*t21[ix][iz];

          // S22 component
          // s_xx=c11*dvxdx
          sxx_fv[ix][iz]+= dtx*c11[ix][iz]*t21[ix][iz];

          // S12 or S21 component
          // s_xz=c15*dvxdx
          sxz_fv[ix][iz]+= dtx*c15[ix][iz]*t21[ix][iz];

            }
          }

      memset(t21[0],0,(nzpad+2)*(nxpad+1)*sizeof(float));

      // compute DV1DE1 and update components
      D1_CC(t21,vz_vv,nzpad+1,nxpad+1);

      for(ix=0; ix<nxpad+1; ix++){
        for(iz=0; iz<nzpad+2; iz++){

          // S11 component
          // s_zz=c33*dvzdz
          szz_fv[ix][iz]+= dtz*c33[ix][iz]*t21[ix][iz];

          // S22 component
          // s_xx=c13*dvzdz
          sxx_fv[ix][iz]+= dtz*c13[ix][iz]*t21[ix][iz];

          // S12 or S21 component
          // sxz = c35*dvzdz
          sxz_fv[ix][iz]+= dtz*c35[ix][iz]*t21[ix][iz];

        }
      }


      memset(t21[0],0,(nzpad+2)*(nxpad+1)*sizeof(float));
      // compute  DV2DE1 and update components

      D1_CC(t21,vx_vv,nzpad+1,nxpad+1);

      for(ix=0; ix<nxpad+1; ix++){
        for(iz=0; iz<nzpad+2; iz++){

          // S11 component
          // s_zz = c35*dvxdz
          szz_fv[ix][iz]+= dtz*c35[ix][iz]*t21[ix][iz];

          // S22 component
          // s_xx=c15*dvxdz
          sxx_fv[ix][iz]+= dtz*c15[ix][iz]*t21[ix][iz];

          // S12 or S21 component
          // s_zx=c55*dvxdz (dv2de1 part)
          sxz_fv[ix][iz]+= dtz*c55[ix][iz]*t21[ix][iz];

          }
        }

    }


/****************************************
/ UPDATE VELOCITY FOR [f,f] GRIDS - CART
***************************************/

void update_velocity_ff_PML(float **vx_ff, float ** vz_ff, float **sxx_fv, float **sxz_fv, float **sxz_vf, \
                      float **szz_vf,float **rho, float **t22, float dtx, float dtz, int nzpad, int nxpad)
{
  int iz,ix;

  memset(t22[0],0,(nzpad+2)*(nxpad+2)*sizeof(float));

  // compute ds22de2 and update components

  D2_CC(t22,sxx_fv,nzpad+2,nxpad+1);

  for(ix=0; ix<nxpad+2; ix++){
    for(iz=0; iz<nzpad+2; iz++){
      vx_ff[ix][iz]+= dtx*t22[ix][iz]/rho[ix][iz];
      }
    }


  memset(t22[0],0,(nzpad+2)*(nxpad+2)*sizeof(float));
  // compute ds12de1 and update components

  D1_CC(t22,sxz_vf,nzpad+1,nxpad+2);

  for(ix=0; ix<nxpad+2; ix++){
    for(iz=0; iz<nxpad+2; iz++){
      vx_ff[ix][iz]+= dtz*t22[ix][iz]/rho[ix][iz];
      }
    }


  memset(t22[0],0,(nzpad+2)*(nxpad+2)*sizeof(float));
  // compute ds11de1 and update components

  D1_CC(t22,szz_vf,nzpad+1,nxpad+2);

  for(ix=0; ix<nxpad+2; ix++){
    for(iz=0; iz<nzpad+2; iz++){

      vz_ff[ix][iz]+= dtz*t22[ix][iz]/rho[ix][iz];

      }
    }

  memset(t22[0],0,(nzpad+2)*(nxpad+2)*sizeof(float));

  // compute ds21de2 and update components

  D2_CC(t22,sxz_fv,nzpad+2,nxpad+1);

  for(ix=0; ix<nxpad+2; ix++){
    for(iz=0; iz<nzpad+2; iz++){
      vz_ff[ix][iz]+= dtx*t22[ix][iz]/rho[ix][iz];
      }
    }

}



/****************************************
/ UPDATE STRESSES FOR [v,v] GRIDS - CART
***************************************/


void update_velocity_vv_PML(float **vx_vv, float **vz_vv, float **szz_fv, float **sxz_fv, float **sxz_vf,
                            float **sxx_vf, float **rho, float **t11, float dtx, float dtz, int nzpad, int nxpad)
{
int iz,ix;


memset(t11[0],0,(nzpad+1)*(nxpad+1)*sizeof(float));


// compute ds11de1 and update components

G1_CC(t11,szz_fv,nzpad+1,nxpad+1);

for(ix=0; ix<nxpad+1; ix++){
  for(iz=0; iz<nzpad+1; iz++){
    vz_vv[ix][iz]+= dtz*t11[ix][iz]/rho[ix][iz];
    }
  }

memset(t11[0],0,(nzpad+1)*(nxpad+1)*sizeof(float));

// compute ds12de2 and update components

G2_CC(t11,sxz_vf,nzpad+1,nxpad+1);

for(ix=0; ix<nxpad+1; ix++){
  for(iz=0; iz<nzpad+1; iz++){
    vz_vv[ix][iz]+= dtx*t11[ix][iz]/rho[ix][iz];
    }
  }

memset(t11[0],0,(nzpad+1)*(nxpad+1)*sizeof(float));

// compute ds22de2 and update components

G2_CC(t11,sxx_vf,nzpad+1,nxpad+1);

for(ix=0; ix<nxpad+1; ix++){
  for(iz=0; iz<nzpad+1; iz++){
    vx_vv[ix][iz]+= dtx*t11[ix][iz]/rho[ix][iz];
    }
  }

memset(t11[0],0,(nzpad+1)*(nxpad+1)*sizeof(float));

// compute ds12de1 and update components

G1_CC(t11,sxz_fv,nzpad+1,nxpad+1);

for(ix=0; ix<nxpad+1; ix++){
  for(iz=0; iz<nzpad+1; iz++){
    vx_vv[ix][iz]+= dtz*t11[ix][iz]/rho[ix][iz];
    }
  }

}
