/* Convolutional Perfectly Matched Layer */
/* Reference: An unsplit convolutional perfectly matched layer improved
at grazing incidence for the seismic wave equation : Komatitsch and Martin 2007 */

#include "math.h"

void pml_coeff(float *a_x, float *a_x_half, float *b_x, float *b_x_half, float *K_x, float *K_x_half,
               float *a_z, float *a_z_half, float *b_z, float *b_z_half, float *K_z, float *K_z_half,
               float damping, float fpml, float dx, float dz, float dt, int npml, int nx, int nz)
{

  const float npower = 2.0; /*  power to compute d0 profile */
  const float k_max_PML = 1.0; /* from Gedney page 8.11 */
  const float PI = 3.141592653589793;
  const float alpha_max_PML = 2.0 * PI * (fpml/2.0); /* from festa and Vilotte */

  float thickness_PML_x, thickness_PML_z, xoriginleft, xoriginright, zoriginbottom, zorigintop;
  float Rcoef , d0_x, d0_z, xval, zval, abscissa, abscissa_norm;
  float d_x, d_x_half, d_z, d_z_half;
  float alpha_prime_x, alpha_prime_x_half, alpha_prime_z, alpha_prime_z_half;


  /* thickness of the PML layer in meters */
  thickness_PML_x = npml*dx;
  thickness_PML_z = npml*dz;

  /* reflection coefficient (INRIA report section 6.1) */
  Rcoef = 0.001;

	/* local variables */
	int i, h;

  /* define profile of absorption in PML region */

  /* compute d0 from INRIA report section 6.1 */
  d0_x = - (npower + 1) * damping * log(Rcoef) / (2.0 * thickness_PML_x);
  d0_z = - (npower + 1) * damping * log(Rcoef) / (2.0 * thickness_PML_z);

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

    xval = dx*i;

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
     if(abs(d_x) > 1.0e-6){ a_x[i] = d_x * (b_x[i] - 1.0) / (K_x[i] * (d_x[i] + K_x[i] * alpha_prime_x));}
     if(abs(d_x_half) > 1.0e-6){ a_x_half[i] = d_x_half * (b_x_half[i] - 1.0) / (K_x_half[i] * (d_x_half + K_x_half[i] * alpha_prime_x_half));}

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

	  xval = dx*i;

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
      if(abs(d_x) > 1.0e-6){ a_x[h] = d_x * (b_x[h] - 1.0) / (K_x[h] * (d_x + K_x[h] * alpha_prime_x));}
      if(abs(d_x_half) > 1.0e-6){ a_x_half[h] = d_x_half * (b_x_half[h] - 1.0) / (K_x_half[h] * (d_x_half + K_x_half[h] * alpha_prime_x_half));}

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
    if(abs(d_z) > 1.0e-6){ a_z[i] = d_z * (b_z[i] - 1.0) / (K_z[i] * (d_z + K_z[i] * alpha_prime_z));}
    if(abs(d_z_half) > 1.0e-6){ a_z_half[i] = d_z_half * (b_z_half[i] - 1.0) / (K_z_half[i] * (d_z_half + K_z_half[i] * alpha_prime_z_half));}

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
    abscissa = zval - zorigintop;

    if(abscissa >= 0.0){

      abscissa_norm = abscissa / thickness_PML_z;
      d_z = d0_z * pow(abscissa_norm,npower);

      /* this taken from Gedney page 8.2 */
      K_z[h] = 1.0 + (k_max_PML - 1.0) * pow(abscissa_norm,npower);
      alpha_prime_z = alpha_max_PML * (1.0 - abscissa_norm);

        }

    /* define damping profile at half the grid points */
    abscissa = zval + dz/2.0 - zorigintop;

    if(abscissa_in_PML >= 0.0){

        abscissa_norm = abscissa / thickness_PML_z;
        d_z_half = d0_z * pow(abscissa_norm,npower);

        /* this taken from Gedney page 8.2 */
        K_z_half[h] = 1.0 + (k_max_PML - 1.0) * pow(abscissa_norm,npower);
        alpha_prime_z_half = alpha_max_PML * (1.0 - abscissa_norm);

          }

     b_z[h] = exp(- (d_z / K_z[h] + alpha_prime_z) * dt);
     b_z_half[h] = exp(- (d_z_half / K_z_half[h] + alpha_prime_z_half) * dt);

     /* avoid division by zero outside the PML */
     if(abs(d_z) > 1.0e-6){ a_z[h] = d_z * (b_z[h] - 1.0) / (K_z[h] * (d_z + K_z[h] * alpha_prime_z));}
     if(abs(d_z_half) > 1.0e-6){ a_z_half[h] = d_z_half * (b_z_half[h] - 1.0) / (K_z_half[h] * (d_z_half + K_z_half[h] * alpha_prime_z_half));}

       } /* end of top boundary */

}

/* PML FSG update */
void update_PML_vv(float **wav, float **conv, float *K_x, float *K_z,  float *b_x, float *b_z, \
                   float *a_x, float *a_z, int npml)
{
  int ix,iz;



}


void update_PML_ff(float **wav, float **conv, float *K_x_half, float *K_z_half, float *b_x_half, float *b_z_half, \
                   float *a_x_half, float *a_z_half, int npml)
{
  int ix,iz;


}


void update_PML_fv(float **wav, float **conv, float *K_x, float *K_z, float *b_x, float *a_x, float *a_z, \
                   float *K_x_half, float *K_z_half, float *b_x_half, float *b_z_half, float *a_x_half, \
                   float *a_z_half, int npml)
{
  int ix,iz;


}

void update_PML_vf(float **wav, float **conv, float *K_x, float *K_z, float *b_x, float *a_x, float *a_z, \
                   float *K_x_half, float *K_z_half, float *b_x_half, float *b_z_half, float *a_x_half, \
                   float *a_z_half, int npml)
{
  int ix,iz;


}
