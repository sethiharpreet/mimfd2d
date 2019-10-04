#include <rsf.h>
#include <stdio.h>
#include <stdlib.h>
#include "fdutil.h"




/*-------------------General Utilities ---------------------------*/

void expand2d(float** in, float** out, fdm2d fdm)
/* expand domain considering the free surface */
{
    int iz,ix,ixp;

    for	(ix=0;ix<fdm->nx;ix++) {
		    for (iz=0;iz<fdm->nz;iz++) {
	    	    out[fdm->nb+ix][iz] = in[ix][iz];
		    }
    }

    for (ix=0; ix<fdm->nxpad+2; ix++) {
		    for (iz=0; iz<fdm->nb+2; iz++) {
	    	    out[ix][fdm->nzpad-iz+1] = out[ix][fdm->nzpad-fdm->nb-1];
		        }
    }

    for (ix=0; ix<fdm->nb+2;    ix++) {
		    for (iz=0; iz<fdm->nzpad+2; iz++) {
            ixp = (ix<fdm->nb) ? ix : fdm->nb-1;
		        out[ixp][iz] = out[fdm->nb][iz];
	    	    out[fdm->nxpad-ix+1][iz] = out[fdm->nxpad-fdm->nb-1][iz];
		}
    }
}

void expand2d_bottom(float **in, float **out, fdm2d fdm)
/*< expand domain with free surface at top >*/
{
    int iz,ix,ixp;

    for(ix=0; ix<fdm->nx; ix++){
	     for(iz=0; iz<fdm->nz; iz++){
	        out[fdm->nb+ix][iz] = in[ix][iz];
	       }
       }

       for (ix=0; ix<fdm->nxpad+2; ix++) {
   		    for (iz=0; iz<fdm->nb+2; iz++) {
   	    	    out[ix][fdm->nzpad-iz+1] = out[ix][fdm->nzpad-fdm->nb-1];
   		        }
       }

       for (ix=0; ix<fdm->nb+2; ix++) {
   		    for (iz=0; iz<fdm->nzpad+2; iz++) {
               ixp = (ix<fdm->nb) ? ix : fdm->nb-1;
   		        out[ixp][iz] = out[fdm->nb][iz];
   	    	    out[fdm->nxpad-ix+1][iz] = out[fdm->nxpad-fdm->nb-1][iz];
   		       }
       }
}

void expand2d_top(float **in, float **out, fdm2d fdm )
/*< expand domain with free surface at top >*/
{
    int iz,ix,ixp;

    for(ix=0; ix<fdm->nx; ix++){
	     for(iz=0; iz<fdm->nz; iz++){
	        out[fdm->nb+ix][iz] = in[ix][iz];
	       }
       }

    for (ix=0; ix<fdm->nxpad+2; ix++) {
   	    out[ix][fdm->nzpad+1] = out[ix][fdm->nzpad-1];
        out[ix][fdm->nzpad]   = out[ix][fdm->nzpad-1];
       }

    for (ix=0; ix<fdm->nb+2;    ix++) {
   		 for (iz=0; iz<fdm->nzpad+2; iz++) {
          ixp = (ix<fdm->nb) ? ix : fdm->nb-1;
   		    out[ixp][iz] = out[fdm->nb][iz];
   	    	out[fdm->nxpad-ix+1][iz] = out[fdm->nxpad-fdm->nb-1][iz];
   		   }
       }

 }

void window2d(float **in, float **out, fdm2d fdm)
/*< extract domain asssuming free surface at top >*/
{
  int ix,iz;

  for (ix=0;ix<fdm->nx;ix++) {
    for (iz=0;iz<fdm->nz;iz++) {
      out[ix][iz] = in[fdm->nb+ix][iz];
    }
  }

}
/*-------------------- Boundary conditions -------------------------------*/

void sponge_ex(float **u, int nxpad, int nzpad, int nb)
/*< apply absorbing boundary condition >*/
{
	int ix,iz,ib,ibx,ibz;
	float w;

  for(ib=0; ib<nb; ib++) {
		float tmp=ib/(sqrt(2.0)*4*nb);;
		w=expf(-tmp*tmp);
		ibz = nzpad-ib-1;
		for(ix=0; ix<nxpad; ix++) {
			u[ix][ibz] *= w; /* bottom sponge */
			}

		ibx = nxpad-ib-1;
		for(iz=0; iz<nzpad; iz++) {
			u[ib ][iz] *= w; /*   left sponge */
			u[ibx][iz] *= w; /*  right sponge */
			}
    }
}


/*-------------------- Mimetic Operations -----------------------------*/

/****************************************
* Divergence on Dimension 1 (Z)
*****************************************/

void D1_CC(float **d1, float **dd, int nz, int nx)
/***************************************
    d1 : Output after Divergence operator
    dd : vv or ff variables
    nz : nzpad+1 for vv grid or nzpad+2 for ff grid
    nx : nxpad+1 for vv grid or nxpad+2 for ff grid
    k  : order
****************************************/
{
  int ix,iz;
  const float D[14]={ -11.0/12.0, 17.0/24.0, 3.0/8.0, -5.0/24.0, 1.0/24.0, \
                      -1.0/24.0,  5.0/24.0, -3.0/8.0, -17.0/24.0, 11.0/12.0, \
                       1.0/24.0, -9.0/8.0,   9.0/8.0, -1.0/24.0 };

  for(ix=0; ix<nx; ix++){
    for(iz=1; iz<nz-1; iz++){
      if(iz==1)         d1[ix][iz]= D[0]*dd[ix][0]+ D[1]*dd[ix][1]+ D[2]*dd[ix][2]+ \
                                    D[3]*dd[ix][3]+ D[4]*dd[ix][4];

      else if(iz==nz-2) d1[ix][iz]= D[5]*dd[ix][nz-5]+ D[6]*dd[ix][nz-4]+ \
                                    D[7]*dd[ix][nz-3]+ D[8]*dd[ix][nz-2]+ D[9]*dd[ix][nz-1];

      else              d1[ix][iz]= D[10]*dd[ix][iz-2]+ D[11]*dd[ix][iz-1]+ D[12]*dd[ix][iz]+  \
                                    D[13]*dd[ix][iz+1];
                  }
                }

  }


/****************************************
* Divergence on Dimension 2 (X)
*****************************************/
void D2_CC(float **d2, float **dd, int nz, int nx)
/***************************************
    d2 : Output after Divergence operator
    dd : vv or ff variables
    nz : nzpad+1 for vv grid or nzpad+2 for ff grid
    nx : nxpad+1 for vv grid or nxpad+2 for ff grid
    s  : starting index
    k  : order
****************************************/
{
  int ix,iz;
  const float D[14]={ -11.0/12.0, 17.0/24.0, 3.0/8.0, -5.0/24.0, 1.0/24.0, \
                      -1.0/24.0,  5.0/24.0,  -3.0/8.0, -17.0/24.0, 11.0/12.0, \
                       1.0/24.0, -9.0/8.0,   9.0/8.0,  -1.0/24.0 };

  for(ix=1; ix<nx-1; ix++){
    for(iz=0; iz<nz; iz++){

      if(ix==1)         d2[ix][iz]= D[0]*dd[0][iz]+ D[1]*dd[1][iz]+ D[2]*dd[2][iz]+ \
                                    D[3]*dd[3][iz]+ D[4]*dd[4][iz];

      else if(ix==nx-2) d2[ix][iz]= D[5]*dd[nx-5][iz] +D[6]*dd[nx-4][iz]+ \
                                    D[7]*dd[nx-3][iz]+ D[8]*dd[nx-2][iz]+ D[9]*dd[nx-1][iz];

      else              d2[ix][iz]= D[10]*dd[ix-2][iz]+ D[11]*dd[ix-1][iz]+ D[12]*dd[ix][iz]+  \
                                    D[13]*dd[ix+1][iz];
                   }
                 }
    }

/****************************************
* Gradient on Dimension 1 (Z)
*****************************************/

void G1_CC(float **g1,float **gg, int nz, int nx)
/***************************************
    g1 : Output after Gradient operator along z axis
    gg : vv or ff variables
    nz : nzpad+1 for vv grid or nzpad+2 for ff grid
    nx : nxpad+1 for vv grid or nxpad+2 for ff grid
    k  : order
****************************************/
{

  int ix,iz;

  const float G[24] = {-352.0/105.0, 35.0/8.0, -35.0/24.0, 21.0/40.0, -5.0/56.0, \
                        16.0/105.0, -31.0/24.0, 29.0/24.0, -3.0/40.0, 1.0/168.0, \
                        -1.0/168.0,  3.0/40.0, -29.0/24.0, 31.0/24.0, -16.0/105.0, \
                        5.0/56.0,   -21.0/40.0, 35.0/24.0, -35.0/8.0, 352.0/105.0, \
                        1.0/24.0,   -9.0/8.0,   9.0/8.0,  -1.0/24.0
                      };


  for(ix=0; ix<nx; ix++){
    for(iz=0; iz<nz; iz++){

      if(iz==0)         g1[ix][iz]= G[0]*gg[ix][0]+ G[1]*gg[ix][1]+ G[2]*gg[ix][2]+ \
                                    G[3]*gg[ix][3]+ G[4]*gg[ix][4];

      else if(iz==1)    g1[ix][iz]= G[5]*gg[ix][0]+ G[6]*gg[ix][1]+ G[7]*gg[ix][2]+ \
                                    G[8]*gg[ix][3]+ G[9]*gg[ix][4];

      else if(iz==nz-2) g1[ix][iz]= G[10]*gg[ix][nz-5]+ G[11]*gg[ix][nz-4]+ G[12]*gg[ix][nz-3]+ \
                                    G[13]*gg[ix][nz-2]+ G[14]*gg[ix][nz-1];

      else if(iz==nz-1) g1[ix][iz]= G[15]*gg[ix][nz-5]+ G[16]*gg[ix][nz-4]+G[17]*gg[ix][nz-3]+ \
                                    G[18]*gg[ix][nz-2]+ G[19]*gg[ix][nz-1];

      else              g1[ix][iz]= G[20]*gg[ix][iz-1]+ G[21]*gg[ix][iz]+ G[22]*gg[ix][iz+1]+  \
                                    G[23]*gg[ix][iz+2];
            }
       }

}

/****************************************
* Gradient on Dimension 2 (X)
*****************************************/

void G2_CC(float **g2,float **gg, int nz, int nx)
/***************************************
    g2 : Output after Gradient operator along x axis
    gg : vv or ff variables
    nz : nzpad+1 for vv grid or nzpad+2 for ff grid
    nx : nxpad+1 for vv grid or nxpad+2 for ff grid
    k  : order
****************************************/
{
  int iz,ix;

  // Corbino coefficients
  const float G[24] = {-352.0/105.0, 35.0/8.0, -35.0/24.0, 21.0/40.0, -5.0/56.0,   \
                        16.0/105.0, -31.0/24.0, 29.0/24.0, -3.0/40.0, 1.0/168.0,   \
                        -1.0/168.0,  3.0/40.0, -29.0/24.0, 31.0/24.0, -16.0/105.0, \
                        5.0/56.0,   -21.0/40.0, 35.0/24.0, -35.0/8.0, 352.0/105.0, \
                        1.0/24.0,   -9.0/8.0,   9.0/8.0,  -1.0/24.0
                      };

  for(ix=0; ix<nx; ix++){
    for(iz=0; iz<nz; iz++){

      if(ix==0)         g2[ix][iz]= G[0]*gg[0][iz]+ G[1]*gg[1][iz]+ G[2]*gg[2][iz]+ \
                                    G[3]*gg[3][iz]+ G[4]*gg[4][iz];

      else if(ix==1)    g2[ix][iz]= G[5]*gg[0][iz]+ G[6]*gg[1][iz]+ G[7]*gg[2][iz]+ \
                                    G[8]*gg[3][iz]+ G[9]*gg[4][iz];

      else if(ix==nx-2) g2[ix][iz]= G[10]*gg[nx-5][iz]+ G[11]*gg[nx-4][iz]+ G[12]*gg[nx-3][iz]+ \
                                    G[13]*gg[nx-2][iz]+ G[14]*gg[nx-1][iz];

      else if(ix==nx-1) g2[ix][iz]= G[15]*gg[nx-5][iz]+ G[16]*gg[nx-4][iz]+ G[17]*gg[nx-3][iz]+ \
                                    G[18]*gg[nx-2][iz]+ G[19]*gg[nx-1][iz];

      else              g2[ix][iz]= G[20]*gg[ix-1][iz]+ G[21]*gg[ix][iz]+ G[22]*gg[ix+1][iz]+  \
                                    G[23]*gg[ix+2][iz];

                   }
                 }

}


/***************************************************************
/ Acoustic functions
/ P-wave velocity and density
/
***************************************************************/


/* Pressure Update */

void update_pressure_fv(float **p_fv, float **vx_ff, float **vz_vv, float **t21, float **ro, float **vp, float dtx, float dtz, int nzpad, int nxpad)
/* vx[f,f] -------> p[f,v]   */
/* vz[v,v] -------> p[f,v]   */
{
  int iz,ix;

  memset(t21[0],0,(nzpad+2)*(nxpad+1)*sizeof(float));

  D1_CC(t21,vz_vv,nzpad+1,nxpad+1);

  for (ix=0; ix<nxpad+1; ix++){
    for(iz=0; iz<nzpad+2; iz++){
      p_fv[ix][iz]-= dtz*ro[ix][iz]*vp[ix][iz]*vp[ix][iz]*t21[ix][iz];
    }
  }

  memset(t21[0],0,(nzpad+2)*(nxpad+1)*sizeof(float));

  G2_CC(t21,vx_ff,nzpad+2,nxpad+1);

  for (ix=0; ix<nxpad+1; ix++){
    for(iz=0; iz<nzpad+2; iz++){
      p_fv[ix][iz]-= dtx*ro[ix][iz]*vp[ix][iz]*vp[ix][iz]*t21[ix][iz];
    }
  }

}


void update_pressure_vf(float **p_vf, float **vz_ff, float **vx_vv, float **t12, float **ro, float **vp, float dtx, float dtz, int nzpad, int nxpad)
/* vx[f,f] -------> p[f,v]   */
/* vz[v,v] -------> p[f,v]   */
{
  int iz,ix;

  memset(t12[0],0,(nzpad+1)*(nxpad+2)*sizeof(float));

  G1_CC(t12,vz_ff,nzpad+1,nxpad+2);

  for (ix=0; ix<nxpad+2; ix++){
    for(iz=0; iz<nzpad+1; iz++){
      p_vf[ix][iz]-= dtz*ro[ix][iz]*vp[ix][iz]*vp[ix][iz]*t12[ix][iz];
    }
  }

  memset(t12[0],0,(nzpad+1)*(nxpad+2)*sizeof(float));

  D2_CC(t12,vx_vv,nzpad+1,nxpad+1);

  for (ix=0; ix<nxpad+2; ix++){
    for(iz=0; iz<nzpad+1; iz++){
      p_vf[ix][iz]-= dtx*ro[ix][iz]*vp[ix][iz]*vp[ix][iz]*t12[ix][iz];
    }
  }

}

/* Acoustic velocity Update */

void update_ac_velocity_ff(float **vz_ff, float **vx_ff, float ** p_vf, float **p_fv, float **t22, float **ro, float dtx, float dtz, int nzpad, int nxpad)
/* p[f,v] -------> vx[f,f]   */
/* p[v,f] -------> vz[f,f]   */
{

  int iz,ix;

  memset(t22[0],0,(nzpad+2)*(nxpad+2)*sizeof(float));

  D1_CC(t22,p_vf,nzpad+1,nxpad+2);

  for (ix=0; ix<nxpad+2; ix++){
    for(iz=0; iz<nzpad+2; iz++){
      vz_ff[ix][iz]-= dtz*t22[ix][iz]/ro[ix][iz];
    }
  }

  memset(t22[0],0,(nzpad+2)*(nxpad+2)*sizeof(float));

  D2_CC(t22,p_fv,nzpad+2,nxpad+1);

  for (ix=0; ix<nxpad+2; ix++){
    for(iz=0; iz<nzpad+2; iz++){
      vx_ff[ix][iz]-= dtx*t22[ix][iz]/ro[ix][iz];
    }
  }

}



void update_ac_velocity_vv(float **vz_vv, float **vx_vv, float **p_vf, float **p_fv, float **t11, float **ro, float dtx, float dtz, int nzpad, int nxpad)
/* p[f,v] -------> vz[v,v]   */
/* p[v,f] -------> vx[v,v]   */
{

  int iz,ix;

  memset(t11[0],0,(nzpad+1)*(nxpad+1)*sizeof(float));

  G1_CC(t11,p_fv,nzpad+1,nxpad+1);

  for (ix=0; ix<nxpad+1; ix++){
    for(iz=0; iz<nzpad+1; iz++){
      vz_vv[ix][iz]-= dtz*t11[ix][iz]/ro[ix][iz];
    }
  }

  memset(t11[0],0,(nzpad+1)*(nxpad+1)*sizeof(float));

  G2_CC(t11,p_vf,nzpad+1,nxpad+1);

  for (ix=0; ix<nxpad+1; ix++){
    for(iz=0; iz<nzpad+1; iz++){
      vx_vv[ix][iz]-= dtx*t11[ix][iz]/ro[ix][iz];
    }
  }
}

/***************************************************************
/ Elastic (Anisotropic) orthorombic footprint - 4 coefficients
/ c11 c13
/   .   c33
/           c55
/  VTI
***************************************************************/


/****************************************
/ UPDATE STRESSES FOR [v,f] GRIDS - CART
*****************************************/


void update_stress_vf(float **vz_ff,  float **vx_ff,  float **vz_vv,  float **vx_vv,          \
                      float **sxx_vf, float **sxz_vf, float **szz_vf, float **rho,            \
                      float **c11,    float **c13,    float **c33,    float **c55,            \
                      float **t12,    float dtx,      float dtz,      int nzpad,  int nxpad)
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

      // S12 or S21 component
      // s_zx=c55*tzx (dv1de2 part)
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
      // s_zz=c13*txx
      szz_vf[ix][iz]+= dtx*c13[ix][iz]*t12[ix][iz];

      // S22 component
      // s_xx=c11*txx
      sxx_vf[ix][iz]+= dtx*c11[ix][iz]*t12[ix][iz];

        }
      }

  memset(t12[0],0,(nzpad+1)*(nxpad+2)*sizeof(float));


  // compute DV1DE1 and update components
  G1_CC(t12,vz_ff,nzpad+1,nxpad+2);

  for(ix=0; ix<nxpad+2; ix++){
    for(iz=0; iz<nzpad+1; iz++){

      // S11 component
      // s_zz=c33*tzz
      szz_vf[ix][iz]+= dtz*c33[ix][iz]*t12[ix][iz];

      // S22 component
      // s_xx=c13*tzz
      sxx_vf[ix][iz]+= dtz*c13[ix][iz]*t12[ix][iz];
    }
  }

  // compute DV2DE1 and update components

  memset(t12[0],0,(nzpad+1)*(nxpad+2)*sizeof(float));

  G1_CC(t12,vx_ff,nzpad+1,nxpad+2);

  for(ix=0; ix<nxpad+2; ix++){
    for(iz=0; iz<nzpad+1; iz++){

      // S12 or S21 component
      // s_zx=c55*tzx (dv2de1 part)
      sxz_vf[ix][iz]+= dtz*c55[ix][iz]*t12[ix][iz];

      }
    }

  }


void update_stress_fv(float **vz_ff,  float **vx_ff,  float **vz_vv,  float **vx_vv,          \
                      float **sxx_fv, float **sxz_fv, float **szz_fv, float **rho,            \
                      float **c11,    float **c13,    float **c33,    float **c55,            \
                      float **t21,    float dtx,      float dtz,      int nzpad, int nxpad)
{
    int iz,ix;

    memset(t21[0],0,(nzpad+2)*(nxpad+1)*sizeof(float));

    //compute DV1DE2 and update components

    //Gradient
    G2_CC(t21,vz_ff,nzpad+2,nxpad+1);

    for(ix=0; ix<nxpad+1; ix++){
      for(iz=0; iz<nzpad+2; iz++){

        // S12 or S21 component
        // s_zx=c55*tzx (dv1de2 part)
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
        // s_zz=c13*txx
        szz_fv[ix][iz]+= dtx*c13[ix][iz]*t21[ix][iz];

        // S22 component
        // s_xx=c11*txx
        sxx_fv[ix][iz]+= dtx*c11[ix][iz]*t21[ix][iz];

          }
        }

    memset(t21[0],0,(nzpad+2)*(nxpad+1)*sizeof(float));

    // compute DV1DE1 and update components
    D1_CC(t21,vz_vv,nzpad+1,nxpad+1);

    for(ix=0; ix<nxpad+1; ix++){
      for(iz=0; iz<nzpad+2; iz++){

        // S11 component
        // s_zz=c33*tzz
        szz_fv[ix][iz]+= dtz*c33[ix][iz]*t21[ix][iz];

        // S22 component
        // s_xx=c13*tzz
        sxx_fv[ix][iz]+= dtz*c13[ix][iz]*t21[ix][iz];

      }
    }


    memset(t21[0],0,(nzpad+2)*(nxpad+1)*sizeof(float));
    // compute  DV2DE1 and update components

    D1_CC(t21,vx_vv,nzpad+1,nxpad+1);

    for(ix=0; ix<nxpad+1; ix++){
      for(iz=0; iz<nzpad+2; iz++){

        // S12 or S21 component
        // s_zx=c55*tzx (dv2de1 part)
        sxz_fv[ix][iz]+= dtz*c55[ix][iz]*t21[ix][iz];

        }
      }

  }

  /***************************************************************
  / Elastic (Anisotropic) orthorombic footprint - 4 coefficients
  / c11 c13   c15
  /   .   c33 c35
  /           c55
  /  TTI
  ***************************************************************/


  /****************************************
  / UPDATE STRESSES FOR [v,f] GRIDS - CART : TTI
  *****************************************/


  void update_stress_vf_tti(float **vz_ff,  float **vx_ff,  float **vz_vv,  float **vx_vv,          \
                            float **sxx_vf, float **sxz_vf, float **szz_vf, float **rho,            \
                            float **c11,    float **c13,    float **c33,    float **c55,            \
                            float **c15,    float **c35,    float **t12,    float dtx, float dtz,   \
                            int nzpad,      int nxpad)
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

    void update_stress_fv_tti(float **vz_ff,  float **vx_ff,  float **vz_vv,  float **vx_vv,          \
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

void update_velocity_ff(float **vx_ff, float ** vz_ff, float **sxx_fv, float **sxz_fv, float **sxz_vf, \
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


void update_velocity_vv(float **vx_vv, float **vz_vv, float **szz_fv, float **sxz_fv, float **sxz_vf,
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


/******************************************
/ Sourcce Injection
*******************************************/

void inject_ps_src_cart(float **p_fv,float **p_vf, float*** ww, lint2d cs, int it)
/* Pressure or Stress source*/
{
  int is2 = cs->jx[0];
  int is1 = cs->jz[0];
  // [v,f] grid
  p_vf[is1][is2]      =   p_vf[is1][is2]     + ww[it][0][0]/4.0;
  p_vf[is1+1][is2]    =   p_vf[is1+1][is2]   + ww[it][0][0]/4.0;

  // [f,v] grids
  p_fv[is1][is2+1]    =   p_fv[is1][is2+1]   + ww[it][0][0]/4.0;
  p_fv[is1+1][is2+1]  =   p_fv[is1+1][is2+1] + ww[it][0][0]/4.0;

}

void inject_vel_src_cart(float **v_ff,float **v_vv, float*** ww, lint2d cs, int it)
/* Vz  source*/
{
  int is2 = cs->jx[0];
  int is1 = cs->jz[0];
  // [f,f] grid
  v_ff[is2+1][is1+1]   =   v_ff[is2+1][is1+1] + ww[it][0][0];

  //[v,v] grid
  v_vv[is2][is1]       =   v_vv[is2][  is1]   + ww[it][0][0]/4.0;
  v_vv[is2+1][is1]     =   v_vv[is2+1][is1]   + ww[it][0][0]/4.0;
  v_vv[is2][is1+1]     =   v_vv[is2][is1+1]   + ww[it][0][0]/4.0;
  v_vv[is2+1][is1+1]    =   v_vv[is2+1][is1+1] + ww[it][0][0]/4.0;

}

/*----------------------------------------------------------------------------------------*/

/*****************************************
/ Coupling conditions
*****************************************/

/* Update Mimetic points using coupling conditions */

/* [v] grid coupling */
void apply_vel_cpl_vv(float **vz1_vv, float **vz2_vv, int nz, int nx)
/* Velocity continuous  */
{
  int ix;

  for(ix=0;ix<nx+1;ix++){
    vz1_vv[ix][nz+1]  = 0.5*(vz1_vv[ix][nz+1]+vz2_vv[ix][0]);
    vz2_vv[ix][0]     = vz1_vv[ix][nz+1];
  }
}

void apply_stress_cpl_vf(float **p1_vf, float **szz_vf, float **sxz_vf, int nz, int nx)
/* Stress continuous */
{
  int ix;

  for(ix=0;ix<nx+2;ix++){
    p1_vf[ix][nz+1]  = 0.5*(p1_vf[ix][nz+1]-szz_vf[ix][0]);
    szz_vf[ix][0]  = -p1_vf[ix][nz+1];
    sxz_vf[ix][0]  = 0.0;
  }
}


/* [f] grid coupling : update mimetic points using derivatives*/


void apply_vel_cpl_ff(float **vz1_ff, float **vz2_ff, float **vx1_ff, float **vx2_ff, \
                      float **vz2_vv, float **vx2_vv, float **vx1_vv, \
                      float **ro, float **vp, float **c13, float **c33, float idz, float idx, \
                      int nzl, int nxl)

/* Compute mimetic velocity from derivatives */
{

  int ix, nz, nx;
  float k1,k2,k3,k4;

  const float G[24] = {-352.0/105.0, 35.0/8.0, -35.0/24.0, 21.0/40.0, -5.0/56.0,   \
                        16.0/105.0, -31.0/24.0, 29.0/24.0, -3.0/40.0, 1.0/168.0,   \
                        -1.0/168.0,  3.0/40.0, -29.0/24.0, 31.0/24.0, -16.0/105.0, \
                        5.0/56.0,   -21.0/40.0, 35.0/24.0, -35.0/8.0, 352.0/105.0, \
                        1.0/24.0,   -9.0/8.0,   9.0/8.0,  -1.0/24.0 };

  const float D[14]= {-11.0/12.0, 17.0/24.0, 3.0/8.0, -5.0/24.0, 1.0/24.0, \
                      -1.0/24.0,  5.0/24.0,  -3.0/8.0, -17.0/24.0, 11.0/12.0, \
                       1.0/24.0, -9.0/8.0,   9.0/8.0,  -1.0/24.0 };


  /* Update vx_ff */
  /* compute D2(vz2_vv) */
  /* Loop unrolled */

  nz = nzl+1;
  nx = nxl+1;

  vx2_ff[0][0]    = 0.0;                                                        //ix=0
  vx2_ff[1][0]    = idx*(D[0]*vz2_vv[0][0]+ D[1]*vz2_vv[1][0]+ D[2]*vz2_vv[2][0]+ \
                         D[3]*vz2_vv[3][0]+ D[4]*vz2_vv[4][0]);                       //ix=1
  vx2_ff[nx-2][0] = idx*(D[5]*vz2_vv[nx-5][0] +D[6]*vz2_vv[nx-4][0]+ \
                         D[7]*vz2_vv[nx-3][0]+ D[8]*vz2_vv[nx-2][0]+ D[9]*vz2_vv[nx-1][0]); //ix=nx-2
  vx2_ff[nx-1][0] = 0.0;                                                              //ix=nx-1

  for(ix=2; ix<nx-2; ix++){
    vx2_ff[ix][0] = idx*(D[10]*vz2_vv[ix-2][0]+ D[11]*vz2_vv[ix-1][0]+ D[12]*vz2_vv[ix][0]+  \
                         D[13]*vz2_vv[ix+1][0]);
      }

  /* compute G1(vxn2_ff) and update mimetic point vx2_ff */

  nz = nzl+2;
  nx = nxl+2;

  for(ix=0; ix<nx; ix++){

    vx2_ff[ix][0]-= idz*(G[1]*vx2_ff[ix][1]+ G[2]*vx2_ff[ix][2]+ G[3]*vx2_ff[ix][3]+ G[4]*vx2_ff[ix][4]);
    // { D2(vz2_ff)-G1(vx2_ff) } /Gz[0]
    vx2_ff[ix][0]/= idz*G[0];

     }

  /* Update vz_ff */

  nz = nzl+1;
  nx = nxl+1;

  /* compute D2(vx1_vv) */
  k1=ro[1][nz-1]*vp[1][nz-1]*vp[1][nz-1];
  k2=ro[nx-2][nz-1]*vp[nx-2][nz-1]*vp[nx-2][nz-1];

  vz2_ff[0][0]    = 0.0;
  vz2_ff[1][0]    = idx*k1*(D[0]*vx1_vv[0][nz-1] +D[1]*vx1_vv[1][nz-1]+ D[2]*vx1_vv[2][nz-1]+ \
                            D[3]*vx1_vv[3][nz-1]+ D[4]*vx1_vv[4][nz-1]);
  vz2_ff[nx-2][0] = idx*k2*(D[5]*vx1_vv[nx-5][nz-1] +D[6]*vx1_vv[nx-4][nz-1]+ \
                            D[7]*vx1_vv[nx-3][nz-1] +D[8]*vx1_vv[nx-2][nz-1]+ D[9]*vx1_vv[nx-1][nz-1]);
  vz2_ff[nx-1][0] = 0.0;

  for(ix=2; ix<nx-2; ix++){
    k3 = ro[ix][nz-1]*vp[ix][nz-1]*vp[ix][nz-1];

    vz2_ff[ix][0] = idx*k3*(D[10]*vx2_vv[ix-2][0]+ D[11]*vx2_vv[ix-1][0]+ D[12]*vx2_vv[ix][0]+  \
                            D[13]*vx2_vv[ix+1][0]);
      }

  nz = nzl+1;
  nx = nxl+1;

  /* compute D2(vx2_vv) */
  vz2_ff[0][0]    -= 0.0;
  vz2_ff[1][0]    -= idx*c13[1][0]*(D[0]*vx2_vv[0][0]+ D[1]*vx2_vv[1][0]+ D[2]*vx2_vv[2][0]+ \
                                    D[3]*vx2_vv[3][0]+ D[4]*vx2_vv[4][0]);

  vz2_ff[nx-2][0] -= idx*c13[nx-2][0]*(D[5]*vx2_vv[nx-5][0] +D[6]*vx2_vv[nx-4][0]+ \
                                       D[7]*vx2_vv[nx-3][0] +D[8]*vx2_vv[nx-2][0]+ D[9]*vx2_vv[nx-1][0]);

  vz2_ff[nx-1][0] -= 0.0;

  for(ix=2; ix<nx-2; ix++){
    vz2_ff[ix][0]  -= idx*c13[ix][0]*(D[10]*vx2_vv[ix-2][0]+ D[11]*vx2_vv[ix-1][0]+ D[12]*vx2_vv[ix][0]+  \
                                      D[13]*vx2_vv[ix+1][0]);
      }

  nz = nzl+2;
  nx = nxl+2;

  for(ix=0; ix<nx; ix++){

    k1=ro[ix][nz-1]*vp[ix][nz-1]*vp[ix][nz-1];

    /* compute G1(vz1_ff) */
    vz2_ff[ix][0] += idz*k1*(G[15]*vz1_ff[ix][nz-5]+ G[16]*vz1_ff[ix][nz-4]+G[17]*vz1_ff[ix][nz-3]+ \
                             G[18]*vz1_ff[ix][nz-2]); //G[19]*gg[ix][nz-1];

    /* compute G1(vz2_ff)*/
    vz2_ff[ix][0] -= idz*c33[ix][0]*(G[1]*vz2_ff[ix][1]+ G[2]*vz2_ff[ix][2]+ G[3]*vz2_ff[ix][3]+ G[4]*vz2_ff[ix][4]);

    vz2_ff[ix][0] /= idz*(c33[ix][0]*G[0] - k1*G[19]);

    /* Set vz1_ff = vz2_ff */
    vz1_ff[ix][nz-1] = vz2_ff[ix][0];

      }

    }


void apply_stress_cpl_fv(float **p_fv, float **szz_fv, float **sxz_fv, float **sxz_vf, \
                         float **ro1, float **ro2, float idz, float idx, int nzl, int nxl)
/* Stress continuous */
{
  int ix,nz,nx;

  const float G[24] = {-352.0/105.0, 35.0/8.0, -35.0/24.0, 21.0/40.0, -5.0/56.0,   \
                        16.0/105.0, -31.0/24.0, 29.0/24.0, -3.0/40.0, 1.0/168.0,   \
                        -1.0/168.0,  3.0/40.0, -29.0/24.0, 31.0/24.0, -16.0/105.0, \
                        5.0/56.0,   -21.0/40.0, 35.0/24.0, -35.0/8.0, 352.0/105.0, \
                        1.0/24.0,   -9.0/8.0,   9.0/8.0,  -1.0/24.0 };

  const float D[14]= {-11.0/12.0, 17.0/24.0, 3.0/8.0, -5.0/24.0, 1.0/24.0, \
                      -1.0/24.0,  5.0/24.0,  -3.0/8.0, -17.0/24.0, 11.0/12.0, \
                       1.0/24.0, -9.0/8.0,   9.0/8.0,  -1.0/24.0 };


  nz = nzl+1;
  nx = nxl+1;

    //G2(sxz_vf)
  szz_fv[0][0]    = -idx*ro1[0][nz-1]*(G[0]*sxz_vf[0][0]+ G[1]*sxz_vf[1][0]+ G[2]*sxz_vf[2][0]+ \
                                   G[3]*sxz_vf[3][0]+ G[4]*sxz_vf[4][0]);
  szz_fv[1][0]    = -idx*ro1[1][nz-1]*(G[5]*sxz_vf[0][0]+ G[6]*sxz_vf[1][0]+ G[7]*sxz_vf[2][0]+ \
                                   G[8]*sxz_vf[3][0]+ G[9]*sxz_vf[4][0]);

  szz_fv[nx-2][0] = -idx*ro1[nx-2][nz-1]*(G[10]*sxz_vf[nx-5][0]+ G[11]*sxz_vf[nx-4][0]+ G[12]*sxz_vf[nx-3][0]+ \
                                          G[13]*sxz_vf[nx-2][0]+G[14]*sxz_vf[nx-1][0]);

  szz_fv[nx-1][0] = -idx*ro1[nx-1][nz-1]*(G[15]*sxz_vf[nx-5][0]+ G[16]*sxz_vf[nx-4][0]+ G[17]*sxz_vf[nx-3][0]+ \
                                      G[18]*sxz_vf[nx-2][0]+ G[19]*sxz_vf[nx-1][0]);
  for(ix=2; ix<nx-2; ix++){
    szz_fv[ix][0] = -idx*ro1[ix][nz-1]*(G[20]*sxz_vf[ix-1][0]+ G[21]*sxz_vf[ix][0]+ G[22]*sxz_vf[ix+1][0]+  \
                                        G[23]*sxz_vf[ix+2][0]);
  }


  for(ix=0; ix<nx; ix++){

    //G1(p_fv)
    szz_fv[ix][0] -= idz*ro2[ix][0]*(G[15]*p_fv[ix][nz-5]+ G[16]*p_fv[ix][nz-4]+G[17]*p_fv[ix][nz-3]+ \
                                     G[18]*p_fv[ix][nz-2]); //G[19]*p_fv[ix][nz-1];
    //G1(szz_fv)
    szz_fv[ix][0] -= idz*ro1[ix][nz-1]*(G[1]*szz_fv[ix][1]+ G[2]*szz_fv[ix][2]+ G[3]*szz_fv[ix][3]+ G[4]*szz_fv[ix][4]);

    // divide by (rho1[0,0]*Gz[0,0]-rho2[-1,-1]*Gz1[-1,-1])

    szz_fv[ix][0] /= idz*(ro1[ix][nz-1]*G[0] - ro2[ix][0]*G[19]);

    p_fv[ix][nz-1] = -szz_fv[ix][0];

    sxz_fv[ix][0] = 0.0;

  }

}

/*-----------------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------------*/

/*****************************************
/ Free surface mimetic point updates
*****************************************/

/* acoustic free surface */
s




/* elastic anisotropic free surface */



/*-----------------------------------------------------------------------------------------------*/

/* Wavefield manipulation */

void extract_wfld_fsg(float** p, float** p_ff, float** p_vv, int nz, int nx)
/* extract wavefield from [f,f] and [v,v] */
{
    int ix, iz;

    for (ix=0; ix<nx; ix++){
        for (iz=0; iz<nz; iz++){
            p[ix][iz] = p_vv[ix][iz]*0.5 + (p_ff[ix][iz] + p_ff[ix][iz+1] + p_ff[ix+1][iz] + p_ff[ix+1][iz+1])*0.125;
        }
      }
  }


void combine_wfld(float **wout, float** win1, float** win2, int nz1, int nz2, int nx)
 {
   int ix,iz,nzf;
   nzf = (nz1>nz2) ? nz1 : nz2;
   for (ix=0; ix<nx; ix++){
       for (iz=0; iz<nzf; iz++){
           if(iz<nz1) wout[ix][iz]     = win1[ix][iz];
           if(iz<nz2) wout[ix][iz+nz1] = win2[ix][iz];
       }
     }
 }


/* cut a rectangular wavefield subset considering the free surface */
void cut2dfree(float** in, float** out, fdm2d fdm, sf_axis cz, sf_axis cx){
    int iz,ix;
    int fz,fx;

    fz = (floor)((sf_o(cz)-fdm->ozpad)/fdm->dz);
    fx = (floor)((sf_o(cx)-fdm->oxpad)/fdm->dx);
    //sf_warning("fx: %i fz: %i",fx,fz);
    for (ix=0;ix<fdm->nx;ix++) {
	    for (iz=0;iz<fdm->nz;iz++) {
    	    out[ix][iz] = in[fx+ix][fz+iz];
    	}
    }
}
