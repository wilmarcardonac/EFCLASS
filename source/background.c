/** @file background.c Documented background module
 *
 * Julien Lesgourgues, 17.04.2011
 * routines related to ncdm written by T. Tram in 2011
 *
 * Deals with the cosmological background evolution.
 * This module has two purposes:
 *
 * - at the beginning, to initialize the background, i.e. to integrate
 *    the background equations, and store all background quantities
 *    as a function of conformal time inside an interpolation table.
 *
 * - to provide routines which allow other modules to evaluate any
 *    background quantity for a given value of the conformal time (by
 *    interpolating within the interpolation table), or to find the
 *    correspondance between redhsift and conformal time.
 *
 *
 * The overall logic in this module is the following:
 *
 * 1. most background parameters that we will call {A}
 * (e.g. rho_gamma, ..) can be expressed as simple analytical
 * functions of a few variables that we will call {B} (in simplest
 * models, of the scale factor 'a'; in extended cosmologies, of 'a'
 * plus e.g. (phi, phidot) for quintessence, or some temperature for
 * exotic particles, etc...).
 *
 * 2. in turn, quantitites {B} can be found as a function of conformal
 * time by integrating the background equations.
 *
 * 3. some other quantitites that we will call {C} (like e.g. the
 * sound horizon or proper time) also require an integration with
 * respect to time, that cannot be infered analytically from
 * parameters {B}.
 *
 * So, we define the following routines:
 *
 * - background_functions() returns all background
 *    quantitites {A} as a function of quantitites {B}.
 *
 * - background_solve() integrates the quantities {B} and {C} with
 *    respect to conformal time; this integration requires many calls
 *    to background_functions().
 *
 * - the result is stored in the form of a big table in the background
 *    structure. There is one column for conformal time 'tau'; one or
 *    more for quantitites {B}; then several columns for quantities {A}
 *    and {C}.
 *
 * Later in the code, if we know the variables {B} and need some
 * quantity {A}, the quickest and most procise way is to call directly
 * background_functions() (for instance, in simple models, if we want
 * H at a given value of the scale factor). If we know 'tau' and want
 * any other qunatity, we can call background_at_tau(), which
 * interpolates in the table and returns all values. Finally it can be
 * useful to get 'tau' for a given redshift 'z': this can be done with
 * background_tau_of_z(). So if we are somewhere in the code, knowing
 * z and willing to get background quantitites, we should call first
 * background_tau_of_z() and then background_at_tau().
 *
 *
 * In order to save time, background_at_tau() can be called in three
 * modes: short_info, normal_info, long_info (returning only essential
 * quantities, or useful quantitites, or rarely useful
 * quantities). Each line in the interpolation table is a vector which
 * first few elements correspond to the short_info format; a larger
 * fraction contribute to the normal format; and the full vector
 * corresponds to the long format. The guideline is that short_info
 * returns only geometric quantitites like a, H, H'; normal format
 * returns quantities strictly needed at each step in the integration
 * of perturbations; long_info returns quantitites needed only
 * occasionally.
 *
 * In summary, the following functions can be called from other modules:
 *
 * -# background_init() at the beginning
 * -# background_at_tau(), background_tau_of_z() at any later time
 * -# background_free() at the end, when no more calls to the previous functions are needed
 */

#include "background.h"
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_bessel.h>

/**
 * Background quantities at given conformal time tau.
 *
 * Evaluates all background quantities at a given value of
 * conformal time by reading the pre-computed table ant interpolating.
 *
 * @param pba           Input: pointer to background structure (containing pre-computed table)
 * @param tau           Input: value of conformal time
 * @param return_format Input: format of output vector (short, normal, long)
 * @param intermode     Input: interpolation mode (normal or closeby)
 * @param last_index    Input/Ouput: index of the previous/current point in the interpolation array (input only for closeby mode, output for both)
 * @param pvecback      Output: vector (assumed to be already allocated)
 * @return the error status
 */

int background_at_tau(
                      struct background *pba,
                      double tau,
                      short return_format,
                      short intermode,
                      int * last_index,
                      double * pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size comptible with return_format) */
                      ) {

  /** Summary: */

  /** - define local variables */

  /* size of output vector, controlled by input parameter return_format */
  int pvecback_size;

  /** - check that tau is in the pre-computed range */

  class_test(tau < pba->tau_table[0],
             pba->error_message,
             "out of range: tau=%e < tau_min=%e, you should decrease the precision parameter a_ini_over_a_today_default\n",tau,pba->tau_table[0]);

  class_test(tau > pba->tau_table[pba->bt_size-1],
             pba->error_message,
             "out of range: tau=%e > tau_max=%e\n",tau,pba->tau_table[pba->bt_size-1]);

  /** - deduce length of returned vector from format mode */

  if (return_format == pba->normal_info) {
    pvecback_size=pba->bg_size_normal;
  }
  else {
    if (return_format == pba->short_info) {
      pvecback_size=pba->bg_size_short;
    }
    else {
      pvecback_size=pba->bg_size;
    }
  }

  /** - interpolate from pre-computed table with array_interpolate()
      or array_interpolate_growing_closeby() (depending on
      interpolation mode) */

  if (intermode == pba->inter_normal) {
    class_call(array_interpolate_spline(
                                        pba->tau_table,
                                        pba->bt_size,
                                        pba->background_table,
                                        pba->d2background_dtau2_table,
                                        pba->bg_size,
                                        tau,
                                        last_index,
                                        pvecback,
                                        pvecback_size,
                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }
  if (intermode == pba->inter_closeby) {
    class_call(array_interpolate_spline_growing_closeby(
                                                        pba->tau_table,
                                                        pba->bt_size,
                                                        pba->background_table,
                                                        pba->d2background_dtau2_table,
                                                        pba->bg_size,
                                                        tau,
                                                        last_index,
                                                        pvecback,
                                                        pvecback_size,
                                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }

  return _SUCCESS_;
}

/**
 * Conformal time at given redhsift.
 *
 * Returns tau(z) by interpolation from pre-computed table.
 *
 * @param pba Input: pointer to background structure
 * @param z   Input: redshift
 * @param tau Output: conformal time
 * @return the error status
 */

int background_tau_of_z(
                        struct background *pba,
                        double z,
                        double * tau
                        ) {

  /** Summary: */

  /** - define local variables */

  /* necessary for calling array_interpolate(), but never used */
  int last_index;

  /** - check that \f$ z \f$ is in the pre-computed range */
  class_test(z < pba->z_table[pba->bt_size-1],
             pba->error_message,
             "out of range: z=%e < z_min=%e\n",z,pba->z_table[pba->bt_size-1]);

  class_test(z > pba->z_table[0],
             pba->error_message,
             "out of range: a=%e > a_max=%e\n",z,pba->z_table[0]);

  /** - interpolate from pre-computed table with array_interpolate() */
  class_call(array_interpolate_spline(
                                      pba->z_table,
                                      pba->bt_size,
                                      pba->tau_table,
                                      pba->d2tau_dz2_table,
                                      1,
                                      z,
                                      &last_index,
                                      tau,
                                      1,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  return _SUCCESS_;
}

/**
 * Initialize the background structure, and in particular the
 * background interpolation table.
 *
 * @param ppr Input : pointer to precision structure
 * @param pba Input/Output : pointer to initialized background structure
 * @return the error status
 */

int background_init(
                    struct precision * ppr,
                    struct background * pba
                    ) {

  /** Summary: */

  /** - local variables : */
  int n_ncdm;
  double rho_ncdm_rel,rho_nu_rel;
  int filenum=0;

  /** - in verbose mode, provide some information */
  if (pba->background_verbose > 0) {
    printf("Running CLASS version %s\n",_VERSION_);
    printf("Computing background\n");

    /* below we want to inform the user about ncdm species*/
    if (pba->N_ncdm > 0) {

      /* loop over ncdm species */
      for (n_ncdm=0;n_ncdm<pba->N_ncdm; n_ncdm++) {

        /* inform if p-s-d read in files */
        if (pba->got_files[n_ncdm] == _TRUE_) {
          printf(" -> ncdm species i=%d read from file %s\n",n_ncdm+1,pba->ncdm_psd_files+filenum*_ARGUMENT_LENGTH_MAX_);
          filenum++;
        }

        /* call this function to get rho_ncdm */
        background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm],
                                pba->w_ncdm_bg[n_ncdm],
                                pba->q_size_ncdm_bg[n_ncdm],
                                0.,
                                pba->factor_ncdm[n_ncdm],
                                0.,
                                NULL,
                                &rho_ncdm_rel,
                                NULL,
                                NULL,
                                NULL);

        /* inform user of the contribution of each species to
           radiation density (in relativistic limit): should be
           between 1.01 and 1.02 for each active neutrino species;
           evaluated as rho_ncdm/rho_nu_rel where rho_nu_rel is the
           density of one neutrino in the instantaneous decoupling
           limit, i.e. assuming T_nu=(4/11)^1/3 T_gamma (this comes
           from the definition of N_eff) */
        rho_nu_rel = 56.0/45.0*pow(_PI_,6)*pow(4.0/11.0,4.0/3.0)*_G_/pow(_h_P_,3)/pow(_c_,7)*
          pow(_Mpc_over_m_,2)*pow(pba->T_cmb*_k_B_,4);

        printf(" -> ncdm species i=%d sampled with %d (resp. %d) points for purpose of background (resp. perturbation) integration. In the relativistic limit it gives N_eff = %g\n",
               n_ncdm+1,
               pba->q_size_ncdm_bg[n_ncdm],
               pba->q_size_ncdm[n_ncdm],
               rho_ncdm_rel/rho_nu_rel);
      }
    }
  }

  /** - assign values to all indices in vectors of background quantities with background_indices()*/
  class_call(background_indices(pba),
             pba->error_message,
             pba->error_message);

  /** - control that cosmological parameter values make sense */

  /* H0 in Mpc^{-1} */
  class_test((pba->H0 < _H0_SMALL_)||(pba->H0 > _H0_BIG_),
             pba->error_message,
             "H0=%g out of bounds (%g<H0<%g) \n",pba->H0,_H0_SMALL_,_H0_BIG_);

  class_test(fabs(pba->h * 1.e5 / _c_  / pba->H0 -1.)>ppr->smallest_allowed_variation,
             pba->error_message,
             "inconsistency between Hubble and reduced Hubble parameters: you have H0=%f/Mpc=%fkm/s/Mpc, but h=%f",pba->H0,pba->H0/1.e5* _c_,pba->h);

  /* T_cmb in K */
  class_test((pba->T_cmb < _TCMB_SMALL_)||(pba->T_cmb > _TCMB_BIG_),
             pba->error_message,
             "T_cmb=%g out of bounds (%g<T_cmb<%g)",pba->T_cmb,_TCMB_SMALL_,_TCMB_BIG_);

  /* H0 in Mpc^{-1} */
  class_test((pba->Omega0_k < _OMEGAK_SMALL_)||(pba->Omega0_k > _OMEGAK_BIG_),
             pba->error_message,
             "Omegak = %g out of bounds (%g<Omegak<%g) \n",pba->Omega0_k,_OMEGAK_SMALL_,_OMEGAK_BIG_);

  /* fluid equation of state */
  if (pba->has_fld == _TRUE_) {
    class_test(pba->w0_fld+pba->wa_fld>=1./3.,
               pba->error_message,
               "Your choice for w0_fld+wa_fld=%g is suspicious, there would not be radiation domination at early times\n",
               pba->w0_fld+pba->wa_fld);
  }

  /* in verbose mode, inform the user about the value of the ncdm
     masses in eV and about the ratio [m/omega_ncdm] in eV (the usual
     93 point something)*/
  if ((pba->background_verbose > 0) && (pba->has_ncdm == _TRUE_)) {
    for (n_ncdm=0; n_ncdm < pba->N_ncdm; n_ncdm++) {
      printf(" -> non-cold dark matter species with i=%d has m_i = %e eV (so m_i / omega_i =%e eV)\n",
             n_ncdm+1,
             pba->m_ncdm_in_eV[n_ncdm],
             pba->m_ncdm_in_eV[n_ncdm]/pba->Omega0_ncdm[n_ncdm]/pba->h/pba->h);
    }
  }

  /* check other quantities which would lead to segmentation fault if zero */
  class_test(pba->a_today <= 0,
             pba->error_message,
             "input a_today = %e instead of strictly positive",pba->a_today);

  class_test(_Gyr_over_Mpc_ <= 0,
             pba->error_message,
             "_Gyr_over_Mpc = %e instead of strictly positive",_Gyr_over_Mpc_);

  /** - this function integrates the background over time, allocates
      and fills the background table */
  class_call(background_solve(ppr,pba),
             pba->error_message,
             pba->error_message);

  return _SUCCESS_;

}

/**
 * Free all memory space allocated by background_init().
 *
 *
 * @param pba Input : pointer to background structure (to be freed)
 * @return the error status
 */

int background_free(
                    struct background *pba
                    ) {
  int k;

  free(pba->tau_table);
  free(pba->z_table);
  free(pba->d2tau_dz2_table);
  free(pba->background_table);
  free(pba->d2background_dtau2_table);

  if (pba->has_ncdm == _TRUE_) {
    for(k=0; k<pba->N_ncdm; k++){
      free(pba->q_ncdm[k]);
      free(pba->w_ncdm[k]);
      free(pba->q_ncdm_bg[k]);
      free(pba->w_ncdm_bg[k]);
      free(pba->dlnf0_dlnq_ncdm[k]);
    }
    free(pba->q_ncdm);
    free(pba->w_ncdm);
    free(pba->q_ncdm_bg);
    free(pba->w_ncdm_bg);
    free(pba->dlnf0_dlnq_ncdm);
    free(pba->q_size_ncdm);
    free(pba->q_size_ncdm_bg);
    free(pba->M_ncdm);
    free(pba->T_ncdm);
    free(pba->ksi_ncdm);
    free(pba->deg_ncdm);
    free(pba->Omega0_ncdm);
    free(pba->m_ncdm_in_eV);
    free(pba->factor_ncdm);
    if(pba->got_files!=NULL) free(pba->got_files);
    if(pba->ncdm_psd_files!=NULL)  free(pba->ncdm_psd_files);
    if(pba->ncdm_psd_parameters!=NULL)  free(pba->ncdm_psd_parameters);
  }

  return _SUCCESS_;
}

/**
 * Assign value to each relevant index in vectors of background quantities.
 *
 * @param pba Input : pointer to background structure
 * @return the error status
 */

int background_indices(
                       struct background *pba
                       ) {

  /** Summary: */

  /** - define local variables */

  /* a running index for the vector of background quantities */
  int index_bg;
  /* a running index for the vector of background quantities to be integrated */
  int index_bi;

  /** - initialize all flags: which species are present? */

  pba->has_cdm = _FALSE_;
  pba->has_ncdm = _FALSE_;
  pba->has_lambda = _FALSE_;
  pba->has_fld = _FALSE_;
  pba->has_ur = _FALSE_;
  pba->has_curvature = _FALSE_;
  pba->has_fR = _FALSE_;

  if (pba->Omega0_cdm != 0.)
    {
      pba->has_cdm = _TRUE_;
      printf("Omega0_cdm IS = %10e\n",pba->Omega0_cdm);
    }

  if (pba->Omega0_ncdm_tot != 0.)
    pba->has_ncdm = _TRUE_;

  if (pba->Omega0_lambda != 0.)
    {
      pba->has_lambda = _TRUE_;
      printf("Omega0_lambda IS = %.5e\n",pba->Omega0_lambda);
      if ( (pba->b_pi != 0.) && (pba->bhs != 0.) )
	{
	  printf("ONE CAN ONLY GIVE b_pi OR bhs \n");
	  exit(1);
	}
      else
	{
	  if (pba->b_pi != 0.)
	    {
	      if ( (pba->c0_des != 0.) || (pba->j0_des != 0.) )
		{
		  pba->has_fR = _TRUE_;
		  printf("RUNNING HORNDESKI DESIGNER MODEL WITH w=-1 \n");
		}
	      else
		{
		  pba->has_fR = _TRUE_;
		  printf("RUNNING f(R) DESIGNER MODEL WITH w=-1 \n");
		}
	    }
	  if (pba->bhs != 0.)
	    {
	      pba->has_fR = _TRUE_;
	      printf("RUNNING HU-SAWICKI MODEL WITH SAVVAS PARAMETRISATION \n");
	    }
	}
    }
  if (pba->Omega0_fld != 0.)
    {
      pba->has_fld = _TRUE_;
      printf("Omega0_fld IS = %.5e\n",pba->Omega0_fld);
    }
  if (pba->Omega0_ur != 0.)
    pba->has_ur = _TRUE_;

  if (pba->sgnK != 0)
    pba->has_curvature = _TRUE_;

  /** - intialization of all indices */

  index_bg=0;

  /* index for scale factor */
  class_define_index(pba->index_bg_a,_TRUE_,index_bg,1);

  /* - indices for H and its conformal-time-derivative */
  class_define_index(pba->index_bg_H,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_H_prime,_TRUE_,index_bg,1);

  /* - end of indices in the short vector of background values */
  pba->bg_size_short = index_bg;

  /* - index for rho_g (photon density) */
  class_define_index(pba->index_bg_rho_g,_TRUE_,index_bg,1);

  /* - index for rho_b (baryon density) */
  class_define_index(pba->index_bg_rho_b,_TRUE_,index_bg,1);

  /* - index for rho_cdm */
  class_define_index(pba->index_bg_rho_cdm,pba->has_cdm,index_bg,1);

  /* - indices for ncdm. We only define the indices for ncdm1
     (density, pressure, pseudo-pressure), the other ncdm indices
     are contiguous */
  class_define_index(pba->index_bg_rho_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);
  class_define_index(pba->index_bg_p_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);
  class_define_index(pba->index_bg_pseudo_p_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);

  /* - index for Lambda */
  class_define_index(pba->index_bg_rho_lambda,pba->has_lambda,index_bg,1);

  /* - index for fluid */
  //  class_define_index(pba->index_bg_rho_fld,pba->has_fld,index_bg,1);
  class_define_index(pba->index_bg_w_fR,pba->has_fR,index_bg,1);
  //class_define_index(pba->index_bg_w_prime_fR,pba->has_fR,index_bg,1);
  //class_define_index(pba->index_bg_conformal_H_fR,pba->has_fR,index_bg,1);
  //class_define_index(pba->index_bg_derivative_conformal_H_fR,pba->has_fR,index_bg,1);
  //class_define_index(pba->index_bg_f_fR,pba->has_fR,index_bg,1);
  class_define_index(pba->index_bg_F_fR,pba->has_fR,index_bg,1);
  class_define_index(pba->index_bg_Fprime_fR,pba->has_fR,index_bg,1);
  //class_define_index(pba->index_bg_Fdoubleprime_fR,pba->has_fR,index_bg,1);
  class_define_index(pba->index_bg_FR_fR,pba->has_fR,index_bg,1);
  //class_define_index(pba->index_bg_FRprime_fR,pba->has_fR,index_bg,1);
  //class_define_index(pba->index_bg_FRdoubleprime_fR,pba->has_fR,index_bg,1);

  /* - index for ultra-relativistic neutrinos/species */
  class_define_index(pba->index_bg_rho_ur,pba->has_ur,index_bg,1);

  /* - index for Omega_r (relativistic density fraction) */
  class_define_index(pba->index_bg_Omega_r,_TRUE_,index_bg,1);

  /* - put here additional ingredients that you want to appear in the
     normal vector */
  /*    */
  /*    */

  /* - end of indices in the normal vector of background values */
  pba->bg_size_normal = index_bg;

  /* - indices in the long version : */

  /* -> critical density */
  class_define_index(pba->index_bg_rho_crit,_TRUE_,index_bg,1);

  /* - index for Omega_m (non-relativistic density fraction) */
  class_define_index(pba->index_bg_Omega_m,_TRUE_,index_bg,1);

  /* -> conformal distance */
  class_define_index(pba->index_bg_conf_distance,_TRUE_,index_bg,1);

  /* -> angular diameter distance */
  class_define_index(pba->index_bg_ang_distance,_TRUE_,index_bg,1);

  /* -> luminosity distance */
  class_define_index(pba->index_bg_lum_distance,_TRUE_,index_bg,1);

  /* -> proper time (for age of the Universe) */
  class_define_index(pba->index_bg_time,_TRUE_,index_bg,1);

  /* -> conformal sound horizon */
  class_define_index(pba->index_bg_rs,_TRUE_,index_bg,1);

  /* -> density growth factor in dust universe */
  class_define_index(pba->index_bg_D,_TRUE_,index_bg,1);

  /* -> velocity growth factor in dust universe */
  class_define_index(pba->index_bg_f,_TRUE_,index_bg,1);

  /* -> put here additional quantities describing background */
  /*    */
  /*    */

  /* -> end of indices in the long vector of background values */
  pba->bg_size = index_bg;

  /* - now, indices in vector of variables to integrate */

  index_bi=0;

  /* -> scale factor */
  class_define_index(pba->index_bi_a,_TRUE_,index_bi,1);

  /* -> proper time (for age of the Universe) */
  class_define_index(pba->index_bi_time,_TRUE_,index_bi,1);

  /* -> sound horizon */
  class_define_index(pba->index_bi_rs,_TRUE_,index_bi,1);

  /* -> integral for growth factor */
  class_define_index(pba->index_bi_growth,_TRUE_,index_bi,1);

  /* -> index for conformal time in vector of variables to integrate */
  class_define_index(pba->index_bi_tau,_TRUE_,index_bi,1);

  /* -> end of indices in the vector of variables to integrate */
  pba->bi_size = index_bi;

  /* index_bi_tau must be the last index, because tau is part of this vector for the purpose of being stored, */
  /* but it is not a quantity to be integrated (since integration is over tau itself) */
  class_test(pba->index_bi_tau != index_bi-1,
             pba->error_message,
             "background integration requires index_bi_tau to be the last of all index_bi's");

  /* flags for calling the interpolation routine */

  pba->short_info=0;
  pba->normal_info=1;
  pba->long_info=2;

  pba->inter_normal=0;
  pba->inter_closeby=1;

  return _SUCCESS_;

}

/**
 * Background quantities at given a.
 *
 * Function evaluating all background quantities which can be computed
 * analytically as a function of parameters like the scale factor 'a'
 * (see discussion at the beginnign of this file). In extended
 * comsological models, we would include other input parameters than
 * just 'a', e.g. (phi, phidot) for quintessence, some temperature of
 * exotic relics, etc...
 *
 * @param pba           Input: pointer to background structure
 * @param a             Input: value of scale factor
 * @param return_format Input: format of output vector
 * @param pvecback      Output: vector of background quantities (assmued to be already allocated)
 * @return the error status
 */

int background_functions(
                         struct background *pba,
                         double a, /* in extended models there could be more than one argument: phi, phidot of quintessence; temperature of some particles; etc. */
                         short return_format,
                         double * pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size comptible with return_format) */
                         ) {

  /** Summary: */

  /** - define local variables */

  /* total density */
  double rho_tot;
  /* total pressure */
  double p_tot;
  /* total relativistic density */
  double rho_r;
  /* total non-relativistic density */
  double rho_m;
  /* scale factor relative to scale factor today */
  double a_rel;
  /* background ncdm quantities */
  double rho_ncdm,p_ncdm,pseudo_p_ncdm;
  /* index for n_ncdm species */
  int n_ncdm;
  /* MATTER PARAMETER DENSITY */
  double om0,or0;
  /* COSMOLOGICAL CONSTANT */
  //double Lambda;
  /* CONSTANTS f(R) */
  double des_c0;//,a2,b,c,a3;
  //double R0,x0,y0,alpha_fR;
  //  double cH,cH_prime,cH0,cH0_prime,cH_twoprime,cH_threeprime,R_prime,R_twoprime,F_prime,F_doubleprime;
  //double F_prime;
  double fR0HS,HSb;
  /* VARIABLES f(R) */
  //double x,y,R;


  /** - initialize local variables */
  rho_tot = 0.;
  p_tot = 0.;
  rho_r=0.;
  rho_m=0.;
  a_rel = a / pba->a_today;
  om0 = pba->Omega0_b+pba->Omega0_cdm;
  or0 = pba->Omega0_g;
  //Lambda = 3.*pow(pba->H0,2)*(1. - Omega0_M);
  des_c0 = (-7. + sqrt(73.))/12.;
  fR0HS = pba->bhs;//-0.1 ;//59285033891483;
  HSb = -(fR0HS/(pow(-1+om0+or0,2)/pow(-4+3*om0+4*or0,2)+sqrt((pow(-1+om0+or0,3)*((-1+om0+or0)*pow(-4+3*om0+4*or0,3)-4*fR0HS*(135*pow(om0,4)+128*pow(-1+or0,4)+96*om0*pow(-1+or0,2)*(-5+6*or0)+36*pow(om0,3)*(-11+16*or0)+12*pow(om0,2)*(-1+or0)*(-51+74*or0))))/pow(-4+3*om0+4*or0,7))));

  //-(pba->bhs*pow(4.-3.*om0,2.))/(2.*pow(1.-om0,2.));

  //a2 = 1. + b2;
  //a3 = 2. + b2;
  //b = 1.5 + b2; 
  //c = 13./6. + 2.*b2;


  class_test(a_rel <= 0.,
             pba->error_message,
             "a = %e instead of strictly positive",a_rel);

  /** - pass value of a to output */
  pvecback[pba->index_bg_a] = a;

  /** - compute each component's density and pressure */

  /* photons */
  pvecback[pba->index_bg_rho_g] = pba->Omega0_g * pow(pba->H0,2) / pow(a_rel,4);
  rho_tot += pvecback[pba->index_bg_rho_g];
  p_tot += (1./3.) * pvecback[pba->index_bg_rho_g];
  rho_r += pvecback[pba->index_bg_rho_g];

  /* baryons */
  pvecback[pba->index_bg_rho_b] = pba->Omega0_b * pow(pba->H0,2) / pow(a_rel,3);
  rho_tot += pvecback[pba->index_bg_rho_b];
  p_tot += 0;
  rho_m += pvecback[pba->index_bg_rho_b];

  /* cdm */
  if (pba->has_cdm == _TRUE_) {
    pvecback[pba->index_bg_rho_cdm] = pba->Omega0_cdm * pow(pba->H0,2) / pow(a_rel,3);
    rho_tot += pvecback[pba->index_bg_rho_cdm];
    p_tot += 0.;
    rho_m += pvecback[pba->index_bg_rho_cdm];
  }

  /* ncdm */
  if (pba->has_ncdm == _TRUE_) {

    /* Loop over species: */
    for(n_ncdm=0; n_ncdm<pba->N_ncdm; n_ncdm++){

      /* function returning background ncdm[n_ncdm] quantities (only
         those for which non-NULL pointers are passed) */
      class_call(background_ncdm_momenta(
                                         pba->q_ncdm_bg[n_ncdm],
                                         pba->w_ncdm_bg[n_ncdm],
                                         pba->q_size_ncdm_bg[n_ncdm],
                                         pba->M_ncdm[n_ncdm],
                                         pba->factor_ncdm[n_ncdm],
                                         1./a_rel-1.,
                                         NULL,
                                         &rho_ncdm,
                                         &p_ncdm,
                                         NULL,
                                         &pseudo_p_ncdm),
                 pba->error_message,
                 pba->error_message);

      pvecback[pba->index_bg_rho_ncdm1+n_ncdm] = rho_ncdm;
      rho_tot += rho_ncdm;
      pvecback[pba->index_bg_p_ncdm1+n_ncdm] = p_ncdm;
      p_tot += p_ncdm;
      pvecback[pba->index_bg_pseudo_p_ncdm1+n_ncdm] = pseudo_p_ncdm;

      /* (3 p_ncdm1) is the "relativistic" contrinution to rho_ncdm1 */
      rho_r += 3.* p_ncdm;

      /* (rho_ncdm1 - 3 p_ncdm1) is the "non-relativistic" contribution
         to rho_ncdm1 */
      rho_m += rho_ncdm - 3.* p_ncdm;
    }
  }

  /* Lambda */
  if (pba->has_lambda == _TRUE_)
    {
      if (pba->has_fR == _TRUE_)
	{
	  if (pba->bhs != 0.)
	    {
	      /* w(a) */
	      pvecback[pba->index_bg_w_fR] = -1-(4*pow(a,2)*HSb*(-1+om0+or0)*(-9*a*pow(om0,3)-7*pow(om0,2)*or0-63*pow(a,4)*pow(om0,2)*(-1+om0+or0)-88*pow(a,3)*om0*or0*(-1+om0+or0)+72*pow(a,7)*om0*pow(-1+om0+or0,2)+32*pow(a,6)*or0*pow(-1+om0+or0,2)))/(3.*pow(om0-4*pow(a,3)*(-1+om0+or0),4))-(2*pow(a,4)*pow(HSb,2)*(1-om0-or0)*(-1+om0+or0)*(327*pow(a,2)*pow(om0,7)+520*a*pow(om0,6)*or0+196*pow(om0,5)*pow(or0,2)-19944*pow(a,5)*pow(om0,6)*(-1+om0+or0)-32960*pow(a,4)*pow(om0,5)*or0*(-1+om0+or0)-13328*pow(a,3)*pow(om0,4)*pow(or0,2)*(-1+om0+or0)-195084*pow(a,8)*pow(om0,5)*pow(-1+om0+or0,2)-447296*pow(a,7)*pow(om0,4)*or0*pow(-1+om0+or0,2)-244352*pow(a,6)*pow(om0,3)*pow(or0,2)*pow(-1+om0+or0,2)+55344*pow(a,11)*pow(om0,4)*pow(-1+om0+or0,3)-417664*pow(a,10)*pow(om0,3)*or0*pow(-1+om0+or0,3)-503296*pow(a,9)*pow(om0,2)*pow(or0,2)*pow(-1+om0+or0,3)+225792*pow(a,14)*pow(om0,3)*pow(-1+om0+or0,4)+363008*pow(a,13)*pow(om0,2)*or0*pow(-1+om0+or0,4)-243712*pow(a,12)*om0*pow(or0,2)*pow(-1+om0+or0,4)-239616*pow(a,17)*pow(om0,2)*pow(-1+om0+or0,5)-53248*pow(a,16)*om0*or0*pow(-1+om0+or0,5)+8192*pow(a,15)*pow(or0,2)*pow(-1+om0+or0,5)+122880*pow(a,20)*om0*pow(-1+om0+or0,6)+49152*pow(a,19)*or0*pow(-1+om0+or0,6)))/(3.*pow(om0-4*pow(a,3)*(-1+om0+or0),9));

	      /*-1.-(12.*pow(a,3.)*pow(pba->a_today,3.)*HSb*(-1. + om0)*om0*(pow(a,3.)*(-1. + om0) - pow(pba->a_today,3.)*om0)*(8.*pow(a,3.)*(-1. + om0) + pow(pba->a_today,3.)*om0))/pow(-4.*pow(a,3.)*(-1. + om0) + pow(pba->a_today,3.)*om0,4.) + (2.*pow(a,6.)*pow(pba->a_today,3.)*pow(HSb,2.)*pow(-1. + om0,2.)*om0*(-40960.*pow(a,18.)*pow(-1. + om0,6.) +79872.*pow(a,15.)*pow(pba->a_today,3.)*pow(-1. + om0,5.)*om0 -75264.*pow(a,12.)*pow(pba->a_today,6.)*pow(-1. + om0,4.)*pow(om0,2.) -18448.*pow(a,9.)*pow(pba->a_today,9.)*pow(-1. + om0,3.)*pow(om0,3.) +65028*pow(a,6.)*pow(pba->a_today,12.)*pow(-1. + om0,2.)*pow(om0,4.) +6648.*pow(a,3.)*pow(pba->a_today,15.)*(-1. + om0)*pow(om0,5.) - 109.*pow(pba->a_today,1.)*pow(om0,6.)))/pow(4.*pow(a,3.)*(-1. + om0) -pow(pba->a_today,3.)*om0,9.); */

	      //-1 - 12*HSb*(-1 + om0)*om0*pow(a,3)*(-om0 + (-1 + om0)*pow(a,3))*(om0 + 8*(-1 + om0)*pow(a,3))*pow(om0 - 4*(-1 + om0)*pow(a,3),-4) - 
	      //2*om0*pow(a,6)*pow(HSb,2)*pow(-1 + om0,2)*(-79872*om0*pow(a,15)*pow(-1 + om0,5) + 40960*pow(a,18)*pow(-1 + om0,6) + 75264*pow(a,12)*pow(-1 + om0,4)*pow(om0,2) + 18448*pow(a,9)*pow(-1 + om0,3)*pow(om0,3) - 65028*pow(a,6)*pow(-1 + om0,2)*pow(om0,4) - 6648*(-1 + om0)*pow(a,3)*pow(om0,5) + 109*pow(om0,6))*pow(-om0 + 4*(-1 + om0)*pow(a,3),-9);
	  /* w_prime(a) */
	    /*	    pvecback[pba->index_bg_w_prime_fR] = 12*HSb*(-1 + om0)*om0*pow(a,2)*(18*HSb*pow(a,6)*pow(-1 + om0,2)*(-79872*om0*pow(a,15)*pow(-1 + om0,5) + 40960*pow(a,18)*pow(-1 + om0,6) + 75264*pow(a,12)*pow(-1 + om0,4)*pow(om0,2) + 18448*pow(a,9)*pow(-1 + om0,3)*pow(om0,3) - 65028*pow(a,6)*pow(-1 + om0,2)*pow(om0,4) - 6648*(-1 + om0)*pow(a,3)*pow(om0,5) + 109*pow(om0,6))*pow(om0 - 4*(-1 + om0)*pow(a,3),-10) - 48*(-1 + om0)*pow(a,3)*(-om0 + (-1 + om0)*pow(a,3))*(om0 + 8*(-1 + om0)*pow(a,3))*pow(om0 - 4*(-1 + om0)*pow(a,3),-5) - 24*(-1 + om0)*pow(a,3)*(-om0 + (-1 + om0)*pow(a,3))*pow(om0 - 4*(-1 + om0)*pow(a,3),-4) - 3*(-1 + om0)*pow(a,3)*(om0 + 8*(-1 + om0)*pow(a,3))*pow(om0 - 4*(-1 + om0)*pow(a,3),-4) - 3*(-om0 + (-1 + om0)*pow(a,3))*(om0 + 8*(-1 + om0)*pow(a,3))*pow(om0 - 4*(-1 + om0)*pow(a,3),-4) - 12*HSb*pow(a,6)*pow(-1 + om0,2)*(-16640*om0*pow(a,12)*pow(-1 + om0,4) + 10240*pow(a,15)*pow(-1 + om0,5) + 12544*pow(a,9)*pow(-1 + om0,3)*pow(om0,2) + 2306*pow(a,6)*pow(-1 + om0,2)*pow(om0,3) - 5419*(-1 + om0)*pow(a,3)*pow(om0,4) - 277*pow(om0,5))*pow(-om0 + 4*(-1 + om0)*pow(a,3),-9) - HSb*(-1 + om0)*pow(a,3)*(-79872*om0*pow(a,15)*pow(-1 + om0,5) + 40960*pow(a,18)*pow(-1 + om0,6) + 75264*pow(a,12)*pow(-1 + om0,4)*pow(om0,2) + 18448*pow(a,9)*pow(-1 + om0,3)*pow(om0,3) - 65028*pow(a,6)*pow(-1 + om0,2)*pow(om0,4) - 6648*(-1 + om0)*pow(a,3)*pow(om0,5) + 109*pow(om0,6))*pow(-om0 + 4*(-1 + om0)*pow(a,3),-9)); */

	      pvecback[pba->index_bg_rho_lambda] = ((1-om0-or0+(2*pow(a,2)*HSb*pow(-1+om0+or0,2)*(-6*a*pow(om0,2)-7*om0*or0+3*pow(a,4)*om0*(-1+om0+or0)+4*pow(a,3)*or0*(-1+om0+or0)+12*pow(a,7)*pow(-1+om0+or0,2)))/pow(-om0+4*pow(a,3)*(-1+om0+or0),3)-(pow(a,5)*pow(HSb,2)*pow(-1+om0+or0,3)*(-37*a*pow(om0,6)-40*pow(om0,5)*or0+4656*pow(a,4)*pow(om0,5)*(-1+om0+or0)+8692*pow(a,3)*pow(om0,4)*or0*(-1+om0+or0)+4032*pow(a,2)*pow(om0,3)*pow(or0,2)*(-1+om0+or0)+7452*pow(a,7)*pow(om0,4)*pow(-1+om0+or0,2)+25728*pow(a,6)*pow(om0,3)*or0*pow(-1+om0+or0,2)+17856*pow(a,5)*pow(om0,2)*pow(or0,2)*pow(-1+om0+or0,2)-25408*pow(a,10)*pow(om0,3)*pow(-1+om0+or0,3)-22016*pow(a,9)*pow(om0,2)*or0*pow(-1+om0+or0,3)+9216*pow(a,8)*om0*pow(or0,2)*pow(-1+om0+or0,3)+22848*pow(a,13)*pow(om0,2)*pow(-1+om0+or0,4)+2048*pow(a,12)*om0*or0*pow(-1+om0+or0,4)-9216*pow(a,16)*om0*pow(-1+om0+or0,5)-3072*pow(a,15)*or0*pow(-1+om0+or0,5)-1024*pow(a,19)*pow(-1+om0+or0,6)))/pow(om0-4*pow(a,3)*(-1+om0+or0),8))/(1-om0-or0))*pow(pba->H0,2)*pba->Omega0_lambda;

		/*exp((3.*HSb*om0*((64. + 21.*(-4. + om0)*om0)/pow(-4. + 3.*om0,3.) +(-64.*pow(a,6.)*pow(pba->a_today,3.)*pow(-1. + om0,2.) +44.*pow(a,3.)*pow(pba->a_today,6.)*(-1 + om0)*om0 - pow(pba->a_today,9.)*pow(om0,2.))/pow(4.*pow(a,3.)*(-1. + om0) - pow(pba->a_today,3.)*om0,3.)))/8. +(pow(pba->a_today,3.)*pow(HSb,2.)*om0*((-2621440.*pow(a,21.)*pow(-1. + om0,7.) + 4849664.*pow(a,18.)*pow(pba->a_today,3.)*pow(-1. + om0,6.)*om0 -4030464.*pow(a,15.)*pow(pba->a_today,6.)*pow(-1. + om0,5.)*pow(om0,2.) +964352.*pow(a,12.)*pow(pba->a_today,9.)*pow(-1. + om0,4.)*pow(om0,3.) +639488.*pow(a,9.)*pow(pba->a_today,12.)*pow(-1. + om0,3.)*pow(om0,4.) -9024.*pow(a,6.)*pow(pba->a_today,15.)*pow(-1.+om0,2.)*pow(om0,5.) -352.*pow(a,3.)*pow(pba->a_today,18.)*(-1. + om0)*pow(om0,6.) + 11.*pow(pba->a_today,21.)*pow(om0,7.))/pow(-4.*pow(a,3.)*(-1. + om0) + pow(pba->a_today,3.)*om0,8) +(-2621440. + om0*(13500416. + 3.*om0*(-9994240. + 3.*om0*(4243712. +9.*om0*(-377344. + 3.*om0*(61376. + om0*(-15008. + 855.*om0)))))))/(pow(pba->a_today,3.)*pow(4. - 3.*om0,8.))))/128.)*pow(pba->H0,2)*pba->Omega0_lambda; */

	      /*-3*(-1 + om0)*pow(pba->H0,2) + 3*pow(a,6)*pow(HSb,2)*pow(pba->H0,2)*pow(-1 + om0,3)*(9216*om0*pow(a,15)*pow(-1 + om0,5) + 1024*pow(a,18)*pow(-1 + om0,6) - 22848*pow(a,12)*pow(-1 + om0,4)*pow(om0,2) + 25408*pow(a,9)*pow(-1 + om0,3)*pow(om0,3) - 7452*pow(a,6)*pow(-1 + om0,2)*pow(om0,4) - 4656*(-1 + om0)*pow(a,3)*pow(om0,5) + 37*pow(om0,6))*pow(om0 - 4*(-1 + om0)*pow(a,3),-8) - 18*HSb*pow(a,3)*pow(pba->H0,2)*pow(-1 + om0,2)*((-1 + om0)*om0*pow(a,3) + 4*pow(a,6)*pow(-1 + om0,2) - 2*pow(om0,2))*pow(om0 - 4*(-1 + om0)*pow(a,3),-3);*/
	      rho_tot += pvecback[pba->index_bg_rho_lambda];
	      p_tot +=  pvecback[pba->index_bg_w_fR]*pvecback[pba->index_bg_rho_lambda];
	    }
	  if (pba->b_pi != 0.)
	    {
	      pvecback[pba->index_bg_rho_lambda] = pba->Omega0_lambda * pow(pba->H0,2);
	      rho_tot += pvecback[pba->index_bg_rho_lambda];
	      p_tot -= pvecback[pba->index_bg_rho_lambda];
	    }
	}
      else
	{
	  pvecback[pba->index_bg_rho_lambda] = pba->Omega0_lambda * pow(pba->H0,2);
	  rho_tot += pvecback[pba->index_bg_rho_lambda];
	  p_tot -= pvecback[pba->index_bg_rho_lambda];
	}
    }

  if (pba->has_fld == _TRUE_) 
    {
      /* fluid with w=w0+wa(1-a/a0) and constant cs2 */
      pvecback[pba->index_bg_rho_fld] = pba->Omega0_fld * pow(pba->H0,2)
	/ pow(a_rel,3.*(1.+pba->w0_fld+pba->wa_fld))
	* exp(3.*pba->wa_fld*(a_rel-1.));
      rho_tot += pvecback[pba->index_bg_rho_fld];
      p_tot += (pba->w0_fld+pba->wa_fld*(1.-a_rel)) * pvecback[pba->index_bg_rho_fld];
    }

  /* relativistic neutrinos (and all relativistic relics) */
  if (pba->has_ur == _TRUE_) {
    pvecback[pba->index_bg_rho_ur] = pba->Omega0_ur * pow(pba->H0,2) / pow(a_rel,4);
    rho_tot += pvecback[pba->index_bg_rho_ur];
    p_tot += (1./3.) * pvecback[pba->index_bg_rho_ur];
    rho_r += pvecback[pba->index_bg_rho_ur];
  }

  /** - compute expansion rate H from Friedmann equation: this is the
      unique place where the Friedmann equation is assumed. Remember
      that densities are all expressed in units of [3c^2/8piG], ie
      rho_class = [8 pi G rho_physical / 3 c^2] */
  pvecback[pba->index_bg_H] = sqrt(rho_tot-pba->K/a/a);

  /** - compute derivative of H with respect to conformal time */
  pvecback[pba->index_bg_H_prime] = - (3./2.) * (rho_tot + p_tot) * a + pba->K/a;

  if (pba->has_fR == _TRUE_) 
    {
      if (pba->b_pi != 0.)
	{
	  pvecback[pba->index_bg_Fprime_fR] = 0. ; 
	    /*((pow(a,2.+3.*des_c0)*(7.+6.*des_c0)*fR0HS*pow(pow(a,3.)*(1.-om0)+om0,-1.-des_c0)*gsl_sf_hyperg_2F1(des_c0,1.5+des_c0,13./6.+2.*des_c0,(pow(a,3.)*(1.-om0))/(pow(a,3.)*(1.-om0)+om0)))/(2.*gsl_sf_hyperg_2F1(1.+des_c0,1.5+des_c0,13./6.+2.*des_c0,1.-om0))+(pow(a,2.+3.*des_c0)*fR0HS*(2.*pow(a,3.)*(1.-om0)-om0)*pow(pow(a,3.)*(1.-om0)+om0,-2.-des_c0)*gsl_sf_hyperg_2F1(1.+des_c0,1.5+des_c0,13./6.+2.*des_c0,(pow(a,3.)*(1.-om0))/(pow(a,3.)*(1.-om0)+om0)))/(2.*gsl_sf_hyperg_2F1(1.+des_c0,1.5+des_c0,13./6+2.*des_c0,1.-om0)));*/
	  pvecback[pba->index_bg_F_fR] = 0. ; /*1 + (pow(a,3.)*fR0HS*pow((pow(a,3.)*(-1 + om0))/(pow(a,3.)*(-1. + om0) - om0),des_c0)*gsl_sf_hyperg_2F1(1. + des_c0,1.5 + des_c0,13./6. +2.*des_c0,(pow(a,3.)*(-1. + om0))/(pow(a,3.)*(-1. + om0) - om0)))/(pow(1. - om0,des_c0)*(-(pow(a,3.)*(-1. + om0)) + om0)*gsl_sf_hyperg_2F1(1. + des_c0,1.5 + des_c0,13./6. + 2.*des_c0,1. - om0)); */
	  pvecback[pba->index_bg_FR_fR] = 0. ;/*(pow(3.,2. + des_c0)*fR0HS*pow((pow(a,3.)*(-1. + om0))/(3.*pow(a,3.)*(-1. + om0) - 3.*om0),des_c0)*pow(-(pow(a,3.)*(-1. + om0)) + om0,2.)*((7. + 6.*des_c0)*(pow(a,3.)*(-1. + om0) - om0)*gsl_sf_hyperg_2F1(des_c0,1.5 + des_c0,13./6. + 2.*des_c0,(pow(a,3.)*(-1. + om0))/(pow(a,3.)*(-1. + om0) - om0)) +(2.*pow(a,3.)*(-1. + om0) + om0)*gsl_sf_hyperg_2F1(1. + des_c0,1.5 + des_c0,13./6. + 2.*des_c0, (pow(a,3.)*(-1. + om0))/(pow(a,3.)*(-1. + om0) - om0))))/(2.*pow(a,6.)*pow(1. - om0,des_c0)*om0*pow(3. + 3.*(-1. + pow(a,-3.))*om0,4.)*gsl_sf_hyperg_2F1(1. + des_c0,1.5 + des_c0,13./6. + 2.*des_c0,1. - om0)); */
	}
      else
	{
	  pvecback[pba->index_bg_FR_fR] = (-4*pow(a,9)*HSb*pow(-1+om0+or0,2))/(3.*pow(-om0+4*pow(a,3)*(-1+om0+or0),3))-(4*pow(a,12)*pow(HSb,2)*pow(-1+om0+or0,3)*(-pow(om0,4)+100*pow(a,3)*pow(om0,3)*(-1+om0+or0)+84*pow(a,2)*pow(om0,2)*or0*(-1+om0+or0)-60*pow(a,6)*pow(om0,2)*pow(-1+om0+or0,2)+96*pow(a,5)*om0*or0*pow(-1+om0+or0,2)-32*pow(a,9)*om0*pow(-1+om0+or0,3)+128*pow(a,12)*pow(-1+om0+or0,4)))/pow(om0-4*pow(a,3)*(-1+om0+or0),8);

	    /*(4.*pow(a,9.)*HSb*pow(-1. + om0,2.))/(3.*pow(-4.*pow(a,3.)*(-1. + om0) + om0,3.)) -(4.*pow(a,12.)*pow(HSb,2.)*pow(-1. + om0,3.)*(128.*pow(a,12.)*pow(-1. + om0,4.) -32.*pow(a,9.)*pow(-1. + om0,3.)*om0 -60.*pow(a,6.)*pow(-1. + om0,2.)*pow(om0,2.) +100.*pow(a,3.)*(-1. + om0)*pow(om0,3.) -pow(om0,4.)))/pow(-4.*pow(a,3.)*(-1. + om0) + om0,8.); */

	  pvecback[pba->index_bg_F_fR] = 1-(2*pow(a,6)*HSb*pow(-1+om0+or0,2))/pow(om0-4*pow(a,3)*(-1+om0+or0),2)-(4*pow(a,9)*pow(HSb,2)*pow(-1+om0+or0,3)*(-pow(om0,4)+100*pow(a,3)*pow(om0,3)*(-1+om0+or0)+84*pow(a,2)*pow(om0,2)*or0*(-1+om0+or0)-60*pow(a,6)*pow(om0,2)*pow(-1+om0+or0,2)+96*pow(a,5)*om0*or0*pow(-1+om0+or0,2)-32*pow(a,9)*om0*pow(-1+om0+or0,3)+128*pow(a,12)*pow(-1+om0+or0,4)))/pow(-om0+4*pow(a,3)*(-1+om0+or0),7);

	    /*1.-(2.*pow(a,6.)*HSb*pow(-1.+om0,2.))/pow(-4.*pow(a,3.)*(-1.+om0)+om0,2.) -(4.*pow(a,9.)*pow(HSb,2.)*pow(-1.+om0,3.)*(128.*pow(a,12.)*pow(-1.+om0,4.)-32.*pow(a,9.)*pow(-1. +om0,3.)*om0-60.*pow(a,6.)*pow(-1.+om0,2.)*pow(om0,2.)+100.*pow(a,3.)*(-1.+om0)*pow(om0,3.) -pow(om0,4.)))/pow(4.*pow(a,3.)*(-1.+om0)-om0,7.); */

	  pvecback[pba->index_bg_Fprime_fR] = (12*pow(a,5)*HSb*om0*pow(-1+om0+or0,2))/pow(-om0+4*pow(a,3)*(-1+om0+or0),3)+(12*pow(a,8)*pow(HSb,2)*om0*pow(-1+om0+or0,3)*(-3*pow(om0,4)+384*pow(a,3)*pow(om0,3)*(-1+om0+or0)+308*pow(a,2)*pow(om0,2)*or0*(-1+om0+or0)+900*pow(a,6)*pow(om0,2)*pow(-1+om0+or0,2)+1568*pow(a,5)*om0*or0*pow(-1+om0+or0,2)-672*pow(a,9)*om0*pow(-1+om0+or0,3)+896*pow(a,8)*or0*pow(-1+om0+or0,3)+768*pow(a,12)*pow(-1+om0+or0,4)))/pow(om0-4*pow(a,3)*(-1+om0+or0),8);

	    /*((12.*pow(a,5.)*HSb*pow(-1. + om0,2.)*om0)/pow(4.*pow(a,3.)*(-1. + om0) - om0,3.) +(36.*pow(a,8.)*pow(HSb,2.)*pow(-1. + om0,3.)*om0*(256.*pow(a,12.)*pow(-1. + om0,4.)- 24.*pow(a,9.)*pow(-1. + om0,3.)*om0 +300.*pow(a,6.)*pow(-1. + om0,2.)*pow(om0,2.) +128.*pow(a,3.)*(-1. + om0)*pow(om0,3.) -pow(om0,4.)))/pow(-4.*pow(a,3.)*(-1. + om0) + om0,8.));*/

	  //printf("FR = %.5e F = %.5e Fprime = %.5e HSb = %.5e bhs = %.5e\n",pvecback[pba->index_bg_FR_fR],pvecback[pba->index_bg_F_fR],pvecback[pba->index_bg_Fprime_fR],HSb,pba->bhs);
	  //exit(1);
	}
    }

  /** - compute relativistic density to total density ratio */
  pvecback[pba->index_bg_Omega_r] = rho_r / rho_tot;

  /** - compute other quantities in the exhaustive, redundent format: */
  if (return_format == pba->long_info) {

    /** - compute critical density */
    pvecback[pba->index_bg_rho_crit] = rho_tot-pba->K/a/a;
    class_test(pvecback[pba->index_bg_rho_crit] <= 0.,
               pba->error_message,
               "rho_crit = %e instead of strictly positive",pvecback[pba->index_bg_rho_crit]);

    /** - compute Omega_m */
    pvecback[pba->index_bg_Omega_m] = rho_m / rho_tot;

    /* one can put other variables here */
    /*  */
    /*  */

  }

  return _SUCCESS_;

}

/**
 * This is the routine where the distribution function f0(q) of each
 * ncdm species is specified (it is the only place to modify if you
 * need a partlar f0(q))
 *
 * @param pbadist Input:  structure containing all parameters defining f0(q)
 * @param q       Input:  momentum
 * @param f0      Output: phase-space distribution
 */

int background_ncdm_distribution(
                                 void * pbadist,
                                 double q,
                                 double * f0
                                 ) {
  struct background * pba;
  struct background_parameters_for_distributions * pbadist_local;
  int n_ncdm,lastidx;
  double ksi;
  double qlast,dqlast,f0last,df0last;
  double *param;
  /** Variables corresponing to entries in param: */
  //double square_s12,square_s23,square_s13;
  //double mixing_matrix[3][3];
  //int i;

  /** - extract from the input structure pbadist all the relevant information */
  pbadist_local = pbadist;          /* restore actual format of pbadist */
  pba = pbadist_local->pba;         /* extract the background structure from it */
  param = pba->ncdm_psd_parameters; /* extract the optional parameter list from it */
  n_ncdm = pbadist_local->n_ncdm;   /* extract index of ncdm species under consideration */
  ksi = pba->ksi_ncdm[n_ncdm];      /* extract chemical potential */

  /** - shall we interpolate in file, or shall we use analytical formula below? */

  /** -> deal first with the case of interpolating in files */
  if (pba->got_files[n_ncdm]==_TRUE_) {

    lastidx = pbadist_local->tablesize-1;
    if(q<pbadist_local->q[0]){
      //Handle q->0 case:
      *f0 = pbadist_local->f0[0];
    }
    else if(q>pbadist_local->q[lastidx]){
      //Handle q>qmax case (ensure continuous and derivable function with Boltzmann tail):
      qlast=pbadist_local->q[lastidx];
      f0last=pbadist_local->f0[lastidx];
      dqlast=qlast - pbadist_local->q[lastidx-1];
      df0last=f0last - pbadist_local->f0[lastidx-1];

      *f0 = f0last*exp(-(qlast-q)*df0last/f0last/dqlast);
    }
    else{
      //Do interpolation:
      class_call(array_interpolate_spline(
                                          pbadist_local->q,
                                          pbadist_local->tablesize,
                                          pbadist_local->f0,
                                          pbadist_local->d2f0,
                                          1,
                                          q,
                                          &pbadist_local->last_index,
                                          f0,
                                          1,
                                          pba->error_message),
                 pba->error_message,     pba->error_message);
    }
  }

  /** -> deal now with case of reading analytical function */
  else{
    /**
       Enter here your analytic expression(s) for the p.s.d.'s. If
       you need different p.s.d.'s for different species, put each
       p.s.d inside a condition, like for instance: if (n_ncdm==2) {
       *f0=...}.  Remember that n_ncdm = 0 refers to the first
       species.
    */

    /**************************************************/
    /*    FERMI-DIRAC INCLUDING CHEMICAL POTENTIALS   */
    /**************************************************/

    *f0 = 1.0/pow(2*_PI_,3)*(1./(exp(q-ksi)+1.) +1./(exp(q+ksi)+1.));

    /**************************************************/

    /** This form is only appropriate for approximate studies, since in
        reality the chemical potential are associated with flavor
        eigenstates, not mass eigenstates. It is easy to take this into
        account by introducing the mixing angles. In the later part
        (not read by the code) we illustrate how to do this */

    if (_FALSE_) {

      /* We must use the list of extra parameters read in input, stored in the
         ncdm_psd_parameter list, extracted above from the structure
         and now called param[..] */

      /* check that this list has been read */
      class_test(param == NULL,
                 pba->error_message,
                 "Analytic expression wants to use 'ncdm_psd_parameters', but they have not been entered!");

      /* extract values from the list (in this example, mixing angles) */
      double square_s12=param[0];
      double square_s23=param[1];
      double square_s13=param[2];

      /* infer mixing matrix */
      double mixing_matrix[3][3];
      int i;

      mixing_matrix[0][0]=pow(fabs(sqrt((1-square_s12)*(1-square_s13))),2);
      mixing_matrix[0][1]=pow(fabs(sqrt(square_s12*(1-square_s13))),2);
      mixing_matrix[0][2]=fabs(square_s13);
      mixing_matrix[1][0]=pow(fabs(sqrt((1-square_s12)*square_s13*square_s23)+sqrt(square_s12*(1-square_s23))),2);
      mixing_matrix[1][1]=pow(fabs(sqrt(square_s12*square_s23*square_s13)-sqrt((1-square_s12)*(1-square_s23))),2);
      mixing_matrix[1][2]=pow(fabs(sqrt(square_s23*(1-square_s13))),2);
      mixing_matrix[2][0]=pow(fabs(sqrt(square_s12*square_s23)-sqrt((1-square_s12)*square_s13*(1-square_s23))),2);
      mixing_matrix[2][1]=pow(sqrt((1-square_s12)*square_s23)+sqrt(square_s12*square_s13*(1-square_s23)),2);
      mixing_matrix[2][2]=pow(fabs(sqrt((1-square_s13)*(1-square_s23))),2);

      /* loop over flavor eigenstates and compute psd of mass eigenstates */
      *f0=0.0;
      for(i=0;i<3;i++){

    	*f0 += mixing_matrix[i][n_ncdm]*1.0/pow(2*_PI_,3)*(1./(exp(q-pba->ksi_ncdm[i])+1.) +1./(exp(q+pba->ksi_ncdm[i])+1.));

      }
    } /* end of region not used, but shown as an example */
  }

  return _SUCCESS_;
}

/**
 * This function is only used for the purpose of finding optimal
 * quadrature weigths. The logic is: if we can convolve accurately
 * f0(q) with this function, then we can convolve it accuractely with
 * any other relevant function.
 *
 * @param pbadist Input:  structure containing all parameters defining f0(q)
 * @param q       Input:  momentum
 * @param f0      Output: phase-space distribution
 */

int background_ncdm_test_function(
                                  void * pbadist,
                                  double q,
                                  double * test
                                  ) {

  double c = 2.0/(3.0*_zeta3_);
  double d = 120.0/(7.0*pow(_PI_,4));
  double e = 2.0/(45.0*_zeta5_);

  /** Using a + bq creates problems for otherwise acceptable distributions
      which diverges as 1/r or 1/r^2 for r->0 */
  *test = pow(2.0*_PI_,3)/6.0*(c*q*q-d*q*q*q-e*q*q*q*q);

  return _SUCCESS_;
}

/**
 * This function finds optimal quadrature weights for each ncdm
 * species
 *
 * @param ppr Input: precision structure
 * @param pba Input/Output: background structure
 */

int background_ncdm_init(
                         struct precision *ppr,
                         struct background *pba
                         ) {

  int index_q, k,tolexp,row,status,filenum;
  double f0m2,f0m1,f0,f0p1,f0p2,dq,q,df0dq,tmp1,tmp2;
  struct background_parameters_for_distributions pbadist;
  FILE *psdfile;

  pbadist.pba = pba;

  /* Allocate pointer arrays: */
  class_alloc(pba->q_ncdm, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->w_ncdm, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->q_ncdm_bg, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->w_ncdm_bg, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->dlnf0_dlnq_ncdm, sizeof(double*)*pba->N_ncdm,pba->error_message);

  /* Allocate pointers: */
  class_alloc(pba->q_size_ncdm,sizeof(int)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->q_size_ncdm_bg,sizeof(int)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->factor_ncdm,sizeof(double)*pba->N_ncdm,pba->error_message);

  for(k=0, filenum=0; k<pba->N_ncdm; k++){
    pbadist.n_ncdm = k;
    pbadist.q = NULL;
    pbadist.tablesize = 0;
    /*Do we need to read in a file to interpolate the distribution function? */
    if ((pba->got_files!=NULL)&&(pba->got_files[k]==_TRUE_)){
      psdfile = fopen(pba->ncdm_psd_files+filenum*_ARGUMENT_LENGTH_MAX_,"r");
      class_test(psdfile == NULL,pba->error_message,
                 "Could not open file %s!",pba->ncdm_psd_files+filenum*_ARGUMENT_LENGTH_MAX_);
      // Find size of table:
      for (row=0,status=2; status==2; row++){
        status = fscanf(psdfile,"%lf %lf",&tmp1,&tmp2);
      }
      rewind(psdfile);
      pbadist.tablesize = row-1;

      /*Allocate room for interpolation table: */
      class_alloc(pbadist.q,sizeof(double)*pbadist.tablesize,pba->error_message);
      class_alloc(pbadist.f0,sizeof(double)*pbadist.tablesize,pba->error_message);
      class_alloc(pbadist.d2f0,sizeof(double)*pbadist.tablesize,pba->error_message);
      for (row=0; row<pbadist.tablesize; row++){
        status = fscanf(psdfile,"%lf %lf",
                        &pbadist.q[row],&pbadist.f0[row]);
        //		printf("(q,f0) = (%g,%g)\n",pbadist.q[row],pbadist.f0[row]);
      }
      fclose(psdfile);
      /* Call spline interpolation: */
      class_call(array_spline_table_lines(pbadist.q,
                                          pbadist.tablesize,
                                          pbadist.f0,
                                          1,
                                          pbadist.d2f0,
                                          _SPLINE_EST_DERIV_,
                                          pba->error_message),
                 pba->error_message,
                 pba->error_message);
      filenum++;
    }

    /* Handle perturbation qsampling: */
    class_alloc(pba->q_ncdm[k],_QUADRATURE_MAX_*sizeof(double),pba->error_message);
    class_alloc(pba->w_ncdm[k],_QUADRATURE_MAX_*sizeof(double),pba->error_message);

    class_call(get_qsampling(pba->q_ncdm[k],
                             pba->w_ncdm[k],
                             &(pba->q_size_ncdm[k]),
                             _QUADRATURE_MAX_,
                             ppr->tol_ncdm,
                             pbadist.q,
                             pbadist.tablesize,
                             background_ncdm_test_function,
                             background_ncdm_distribution,
                             &pbadist,
                             pba->error_message),
               pba->error_message,
               pba->error_message);
    pba->q_ncdm[k]=realloc(pba->q_ncdm[k],pba->q_size_ncdm[k]*sizeof(double));
    pba->w_ncdm[k]=realloc(pba->w_ncdm[k],pba->q_size_ncdm[k]*sizeof(double));


    if (pba->background_verbose > 0)
      printf("ncdm species i=%d sampled with %d points for purpose of perturbation integration\n",
             k+1,
             pba->q_size_ncdm[k]);

    /* Handle background q_sampling: */
    class_alloc(pba->q_ncdm_bg[k],_QUADRATURE_MAX_BG_*sizeof(double),pba->error_message);
    class_alloc(pba->w_ncdm_bg[k],_QUADRATURE_MAX_BG_*sizeof(double),pba->error_message);

    class_call(get_qsampling(pba->q_ncdm_bg[k],
                             pba->w_ncdm_bg[k],
                             &(pba->q_size_ncdm_bg[k]),
                             _QUADRATURE_MAX_BG_,
                             ppr->tol_ncdm_bg,
                             pbadist.q,
                             pbadist.tablesize,
                             background_ncdm_test_function,
                             background_ncdm_distribution,
                             &pbadist,
                             pba->error_message),
               pba->error_message,
               pba->error_message);


    pba->q_ncdm_bg[k]=realloc(pba->q_ncdm_bg[k],pba->q_size_ncdm_bg[k]*sizeof(double));
    pba->w_ncdm_bg[k]=realloc(pba->w_ncdm_bg[k],pba->q_size_ncdm_bg[k]*sizeof(double));

    /** - in verbose mode, inform user of number of sampled momenta
        for background quantities */
    if (pba->background_verbose > 0)
      printf("ncdm species i=%d sampled with %d points for purpose of background integration\n",
             k+1,
             pba->q_size_ncdm_bg[k]);

    class_alloc(pba->dlnf0_dlnq_ncdm[k],
                pba->q_size_ncdm[k]*sizeof(double),
                pba->error_message);


    for (index_q=0; index_q<pba->q_size_ncdm[k]; index_q++) {
      q = pba->q_ncdm[k][index_q];
      class_call(background_ncdm_distribution(&pbadist,q,&f0),
                 pba->error_message,pba->error_message);

      //Loop to find appropriate dq:
      for(tolexp=_PSD_DERIVATIVE_EXP_MIN_; tolexp<_PSD_DERIVATIVE_EXP_MAX_; tolexp++){

        if (index_q == 0){
          dq = MIN((0.5-ppr->smallest_allowed_variation)*q,2*exp(tolexp)*(pba->q_ncdm[k][index_q+1]-q));
        }
        else if (index_q == pba->q_size_ncdm[k]-1){
          dq = exp(tolexp)*2.0*(pba->q_ncdm[k][index_q]-pba->q_ncdm[k][index_q-1]);
        }
        else{
          dq = exp(tolexp)*(pba->q_ncdm[k][index_q+1]-pba->q_ncdm[k][index_q-1]);
        }

        class_call(background_ncdm_distribution(&pbadist,q-2*dq,&f0m2),
                   pba->error_message,pba->error_message);
        class_call(background_ncdm_distribution(&pbadist,q+2*dq,&f0p2),
                   pba->error_message,pba->error_message);

        if (fabs((f0p2-f0m2)/f0)>sqrt(ppr->smallest_allowed_variation)) break;
      }

      class_call(background_ncdm_distribution(&pbadist,q-dq,&f0m1),
                 pba->error_message,pba->error_message);
      class_call(background_ncdm_distribution(&pbadist,q+dq,&f0p1),
                 pba->error_message,pba->error_message);
      //5 point estimate of the derivative:
      df0dq = (+f0m2-8*f0m1+8*f0p1-f0p2)/12.0/dq;
      //printf("df0dq[%g] = %g. dlf=%g ?= %g. f0 =%g.\n",q,df0dq,q/f0*df0dq,
      //Avoid underflow in extreme tail:
      if (fabs(f0)==0.)
        pba->dlnf0_dlnq_ncdm[k][index_q] = -q; /* valid for whatever f0 with exponential tail in exp(-q) */
      else
        pba->dlnf0_dlnq_ncdm[k][index_q] = q/f0*df0dq;
    }

    pba->factor_ncdm[k]=pba->deg_ncdm[k]*4*_PI_*pow(pba->T_cmb*pba->T_ncdm[k]*_k_B_,4)*8*_PI_*_G_
      /3./pow(_h_P_/2./_PI_,3)/pow(_c_,7)*_Mpc_over_m_*_Mpc_over_m_;

    /* If allocated, deallocate interpolation table:  */
    if ((pba->got_files!=NULL)&&(pba->got_files[k]==_TRUE_)){
      free(pbadist.q);
      free(pbadist.f0);
      free(pbadist.d2f0);
    }
  }


  return _SUCCESS_;
}

/**
 * For a given ncdm sepcies: given the quadrature weights, the mass
 * and the redshift, find background quantities by a quick weighted
 * sum over.  Input parameters passed as NULL pointers are not
 * evaluated for speed-up
 *
 * @param qvec     Input: smapled momenta
 * @param wvec     Input: quadrature weigths
 * @param qsize    Input: number of momenta/weigths
 * @param M        Input: mass
 * @param factor   Input: normalization factor for the p.s.d.
 * @param z        Input: redhsift
 * @param n        Output: number density
 * @param rho      Output: energy density
 * @param p        Output: pressure
 * @param drho_dM  Output: derivative used in next function
 * @param pseudo_p Ouput: pseudo-pressure used in perturbation module for fluid approx
 *
 */

int background_ncdm_momenta(
                            /* Only calculate for non-NULL pointers: */
                            double * qvec,
                            double * wvec,
                            int qsize,
                            double M,
                            double factor,
                            double z,
                            double * n,
                            double * rho, // density
                            double * p,   // pressure
                            double * drho_dM,  // d rho / d M used in next function
                            double * pseudo_p  // pseudo-p used in ncdm fluid approx
                            ) {

  int index_q;
  double epsilon;
  double q2;
  double factor2;
  double density_neutrino_at_a() ;
  double pressure_neutrino_at_a();
  double x,y;

  /** - rescale normalization at given redshift */
  factor2 = factor*pow(1+z,4);

  /** - initialize quantities */
  if (n!=NULL) *n = 0.;
  if (rho!=NULL) *rho = 0.;
  if (p!=NULL) *p = 0.;
  if (drho_dM!=NULL) *drho_dM = 0.;
  if (pseudo_p!=NULL) *pseudo_p = 0.;

  /** - loop over momenta */
  for (index_q=0; index_q<qsize; index_q++) {

    /* squared momentum */
    q2 = qvec[index_q]*qvec[index_q];

    /* energy */
    epsilon = sqrt(q2+M*M/(1.+z)/(1.+z));

    /* integrand of the various quantities */
    if (n!=NULL) *n += q2*wvec[index_q];
    if (rho!=NULL) *rho += q2*epsilon*wvec[index_q];
    if (p!=NULL) *p += q2*q2/3./epsilon*wvec[index_q];
    if (drho_dM!=NULL) *drho_dM += q2*M/(1.+z)/(1.+z)/epsilon*wvec[index_q];
    if (pseudo_p!=NULL) *pseudo_p += pow(q2/epsilon,3)/3.0*wvec[index_q];
  }

  /** - ajust normalization */
  if (n!=NULL) *n *= factor2*(1.+z);
  if (rho!=NULL) *rho *= factor2;
  if (p!=NULL) *p *= factor2;
  if (drho_dM!=NULL) *drho_dM *= factor2;
  if (pseudo_p!=NULL) *pseudo_p *=factor2;

  x = density_neutrino_at_a(1./(1.+z),M);
  y = pressure_neutrino_at_a(1./(1.+z),M);
  printf("rho_nu = %1.12e p_nu = %1.12e at a = %1.12e \n",x,y,1./(1.+z));
  //printf("p_nu = %e at a = %e \n",y,1./(1.+z));
  exit(1);
  return _SUCCESS_;
}

/* MODIFICATIONS INTRODUCED TO AVOID NUMERICAL INTEGRATION WITH WEIGHTS */

double density_neutrino_at_a(double a,double m)
{
  double result;
  //double z_nu;
  int index_sum;
  int index_sum_max = 5;
  double Tnu_0;
  double Bessel0,Bessel1;
  double Struve1, Struve0;
  double struve();
  double param0,param1;
  double g,x,y ; 

  Tnu_0 = 2.726*pow(4./11.,1./3.)*_k_B_/1.6e-19; // T_ncdm in eV
  y = 1./Tnu_0;
  result = 0.;
  param0 = 0.;
  param1 = 1.;
  g = 2.;

  for (index_sum=1;index_sum<index_sum_max;index_sum++)
    {
      //z_nu = pow(index_sum*a*m/2./Tnu_0,2);
      x = 1./Tnu_0; //2.131330993153e6 ;//index_sum*a*m/Tnu_0;
      Bessel0 = gsl_sf_bessel_Y0(1./Tnu_0);
      Bessel1 = gsl_sf_bessel_Y1(x);
      Struve0 = struve(param0,x);
      Struve1 = struve(param1,x);
      printf("Y0 = %1.12e  Y1 = %1.12e H0 = %1.12e  H1 = %1.12e at x = %1.12e \n",Bessel0,Bessel1,Struve0,Struve1,y);
      result += pow(-1.,index_sum)/pow(x,3)*( -pow(x,2) + 3.*_PI_*x*( Struve0 - Bessel0 )/2. + _PI_*(3. - pow(x,2)/2.)*(Bessel1 - Struve1) );
    }
  result = 8*_PI_*_G_/3.*g*pow(m,4)*result/2./pow(_PI_,2);
  
  return result;
}

double pressure_neutrino_at_a(double a,double m)
{
  double result;
  //double z_nu;
  int index_sum;
  int index_sum_max = 10000;
  double Tnu_0;
  double Bessel2;
  double Struve1, Struve0;
  double struve();
  //double yn();
  double param0,param1;
  double g,x; 
  int n = 2;

  Tnu_0 = 2.726*pow(4./11.,1./3.)*_k_B_/1.6e-19; // T_ncdm in eV
  result = 0.;
  param0 = 0.;
  param1 = 1.;
  g = 2.;

  for (index_sum=1;index_sum<index_sum_max;index_sum++)
    {
      //z_nu = pow(index_sum*a*m/2./Tnu_0,2);
      x = index_sum*a*m/Tnu_0 ; 
      //  y = yn(n,x);
      Bessel2 = gsl_sf_bessel_Yn(n,x);
      Struve0 = struve(param0,x);
      Struve1 = struve(param1,x);
      //printf("Y2 = %1.12e  H0 = %e H1 = %e at x = %1.12e \n",Bessel2,Struve0,Struve1,x);
      result += pow(-1.,index_sum+1)/pow(x,3)*( pow(x,2) + _PI_*( pow(x,2) - 3.)*( x*(Bessel2 + Struve0)/2. - Struve1 ));
    }
  result = 8*_PI_*_G_/3.*g*pow(m,4)*8.*result/48./pow(_PI_,2);
  
  return result;
}

/* MODIFICATIONS END */

/**
 * When the user passed in input the density fraction Omeha_ncdm or
 * omega_ncdm but not the mass, infer the mass with Newton iteration method.
 *
 * @param ppr    Input: precision structure
 * @param pba    Input/Output: background structure
 * @param n_ncdm Input: index of ncdm species
 */

int background_ncdm_M_from_Omega(
                                 struct precision *ppr,
                                 struct background *pba,
                                 int n_ncdm
                                 ) {
  double rho0,rho,n,M,deltaM,drhodM;
  int iter,maxiter=50;

  rho0 = pba->H0*pba->H0*pba->Omega0_ncdm[n_ncdm]; /*Remember that rho is defined such that H^2=sum(rho_i) */
  M = 0.0;

  background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm],
                          pba->w_ncdm_bg[n_ncdm],
                          pba->q_size_ncdm_bg[n_ncdm],
                          M,
                          pba->factor_ncdm[n_ncdm],
                          0.,
                          &n,
                          &rho,
                          NULL,
                          NULL,
                          NULL);

  /* Is the value of Omega less than a massless species?*/
  class_test(rho0<rho,pba->error_message,
             "The value of Omega for the %dth species, %g, is less than for a massless species! It should be atleast %g. Check your input.",
             n_ncdm,pba->Omega0_ncdm[n_ncdm],pba->Omega0_ncdm[n_ncdm]*rho/rho0);

  /* In the strict NR limit we have rho = n*(M) today, giving a zero'th order guess: */
  M = rho0/n; /* This is our guess for M. */
  for (iter=1; iter<=maxiter; iter++){

    /* Newton iteration. First get relevant quantities at M: */
    background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm],
                            pba->w_ncdm_bg[n_ncdm],
                            pba->q_size_ncdm_bg[n_ncdm],
                            M,
                            pba->factor_ncdm[n_ncdm],
                            0.,
                            NULL,
                            &rho,
                            NULL,
                            &drhodM,
                            NULL);

    deltaM = (rho0-rho)/drhodM; /* By definition of the derivative */
    if ((M+deltaM)<0.0) deltaM = -M/2.0; /* Avoid overshooting to negative M value. */
    M += deltaM; /* Update value of M.. */
    if (fabs(deltaM/M)<ppr->tol_M_ncdm){
      /* Accuracy reached.. */
      pba->M_ncdm[n_ncdm] = M;
      break;
    }
  }
  class_test(iter>=maxiter,pba->error_message,
             "Newton iteration could not converge on a mass for some reason.");
  return _SUCCESS_;
}

/**
 *  This function integrates the background over time, allocates and
 *  fills the background table
 *
 * @param ppr Input: precision structure
 * @param pba Input/Output: background structure
 */

int background_solve(
                     struct precision *ppr,
                     struct background *pba
                     ) {

  /** Summary: */

  /** - define local variables */

  /* contains all quantities relevant for the integration algorithm */
  struct generic_integrator_workspace gi;
  /* parameters and workspace for the background_derivs function */
  struct background_parameters_and_workspace bpaw;
  /* a growing table (since the number of time steps is not known a priori) */
  growTable gTable;
  /* needed for growing table */
  double * pData;
  /* needed for growing table */
  void * memcopy_result;
  /* initial conformal time */
  double tau_start;
  /* final conformal time */
  double tau_end;
  /* an index running over bi indices */
  int i;
  /* vector of quantities to be integrated */
  double * pvecback_integration;
  /* vector of all background quantities */
  double * pvecback;
  /* necessary for calling array_interpolate(), but never used */
  int last_index=0;
  /* comoving radius coordinate in Mpc (equal to conformal distance in flat case) */
  double comoving_radius=0.;

  bpaw.pba = pba;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);
  bpaw.pvecback = pvecback;

  /** - allocate vector of quantities to be integrated */
  class_alloc(pvecback_integration,pba->bi_size*sizeof(double),pba->error_message);

  /** - initialize generic integrator with initialize_generic_integrator() */

  /* Size of vector to integrate is (pba->bi_size-1) rather than
   * (pba->bi_size), since tau is not integrated.
   */
  class_call(initialize_generic_integrator((pba->bi_size-1),&gi),
             gi.error_message,
             pba->error_message);

  /** - impose initial conditions with background_initial_conditions() */
  class_call(background_initial_conditions(ppr,pba,pvecback,pvecback_integration),
             pba->error_message,
             pba->error_message);

  /* here tau_end is in fact the initial time (in the next loop
     tau_start = tau_end) */
  tau_end=pvecback_integration[pba->index_bi_tau];

  /** - create a growTable with gt_init() */
  class_call(gt_init(&gTable),
             gTable.error_message,
             pba->error_message);

  /* initialize the counter for the number of steps */
  pba->bt_size=0;

  /** - loop over integration steps : call background_functions(), find step size, save data in growTable with gt_add(), perform one step with generic_integrator(), store new value of tau */

  while (pvecback_integration[pba->index_bi_a] < pba->a_today) {

    tau_start = tau_end;

    /* -> find step size (trying to adjust the last step as close as possible to the one needed to reach a=a_today; need not be exact, difference corrected later) */
    class_call(background_functions(pba,pvecback_integration[pba->index_bi_a], pba->short_info, pvecback),
               pba->error_message,
               pba->error_message);

    if ((pvecback_integration[pba->index_bi_a]*(1.+ppr->back_integration_stepsize)) < pba->a_today) {
      tau_end = tau_start + ppr->back_integration_stepsize / (pvecback_integration[pba->index_bi_a]*pvecback[pba->index_bg_H]);
      /* no possible segmentation fault here: non-zeroness of "a" has been checked in background_functions() */
    }
    else {
      tau_end = tau_start + (pba->a_today/pvecback_integration[pba->index_bi_a]-1.) / (pvecback_integration[pba->index_bi_a]*pvecback[pba->index_bg_H]);
      /* no possible segmentation fault here: non-zeroness of "a" has been checked in background_functions() */
    }

    class_test((tau_end-tau_start)/tau_start < ppr->smallest_allowed_variation,
               pba->error_message,
               "integration step: relative change in time =%e < machine precision : leads either to numerical error or infinite loop",(tau_end-tau_start)/tau_start);

    /* -> save data in growTable */
    class_call(gt_add(&gTable,_GT_END_,(void *) pvecback_integration,sizeof(double)*pba->bi_size),
               gTable.error_message,
               pba->error_message);
    pba->bt_size++;

    /* -> perform one step */
    class_call(generic_integrator(background_derivs,
                                  tau_start,
                                  tau_end,
                                  pvecback_integration,
                                  &bpaw,
                                  ppr->tol_background_integration,
                                  ppr->smallest_allowed_variation,
                                  &gi),
               gi.error_message,
               pba->error_message);

    /* -> store value of tau */
    pvecback_integration[pba->index_bi_tau]=tau_end;

  }

  /** - save last data in growTable with gt_add() */
  class_call(gt_add(&gTable,_GT_END_,(void *) pvecback_integration,sizeof(double)*pba->bi_size),
             gTable.error_message,
             pba->error_message);
  pba->bt_size++;


  /* integration finished */

  /** - clean up generic integrator with cleanup_generic_integrator() */
  class_call(cleanup_generic_integrator(&gi),
             gi.error_message,
             pba->error_message);

  /** - retrieve data stored in the growTable with gt_getPtr() */
  class_call(gt_getPtr(&gTable,(void**)&pData),
             gTable.error_message,
             pba->error_message);

  /** - interpolate to get quantities precisely today with array_interpolate() */
  class_call(array_interpolate(
                               pData,
                               pba->bi_size,
                               pba->bt_size,
                               pba->index_bi_a,
                               pba->a_today,
                               &last_index,
                               pvecback_integration,
                               pba->bi_size,
                               pba->error_message),
             pba->error_message,
             pba->error_message);

  /* substitute last line with quantities today */
  for (i=0; i<pba->bi_size; i++)
    pData[(pba->bt_size-1)*pba->bi_size+i]=pvecback_integration[i];

  /** - deduce age of the Universe */
  /* -> age in Gyears */
  pba->age = pvecback_integration[pba->index_bi_time]/_Gyr_over_Mpc_;
  /* -> conformal age in Mpc */
  pba->conformal_age = pvecback_integration[pba->index_bi_tau];

  /** - allocate background tables */
  class_alloc(pba->tau_table,pba->bt_size * sizeof(double),pba->error_message);

  class_alloc(pba->z_table,pba->bt_size * sizeof(double),pba->error_message);

  class_alloc(pba->d2tau_dz2_table,pba->bt_size * sizeof(double),pba->error_message);

  class_alloc(pba->background_table,pba->bt_size * pba->bg_size * sizeof(double),pba->error_message);

  class_alloc(pba->d2background_dtau2_table,pba->bt_size * pba->bg_size * sizeof(double),pba->error_message);

  /** - In a loop over lines, fill background table using the result of the integration plus background_functions() */
  for (i=0; i < pba->bt_size; i++) {

    /* -> establish correspondance between the integrated variable and the bg variables */

    pba->tau_table[i] = pData[i*pba->bi_size+pba->index_bi_tau];

    class_test(pData[i*pba->bi_size+pba->index_bi_a] <= 0.,
               pba->error_message,
               "a = %e instead of strictly positiv",pData[i*pba->bi_size+pba->index_bi_a]);

    pba->z_table[i] = pba->a_today/pData[i*pba->bi_size+pba->index_bi_a]-1.;

    pvecback[pba->index_bg_time] = pData[i*pba->bi_size+pba->index_bi_time];
    pvecback[pba->index_bg_conf_distance] = pba->conformal_age - pData[i*pba->bi_size+pba->index_bi_tau];

    if (pba->sgnK == 0) comoving_radius = pvecback[pba->index_bg_conf_distance];
    else if (pba->sgnK == 1) comoving_radius = sin(sqrt(pba->K)*pvecback[pba->index_bg_conf_distance])/sqrt(pba->K);
    else if (pba->sgnK == -1) comoving_radius = sinh(sqrt(-pba->K)*pvecback[pba->index_bg_conf_distance])/sqrt(-pba->K);

    pvecback[pba->index_bg_ang_distance] = pba->a_today*comoving_radius/(1.+pba->z_table[i]);
    pvecback[pba->index_bg_lum_distance] = pba->a_today*comoving_radius*(1.+pba->z_table[i]);
    pvecback[pba->index_bg_rs] = pData[i*pba->bi_size+pba->index_bi_rs];

    /* -> compute all other quantities depending only on a*/
    class_call(background_functions(pba,pData[i*pba->bi_size+pba->index_bi_a], pba->long_info, pvecback),
               pba->error_message,
               pba->error_message);

    /* -> compute growth functions (valid in dust universe) */

    /* D = H \int [da/(aH)^3] = H \int [dtau/(aH^2)] = H * growth */
    pvecback[pba->index_bg_D] = pvecback[pba->index_bg_H]*pData[i*pba->bi_size+pba->index_bi_growth];

    /* f = [dlnD]/[dln a] = 1/(aH) [dlnD]/[dtau] = H'/(aH^2) + 1/(a^2 H^3 growth) */
    pvecback[pba->index_bg_f] = pvecback[pba->index_bg_H_prime]/pvecback[pba->index_bg_a]/pvecback[pba->index_bg_H]/pvecback[pba->index_bg_H] + 1./(pvecback[pba->index_bg_a]*pvecback[pba->index_bg_a]*pvecback[pba->index_bg_H]*pvecback[pba->index_bg_H]*pvecback[pba->index_bg_H]*pData[i*pba->bi_size+pba->index_bi_growth]);

    /* -> write in the table */
    memcopy_result = memcpy(pba->background_table + i*pba->bg_size,pvecback,pba->bg_size*sizeof(double));

    class_test(memcopy_result != pba->background_table + i*pba->bg_size,
               pba->error_message,
               "cannot copy data back to pba->background_table");
  }

  /** - free the growTable with gt_free() */

  class_call(gt_free(&gTable),
             gTable.error_message,
             pba->error_message);

  /** - fill tables of second derivatives (in view of spline interpolation) */
  class_call(array_spline_table_lines(pba->z_table,
                                      pba->bt_size,
                                      pba->tau_table,
                                      1,
                                      pba->d2tau_dz2_table,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  class_call(array_spline_table_lines(pba->tau_table,
                                      pba->bt_size,
                                      pba->background_table,
                                      pba->bg_size,
                                      pba->d2background_dtau2_table,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  /** - compute remaining "related parameters" */

  /** -> so-called "effective neutrino number", computed at earliest
      time in interpolation table. This should be seen as a
      definition: Neff is the equivalent number of
      instantaneously-decoupled neutrinos accounting for the
      radiation density, beyond photons */
  pba->Neff = (pba->background_table[pba->index_bg_Omega_r]
               *pba->background_table[pba->index_bg_rho_crit]
               -pba->background_table[pba->index_bg_rho_g])
    /(7./8.*pow(4./11.,4./3.)*pba->background_table[pba->index_bg_rho_g]);

  /** - done */
  if (pba->background_verbose > 0) {
    printf(" -> age = %f Gyr\n",pba->age);
    printf(" -> conformal age = %f Mpc\n",pba->conformal_age);
  }

  free(pvecback);
  free(pvecback_integration);

  return _SUCCESS_;

}

/**
 * Assign initial values to background integrated variables.
 *
 * @param ppr                  Input : pointer to precision structure
 * @param pba                  Input : pointer to background structure
 * @param pvecback             Input : vector of background quantitites used as workspace
 * @param pvecback_integration Output : vector of background quantitites to be integrated, returned with proper initial values
 * @return the error status
 */

int background_initial_conditions(
                                  struct precision *ppr,
                                  struct background *pba,
                                  double * pvecback, /* vector with argument pvecback[index_bg] (must be already allocated, normal format is sufficient) */
                                  double * pvecback_integration /* vector with argument pvecback_integration[index_bi] (must be already allocated with size pba->bi_size) */
                                  ) {

  /** Summary: */

  /** - define local variables */

  /* scale factor */
  double a;

  int counter,is_early_enough,n_ncdm;

  /** - fix initial value of \f$ a \f$ */
  a = ppr->a_ini_over_a_today_default * pba->a_today / _SCALE_BACK_;

  /* test the validity of this choice (e.g.: are massive neutrinos
     relativistic? etc.) */

  for (counter=0, is_early_enough = _FALSE_;
       (counter < _MAX_IT_) && (is_early_enough == _FALSE_ );
       counter++) {

    a *= _SCALE_BACK_;

    /** - compute initial H with background_functions() */
    class_call(background_functions(pba,a, pba->normal_info, pvecback),
               pba->error_message,
               pba->error_message);

    is_early_enough = _TRUE_;

    /* check that we are deep inside radiation domination, in order to
       use the approximation \f$ t=1/(2H) \f$ below */
    if (fabs(pvecback[pba->index_bg_Omega_r]-1.) > ppr->tol_initial_Omega_r)
      is_early_enough = _FALSE_;

    /* check that non-cold relics are relativistic, in view of
       integrating their perturbations starting from correct initial
       conditions */
    if (pba->has_ncdm == _TRUE_) {

      for (n_ncdm=0; n_ncdm<pba->N_ncdm; n_ncdm++) {

        if (fabs(pvecback[pba->index_bg_p_ncdm1+n_ncdm]/pvecback[pba->index_bg_rho_ncdm1+n_ncdm]-1./3.) > ppr->tol_ncdm_initial_w)
          is_early_enough = _FALSE_;
      }
    }
  }

  pvecback_integration[pba->index_bi_a] = a;

  /** - compute initial proper time, assuming radiation-dominated
      universe since Big Bang and therefore \f$ t=1/(2H) \f$ (good
      approximation for most purposes) */

  class_test(pvecback[pba->index_bg_H] <= 0.,
             pba->error_message,
             "H = %e instead of strictly positive",pvecback[pba->index_bg_H]);

  pvecback_integration[pba->index_bi_time] = 1./(2.* pvecback[pba->index_bg_H]);

  /** - compute initial conformal time, assuming radiation-dominated
      universe since Big Bang and therefore \f$ \tau=1/(aH) \f$
      (good approximation for most purposes) */
  pvecback_integration[pba->index_bi_tau] = 1./(a * pvecback[pba->index_bg_H]);

  /** - compute initial sound horizon, assuming c_s=1/sqrt(3) initially */
  pvecback_integration[pba->index_bi_rs] = pvecback_integration[pba->index_bi_tau]/sqrt(3.);

  /** - compute initial value of the integral over dtau/(aH^2),
      assumed to be proportional to a^4 during RD, but with arbitrary
      normalization */
  pvecback_integration[pba->index_bi_growth] = 1./(4.*a*a*pvecback[pba->index_bg_H]*pvecback[pba->index_bg_H]*pvecback[pba->index_bg_H]);

  return _SUCCESS_;

}

/**
 * Subroutine evaluating the derivative with respect to conformal time
 * of quantities which are integrated (a, t, etc).
 *
 * This is one of the few functions in the code which are passed to
 * the generic_integrator() routine.  Since generic_integrator()
 * should work with functions passed from various modules, the format
 * of the arguments is a bit special:
 *
 * - fixed input parameters and wokspaces are passed through a generic
 * pointer. Here, this is just a pointer to the background structure
 * and to a background vector, but generic_integrator() doesn't know
 * its fine structure.
 *
 * - the error management is a bit special: errors are not written as
 * usual to pba->error_message, but to a generic error_message passed
 * in the list of arguments.
 *
 * @param tau                      Input : conformal time
 * @param y                        Input : vector of variable
 * @param dy                       Output : its derivative (already allocated)
 * @param parameters_and_workspace Input: pointer to fixed parameters (e.g. indices)
 * @param error_message            Output : error message
 */
int background_derivs(
                      double tau,
                      double* y, /* vector with argument y[index_bi] (must be already allocated with size pba->bi_size) */
                      double* dy, /* vector with argument dy[index_bi]
                                     (must be already allocated with
                                     size pba->bi_size) */
                      void * parameters_and_workspace,
                      ErrorMsg error_message
                      ) {

  /** Summary: */

  /** - define local variables */

  struct background_parameters_and_workspace * pbpaw;
  struct background * pba;
  double * pvecback;

  pbpaw = parameters_and_workspace;
  pba =  pbpaw->pba;
  pvecback = pbpaw->pvecback;

  /** - Calculates functions of /f$ a /f$ with background_functions() */
  class_call(background_functions((struct background *)pba,y[pba->index_bi_a], pba->normal_info, pvecback),
             pba->error_message,
             error_message);

  /** - calculate /f$ a'=a^2 H /f$ */
  dy[pba->index_bi_a] = y[pba->index_bi_a] * y[pba->index_bi_a] * pvecback[pba->index_bg_H];

  /** - calculate /f$ t' = a /f$ */
  dy[pba->index_bi_time] = y[pba->index_bi_a];

  class_test(pvecback[pba->index_bg_rho_g] <= 0.,
             error_message,
             "rho_g = %e instead of strictly positive",pvecback[pba->index_bg_rho_g]);

  /** - calculate rs' = c_s */
  dy[pba->index_bi_rs] = 1./sqrt(3.*(1.+3.*pvecback[pba->index_bg_rho_b]/4./pvecback[pba->index_bg_rho_g]))*sqrt(1.-pba->K*y[pba->index_bi_rs]*y[pba->index_bi_rs]); // TBC: curvature correction

  /** calculate growth' = 1/(aH^2) */
  dy[pba->index_bi_growth] = 1./(y[pba->index_bi_a] * pvecback[pba->index_bg_H] * pvecback[pba->index_bg_H]);

  return _SUCCESS_;

}
