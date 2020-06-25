Changes in the code are presented here.
Every change is marked off with /** CHANGE: (a possible explanation) **/ in the different files of the code.


    -We use the fluid scheme for dark energy, then inside the .ini file is should be "use_ppf = no" to force this scheme. The suppression effect appears when w_DE>-1 and positive coupling constant of the model, for values of w_DE<-1 a negative value of the coupling constant of the model should be used to avoid instabilities, but no suppression effect take place.

    -The mymodel.ini file is full of all possible options for the modified version of this code (copying the idea of the explanatory.ini file of CLASS), so you can use whatever option you want by creating a smaller .ini file. Remember to always use "Omega_Lambda = 0" and "use_ppf = no" so the fluid scheme modified is used, and finally, the modification is done in the Newtonian Gauge, then "gauge = newtonian". The same change can be done using the Synchronous gauge, detailed equations for doing it in arXiv:1902.05532v2. The model has tricky details in the \theta_DM equations for the Synchronous gauge related to the residual gauge freedom. It is explained in arXiv:1902.05532v2 but you can also contact us for more information, if needed.

 

    -The new parameter denoted in arXiv:1902.05532v2 as \alpha* is defined in the code as "alpha_model". Because of that it is defined in:

         + Inside input.c inside "char * const target_namestrings[]" (lines: 230-232), defined as a new parameter as "class_read_double("alpha_model",pba->alpha_model)"; (lines: 608-609) and inside "int input_default_params" for the default value as pba->alpha_model=0.; (lines: 2997-2998)


         +Inside  background.h it is defined the new parameter as double alpha_model;  (lines: 162-163).
       


    -Inside perturbations.c file, the module solving the equations we modify the Newtonian equations for \theta_DM and \theta_DE following arXiv:1902.05532v2.  Then in the Newtonian gauge we modified:
    The changes are for Dark Matter (lines: 7271-7272)

         + Old equation: dy[pv->index_pt_theta_cdm] = - a_prime_over_a*y[pv->index_pt_theta_cdm] + metric_euler;

         + New equation: dy[pv->index_pt_theta_cdm] = - a_prime_over_a*y[pv->index_pt_theta_cdm] + metric_euler + (pba->alpha_model*pow(a,4)*pba->H0/(pba->Omega0_cdm))*(y[pv->index_pt_theta_fld]-y[pv->index_pt_theta_cdm]);

         The coupling term (C.T) in the dark matter equation is following the notation in arXiv:1902.05532v2
         + C.T= (\alpha · a)/(\rho_DM) (\theta_DE -\theta DM) = (\alpha* · H_0 · a^4) / (\Omega_DM) (\theta_DE -\theta DM)


    The changes are for Dark Energy (lines: 7377-7388):


         + Old equation: dy[pv->index_pt_theta_fld] = -(1.-3.*cs2)*a_prime_over_a*y[pv->index_pt_theta_fld] +cs2*k2/(1.+w_fld)*y[pv->index_pt_delta_fld] +metric_euler;

      
         + New equation: dy[pv->index_pt_theta_fld] = 
          -(1.-3.*cs2)*a_prime_over_a*y[pv->index_pt_theta_fld]
          +cs2*k2/(1.+w_fld)*y[pv->index_pt_delta_fld]
          +metric_euler-(pba->alpha_model*pow(a,4+3*w_fld)*pba->H0/(pba->Omega0_fld*(1+w_fld)))*(y[pv->index_pt_theta_fld]-y[pv->index_pt_theta_cdm]);


         The coupling term (C.T) in the dark energy equation is following the notation in arXiv:1902.05532v2
         + C.T= -(\alpha · a)/(\rho_DM) · (\rho_DM)/((1+w_DE)*\rho_DE) (\theta_DE -\theta DM) = -(\alpha* · H_0 · a^(4+3·w_DE) ) / (\Omega_DE · (1+w_DE)) · (\theta_DE -\theta DM)


Version 1.0 of CLASS DE-DM elastic interaction model, developed in CLASS version v2.7.2. Questions, doubts and/or comments about the code are welcomed at davidfiguer@usal.es. There is another version for the same coupling but between DE and Baryons (see paper arXiv:2004.14661 for details of the model) including the relevant modification of the CLASS schemes RSA and TCA, contact for details or for the code.

Forgive old personal comments (probably in spanish) done during the development of the code and not noticed during the depuration of the code.



