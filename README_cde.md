CLASS: Cosmic Linear Anisotropy Solving System  {#mainpage}
==============================================


This is the code for the coupled dark energy model (cde) (Arxi:) that includes both the background and the linear order perturbations.

The model
---------------------------

Cde is a model that include a coupling between cold dark matter (cdm) and a vector field A_{\mu} which is the dark energy source. 

The model has three free functions, $G_{2}$, $G_{3}$, and $f$, which are polynomials. 
 
$G_{2} = b_{2}X^{p_{2}}$, $G_{3} = b_{3}X^{p_{3}}$, and $f = X^{q}$, where $X = -A_{\mu}A^{\mu}/2$. 

$G_{2}$ and $G_{3}$ are the usual functions of Proca theories, while, $f$ is the coupling function. 

The interaction lagrangian is given by $Q f(X)\rho_{cdm}$, where $Q$ is a constant related to the strength of the coupling. 

To simplify the background equations of motion we use a new variable called $s$, it is given by:  $s = p_{2}/p$,  and $p = 2p_{3}-2p_{2}+1$, such that $p_{3}$ is determined by $s$.
Also, We can simplify the model in a way that $b_{3}$ is not necessary to specify.


The shape of the coupling
-------------------------
The model includes another parameter called beta. We have two ways to introduce the coupling term choosing either beta = 1 or beta = 0. The first one corresponds to arXiv:1907.12216v2 where they consider a lagrangian density for cdm and an additional one for the interaction. While the second one is associated with arXiv:2207.13682v1 where there is only one lagrangian density for cdm and the interaction.  

Background
-----------------------
At background level the objective is to solve a diferential equaiton for $\rho_{DE}$. Then, we need to include a variables to asignes the initial conditions, it is called rho_vf_ini. 

rho_vf_ini is a stimated values for $\rho_{DE}(t_{ini})$. Actually, This vaule is found by a procces denominated in Class as shooting, it is mechanism to find zero of functions.    


To find rho_vf_ini, we found the behavior of $\rho_{DE}$ in the radiation dominated epoch. We obtained $\rho_{DE} \propto a^{4s}$. The proportional factor is rho_vf_ini which is determined by Class (from our estimated starting value) using shooting. 



Perturbations
-----------------------
We include two auxiliary variables $\mathrm{Z}$ and $\mathrm{Q}$, these are the dark energy perturbed sources. The variables Z_ini and Q_ini represent their initial conditionals.  

If you want to run Class using dark energy as a vector field energy density, the gauge has to be 'Newtonian'


Paramereters
-----------------------
The model has got five parameters
_vf_parameters_1_ = $s$
_vf_parameters_2_ = $p_{2}$
_vf_parameters_3_ = $Q$
_vf_parameters_4_ = $q$
_vf_parameters_5_ = beta
and three initial conditions 
_vf_parameters[Z_ini, Q_ini, rho_vf_ini]_.

We fix q in background.c as $q=2p_{3}$ based in arXiv:1907.12216v2, however, if you wish to change it, look for '_q_vf = 2.p3_vf_' in background.c and perturbations.c and switch to
_q_vf  =  pba->vf_parameters_4_ .

The option _vf_tuning_index_ in the .ini files correspond to the parameter in _vf_parameters_ that will be used to find the correct $\Omega_DE$ today using the shooting routine such as in the scalar filed case. _vf_tuning_index_ will be the number of the initial condition _rho_vf_ini_ in _vf_parameters_, namely, 2.  


In _/ini_file_cde_ you will find a variety of .ini files with different values for $Q$, and $p_{2}$. $b_{2}$ must be less than zero to avoid negative values for the density energy of dark energy. We choose $b_{2} = -1$, nevertheless, other values can be used. They just will rescale the initial conditions of $\rho_{DE}$.


To run an example with $s=1.0$, $p_{2} = 1.0$, $Q = 0.2$, and $q=2.0$  type 

_./class ini_file_cde/s_1_Q_02.ini_ 

Keep in mind that not all parameter settings are possible.   


Numerical Solution for cdm
-------------------------------

We have an option in .ini files to solve numerically a differential equation for dark matter:
''_num_sol_cdm_vf = yes or no_''

If it is written yes, Class solves the following equations: $\rho_{cdm}'= -3\rho_{cdm} + \frac{q}{p2} \frac{Q f(X)}{(\beta+Qf(X))}\epsilon_{\phi}\rho_{cdm}$.
Where $\epsilon_{\phi}$ is defined in ArXiv:..., $\beta$ is the parameter beta and the derivative is in e-fold number.

On the other hand, if it is written no, Class makes cdm as $\rho_{cdm} = \rho_{cdm}^{0}(\beta + Qf(X))a^{-3}$, wiht $a$ the scale factor. 


In _ini_file_cde/beta_1/cdm_num_ you will find .ini files to execute Class solving numerically cdm. 


The file _vector_field.ini_ can run the vector field model. In this, you can change freely the parameters of the model. 


Our Modifications to Class
----------------------

We modified the following files: input.h, background.h, perturbations.h, input.c, background.c, perturbations.c. All modifications are indexed by "_vf_" 


input.h -we added a new target in the variable 'target_names' and changed '_NUM_TARGETS_' from 6 to 7.  

background.h -We added all the vector field parameters in the structure background.

perturbation.h -We added all the dark energy perturbed sources and parameters in the structure perturbations.

input.c - The following functions were modified to consider a new case ('Omega_vf') to make 'shooting': _input_shooting_, _input_needs_shooting_for_target_, _input_find_root_, _input_get_guess_, _input_try_unknown_parameters_.
Also, we wrote, in the function _input_read_parameters_species_, lines of code to permit Class to read the parameters related to the vector field. In the function _input_default_params_ we added the default values for model's parameters.

backgroun.c -In background_functions we wrote the dark energy source. 
In _background_indices_ we asigned a space in the background pointer (pba->) for the dark energy sources. 
In _background_initial_conditions_ we fixed the initial conditions for the energy density of the vectot file. 
In _background_derivs_ we solved the derivative of $\rho_{DE}$ and $\rho_{cdm}$.
Some changes were made in the following functions as well, _background_free_input_, _background_indices_, _background_solve_, _background_output_titles_, _background_output_budget_.

perturbation.c -In _perturbations_output_data_ we selected the vector field sourses which will be in the output data.
In _perturbations_output_titles_ we wrote the titles for the output of the vector field sources. 
In _perturbations_indices_ were asigned a space for the new sources in the perturbations pointer (ptt->).
In _perturbations_total_stress_energy_ we added the perturbed energy density and velocity for the vector field to 
In _perturbations_sources_ again we included the vector field sources. 
In _perturbations_print_variables_ we wrote the same as in the previous function.
In _perturbations_derivs_ are writen the evolution equations for coupled vector dark energy.
Some changes were made in the following functions as well, _perturbations_prepare_k_output_, _perturbations_vector_init_, _perturbations_initial_conditions_