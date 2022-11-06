CLASS: Cosmic Linear Anisotropy Solving System  {#mainpage}
==============================================


This is the code for the dark vector energy model that includes both the background and the linear order perturbations.


The model has got three parameters ($Q$ and $s$ and $p$) that were modified to obtain the plots in the article.
We fix $s=0.2$ and varying $Q$ and $p$.

To run an example with $Q$ = 0 and $p = 1$, type 

./class ini_file_cde/s_02/Q_0/p_1.ini 

Look at the ini_file_cde/s_02/  folder to see four different cases of $Q$ to execute the code (these correspond to 0, 0.2, -0.5, -1).

If you want to play with the values of p, you can see inside
cd  ini_file_cde/s_02/Q_0/
cd  ini_file_cde/s_02/Q_02/
cd  ini_file_cde/s_02/Q_05/
cd  ini_file_cde/s_02/Q_1/

Here, you will find .ini files with the following values of $p$: 5, 3, 2, 1.5, 1.1, 1

The outputs will be saved in output_cde/.  