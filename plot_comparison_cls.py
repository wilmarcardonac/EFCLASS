import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py
#import pylab as py

# Loading data files into arrays 

# Fiducial model without anisotropic stress 

Clfid = np.loadtxt('../data/Cl_fiducial_lensing_cl.dat',unpack=True)
Clfidnl = np.loadtxt('../data/Cl_fiducial_no_lensing_cl.dat',unpack=True)

# CYAN: e_pi = 0.2 cs2 = 0.4 f_pi = g_pi = 0

Cl_cyan_l = np.loadtxt('./output/Cl_cyan_lensing_cl.dat',unpack=True)
Cl_cyan_nl = np.loadtxt('./output/Cl_cyan_no_lensing_cl.dat',unpack=True)

# GREEN: cs2 = 3.33334 f_pi = 5 g_pi = 1.E10 e_pi = 0 

Cl_green_l = np.loadtxt('./output/Cl_green_lensing_cl.dat',unpack=True)
Cl_green_nl = np.loadtxt('./output/Cl_green_no_lensing_cl.dat',unpack=True)

# RED: cs2 = -1.19999 f_pi = -1.8 g_pi = 1.E-10 e_pi = 0 

Cl_red_l = np.loadtxt('./output/Cl_red_lensing_cl.dat',unpack=True)
Cl_red_nl = np.loadtxt('./output/Cl_red_no_lensing_cl.dat',unpack=True)

# 1-1

#py.loglog(Clfid[0],Clfid[1],label=r'$\Lambda CDM$')
#py.loglog(Cl[0],Cl[1],label='Including anisotropic stress')
#py.loglog(Clfid[0],abs(Clfid[1]-Clfidnl[1]),label='fiducial with lensing - fiducial without lensing')

py.loglog(Clfid[0],abs(Clfid[1]),label='fiducial with lensing')
py.loglog(Clfidnl[0],abs(Clfidnl[1]),label='fiducial without lensing (w/o)')
py.loglog(Cl_cyan_l[0],abs(Cl_cyan_l[1]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Cl_cyan_nl[0],abs(Cl_cyan_nl[1]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Cl_green_l[0],abs(Cl_green_l[1]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Cl_green_nl[0],abs(Cl_green_nl[1]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Cl_red_l[0],abs(Cl_red_l[1]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Cl_red_nl[0],abs(Cl_red_nl[1]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')

#py.loglog(Clbest[0],abs(Clfid[1]-Clbest[1]),label='fiducial with lensing - bestfit without lensing')

py.xlabel(r'$\ell$')

py.ylabel(r'$\ell*(1+\ell) C_\ell/2\pi$')

#py.xlim(1,3000)
#py.ylim(7.e-12,8.e-10)
py.legend(loc=0)

py.title('Correlation bins 1-1')

py.savefig('./output/correlation_1_1.pdf')
py.close()

py.loglog(Clfid[0],abs((Clfid[1]-Cl_cyan_l[1])/Clfid[1])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Clfidnl[0],abs((Clfidnl[1]-Cl_cyan_nl[1])/Clfidnl[1])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[1]-Cl_green_l[1])/Clfid[1])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[1]-Cl_green_nl[1])/Clfidnl[1])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[1]-Cl_red_l[1])/Clfid[1])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[1]-Cl_red_nl[1])/Clfidnl[1])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')
py.xlabel(r'$\ell$')
py.ylabel(r'$|\Delta C_\ell/C_\ell|[\%]$')
py.legend(loc=0)
py.title('Correlation bins 1-1')
py.savefig('./output/diff_correlation_1_1.pdf')
py.close()

#exit()
# 1-2

#py.loglog(Clfid[0],abs(Clfid[2]-Clfidnl[2]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[2]),label='fiducial without lensing')

#py.loglog(Clbest[0],abs(Clfid[2]-Clbest[2]),label='fiducial with lensing - bestfit without lensing')

py.loglog(Clfid[0],abs(Clfid[2]),label='fiducial with lensing')
py.loglog(Clfidnl[0],abs(Clfidnl[2]),label='fiducial without lensing (w/o)')
py.loglog(Cl_cyan_l[0],abs(Cl_cyan_l[2]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Cl_cyan_nl[0],abs(Cl_cyan_nl[2]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Cl_green_l[0],abs(Cl_green_l[2]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Cl_green_nl[0],abs(Cl_green_nl[2]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Cl_red_l[0],abs(Cl_red_l[2]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Cl_red_nl[0],abs(Cl_red_nl[2]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')


py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
#py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 1-2')

py.savefig('./output/correlation_1_2.pdf')

py.close()

py.loglog(Clfid[0],abs((Clfid[2]-Cl_cyan_l[2])/Clfid[2])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Clfidnl[0],abs((Clfidnl[2]-Cl_cyan_nl[2])/Clfidnl[2])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[2]-Cl_green_l[2])/Clfid[2])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[2]-Cl_green_nl[2])/Clfidnl[2])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[2]-Cl_red_l[2])/Clfid[2])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[2]-Cl_red_nl[2])/Clfidnl[2])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')
py.xlabel(r'$\ell$')
py.ylabel(r'$|\Delta C_\ell/C_\ell|[\%]$')
py.legend(loc=0)
py.title('Correlation bins 1-2')
py.savefig('./output/diff_correlation_1_2.pdf')
py.close()

# 1-3

#py.loglog(Clfid[0],abs(Clfid[3]-Clfidnl[3]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[3]),label='fiducial without lensing')

#py.loglog(Clbest[0],abs(Clfid[3]-Clbest[3]),label='fiducial with lensing - bestfit without lensing')

py.loglog(Clfid[0],abs(Clfid[3]),label='fiducial with lensing')
py.loglog(Clfidnl[0],abs(Clfidnl[3]),label='fiducial without lensing (w/o)')
py.loglog(Cl_cyan_l[0],abs(Cl_cyan_l[3]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Cl_cyan_nl[0],abs(Cl_cyan_nl[3]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Cl_green_l[0],abs(Cl_green_l[3]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Cl_green_nl[0],abs(Cl_green_nl[3]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Cl_red_l[0],abs(Cl_red_l[3]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Cl_red_nl[0],abs(Cl_red_nl[3]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
#py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 1-3')

py.savefig('./output/correlation_1_3.pdf')

py.close()

py.loglog(Clfid[0],abs((Clfid[3]-Cl_cyan_l[3])/Clfid[3])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Clfidnl[0],abs((Clfidnl[3]-Cl_cyan_nl[3])/Clfidnl[3])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[3]-Cl_green_l[3])/Clfid[3])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[3]-Cl_green_nl[3])/Clfidnl[3])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[3]-Cl_red_l[3])/Clfid[3])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[3]-Cl_red_nl[3])/Clfidnl[3])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')
py.xlabel(r'$\ell$')
py.ylabel(r'$|\Delta C_\ell/C_\ell|[\%]$')
py.legend(loc=0)
py.title('Correlation bins 1-3')
py.savefig('./output/diff_correlation_1_3.pdf')
py.close()

# 1-4

#py.loglog(Clfid[0],abs(Clfid[4]-Clfidnl[4]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[4]),label='fiducial without lensing')

#py.loglog(Clbest[0],abs(Clfid[4]-Clbest[4]),label='fiducial with lensing - bestfit without lensing')

py.loglog(Clfid[0],abs(Clfid[4]),label='fiducial with lensing')
py.loglog(Clfidnl[0],abs(Clfidnl[4]),label='fiducial without lensing (w/o)')
py.loglog(Cl_cyan_l[0],abs(Cl_cyan_l[4]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Cl_cyan_nl[0],abs(Cl_cyan_nl[4]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Cl_green_l[0],abs(Cl_green_l[4]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Cl_green_nl[0],abs(Cl_green_nl[4]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Cl_red_l[0],abs(Cl_red_l[4]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Cl_red_nl[0],abs(Cl_red_nl[4]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
#py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 1-4')

py.savefig('./output/correlation_1_4.pdf')

py.close()

py.loglog(Clfid[0],abs((Clfid[4]-Cl_cyan_l[4])/Clfid[4])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Clfidnl[0],abs((Clfidnl[4]-Cl_cyan_nl[4])/Clfidnl[4])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[4]-Cl_green_l[4])/Clfid[4])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[4]-Cl_green_nl[4])/Clfidnl[4])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[4]-Cl_red_l[4])/Clfid[4])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[4]-Cl_red_nl[4])/Clfidnl[4])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')
py.xlabel(r'$\ell$')
py.ylabel(r'$|\Delta C_\ell/C_\ell|[\%]$')
py.legend(loc=0)
py.title('Correlation bins 1-4')
py.savefig('./output/diff_correlation_1_4.pdf')
py.close()

# 1-5

#py.loglog(Clfid[0],abs(Clfid[5]-Clfidnl[5]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[5]),label='fiducial without lensing')

#py.loglog(Clbest[0],abs(Clfid[5]-Clbest[5]),label='fiducial with lensing - bestfit without lensing')

py.loglog(Clfid[0],abs(Clfid[5]),label='fiducial with lensing')
py.loglog(Clfidnl[0],abs(Clfidnl[5]),label='fiducial without lensing (w/o)')
py.loglog(Cl_cyan_l[0],abs(Cl_cyan_l[5]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Cl_cyan_nl[0],abs(Cl_cyan_nl[5]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Cl_green_l[0],abs(Cl_green_l[5]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Cl_green_nl[0],abs(Cl_green_nl[5]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Cl_red_l[0],abs(Cl_red_l[5]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Cl_red_nl[0],abs(Cl_red_nl[5]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
#py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 1-5')

py.savefig('./output/correlation_1_5.pdf')

py.close()

py.loglog(Clfid[0],abs((Clfid[5]-Cl_cyan_l[5])/Clfid[5])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Clfidnl[0],abs((Clfidnl[5]-Cl_cyan_nl[5])/Clfidnl[5])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[5]-Cl_green_l[5])/Clfid[5])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[5]-Cl_green_nl[5])/Clfidnl[5])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[5]-Cl_red_l[5])/Clfid[5])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[5]-Cl_red_nl[5])/Clfidnl[5])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')
py.xlabel(r'$\ell$')
py.ylabel(r'$|\Delta C_\ell/C_\ell|[\%]$')
py.legend(loc=0)
py.title('Correlation bins 1-5')
py.savefig('./output/diff_correlation_1_5.pdf')
py.close()


# 2-2

#py.loglog(Clfid[0],abs(Clfid[6]-Clfidnl[6]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[6]),label='fiducial without lensing')

#py.loglog(Clbest[0],abs(Clfid[6]-Clbest[6]),label='fiducial with lensing - bestfit without lensing')

py.loglog(Clfid[0],abs(Clfid[6]),label='fiducial with lensing')
py.loglog(Clfidnl[0],abs(Clfidnl[6]),label='fiducial without lensing (w/o)')
py.loglog(Cl_cyan_l[0],abs(Cl_cyan_l[6]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Cl_cyan_nl[0],abs(Cl_cyan_nl[6]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Cl_green_l[0],abs(Cl_green_l[6]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Cl_green_nl[0],abs(Cl_green_nl[6]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Cl_red_l[0],abs(Cl_red_l[6]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Cl_red_nl[0],abs(Cl_red_nl[6]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
#py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 2-2')

py.savefig('./output/correlation_2_2.pdf')

py.close()

py.loglog(Clfid[0],abs((Clfid[6]-Cl_cyan_l[6])/Clfid[6])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Clfidnl[0],abs((Clfidnl[6]-Cl_cyan_nl[6])/Clfidnl[6])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[6]-Cl_green_l[6])/Clfid[6])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[6]-Cl_green_nl[6])/Clfidnl[6])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[6]-Cl_red_l[6])/Clfid[6])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[6]-Cl_red_nl[6])/Clfidnl[6])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')
py.xlabel(r'$\ell$')
py.ylabel(r'$|\Delta C_\ell/C_\ell|[\%]$')
py.legend(loc=0)
py.title('Correlation bins 2-2')
py.savefig('./output/diff_correlation_2_2.pdf')
py.close()

# 2-3

#py.loglog(Clfid[0],abs(Clfid[7]-Clfidnl[7]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[7]),label='fiducial without lensing')

#py.loglog(Clbest[0],abs(Clfid[7]-Clbest[7]),label='fiducial with lensing - bestfit without lensing')

py.loglog(Clfid[0],abs(Clfid[7]),label='fiducial with lensing')
py.loglog(Clfidnl[0],abs(Clfidnl[7]),label='fiducial without lensing (w/o)')
py.loglog(Cl_cyan_l[0],abs(Cl_cyan_l[7]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Cl_cyan_nl[0],abs(Cl_cyan_nl[7]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Cl_green_l[0],abs(Cl_green_l[7]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Cl_green_nl[0],abs(Cl_green_nl[7]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Cl_red_l[0],abs(Cl_red_l[7]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Cl_red_nl[0],abs(Cl_red_nl[7]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')

#py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 2-3')

py.savefig('./output/correlation_2_3.pdf')

py.close()

py.loglog(Clfid[0],abs((Clfid[7]-Cl_cyan_l[7])/Clfid[7])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Clfidnl[0],abs((Clfidnl[7]-Cl_cyan_nl[7])/Clfidnl[7])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[7]-Cl_green_l[7])/Clfid[7])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[7]-Cl_green_nl[7])/Clfidnl[7])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[7]-Cl_red_l[7])/Clfid[7])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[7]-Cl_red_nl[7])/Clfidnl[7])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')
py.xlabel(r'$\ell$')
py.ylabel(r'$|\Delta C_\ell/C_\ell|[\%]$')
py.legend(loc=0)
py.title('Correlation bins 2-3')
py.savefig('./output/diff_correlation_2_3.pdf')
py.close()

# 2-4

#py.loglog(Clfid[0],abs(Clfid[8]-Clfidnl[8]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[8]),label='fiducial without lensing')

#py.loglog(Clbest[0],abs(Clfid[8]-Clbest[8]),label='fiducial with lensing - bestfit without lensing')

py.loglog(Clfid[0],abs(Clfid[8]),label='fiducial with lensing')
py.loglog(Clfidnl[0],abs(Clfidnl[8]),label='fiducial without lensing (w/o)')
py.loglog(Cl_cyan_l[0],abs(Cl_cyan_l[8]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Cl_cyan_nl[0],abs(Cl_cyan_nl[8]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Cl_green_l[0],abs(Cl_green_l[8]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Cl_green_nl[0],abs(Cl_green_nl[8]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Cl_red_l[0],abs(Cl_red_l[8]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Cl_red_nl[0],abs(Cl_red_nl[8]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
#py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 2-4')

py.savefig('./output/correlation_2_4.pdf')

py.close()

py.loglog(Clfid[0],abs((Clfid[8]-Cl_cyan_l[8])/Clfid[8])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Clfidnl[0],abs((Clfidnl[8]-Cl_cyan_nl[8])/Clfidnl[8])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[8]-Cl_green_l[8])/Clfid[8])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[8]-Cl_green_nl[8])/Clfidnl[8])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[8]-Cl_red_l[8])/Clfid[8])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[8]-Cl_red_nl[8])/Clfidnl[8])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')
py.xlabel(r'$\ell$')
py.ylabel(r'$|\Delta C_\ell/C_\ell|[\%]$')
py.legend(loc=0)
py.title('Correlation bins 2-4')
py.savefig('./output/diff_correlation_2_4.pdf')
py.close()

# 2-5

#py.loglog(Clfid[0],abs(Clfid[9]-Clfidnl[9]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[9]),label='fiducial without lensing')

#py.loglog(Clbest[0],abs(Clfid[9]-Clbest[9]),label='fiducial with lensing - bestfit without lensing')

py.loglog(Clfid[0],abs(Clfid[9]),label='fiducial with lensing')
py.loglog(Clfidnl[0],abs(Clfidnl[9]),label='fiducial without lensing (w/o)')
py.loglog(Cl_cyan_l[0],abs(Cl_cyan_l[9]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Cl_cyan_nl[0],abs(Cl_cyan_nl[9]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Cl_green_l[0],abs(Cl_green_l[9]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Cl_green_nl[0],abs(Cl_green_nl[9]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Cl_red_l[0],abs(Cl_red_l[9]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Cl_red_nl[0],abs(Cl_red_nl[9]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
#py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 2-5')

py.savefig('./output/correlation_2_5.pdf')

py.close()

py.loglog(Clfid[0],abs((Clfid[9]-Cl_cyan_l[9])/Clfid[9])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Clfidnl[0],abs((Clfidnl[9]-Cl_cyan_nl[9])/Clfidnl[9])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[9]-Cl_green_l[9])/Clfid[9])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[9]-Cl_green_nl[9])/Clfidnl[9])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[9]-Cl_red_l[9])/Clfid[9])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[9]-Cl_red_nl[9])/Clfidnl[9])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')
py.xlabel(r'$\ell$')
py.ylabel(r'$|\Delta C_\ell/C_\ell|[\%]$')
py.legend(loc=0)
py.title('Correlation bins 2-5')
py.savefig('./output/diff_correlation_2_5.pdf')
py.close()

# 3-3

#py.loglog(Clfid[0],abs(Clfid[10]-Clfidnl[10]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[10]),label='fiducial without lensing')

#py.loglog(Clbest[0],abs(Clfid[10]-Clbest[10]),label='fiducial with lensing - bestfit without lensing')

py.loglog(Clfid[0],abs(Clfid[10]),label='fiducial with lensing')
py.loglog(Clfidnl[0],abs(Clfidnl[10]),label='fiducial without lensing (w/o)')
py.loglog(Cl_cyan_l[0],abs(Cl_cyan_l[10]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Cl_cyan_nl[0],abs(Cl_cyan_nl[1]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Cl_green_l[0],abs(Cl_green_l[1]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Cl_green_nl[0],abs(Cl_green_nl[10]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Cl_red_l[0],abs(Cl_red_l[1]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Cl_red_nl[0],abs(Cl_red_nl[1]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
#py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 3-3')

py.savefig('./output/correlation_3_3.pdf')

py.close()

py.loglog(Clfid[0],abs((Clfid[10]-Cl_cyan_l[10])/Clfid[10])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Clfidnl[0],abs((Clfidnl[10]-Cl_cyan_nl[10])/Clfidnl[10])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[10]-Cl_green_l[10])/Clfid[10])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[10]-Cl_green_nl[10])/Clfidnl[10])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[10]-Cl_red_l[10])/Clfid[10])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[10]-Cl_red_nl[10])/Clfidnl[10])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')
py.xlabel(r'$\ell$')
py.ylabel(r'$|\Delta C_\ell/C_\ell|[\%]$')
py.legend(loc=0)
py.title('Correlation bins 3-3')
py.savefig('./output/diff_correlation_3_3.pdf')
py.close()

# 3-4

#py.loglog(Clfid[0],abs(Clfid[11]-Clfidnl[11]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[11]),label='fiducial without lensing')

#py.loglog(Clbest[0],abs(Clfid[11]-Clbest[11]),label='fiducial with lensing - bestfit without lensing')

py.loglog(Clfid[0],abs(Clfid[11]),label='fiducial with lensing')
py.loglog(Clfidnl[0],abs(Clfidnl[11]),label='fiducial without lensing (w/o)')
py.loglog(Cl_cyan_l[0],abs(Cl_cyan_l[11]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Cl_cyan_nl[0],abs(Cl_cyan_nl[11]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Cl_green_l[0],abs(Cl_green_l[11]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Cl_green_nl[0],abs(Cl_green_nl[11]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Cl_red_l[0],abs(Cl_red_l[11]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Cl_red_nl[0],abs(Cl_red_nl[11]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
#py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 3-4')

py.savefig('./output/correlation_3_4.pdf')

py.close()

py.loglog(Clfid[0],abs((Clfid[11]-Cl_cyan_l[11])/Clfid[11])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Clfidnl[0],abs((Clfidnl[11]-Cl_cyan_nl[11])/Clfidnl[11])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[11]-Cl_green_l[11])/Clfid[11])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[11]-Cl_green_nl[11])/Clfidnl[11])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[11]-Cl_red_l[11])/Clfid[11])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[11]-Cl_red_nl[11])/Clfidnl[11])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')
py.xlabel(r'$\ell$')
py.ylabel(r'$|\Delta C_\ell/C_\ell|[\%]$')
py.legend(loc=0)
py.title('Correlation bins 3-4')
py.savefig('./output/diff_correlation_3_4.pdf')
py.close()

# 3-5

#py.loglog(Clfid[0],abs(Clfid[12]-Clfidnl[12]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[12]),label='fiducial without lensing')

#py.loglog(Clbest[0],abs(Clfid[12]-Clbest[12]),label='fiducial with lensing - bestfit without lensing')

py.loglog(Clfid[0],abs(Clfid[12]),label='fiducial with lensing')
py.loglog(Clfidnl[0],abs(Clfidnl[12]),label='fiducial without lensing (w/o)')
py.loglog(Cl_cyan_l[0],abs(Cl_cyan_l[12]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Cl_cyan_nl[0],abs(Cl_cyan_nl[12]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Cl_green_l[0],abs(Cl_green_l[12]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Cl_green_nl[0],abs(Cl_green_nl[12]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Cl_red_l[0],abs(Cl_red_l[12]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Cl_red_nl[0],abs(Cl_red_nl[12]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
#py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 3-5')

py.savefig('./output/correlation_3_5.pdf')

py.close()

py.loglog(Clfid[0],abs((Clfid[12]-Cl_cyan_l[12])/Clfid[12])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Clfidnl[0],abs((Clfidnl[12]-Cl_cyan_nl[12])/Clfidnl[12])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[12]-Cl_green_l[12])/Clfid[12])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[12]-Cl_green_nl[12])/Clfidnl[12])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[12]-Cl_red_l[12])/Clfid[12])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[12]-Cl_red_nl[12])/Clfidnl[12])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')
py.xlabel(r'$\ell$')
py.ylabel(r'$|\Delta C_\ell/C_\ell|[\%]$')
py.legend(loc=0)
py.title('Correlation bins 3-5')
py.savefig('./output/diff_correlation_3_5.pdf')
py.close()

# 4-4

#py.loglog(Clfid[0],abs(Clfid[13]-Clfidnl[13]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[13]),label='fiducial without lensing')

#py.loglog(Clbest[0],abs(Clfid[13]-Clbest[13]),label='fiducial with lensing - bestfit without lensing')

py.loglog(Clfid[0],abs(Clfid[13]),label='fiducial with lensing')
py.loglog(Clfidnl[0],abs(Clfidnl[13]),label='fiducial without lensing (w/o)')
py.loglog(Cl_cyan_l[0],abs(Cl_cyan_l[13]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Cl_cyan_nl[0],abs(Cl_cyan_nl[13]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Cl_green_l[0],abs(Cl_green_l[13]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Cl_green_nl[0],abs(Cl_green_nl[13]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Cl_red_l[0],abs(Cl_red_l[13]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Cl_red_nl[0],abs(Cl_red_nl[13]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
#py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 4-4')

py.savefig('./output/correlation_4_4.pdf')

py.close()

py.loglog(Clfid[0],abs((Clfid[13]-Cl_cyan_l[13])/Clfid[13])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Clfidnl[0],abs((Clfidnl[13]-Cl_cyan_nl[13])/Clfidnl[13])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[13]-Cl_green_l[13])/Clfid[13])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[13]-Cl_green_nl[13])/Clfidnl[13])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[13]-Cl_red_l[13])/Clfid[13])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[13]-Cl_red_nl[13])/Clfidnl[13])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')
py.xlabel(r'$\ell$')
py.ylabel(r'$|\Delta C_\ell/C_\ell|[\%]$')
py.legend(loc=0)
py.title('Correlation bins 4-4')
py.savefig('./output/diff_correlation_4_4.pdf')
py.close()

# 4-5

#py.loglog(Clfid[0],abs(Clfid[14]-Clfidnl[14]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[14]),label='fiducial without lensing')

#py.loglog(Clbest[0],abs(Clfid[14]-Clbest[14]),label='fiducial with lensing - bestfit without lensing')

py.loglog(Clfid[0],abs(Clfid[14]),label='fiducial with lensing')
py.loglog(Clfidnl[0],abs(Clfidnl[14]),label='fiducial without lensing (w/o)')
py.loglog(Cl_cyan_l[0],abs(Cl_cyan_l[14]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Cl_cyan_nl[0],abs(Cl_cyan_nl[14]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Cl_green_l[0],abs(Cl_green_l[14]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Cl_green_nl[0],abs(Cl_green_nl[14]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Cl_red_l[0],abs(Cl_red_l[14]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Cl_red_nl[0],abs(Cl_red_nl[14]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
#py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 4-5')

py.savefig('./output/correlation_4_5.pdf')

py.close()

py.loglog(Clfid[0],abs((Clfid[14]-Cl_cyan_l[14])/Clfid[14])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Clfidnl[0],abs((Clfidnl[14]-Cl_cyan_nl[14])/Clfidnl[14])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[14]-Cl_green_l[14])/Clfid[14])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[14]-Cl_green_nl[14])/Clfidnl[14])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[14]-Cl_red_l[14])/Clfid[14])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[14]-Cl_red_nl[14])/Clfidnl[14])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')
py.xlabel(r'$\ell$')
py.ylabel(r'$|\Delta C_\ell/C_\ell|[\%]$')
py.legend(loc=0)
py.title('Correlation bins 4-5')
py.savefig('./output/diff_correlation_4_5.pdf')
py.close()


# 5-5

#py.loglog(Clfid[0],abs(Clfid[15]-Clfidnl[15]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[15]),label='fiducial without lensing')

#py.loglog(Clbest[0],abs(Clfid[15]-Clbest[15]),label='fiducial with lensing - bestfit without lensing')

py.loglog(Clfid[0],abs(Clfid[15]),label='fiducial with lensing')
py.loglog(Clfidnl[0],abs(Clfidnl[15]),label='fiducial without lensing (w/o)')
py.loglog(Cl_cyan_l[0],abs(Cl_cyan_l[15]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Cl_cyan_nl[0],abs(Cl_cyan_nl[15]),label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Cl_green_l[0],abs(Cl_green_l[15]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Cl_green_nl[0],abs(Cl_green_nl[15]),label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Cl_red_l[0],abs(Cl_red_l[15]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Cl_red_nl[0],abs(Cl_red_nl[15]),label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
#py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 5-5')

py.savefig('./output/correlation_5_5.pdf')

py.close()

py.loglog(Clfid[0],abs((Clfid[15]-Cl_cyan_l[15])/Clfid[15])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$')
py.loglog(Clfidnl[0],abs((Clfidnl[15]-Cl_cyan_nl[15])/Clfidnl[15])*1.e2,label=r'$e_{\pi} = 0.2$ $c_{s}^{2} = 0.4$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[15]-Cl_green_l[15])/Clfid[15])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[15]-Cl_green_nl[15])/Clfidnl[15])*1.e2,label=r'$f_{\pi} = 5$ $c_{s}^{2} = 3.3$ $g_{\pi}=10^{10}$ (w/o)')
py.loglog(Clfid[0],abs((Clfid[15]-Cl_red_l[15])/Clfid[15])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$')
py.loglog(Clfidnl[0],abs((Clfidnl[15]-Cl_red_nl[15])/Clfidnl[15])*1.e2,label=r'$f_{\pi} = -1.8$ $c_{s}^{2} = -1.2$ $g_{\pi}=10^{-10}$ (w/o)')
py.xlabel(r'$\ell$')
py.ylabel(r'$|\Delta C_\ell/C_\ell|[\%]$')
py.legend(loc=0)
py.title('Correlation bins 5-5')
py.savefig('./output/diff_correlation_5_5.pdf')
py.close()

exit()
