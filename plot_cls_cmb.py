import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py
#import pylab as py

################################
# Loading data files into arrays
################################ 

# LCDM

cambunits = (2.726e6)**2

#Cl = np.loadtxt('./output/Cl_designer_lcdm_cl.dat',unpack=True)

#Cl1 = np.loadtxt('./output/Cl_lcdm_cl.dat',unpack=True)

Cl_nu = np.loadtxt('./output/Cl_lcdm_nu_cl.dat',unpack=True)

Cl_nu_CLASS = np.loadtxt('./output/Cl_lcdm_nu_CLASS_cl.dat',unpack=True)

#Cl1 = np.loadtxt('./output/Cl_husawicki_cl.dat',unpack=True)
#Cl2 = np.loadtxt('./output/Cl_husawicki_cl_old.dat',unpack=True)
#mPk = np.loadtxt('./output/Cl_fiducial_lensing_ADE_pk.dat',unpack=True)
#pert = np.loadtxt('./output/Cl_fiducial_lensing_ADE_perturbations_k0_s.dat',unpack=True) 

fig = py.figure()

py.plot(Cl_nu_CLASS[0],(Cl_nu[1]-Cl_nu_CLASS[1])/Cl_nu_CLASS[1]*100.,label=r'CLASS vs. Ours',c='black')

py.xlabel(r'$\ell$')

py.ylabel(r'$\Delta C_\ell/C_\ell [\%]$')

py.xscale('log')

#py.xlim(1,30)

#py.ylim(600,1200)
py.legend(loc=0)

py.savefig('./output/relative_difference_cl.pdf')

py.close()

fig = py.figure()

py.plot(Cl_nu[0],Cl_nu[1]*cambunits,label=r'$\Lambda CDM OURS$',c='blue')
py.plot(Cl_nu_CLASS[0],Cl_nu_CLASS[1]*cambunits,label=r'$\Lambda CDM CLASS$',c='red')

py.xlabel(r'$\ell$')

py.ylabel(r'$\ell (\ell+1)C_\ell/2\pi$')

py.xscale('log')

#py.xlim(1,30)

#py.ylim(600,1200)
py.legend(loc=0)

py.savefig('./output/CMB_cl.pdf')

py.close()

exit()

py.loglog(mPkfid[0],mPkfid[1],label=r'$\Lambda CDM$')
py.loglog(mPk[0],mPk[1],label='Including anisotropic stress')
py.xlabel(r'$k$')
py.ylabel(r'$P(k)$')
py.xlim(1.e-4,1.e0)
py.ylim(1.e1,1.e5)
py.legend(loc=0)
py.savefig('./output/DEA_pk.pdf')
py.close()

py.loglog(pertfid[1],abs(pertfid[15]),label=r'$\Lambda CDM$')
py.loglog(pert[1],abs(pert[15]),label='Including anisotropic stress')
py.xlabel(r'$a$')
py.ylabel(r'$\delta$')
#py.xlim(1.e-4,1.e0)
#py.ylim(1.e1,1.e5)
py.legend(loc=0)
py.savefig('./output/perturbations_evolution.pdf')
py.close()

exit()
# 1-2

py.loglog(Clfid[0],abs(Clfid[2]-Clfidnl[2]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[2]),label='fiducial without lensing')

py.loglog(Clbest[0],abs(Clfid[2]-Clbest[2]),label='fiducial with lensing - bestfit without lensing')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 1-2')

py.savefig('correlation_1_2.pdf')

py.close()

# 1-3

py.loglog(Clfid[0],abs(Clfid[3]-Clfidnl[3]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[3]),label='fiducial without lensing')

py.loglog(Clbest[0],abs(Clfid[3]-Clbest[3]),label='fiducial with lensing - bestfit without lensing')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 1-3')

py.savefig('correlation_1_3.pdf')

py.close()

# 1-4

py.loglog(Clfid[0],abs(Clfid[4]-Clfidnl[4]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[4]),label='fiducial without lensing')

py.loglog(Clbest[0],abs(Clfid[4]-Clbest[4]),label='fiducial with lensing - bestfit without lensing')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 1-4')

py.savefig('correlation_1_4.pdf')

py.close()

# 1-5

py.loglog(Clfid[0],abs(Clfid[5]-Clfidnl[5]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[5]),label='fiducial without lensing')

py.loglog(Clbest[0],abs(Clfid[5]-Clbest[5]),label='fiducial with lensing - bestfit without lensing')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 1-5')

py.savefig('correlation_1_5.pdf')

py.close()

# 2-2

py.loglog(Clfid[0],abs(Clfid[6]-Clfidnl[6]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[6]),label='fiducial without lensing')

py.loglog(Clbest[0],abs(Clfid[6]-Clbest[6]),label='fiducial with lensing - bestfit without lensing')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 2-2')

py.savefig('correlation_2_2.pdf')

py.close()

# 2-3

py.loglog(Clfid[0],abs(Clfid[7]-Clfidnl[7]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[7]),label='fiducial without lensing')

py.loglog(Clbest[0],abs(Clfid[7]-Clbest[7]),label='fiducial with lensing - bestfit without lensing')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')

py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 2-3')

py.savefig('correlation_2_3.pdf')

py.close()

# 2-4

py.loglog(Clfid[0],abs(Clfid[8]-Clfidnl[8]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[8]),label='fiducial without lensing')

py.loglog(Clbest[0],abs(Clfid[8]-Clbest[8]),label='fiducial with lensing - bestfit without lensing')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 2-4')

py.savefig('correlation_2_4.pdf')

py.close()

# 2-5

py.loglog(Clfid[0],abs(Clfid[9]-Clfidnl[9]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[9]),label='fiducial without lensing')

py.loglog(Clbest[0],abs(Clfid[9]-Clbest[9]),label='fiducial with lensing - bestfit without lensing')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 2-5')

py.savefig('correlation_2_5.pdf')

py.close()

# 3-3

py.loglog(Clfid[0],abs(Clfid[10]-Clfidnl[10]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[10]),label='fiducial without lensing')

py.loglog(Clbest[0],abs(Clfid[10]-Clbest[10]),label='fiducial with lensing - bestfit without lensing')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 3-3')

py.savefig('correlation_3_3.pdf')

py.close()

# 3-4

py.loglog(Clfid[0],abs(Clfid[11]-Clfidnl[11]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[11]),label='fiducial without lensing')

py.loglog(Clbest[0],abs(Clfid[11]-Clbest[11]),label='fiducial with lensing - bestfit without lensing')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 3-4')

py.savefig('correlation_3_4.pdf')

py.close()

# 3-5

py.loglog(Clfid[0],abs(Clfid[12]-Clfidnl[12]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[12]),label='fiducial without lensing')

py.loglog(Clbest[0],abs(Clfid[12]-Clbest[12]),label='fiducial with lensing - bestfit without lensing')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 3-5')

py.savefig('correlation_3_5.pdf')

py.close()

# 4-4

py.loglog(Clfid[0],abs(Clfid[13]-Clfidnl[13]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[13]),label='fiducial without lensing')

py.loglog(Clbest[0],abs(Clfid[13]-Clbest[13]),label='fiducial with lensing - bestfit without lensing')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 4-4')

py.savefig('correlation_4_4.pdf')

py.close()

# 4-5

py.loglog(Clfid[0],abs(Clfid[14]-Clfidnl[14]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[14]),label='fiducial without lensing')

py.loglog(Clbest[0],abs(Clfid[14]-Clbest[14]),label='fiducial with lensing - bestfit without lensing')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 4-5')

py.savefig('correlation_4_5.pdf')

py.close()

# 5-5

py.loglog(Clfid[0],abs(Clfid[15]-Clfidnl[15]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[15]),label='fiducial without lensing')

py.loglog(Clbest[0],abs(Clfid[15]-Clbest[15]),label='fiducial with lensing - bestfit without lensing')

py.xlabel(r'$\ell$')

py.ylabel(r'$|\Delta C_\ell|$')
    
py.xlim(1,450)

py.legend(loc=0)

py.title('Correlation bins 5-5')

py.savefig('correlation_5_5.pdf')

py.close()

exit()
