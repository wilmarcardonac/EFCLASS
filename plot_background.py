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

#BACKGROUND
z,propertime,conformaltime,Hubbleoverc,comovingdistance,angdiadist,lumdist,comovsndhrz,rhog,rhob,rhocdm,rhoncdm,pncdm,pseudopncdm,rhofld,rhour,rhocrit = np.loadtxt('./output/Cl_lcdm_nu_background.dat',unpack=True,usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])

z1,propertime1,conformaltime1,Hubbleoverc1,comovingdistance1,angdiadist1,lumdist1,comovsndhrz1,rhog1,rhob1,rhocdm1,rhoncdm1,pncdm1,pseudopncdm1,rhofld1,rhour1,rhocrit1 = np.loadtxt('./output/Cl_lcdm_nu_CLASS_background.dat',unpack=True,usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])

z2,propertime2,conformaltime2,Hubbleoverc2,comovingdistance2,angdiadist2,lumdist2,comovsndhrz2,rhog2,rhob2,rhocdm2,rhoncdm2,pncdm2,pseudopncdm2,rhofld2,rhour2,rhocrit2 = np.loadtxt('./output/Cl_lcdm_nu_test_background.dat',unpack=True,usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])

#a,w,wprime,H,cs2,ceff2,Omegam,OmegaDE,dprho,dm,Vm,dde,Vde,pi,GeffGN,Qeff = np.loadtxt('../output/functions.txt',unpack=True,usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])

fig = py.figure()

#py.loglog(1./(1.+z),Hubbleoverc/Hubbleoverc[-1],color='blue')

py.loglog(1./(1.+z),rhoncdm,color='red',label=r'$\rho_{ncdm}$')

py.loglog(1./(1.+z),pncdm,color='red',label=r'$p_{ncdm}$',ls='dashed')

py.loglog(1./(1.+z),rhog,color='blue',label=r'$\rho_{rad}$')

py.loglog(1./(1.+z),rhocdm,color='green',label=r'$\rho_{cdm}$')

py.loglog(1./(1.+z),rhob,color='black',label=r'$\rho_{b}$')

py.loglog(1./(1.+z),rhofld,color='orange',label=r'$\rho_{\Lambda}$')

py.loglog(1./(1.+z),rhour,color='purple',label=r'$\rho_{ur}$')

py.xlabel(r'$a$',fontsize='large')

py.ylabel(r'$8 \pi G \rho /3$',fontsize='large')

py.legend(loc=0)

py.savefig('./output/rhoncdm.pdf')

py.close()

fig = py.figure()

py.plot(1./(1.+z),pncdm/rhoncdm,color='red',label=r'$w_{ncdm}$')

py.hlines(0.,1.e-14,1.,color='black',label=r'$w_{cdm}$',linestyles='dotted')

py.hlines(1./3.,1.e-14,1.,color='blue',label=r'$w_{rad}$',linestyles='dashed')

py.xlabel(r'$a$',fontsize='large')

py.ylabel(r'$w_{ncdm}$',fontsize='large')

py.xscale('log')

py.ylim(-1.,1.)

py.legend(loc=0)

py.savefig('./output/wncdm.pdf')

py.close()

fig = py.figure()

#py.plot(1./(1.+z),(rhoncdm-rhoncdm1)/rhoncdm1*1.e2,color='red',label=r'$\rho_{ncdm}$')

#py.plot(1./(1.+z),(rhoncdm2-rhoncdm1)/rhoncdm2*1.e2,color='red',label=r'$\rho_{ncdm}$',ls='dotted')

#py.plot(1./(1.+z),(pncdm-pncdm1)/pncdm1*1.e2,color='blue',label=r'$p_{ncdm}$')

#py.plot(1./(1.+z),(pncdm2-pncdm1)/pncdm2*1.e2,color='blue',label=r'$p_{ncdm}$',ls='dashed')

#py.plot(1./(1.+z),(rhour-rhour1)/rhour1*1.e2,color='purple',label=r'$\rho_{ur}$',ls='dotted')

py.plot(1./(1.+z),(pseudopncdm-pseudopncdm1)/pseudopncdm1*1.e2,color='green',label=r'$pseudo_{ncdm}$')

#py.plot(1./(1.+z),(pseudopncdm2-pseudopncdm1)/pseudopncdm2*1.e2,color='green',label=r'$pseudo_{ncdm}$')

py.xlabel(r'$a$',fontsize='large')

py.ylabel(r'$\%$',fontsize='large')

py.xscale('log')

py.legend(loc=0)

py.savefig('./output/difference_rhoncdm.pdf')

py.close()

exit()

# Error file
#Clbest = np.loadtxt('../data/Cl_bestfit_no_lensing_cl.dat',unpack=True)

# Fiducial model without lensing
#Clfidnl = np.loadtxt('../data/Cl_fiducial_no_lensing_cl.dat',unpack=True) 

# 1-1

#py.loglog(Clfid[0],Clfid[1],label=r'$\Lambda CDM$')
py.plot(Cl[0],Cl[1]*cambunits,label=r'$designer$',c='green')
py.plot(Cl1[0],Cl1[1]*cambunits,label=r'$\Lambda CDM$',c='blue')#'r'$fr0 = -0.06$')
#py.plot(Cl1[0],Cl1[1]*cambunits,label=r'Hu-Sawicki')#'r'$fr0 = -0.06$')
#py.plot(Cl2[0],Cl2[1]*cambunits,label=r'Hu-Sawicki old')#'r'$fr0 = -0.06$')
#py.loglog(Clfid[0],abs(Clfid[1]-Clfidnl[1]),label='fiducial with lensing - fiducial without lensing')

#py.loglog(Clfidnl[0],abs(Clfidnl[1]),label='fiducial without lensing')

#py.loglog(Clbest[0],abs(Clfid[1]-Clbest[1]),label='fiducial with lensing - bestfit without lensing')

py.xlabel(r'$\ell$')

py.ylabel(r'$\ell*(1+\ell) C_\ell/2\pi$')

py.xscale('log')

#py.xlim(1,30)

#py.ylim(600,1200)
py.legend(loc=0)

#py.title('Correlation bins 1-1')

py.savefig('./output/CMB_cl.pdf')
#py.savefig('correlation_1_1.pdf')
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
