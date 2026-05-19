import numpy as np
import os as os
from astropy import units as u
from astropy import constants as c
from tqdm import tqdm

au = u.au.to(u.cm)
k_b = c.k_B.cgs.value
m_p = c.m_p.cgs.value
Grav = c.G.cgs.value
Msun = c.M_sun.cgs.value
tgas = 280
aH2 = 2e-15


def stokesnr (swrm) :    
    rdis = swrm.rdis*au
    zdis = swrm.zdis*au
    size = swrm.grain_size
    indens = swrm.indens
    cs = cs_speed(rdis) #sound speed
    omegaK = omega_K(rdis) #keplerian angular frequency
    rho_g = densg(rdis, zdis) #gas density
    lmfp = 2.3*m_p/(rho_g*aH2)
    #Stokes regime
    if (size > 2.25*lmfp) :
        stokesnr = np.sqrt(2.*np.pi) * size**2. * indens * omegaK * aH2 / (9. * cs * 2.3 * m_p)
    #Epstein regime
    else :
        stokesnr = size * indens / (np.sqrt(8./np.pi) * cs * rho_g) * omegaK
    return stokesnr

def omega_K (rdis) :
    keplerfreq = np.sqrt(Grav*Msun/(rdis)**3)
    return keplerfreq
def temperature (rdis,q,t0) :
    temp = t0 * (rdis/au)**(-1*q)
    return temp
def cs_speed (rdis,q,t0) :
    cs0 = np.sqrt(k_b*t0/(2.3*m_p))
    cs = cs0 *(rdis/au)**(-0.5*q)
    return cs
def H_g (rdis,q,t0) :
    omegaK = omega_K(rdis)
    cs = cs_speed(rdis,q,t0)
    Hg = cs/omegaK
    return Hg
def sigma_g (rdis,p,sigmag0,r0=au) :
    sigmagas  = sigmag0*(rdis/r0)**(-1*p)
    return sigmagas
def densg (rdis,zdis,p,q,sigmag0,t0) :
    Hgas = H_g(rdis,q,t0)
    sigmagas = sigma_g(rdis,p,sigmag0)
    exponent = (zdis)**2/(2 * Hgas**2)
    rho_g = sigmagas / np.sqrt(2 * np.pi) / Hgas * np.exp(-1 * exponent) #gas density
    return rho_g
def pressure(rdis,zdis,p,q,sigmag0,t0) :
    rho_g = densg(rdis,zdis,p,q,sigmag0,t0)
    cs = cs_speed(rdis,q,t0)
    pg = rho_g*(cs**2)
    return pg

def dlogP(rdis, zdis):
    x = rdis
    z = zdis
    Deltar = 1.e4
    dlogPg = 0.5 * x * (pressure(x+Deltar, z) - pressure(x-Deltar,z)) / pressure(x,z) / Deltar
    return dlogPg
 
class Diskbuild():
    zh = []
    r = []
    zwalls = []
    sigmag = []
    sigmad = []
    cs = []
    Hg = []
    T= []
    rhog = []
    Pgas = []
    omegaK = []
    Stfrag = []
    Stdrift = []
    Stdf =[]
    adf = []
    afrag = []
    adrift = []
    tmix = []
    tgrowth = []
    agrid = []
    rhod = []
    stokesnr = []
    dtg = []
    vfrag = []
    alpha = []
    matdens = []
    taugrowth = []
    def __init__(self,rbins=100,zbins=100,rmin=1,rmax=100) :
        self.r = np.linspace(rmin,rmax,rbins+1)*au
        self.rcents = 0.5*(self.r[1:]+self.r[:-1])
        self.drdis = self.r[1:]-self.r[:-1]
        self.zh = np.linspace(0,5,zbins)
        self.zwalls = np.zeros(zbins)
        self.rlog = np.logspace(np.log10(rmin),np.log10(rmax),rbins+1)*au
        self.agrid = np.logspace(-4,4,100)
        self.sigmag = np.zeros(rbins)
        self.sigmad = np.zeros(rbins)
        self.cs = np.zeros(rbins)
        self.T = np.zeros(rbins)
        self.Hg = np.zeros(rbins)
        self.omegaK = np.zeros(rbins)
        self.rhog = np.zeros((zbins,rbins))
        self.Pgas = np.zeros((zbins,rbins))
        self.alpha = []
        self.vfrag = []
        self.Stfrag = np.zeros(rbins)
        self.afrag = np.zeros(rbins)
        self.Stdrift = np.zeros(rbins)
        self.Stdf = np.zeros(rbins)
        self.adrift = np.zeros(rbins)
        self.adf = np.zeros(rbins)
        self.taugrowth = np.zeros(rbins)
        self.tgrowth = np.zeros(rbins)
        self.dtg = 0.01
    def build(self,sigmag0=800.,p=1,q=0.5,t0=280.,alpha=1e-3, vfrag=1000.,matdens=1.2):
        self.alpha = alpha
        self.vfrag = vfrag
        self.matdens = matdens
        for i in range(np.size(self.rcents)) :
            rdis = self.rcents[i]
            self.omegaK[i] = omega_K(rdis)
            self.sigmag[i] = sigma_g(rdis,p,sigmag0,au)
            self.sigmad[i] = self.dtg * self.sigmag[i]
            self.cs[i]     = cs_speed(rdis,q,t0)
            self.Hg[i]     = H_g(rdis,q,t0)
            self.T[i]      = temperature(rdis,q,t0)
            
            for j in range(np.size(self.zh)) :
                zdis = self.zh[j]*self.Hg[i]
                self.rhog[j,i] = densg(rdis,zdis,p,q,sigmag0,t0)
                self.Pgas[j,i] = pressure(rdis,zdis,p,q,sigmag0,t0)
            self.taugrowth = 1/self.dtg/self.omegaK[:]
            self.tmix = 1./(1e-3 * self.omegaK[:])
            self.rhog_mid  = self.sigmag/np.sqrt(2*np.pi)/self.Hg
    def vel_vn(self,p,q,sigmag0,t0):
        vn = np.zeros_like(self.Pgas)
        for i in range(np.size(self.rcents)):
            rdis = self.rcents[i]
            for j in range(np.size(self.z)):
                zdis = self.zh[j]*Hg[i]
                
                Pg2 = pressure(rdis+1,zdis,p,q,sigmag0,t0)
                Pg1 = pressure(rdis-1,zdis,p,q,sigmag0,t0)
                vn[j,i] = 0.25 * (Pg2-Pg1) / dens[j,i] / omega[i]
        return vn
    def calc_growthbarriers(self):
        tgrowth = 1/self.dtg/self.omegaK[:]
        self.Stfrag[:]  = 1/3*self.vfrag**2/self.alpha/self.cs[:]**2
        gamma = 11/4
        self.Stdf[:] = self.vfrag*self.rcents[:]*self.omegaK[:]/self.cs[:]**2/0.5/gamma
        self.Stdrift[:] = self.dtg*(self.rcents[:]*self.omegaK[:])**2/self.cs[:]**2/gamma
        self.afrag[:] = self.Stfrag[:] * self.rhog_mid[:] * self.cs[:] * np.sqrt(8. / np.pi) / (self.matdens * self.omegaK[:])
        self.adf[:]   = self.Stdf[:] * self.rhog_mid[:] * self.cs[:] * np.sqrt(8. / np.pi) / (self.matdens * self.omegaK[:])
        self.adrift[:] = self.Stdrift[:] * self.rhog_mid[:] * self.cs[:] * np.sqrt(8./ np.pi) / (self.matdens * self.omegaK[:])
        a1 = np.minimum(self.afrag,self.adrift)
        self.tgrowth = self.taugrowth[:]*np.log(a1[:]/1.e-4)