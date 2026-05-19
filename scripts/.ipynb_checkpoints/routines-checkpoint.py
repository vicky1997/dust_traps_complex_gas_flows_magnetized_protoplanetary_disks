import numpy as np
from astropy import constants as c
from astropy import units as u
from mcdisk import densg, omega_K, temperature, cs_speed, H_g, sigma_g, pressure
import random
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as mticker
from matplotlib import animation, rc
from IPython.display import HTML
from tqdm import tqdm
import matplotlib as mpl
from matplotlib.colors import LogNorm
from mcdust import Swarm, Params, Plots, Simulation
from mcdisk import Diskbuild
from copy import deepcopy
from pathlib import Path
import seaborn as sns
from matplotlib.gridspec import GridSpec
import matplotlib as mpl 
from mcdisk import dlogP
import scipy.interpolate as interpolate
import os as os
import h5py
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['font.size'] = 14
cmap = sns.color_palette("viridis", as_cmap=True)
mpl.rcParams['axes.linewidth'] = 1.
au = u.au.to(u.cm)
year = u.year.to(u.s)
k_b = c.k_B.cgs.value
m_p = c.m_p.cgs.value
Grav = c.G.cgs.value
M_sun = c.M_sun.cgs.value
year = u.year.to(u.s)
AH2 = 2.e-15
p=1
q=0.5

class Swarm:
    def __init__(self, par):
        nswrmarray = np.full((par.ntime),1,dtype=int)
        for iout in tqdm(range(par.ntime)):
            string="%05d"%(iout)
            fname = os.path.join(par.datadir+'swarms-'+string+'.h5')
            with h5py.File(fname,"r") as f:
                dset = f['swarms/swarmsout']
                dset1 = f['times/timesout']
                swarmlist = dset[...]
                timesout = dset1[...]
                f.close()
            #tempidnr = np.array(swarmlist['id_number'])
            #tempz    = np.array(swarmlist['height_above_midplane_[AU]'])
            #tempr    = np.array(swarmlist['cylindrical_radius_[AU]'])
            #print(np.shape(swarmlist['id_number']))
            nswrm = np.size(swarmlist['id_number'])
            nswrmarray[iout] = nswrm
        maxnswrm = np.max(nswrmarray)
        self.idnr = np.zeros((par.ntime,maxnswrm))
        self.mass = np.zeros_like(self.idnr)
        self.rdis = np.zeros_like(self.idnr)
        self.zdis = np.zeros_like(self.idnr)
        self.St = np.zeros_like(self.idnr)
        self.grain_size = np.zeros_like(self.idnr)
        self.indens = np.ones_like(self.idnr)
        self.velr = np.zeros_like(self.idnr)
        self.velz = np.zeros_like(self.idnr)
        self.mswarm = None
        self.snapt = np.zeros(par.ntime)
        self.sigmad = []



    def read_hdf5_new(self, par, nancheck=True, readatt=True):
        snapshottime = np.zeros(par.ntime)
        nswrmarray = np.full((par.ntime),1,dtype=int)
        if(readatt):
            f1 = os.path.join(par.datadir+'swarms-00000.h5')
            with h5py.File(f1,"r") as f:
                a1 = f.attrs['author(s)'].decode('UTF-8')
                a2 = f.attrs['code'].decode('UTF-8')
                a3 = f.attrs['mass_of_swarm']
                a4 = f.attrs['output_number']
                print(str(a2))
                print("Author : " + str(a1))
                self.mswarm = a3
                f.close()
        print("Please cite Drążkowska, Windmark & Dullemond (2013) A&A 556, A37")
        print("Reading data ...")
        nswrmarray = np.full((par.ntime),1,dtype=int)
        for iout in tqdm(range(par.ntime)):
            string="%05d"%(iout)
            fname = os.path.join(par.datadir+'swarms-'+string+'.h5')
            with h5py.File(fname,"r") as f:
                dset = f['swarms/swarmsout']
                dset1 = f['times/timesout']
                swarmlist = dset[...]
                timesout = dset1[...]
                f.close()
            #tempidnr = np.array(swarmlist['id_number'])
            #tempz    = np.array(swarmlist['height_above_midplane_[AU]'])
            #tempr    = np.array(swarmlist['cylindrical_radius_[AU]'])
            #print(np.shape(swarmlist['id_number']))
            nswrm = np.size(swarmlist['id_number'])
            nswrmarray[iout] = nswrm
            self.idnr[iout,:nswrm] = np.array(swarmlist['id_number'])[0,:]
            self.rdis[iout,:nswrm] = np.array(swarmlist['cylindrical_radius_[AU]'])[0,:]
            self.zdis[iout,:nswrm] = np.array(swarmlist['height_above_midplane_[AU]'])[0,:]
            self.mass[iout,:nswrm] = np.array(swarmlist['mass_of_a_particle_[g]'])[0,:]
            self.St[iout,:nswrm]   = np.array(swarmlist['Stokes_number'])[0,:]
            self.velr[iout,:nswrm] = np.array(swarmlist['Radial_velocity_v_r_[cm/s]'])[0,:]
            self.velz[iout,:nswrm] = np.array(swarmlist['Vertical_velocity_v_z_[cm/s]'])[0,:]
            #for i in range(np.size(swarmlist['id_number'])):
            #    self.idnr[iout,i] = np.array(swarmlist['id_number'])[:,i]
            #    self.rdis[iout,i] = np.array(swarmlist['cylindrical_radius_[AU]'])[:,i]
            #    self.zdis[iout,i] = np.array(swarmlist['height_above_midplane_[AU]'])[:,i]
            #    self.mass[iout,i] = np.array(swarmlist['mass_of_a_particle_[g]'])[:,i]
            #    self.St[iout,i]   = np.array(swarmlist['Stokes_number'])[:,i]
            #    self.velr[iout,i] = np.array(swarmlist['Radial_velocity_v_r_[cm/s]'])[:,i]
            #    self.velz[iout,i] = np.array(swarmlist['Vertical_velocity_v_z_[cm/s]'])[:,i]
                

            snapshottime[iout]= timesout
        self.indens = self.indens*par.rho_s
        con = np.zeros_like(self.St)
        con = (0.75/np.pi/self.indens)**(1./3.)
        mass_to_one_third = np.zeros_like(con)
        mass_to_one_third = self.mass[:,:]**(1./3.)
        self.grain_size = np.multiply(mass_to_one_third,con)
        if(nancheck) :    
            if(np.isnan(self.St).any()) :
                warnings.warn("there is NaN in St")
            if(np.isnan(self.mass).any()) :
                warnings.warn("there is NaN in mass")
            if(np.isnan(self.rdis).any()) :
                warnings.warn("there is NaN in r")
            if(np.isnan(self.zdis).any()) :
                warnings.warn("there is NaN in z")
            if(np.isnan(self.velr).any()) :
                warnings.warn("there is NaN in v_r")
            if(np.isnan(self.velz).any()) :
                warnings.warn("there is NaN in v_z")    
            else :
                print('No NaN all good upto now')
        self.snapt = snapshottime
        print("Done!")



def vertical_data(sim,rloc,dr=0.1,time=0.):
    it = sim.swarms.snapt.searchsorted(time)
    rmask = (sim.swarms.rdis[it,:] <= rloc + dr) & (sim.swarms.rdis[it,:] >= rloc-dr)
    rmasked = sim.swarms.rdis[it,rmask]
    zmasked = sim.swarms.zdis[it,rmask]
    amasked = sim.swarms.grain_size[it,rmask]
    Stmasked = sim.swarms.St[it,rmask]

    argsort_z = np.argsort(zmasked)
    z_sort = zmasked[argsort_z]
    r_sort = rmasked[argsort_z]
    a_sort = amasked[argsort_z]
    St_sort = Stmasked[argsort_z]
    zh_sort = z_sort*au/H_g(rdis=rloc*au,q=0.5,t0=sim.pars.tgas)
    
    return a_sort


def extract_v_mid(file,rgrid):
    
    r_gas = np.loadtxt(file,skiprows=1,usecols=0,max_rows=300)
    z_gas = np.loadtxt(file,skiprows=1,usecols=1,max_rows=300)
    v_gas = np.loadtxt(file,skiprows=302,usecols=1,max_rows=90000)
    RR,ZZ = np.meshgrid(r_gas,z_gas)
    z_mask = (ZZ <= 0.03) & (ZZ >= -0.03)
    VV_gas = v_gas.reshape(300,300)
    r_sort= np.argsort(RR[z_mask])
    rr_mask = RR[z_mask]
    VV_mask = VV_gas[z_mask]
    rr_unique,idx_unique =np.unique(rr_mask[r_sort],return_index=True)
    VV_unique = VV_mask[r_sort][idx_unique]
    interp_func = interpolate.interp1d(rr_unique,VV_unique)
    v_gas_mid = interp_func(rgrid)
    
    return v_gas_mid

def sigmad_evol(sim,rbins=128,rmin=6.,rmax=23):
    rwalls = np.linspace(rmin,rmax, rbins+1)*au
    rcents = 0.5*(rwalls[1:]+rwalls[:-1])
    drdis = rwalls[1:]-rwalls[:-1]
    sigma = np.zeros((np.size(sim.swarms.snapt),rbins))
    sigmag = sim.pars.sigmagas*(rcents/au)**(-1)
    sigmag[sigmag<1.e-9] = 1.e-9
    for i in range(np.size(sim.swarms.snapt)):
        rdis = sim.swarms.rdis[i,:]
        rcounts,binsr = np.histogram(rdis,bins=rwalls/au)
        sigma[i,:] = rcounts[:] *sim.swarms.mswarm / (2 * np.pi * rcents*drdis)
    sigmad = {
        'rcents': rcents,
        'rwalls': rwalls,
        'sigma_d': sigma,
        'sigma_g': sigmag,
        'dr':drdis
        
    }
    return sigmad

def St_advection_mid(r,Z,vgas):
    vk = r* omega_K(r)
    #print(vgas/vk)
    #print(np.abs(vgas)/vk)
    Stadvsq = np.maximum(np.abs(vgas)/(Z*vk) - 1,0.)
    return Stadvsq

def St_drift_dustpy(r,Z):
    #gamma = 11./4.
    vK = r*omega_K(r)
    Deltar=1.e4
    dlogPg = 0.5 * r * (pressure(r+Deltar, 0,p=1,q=0.5,sigmag0=800.,t0=200.) - pressure(r-Deltar,0.,p=1,q=0.5,sigmag0=800,t0=200.)) / pressure(r,0.,p=1,q=0.5,sigmag0=800.,t0=200.) / Deltar
    gamma = np.abs(dlogPg)
    cs = cs_speed(r,q=0.5,t0=200.)
    St_drift = (Z/gamma) * (vK/cs)**2
    return St_drift

def St_drift(r,Z):
    hg = H_g(r,q=0.5,t0=200.)
    gamma = 11./4.
    St_drift = (Z/gamma) * (hg/r)**(-2)
    return St_drift

def St_frag_drift(r,vfrag):
    vk = r * omega_K(r)
    cs = cs_speed(r,q=0.5,t0=200.)
    gamma = 11./4.
    St_df = 2 * vfrag * vk/gamma/cs**2
    return St_df

def St_frag(r,vfrag,alpha):
    cs = cs_speed(r,q=0.5,t0=200.)
    Stfrag = 1/3*(vfrag**2)/alpha/(cs**2)
    return Stfrag

def St_frag_dustpy(r,vfrag,alpha):
    cs = cs_speed(r,q=0.5,t0=200.)
    b = vfrag**2 / (alpha * cs**2)
    StFr = 1 / (2 * b) * (3 - np.sqrt(9 - 4 * b**2))
    return StFr

def St_adv_frag(vfrag,vgas, N=0.5):
    Stsq = vfrag/(np.abs(vgas)*(1-(N**2)))
    return np.sqrt(Stsq)

def St_to_a(r,St,sigmag,matdens):
    Hg = H_g(r,q=0.5,t0=200.)
    cs = cs_speed(r,q=0.5,t0=200.)
    omegaK = omega_K(r)
    rhog_mid = sigmag/np.sqrt(2*np.pi)/Hg
    a = St * rhog_mid[:] * cs[:] * np.sqrt(8. / np.pi) / (matdens * omegaK[:])
    return a

def sigma2devol(sim,timesplot=[1e4,2.5e4,5e4,8e4],savefig=False,figname='',simtype='wind'):
    it=-1
    sigMax = 4.
    sigMin = -6.
    levels = np.arange(sigMin, sigMax, 1)
    sigmafloor = 1.e-15
    rbins=64
    nabins = 64
    sigma_sim = sigmad_evol(sim,rbins=rbins)
    awalls = np.logspace(np.log10(1.e-5),np.log10(np.max(1.e1)),nabins+1)
    acents = np.sqrt(awalls[1:]*awalls[:-1])
    #rwalls = np.linspace(sim.pars.minr,sim.pars.maxr,rbins+1)*au
    #rcents = 0.5*(rwalls[1:]+rwalls[:-1])
    #dr = rwalls[1:]-rwalls[:-1]
    rwalls = sigma_sim['rwalls']
    rcents = sigma_sim['rcents']
    dr = sigma_sim['dr']
    dtg = sigma_sim['sigma_d']/sigma_sim['sigma_g']
    dtg_interp = interpolate.interp1d(sigma_sim['rcents']/au,dtg)(rcents/au)
    sigmag = sim.pars.sigmagas*(rcents/au)**(-1)
    if(simtype=='wind'):
        v_gas_mid_cs = extract_v_mid('/scratch/vaikundaraman/bai2017averaged.inp',rcents/au)
        cs_mid = cs_speed(rcents,q=0.5,t0=200)
        v_gas_mid = v_gas_mid_cs*cs_mid
    if(simtype=='steadystate'):
        nu = sim.pars.alpha*cs_speed(rcents,q=0.5,t0=200.)*H_g(rcents,q=0.5,t0=200.)
        v_gas_mid = 3*nu/2/rcents
    
    Stadvsq2d = np.zeros_like(dtg_interp)
    Stadvsq2d_2 = np.zeros_like(dtg_interp)
    Stadvsq2d_2 = St_advection_mid(rcents,dtg_interp,v_gas_mid)
    Stadvsq2d = np.copy(Stadvsq2d_2)
    Stadvsq2d[...,Stadvsq2d[:,:]<0.] = 0.
    
    
    
    St_df = St_frag_drift(rcents, sim.pars.vfrag)
    St_f = St_frag(rcents,sim.pars.vfrag,sim.pars.alpha)
    St_f2 = St_frag_dustpy(rcents,sim.pars.vfrag,sim.pars.alpha)
    Stadvf = St_adv_frag(sim.pars.vfrag,v_gas_mid)
    a_frag = St_to_a(rcents,St_f,sigmag,sim.pars.rho_s)
    a_frag_2 = St_to_a(rcents,St_f2,sigmag,sim.pars.rho_s)
    a_df = St_to_a(rcents,St_df,sigmag,sim.pars.rho_s)
    a_adv_f = St_to_a(rcents,Stadvf,sigmag,sim.pars.rho_s)

    
    
    
    adens = np.zeros((rbins,nabins))   
    fig = plt.figure(figsize=(5,4),dpi=300)
    gs = GridSpec(2, 3, width_ratios=[1, 1, 0.06], height_ratios=[1,1])
    gs.update(left=0.0, right=1, wspace=0.2, hspace=0.1)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1,0])
    ax4 = fig.add_subplot(gs[1,1])
    axc = fig.add_subplot(gs[:,2])
    for axs,_time in zip([ax1,ax2,ax3,ax4],timesplot):
        axs.set_xlim(6,23.5)
        axs.set_ylim(1.e-4,1e1)
        axs.set_yscale('log')
        axs.xaxis.set_tick_params(labelsize='x-small',which='both')
        axs.yaxis.set_tick_params(labelsize='x-small',which='both')
        axs.text(6.5,2,'t = {:.0f} yr.'.format(_time),fontsize='x-small',bbox=dict(facecolor='white',alpha=0.5))
        axs.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
        axs.xaxis.set_minor_formatter(mpl.ticker.ScalarFormatter())
        it = sim.swarms.snapt.searchsorted(_time)
        counts,xedges,yedges=np.histogram2d(sim.swarms.rdis[it,:],sim.swarms.grain_size[it,:],bins=(rwalls/au,awalls))
        for i in range(np.size(rwalls)-1):
            adens[i,:] = counts[i,:]* sim.swarms.mswarm /(2.*np.pi*rcents[i]*dr[i])
        adens[adens==0] = sigmafloor
        pcf = axs.contourf(rcents/au, acents, np.log10(adens.T), levels=levels, extend="both", cmap="magma")
        axs.plot(rcents/au,a_frag,label='fragmentation limit')
        #axs.plot(rcents/au,a_frag_2,label='fragmentation limit 2')
        Std = 0.37*St_drift(rcents,dtg[it,:])
        a_drift = St_to_a(rcents,Std,sigmag,sim.pars.rho_s)
        Std2 = 0.55*St_drift_dustpy(rcents,dtg[it,:])
        ad2 = St_to_a(rcents,Std2,sigmag,sim.pars.rho_s)
        a_adv = St_to_a(rcents,np.sqrt(Stadvsq2d[it,:]),sigmag,sim.pars.rho_s)
        axs.plot(rcents/au,a_drift,label='drift limit')
        #axs.plot(rcents/au,ad2,label='drift limit 2')
        if(simtype=='wind'):
            axs.plot(rcents/au,a_adv,label='advection limit')
        #axs.plot(rcents/au,a_df,label='drift-frag limit')
        #axs.plot(rcents/au,a_adv_f,label='advection-frag')
    cbar = plt.colorbar(pcf, cax=axc)
    cbar.ax.set_ylabel(r"$\log\ \sigma_\mathrm{d}$ [g/cm²]",fontsize='small')
    ax1.set_ylabel('grain size [cm]',fontsize='x-small')
    ax3.set_ylabel('grain size [cm]',fontsize='x-small')
    ax3.set_xlabel('distance from star [au]',fontsize='x-small')
    ax4.set_xlabel('distance from star [au]',fontsize='x-small')
    ax2.yaxis.set_ticklabels([])
    ax1.xaxis.set_ticklabels([])
    ax2.xaxis.set_ticklabels([])
    ax4.yaxis.set_ticklabels([])
    ax3.legend(fontsize='xx-small',framealpha=0.5)
    cbar.ax.tick_params(labelsize='medium')
    if(savefig==True):
        plt.savefig(figname,bbox_inches='tight')
    else:
        plt.show()
        
 

def scatter_time_evol(sim,timesplot=[1e4,2.5e4,5e4,8e4],step=1,scatter=[],vmin=0,vmax=100,savefig=False,figname='',lognorm=False,cbarlabel='',cmap='viridis',size=2):
    fig = plt.figure(figsize=(5,4),dpi=300)
    gs = GridSpec(2, 3, width_ratios=[1, 1, 0.06], height_ratios=[1,1])
    gs.update(left=0.0, right=1, wspace=0.2, hspace=0.1)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1,0])
    ax4 = fig.add_subplot(gs[1,1])
    axc = fig.add_subplot(gs[:,2])
    iternum = 0
    for axs,_time in zip([ax1,ax2,ax3,ax4],timesplot):
        axs.set_xlim(6,20)
        axs.set_ylim(-3.5,3.5)
        #axs.set_xscale('log')
        axs.xaxis.set_tick_params(labelsize='x-small',which='both')
        axs.yaxis.set_tick_params(labelsize='x-small',which='both')
        axs.text(6.5,2.3,'t = {:.0f} yr.'.format(_time),fontsize='small')
        axs.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
        axs.xaxis.set_minor_formatter(mpl.ticker.ScalarFormatter())
        it = sim.swarms.snapt.searchsorted(_time)
        swarmH_g = H_g(sim.swarms.rdis[it,::step]*au,q=0.5,t0=sim.pars.tgas)
        if (lognorm==True):
            pcf = axs.scatter(sim.swarms.rdis[it,::step],sim.swarms.zdis[it,::step],c=scatter[it,::step],norm=LogNorm(vmin=vmin,vmax=vmax),cmap=cmap,s=size)
        else:
            pcf = axs.scatter(sim.swarms.rdis[it,::step],sim.swarms.zdis[it,::step],c=scatter[it,::step],vmin=vmin,vmax=vmax,cmap=cmap,s=size)
        iternum+=1
    #fig.supylabel('height above midplane [au]')
    ax1.set_ylabel('vertical height [au]',fontsize='x-small')
    ax3.set_ylabel('vertical height [au]',fontsize='x-small')
    ax3.set_xlabel('distance from central star [au]',fontsize='x-small')
    ax4.set_xlabel('distance from central star [au]',fontsize='x-small')
    ax2.yaxis.set_ticklabels([])
    ax1.xaxis.set_ticklabels([])
    ax2.xaxis.set_ticklabels([])
    ax4.yaxis.set_ticklabels([])
    cbar = plt.colorbar(pcf, cax=axc,)
    cbar.ax.tick_params(labelsize='medium')
    cbar.ax.set_ylabel(cbarlabel,fontsize='small')
    if(savefig==True):
        plt.savefig(figname,bbox_inches='tight')
    else:
        plt.show()
        
        
def vn(rdis,zdis):
    Pg2 = pressure(rdis+1,zdis,p=1,q=0.5,sigmag0=800.,t0=200.)
    Pg1 = pressure(rdis-1,zdis,p=1,q=0.5,sigmag0=800.,t0=200.)
    dens = densg(rdis,zdis,p=1,q=0.5,sigmag0=800.,t0=200.)
    omega = omega_K(rdis)
    vn = 0.25 * (Pg2-Pg1)/dens/omega
    return vn

def t_grow(r,Z,a1=0.1,a0=1.e-4):
    taugrow = 1/(omega_K(r)*Z)
    tgrow = taugrow*np.log(a1/a0)
    return tgrow/year

def t_advection(r,vgas,St):
    vadv = vgas/(1+St**2)
    tadv = r/np.abs(vadv)
    return tadv/year

def t_drift(r,St):
    v_d = 2*vn(r,0.)/(St + (St**(-1)))
    tdrift = r/np.abs(v_d)
    return tdrift/year

def plot_timescales(sim,simtype='wind',showfig=False,times=[1.e4,1.e4,3.e4,4.e4]):
    rbins=64
    sigma_sim = sigmad_evol(sim,rbins=rbins)
    dtg = sigma_sim['sigma_d']/sigma_sim['sigma_g']
    dtg_interp = interpolate.interp1d(sigma_sim['rcents']/au,dtg)(sigma_sim['rcents']/au)
    
    sigmag = sim.pars.sigmagas*(sigma_sim['rcents']/au)**(-1)
    if(simtype=='wind'):
        v_gas_mid_cs = extract_v_mid('/scratch/vaikundaraman/bai2017averaged.inp',sigma_sim['rcents']/au)
        cs_mid = cs_speed(sigma_sim['rcents'],q=0.5,t0=200)
        v_gas_mid = v_gas_mid_cs*cs_mid
    if(simtype=='steadystate'):
        nu = sim.pars.alpha*cs_speed(sigma_sim['rcents'],q=0.5,t0=200.)*H_g(sigma_sim['rcents'],q=0.5,t0=200.)
        v_gas_mid = 3*nu/2/sigma_sim['rcents']
    
    tgrow2d = np.zeros_like(dtg_interp)
    tadv = t_advection(sigma_sim['rcents'],v_gas_mid,St=np.ones_like(v_gas_mid)*1.e-2)
    a1 = St_to_a(sigma_sim['rcents'],1.e-2,sigmag,sim.pars.rho_s)
    for i in range(np.size(dtg_interp,axis=0)):
        tgrow2d[i,:] = t_grow(sigma_sim['rcents'],dtg_interp[i,:],a1=a1)
    tdrift = t_drift(sigma_sim['rcents'],St=1.e-2)
    if(showfig==True):
        #times = [10000.,25000.,50000.,80000.]
        it = sim.swarms.snapt.searchsorted(times)
        f,ax = plt.subplots(dpi=200)
        ax.plot(sigma_sim['rcents']/au,tdrift[:],'--',label='drift',c='b',linewidth=3)
        #plt.plot(r_interp,tdrift2d[::30,:].T)
        ax.plot(sigma_sim['rcents']/au,tadv,label='advection',c='purple',linewidth=3)
        ax.semilogy(sigma_sim['rcents']/au,tgrow2d[0,:].T,'-.',label='growth',c='k',linewidth=3)
        for its,timeplot in zip(it,times):
            ax.semilogy(sigma_sim['rcents']/au,tgrow2d[its,:].T,'-.',label='t = {:.0f} yr'.format(timeplot),linewidth=3)
        #ax.semilogy(sigma_sim['rcents']/au,tgrow2d[it[1],:].T,'-.',label='t = {:.0f} yr'.format(times[1]),linewidth=3)
        #ax.semilogy(sigma_sim['rcents']/au,tgrow2d[it[0],:].T,'-.',label='t = {:.0f} yr'.format(times[0]),linewidth=3)
        #ax.semilogy(sigma_sim['rcents']/au,tgrow2d[it[1],:].T,'-.',label='t = {:.0f} yr'.format(times[1]),linewidth=3)
        #ax.semilogy(sigma_sim['rcents']/au,tgrow2d[it[2],:].T,'-.',label='t = {:.0f} yr'.format(times[2]),linewidth=3)
        #ax.semilogy(sigma_sim['rcents']/au,tgrow2d[it[3],:].T,'-.',label='t = {:.0f} yr'.format(times[3]),linewidth=3)
        
        ax.legend(loc='upper left')
        ax.set_xlabel('distance from central star [au]')
        ax.set_ylabel('timescale [yr]')
        ax.set_xlim(6.5,22)
    tscales = {
            'growth': tgrow2d,
            'advection': tadv,
            'drift': tdrift
        }
        
    return sigma_sim, tscales

def returneps(sim,time,rbins=50,nzbins=50,binning='log',Stmin=1.e-3):
    if(time==-1):
        it = -1
    else:
        it = sim.swarms.snapt.searchsorted(time)
    zwalls = np.linspace(-4,4,nzbins+1)*au
    zcents = 0.5*(zwalls[1:]+zwalls[:-1])
    dzdis = zwalls[1:] - zwalls[:-1]
    if (binning=='log'):
        rwalls = np.logspace(np.log10(sim.pars.minr),np.log10(sim.pars.maxr),rbins+1)*au
        rcents = np.sqrt(rwalls[1:]*rwalls[:-1])
        drdis = rwalls[1:]-rwalls[:-1]
    else:
        rwalls = np.linspace(sim.pars.minr,sim.pars.maxr,rbins+1)*au
        rcents = 0.5*(rwalls[1:]+rwalls[:-1])
        drdis = rwalls[1:] - rwalls[:-1]
    width=7.
    zdens = np.zeros((sim.pars.ntime,rbins,nzbins))
    zdens1 = np.zeros((rbins,nzbins),dtype='float64')
    pebble_mask = sim.swarms.St[it,:]>= Stmin
    pebbles_r = sim.swarms.rdis[it,pebble_mask]
    pebbles_z = sim.swarms.zdis[it,pebble_mask]
    rho2d = np.zeros((rbins,nzbins))
            
            
    f,ax = plt.subplots(figsize=(width/2., width/2), dpi=300)
    ax.set_xlabel('distance from star [au]')
    ax.set_ylabel('particle mass [g]')
    if (binning=='log'):
        ax.set_xscale('log')
    #ax.set_yscale('log')
    counts,xedges,yedges,im=ax.hist2d(pebbles_r,pebbles_z,bins=(rwalls/au,zwalls/au),norm=LogNorm(),cmap='magma_r')
    #counts,xedges,yedges,im=ax.hist2d(sim.swarms.rdis[it,:],sim.swarms.zdis[it,:],bins=(rwalls/au,zwalls/au),norm=LogNorm(),cmap='magma_r')
    for i in range(rbins):
        for j in range(nzbins):
            rho2d[i,j] = densg(rcents[i],zcents[j],p=1,q=0.5, sigmag0=800.,t0=200.)
            zdens1[i,j] = counts[i,j]* sim.swarms.mswarm /(2.*np.pi*rcents[i]*drdis[i]*dzdis[j])
    f.colorbar(im)
    time = int(sim.swarms.snapt[it])
    ax.set_title('t = {:.1f} yr'.format(sim.swarms.snapt[it]))
    plt.close()
    rho_dict = {
        'rho_d' : zdens1,
        'rho_g' : rho2d,
        'zcents' : zcents,
        'rcents' : rcents
    }
    return rho_dict

def plot_dtgevol(sim,sigma_sim,text,times=[0,10000,20000,30000,40000,50000,60000,70000,80000]):
    f,axs = plt.subplots(dpi=150)
    its = sim.swarms.snapt.searchsorted(times)
    axs.set_xlim(6.2,23.5)
    axs.set_ylim(1.e-3,0.1)
    axs.set_yscale('log')
    axs.text(6.5,0.07,text,fontsize='small')
    axs.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    axs.xaxis.set_minor_formatter(mpl.ticker.ScalarFormatter())
    axs.set_ylabel('dust-to-gas ratio $Z$',fontsize='small')
    axs.set_xlabel('distance from central star [au]',fontsize='small')
    lw=3
    for it,time in zip(its,times):
        axs.semilogy(sigma_sim['rcents']/au,sigma_sim['sigma_d' ][it,:].T/sigma_sim['sigma_g'],'-',label= 't = {:.0f} yr'.format(time),linewidth=lw)
    axs.legend(loc='upper right',framealpha=0.5,fontsize='x-small')    

    
 
    





