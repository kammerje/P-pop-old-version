"""
P_pop, a Python package to simulate exoplanet populations based on Kepler
statistics. Information about the method can be found in this paper:
https://ui.adsabs.harvard.edu/abs/2018A%26A...609A...4K/abstract. This library
is maintained on GitHub at https://github.com/kammerje/P_pop.

===============================================================================
Please cite Kammerer & Quanz 2018 (2018A&A...609A...4K) if you use P_pop for
your research. Thank you!
===============================================================================

This repository makes use of "forecaster"
(https://github.com/chenjj2/forecaster.git, Chen & Kipping 2017,
2017ApJ...834...17C) which must be cloned into the folder called "forecaster".

Author: Jens Kammerer
Version: 3.1.0
Last edited: 30.07.19
"""


# PREAMBLE
#==============================================================================
import matplotlib.pyplot as plt
import numpy as np

import scipy.integrate as si
import scipy.interpolate as interpolate
import scipy.optimize as so

import sys
sys.path.append('/Volumes/OZ 1/_p_pop/forecaster/forecaster')
sys.path.append('/home/kjens/Downloads/p_pop/forecaster/forecaster')
import mr_forecast as mr


# MAIN
#==============================================================================
class system():
    
    def __init__(self,
                 stype='G',
                 dist=10,
                 Rs=1,
                 Ts=5772,
                 Ms=1,
                 model=0,
                 rates='sag13',
                 albedo_model='uniform'):
        """
        Initializes system instance
        """
        
        # PARAMETERS FROM INPUT
        # Host star spectral type
        self.stype = str(stype)
        # Host star distance from Earth
        self.dist = float(dist)
        # Host star radius
        self.Rs = float(Rs)
        # Host star effective temperature
        self.Ts = float(Ts)
        # Host star mass
        self.Ms = float(Ms)
        # Planet occurrence rate model
        self.model = int(model)
        # Planet occurrence rates
        self.rates = str(rates)
        
        # RANDOMLY DRAWN PARAMETERS
        # System orbital inclination
        self.inc = np.arccos(2.*np.random.rand()-1.)
        # System exozodiacal level
        mu = np.log(26.); zeta = 1.2 # clean Sun-like (Ertel et al. 2018)
#        mu = np.log(13.); zeta = 1.5 # clean (Ertel et al. 2018)
        self.zodis = np.random.lognormal(mu, zeta)
        
        # Draw planet radii and planet orbital periods from Kepler statistics
        Stats = stats(stype=self.stype,
                      model=self.model)
        if (self.rates == 'sag13'):
            dummy = Stats.sag13()
        elif (self.rates == 'hsu2019'):
            dummy = Stats.hsu2019()
        elif (self.rates == 'dressing2015'):
            dummy = Stats.dressing2015()
        elif (self.rates == 'fressin2013'):
            dummy = Stats.fressin2013()
        elif (self.rates == 'sag13_dressing2015'):
            dummy = Stats.sag13_dressing2015()
        elif (self.rates == 'hsu2019_dressing2015'):
            dummy = Stats.hsu2019_dressing2015()
        else:
            raise UserWarning('Planet occurrence rates '+self.rates+' are not known')
        
        self.nplanets = len(dummy[0])
        self.planets = []
        for i in range(self.nplanets):
            self.planets += [exoplanet(Rp=dummy[0][i],
                                       Porb=dummy[1][i],
                                       inc=self.inc,
                                       zodis=self.zodis,
                                       stype=self.stype,
                                       dist=self.dist,
                                       Rs=self.Rs,
                                       Ts=self.Ts,
                                       Ms=self.Ms,
                                       albedo_model=albedo_model)]
        
        pass

class exoplanet():
    
    def __init__(self,
                 Rp=1,
                 Porb=365,
                 inc=0,
                 zodis=1,
                 stype='G',
                 dist=10,
                 Rs=1,
                 Ts=5772,
                 Ms=1,
                 albedo_model='uniform'):
        """
        Initializes exoplanet instance
        """
        
        # PARAMETERS FROM INPUT
        # Planet radius
        self.Rp = float(Rp)
        # Planet orbital period
        self.Porb = float(Porb)
        # System orbital inclination
        self.inc = float(inc)
        # System exozodiacal level
        self.zodis = float(zodis)
        # Host star spectral type
        self.stype = str(stype)
        # Host star distance from Earth
        self.dist = float(dist)
        # Host star radius
        self.Rs = float(Rs)
        # Host star effective temperature
        self.Ts = float(Ts)
        # Host star mass
        self.Ms = float(Ms)
        
        # RANDOMLY DRAWN PARAMETERS
        # Planet orbital eccentricity
        self.ecc = 0.
        # Planet Bond albedo
        if (albedo_model == 'uniform'):
            self.Abond = 0.8*np.random.rand()
        elif (albedo_model == 'cahoy2010'):
            self.Abond = 0.8*np.random.rand()
        else:
            raise UserWarning('Albedo model must be either uniform or cahoy2010')
        # Planet geometric albedo
        if (albedo_model == 'uniform'):
            self.AgeomMIR = 0.1*np.random.rand()
            self.AgeomVIS = 0.6*np.random.rand()
        elif (albedo_model == 'cahoy2010'):
            self.AgeomMIR = 0.1*np.random.rand()
            self.AgeomVIS = self.__AgeomVIS_cahoy2010()
        else:
            raise UserWarning('Albedo model must be either uniform or cahoy2010')
        # Planet longitude of ascending node
        self.Omega = 2.*np.pi*np.random.rand()
#        self.Omega = np.pi/2. # quadrature
        # Planet argument of periapsis
        self.omega = 2.*np.pi*np.random.rand()
#        self.omega = 0. # quadrature
        # Planet true anomaly
        self.theta = self.__theta()
#        self.theta = 0. # quadrature
        # Planet mass
#        Mmedian, Mplus, Mminus = mr.Rstat2M(mean=self.Rp, std=1E-10, unit='Earth', sample_size=1e3, grid_size=1e3)
        self.Mp = mr.Rpost2M(radius=[self.Rp], unit='Earth', grid_size=1e3, classify='No')
        
        pass
    
    def __AgeomVIS_cahoy2010(self):
        """
        Computes randomly distributed planet geometric albedo (VIS)
        
        Credit: EXOSIMS (https://github.com/dsavransky/EXOSIMS)
        """
        
        # Table 4 of Cahoy et al. 2010
        sas = [0.8, 2., 5., 10.]
        fes = [1., 3., 10., 30.]
        pts = np.array([[0.322, 0.241, 0.209, 0.142], [0.742, 0.766, 0.728, 0.674], [0.567, 0.506, 0.326, 0.303], [0.386, 0.260, 0.295, 0.279]])
        
        grid_sa, grid_fe = np.meshgrid(sas, fes)
        albedo_grid = np.vstack((grid_sa.flatten(), grid_fe.flatten())).T
        albedo_vals = pts.T.flatten()
        
        # Planet semi-major axis & metallicity
        sa_temp = np.clip(self.a(), 0.8, 10.)
        fe_temp = np.random.uniform(low=1., high=30.)
        
        p = interpolate.griddata(albedo_grid, albedo_vals, (sa_temp, fe_temp), method='cubic')
        
        # Return randomly distributed planet geometric albedo (VIS)
        return p
    
    def __theta(self):
        """
        Computes randomly distributed true anomaly (in radians)
        """
        
        # Draw randomly distributed mean anomaly
        M = 2.*np.pi*np.random.rand()
        
        # Calculate eccentric anomaly
        func = lambda x: x-self.ecc*np.sin(x)-M
        E = so.fsolve(func, M)[0]
        
        # Calculate true anomaly
        func = lambda x: (self.ecc+np.cos(x))/(1.+self.ecc*np.cos(x))-np.cos(E)
        
        # Return randomly distributed true anomaly
        return so.fsolve(func, E)[0] % (2.*np.pi)
    
    def a(self):
        """
        Computes planet semi-major axis (in astronomical units)
        """
        
        # Define constants
        G = 6.67408E-11
        Msun = 1.989E+30
        au = 149597870700.
        
        # Return planet semi-major axis
        return ((G*self.Ms*Msun*(self.Porb*86400.)**2.)/(4.*np.pi**2.))**(1./3.)/au
    
    def alpha(self):
        """
        Computes planet phase angle (in radians)
        """
        
        # Return planet phase angle
        return np.arccos(-np.sin(self.inc)*np.sin(self.omega+self.theta))
    
    def ang_sep(self):
        """
        Computes angular planet host star separation (in arcseconds)
        """
        
        # Return angular planet host star separation
        return self.rp_proj()/self.dist
    
    def ang_sep_max(self):
        """
        Computes maximal angular planet host star separation (in arcseconds)
        """
        
        #
        a_temp = self.a()
        
        def ftheta(x,
                   a_temp):
            
            # Compute projected planet host star separation
            return a_temp*(1.-self.ecc**2.)/(1.+self.ecc*np.cos(x))*np.sqrt(np.cos(self.omega+x)**2.+np.cos(self.inc)**2.*np.sin(self.omega+x)**2.)
        
        def dfdtheta(x,
                     a_temp):
            
            # Compute derivative of projected planet host star separation with respect to the true anomaly
            return a_temp*(1.-self.ecc**2.)*(-1.)*(1.+self.ecc*np.cos(x))**(-2.)*(-self.ecc*np.sin(x))*np.sqrt(np.cos(self.omega+x)**2.+np.cos(self.inc)**2.*np.sin(self.omega+x)**2.)+a_temp*(1.-self.ecc**2.)/(1.+self.ecc*np.cos(x))*0.5*(np.cos(self.omega+x)**2.+np.cos(self.inc)**2.*np.sin(self.omega+x)**2.)**(-0.5)*(2.*np.cos(self.omega+x)*(-np.sin(self.omega+x))+np.cos(self.inc)**2.*2.*np.sin(self.omega+x)*np.cos(self.omega+x))
        
        # Compute roots of derivative of projected planet host star separation with respect to the true anomaly
        xx = np.linspace(0., 2.*np.pi, 8, endpoint=False)
        x0 = so.fsolve(dfdtheta, xx, args=(a_temp))
        f0 = ftheta(x0, a_temp)
        
        # Compute maximal angular planet host star separation
        theta_buffer = self.theta.copy()
        self.theta = x0[np.argmax(f0)] % (2.*np.pi)
        ang_sep_max = self.rp_proj()/self.dist
        self.theta = theta_buffer
        
        # Return maximal angular planet host star separation
        return ang_sep_max
    
    def f(self):
        """
        Computes Lambertian reflectance
        """
        
        #
        alpha = self.alpha()
        
        # Return Lambertian reflectance
        return np.abs((np.sin(alpha)+(np.pi-alpha)*np.cos(alpha))/np.pi)
    
    def Finc(self):
        """
        Computes incoming stellar flux (in solar constants)
        """
        
        # Define constants
        sigma = 5.6704E-8
        Rsun = 695700000.
        au = 149597870700.
        
        # Return incoming stellar flux
        return sigma*self.Ts**4.*(self.Rs*Rsun)**2./(self.rp()*au)**2./1361.
    
    def rp(self):
        """
        Computes planet host star separation (in astronomical units)
        """
        
        # Return planet host star separation
        return self.a()*(1.-self.ecc**2.)/(1.+self.ecc*np.cos(self.theta))
    
    def rp_proj(self):
        """
        Computes projected planet host star separation (in astronomical units)
        """
        
        # Return projected planet host star separation
        return self.rp()*np.sqrt(np.cos(self.omega+self.theta)**2.+np.cos(self.inc)**2.*np.sin(self.omega+self.theta)**2.)
    
    def Tp(self):
        """
        Computes planet equilibrium temperature (in Kelvins)
        """
        
        # Define constants
        Rsun = 695700000.
        au = 149597870700.
        
        # Return planet equilibrium temperature
        return ((self.Rs*Rsun)**2.*(1.-self.Abond)/(4*(self.rp()*au)**2.))**(1./4.)*self.Ts

class stats():
    
    def __init__(self,
                 stype='G',
                 model=0):
        """
        Initializes stats instance
        """
        
        # Host star type
        self.stype = str(stype)
        # Planet occurrence rate model
        self.model = int(model)
        
        pass
    
    def sag13(self,
              Rp_range=[0.5, 16.0],
              Porb_range=[0.5, 500.0]):
        """
        Draws planet radii and planet orbital periods according to SAG13
        statistics
        (https://ui.adsabs.harvard.edu/abs/2018ApJ...856..122K/abstract)
        """
        
        # Create lists for planet radii and planet orbital periods
        Rp_sag13 = []
        Porb_sag13 = []
        
        # Check if host star type suits
        if (True):
            
            # Use baseline, min or max planet occurrence rates
            if (self.model == 0):
                Gamma = [0.38, 0.73]
                alpha = [-0.19, -1.18]
                beta = [0.26, 0.59]
                F0 = [2.42, 0.25]
            elif (self.model == 1):
                Gamma = [0.138, 0.72]
                alpha = [0.277, -1.56]
                beta = [0.204, 0.51]
                F0 = [1.14, 0.14]
            elif (self.model == 2):
                Gamma = [1.06, 0.78]
                alpha = [-0.68, -0.82]
                beta = [0.32, 0.67]
                F0 = [5.60, 0.46]
            else:
                raise UserWarning('Model must be 0, 1 or 2')
            Rp_range=[0.5, 16.0]
            Porb_range=[0.5, 500.0]
            Rbrk = [0., 3.4, np.inf]
            ytod = 365.24
            
#            # Compute total number of planets per star for the given planet
#            # radius and planet orbital period ranges
#            dlogR = 0.001
#            dlogP = 0.001
#            logR = np.arange(np.log(Rp_range[0])+dlogR/2., np.log(Rp_range[1])+dlogR/2., dlogR)
#            logP = np.arange(np.log(Porb_range[0]/ytod)+dlogP/2., np.log(Porb_range[1]/ytod)+dlogP/2., dlogP)
#            R = np.exp(logR)
#            P = np.exp(logP)
#            N = 0.
#            for i in range(2):
#                for j in range(len(logR)):
#                    for k in range(len(logP)):
#                        if (Rbrk[i] <= R[j] < Rbrk[i+1]):
#                            N += Gamma[i]*R[j]**alpha[i]*P[k]**beta[i]
#            N *= dlogR*dlogP
#            import pdb; pdb.set_trace()
            
            #
            CR0 = 1./(Rbrk[1]**alpha[0]/alpha[0]-Rp_range[0]**alpha[0]/alpha[0])
            CP0 = 1./((Porb_range[1]/ytod)**beta[0]/beta[0]-(Porb_range[0]/ytod)**beta[0]/beta[0])
            
            def iCDF_R0(x):
                return (alpha[0]*x/CR0+Rp_range[0]**alpha[0])**(1./alpha[0])
            def iCDF_P0(x):
                return (beta[0]*x/CP0+(Porb_range[0]/ytod)**beta[0])**(1./beta[0])
            
            # Distribute number of planets per star Poissonian
            for k in range(np.random.poisson(F0[0])):
                dummy1 = iCDF_R0(np.random.rand())
                dummy2 = iCDF_P0(np.random.rand())*ytod
                
                # Check if planet radius and planet orbital period lie inside range
                if (float(Rp_range[0]) <= dummy1 <= float(Rp_range[1])):
                    if (float(Porb_range[0]) <= dummy2 <= float(Porb_range[1])):
                        
                        # Fill lists for planet radii and planet orbital periods
                        Rp_sag13 += [dummy1]
                        Porb_sag13 += [dummy2]
            
            CR1 = 1./(Rp_range[1]**alpha[1]/alpha[1]-Rbrk[1]**alpha[1]/alpha[1])
            CP1 = 1./((Porb_range[1]/ytod)**beta[1]/beta[1]-(Porb_range[0]/ytod)**beta[1]/beta[1])
            
            def iCDF_R1(x):
                return (alpha[1]*x/CR1+Rbrk[1]**alpha[1])**(1./alpha[1])
            def iCDF_P1(x):
                return (beta[1]*x/CP1+(Porb_range[0]/ytod)**beta[1])**(1./beta[1])
            
            # Distribute number of planets per star Poissonian
            for k in range(np.random.poisson(F0[1])):
                dummy1 = iCDF_R1(np.random.rand())
                dummy2 = iCDF_P1(np.random.rand())*ytod
                
                # Check if planet radius and planet orbital period lie inside range
                if (float(Rp_range[0]) <= dummy1 <= float(Rp_range[1])):
                    if (float(Porb_range[0]) <= dummy2 <= float(Porb_range[1])):
                        
                        # Fill lists for planet radii and planet orbital periods
                        Rp_sag13 += [dummy1]
                        Porb_sag13 += [dummy2]
        
        # Return lists with planet radii and planet orbital periods
        return Rp_sag13, Porb_sag13
    
    def hsu2019(self,
                Rp_range=[0.5, 16.0],
                Porb_range=[0.5, 500.0]):
        """
        Draws planet radii and planet orbital periods according to Hsu et al.
        2019 statistics
        (https://ui.adsabs.harvard.edu/abs/2019arXiv190201417H/abstract)
        """
        
        # Create lists for planet radii and planet orbital periods
        Rp_hsu2019 = []
        Porb_hsu2019 = []
        
        # Check if host star type suits
        if (True):
            
            # Parameters from Hsu 2019
            bins_Rp = np.log(np.array([0.5, 0.75, 1., 1.25, 1.5, 1.75, 2., 2.5, 3., 4., 6., 8., 12., 16.]))
            bins_Porb = np.log(np.array([0.5, 1., 2., 4., 8., 16., 32., 64., 128., 256., 500.]))
#            rates = np.array([[0.35, 0.31, 2.6, 7.77, 8.33, 0., 0., 0., 0., 0.], [0.206, 0.21, 1.01, 3.75, 5., 8.29, 5.2, 0., 0., 0.], [0.27, 0.136, 0.27, 2.46, 3.77, 4.8, 6., 11.2, 0., 0.], [0.167, 0.244, 0.917, 1.94, 3.24, 2.5, 1.1, 2.23, 0., 0.], [0.079, 0.085, 0.624, 1.37, 2.17, 2.83, 2.1, 2.1, 3.24, 18.], [0.025, 0.0591, 0.24, 0.7, 1.11, 0.58, 1.27, 1.3, 2.2, 17.], [0., 0.034, 0.207, 1.03, 3.31, 5.5, 5.37, 6.48, 5.6, 5.4], [0., 0.035, 0.093, 1.08, 2.52, 3.74, 6.19, 3.72, 7.78, 7.55], [0., 0., 0.091, 0.3, 0.74, 2.21, 2.25, 3.6, 2.9, 5.1], [0.033, 0., 0.0644, 0.14, 0.34, 0.33, 0.374, 1.51, 0.92, 2.6], [0., 0., 0., 0.0594, 0.087, 0.24, 0.61, 0.51, 1.1, 0.], [0., 0.025, 0.092, 0.12, 0.165, 0.21, 0.317, 0.93, 1.8, 2.], [0., 0.031, 0.183, 0.13, 0.091, 0.23, 0.501, 0.453, 0., 0.]])*1E-2
#            err_u = np.array([[0.25, 0.35, 1.3, 2.3, 5.1, 0., 0., 0., 0., 0.], [0.13, 0.17, 0.51, 1.3, 2.1, 5.7, 3.3, 0., 0., 0.], [0.21, 0.12, 0.35, 0.65, 1.9, 2.5, 4., 5.8, 0., 0.], [0.12, 0.1, 0.28, 0.44, 1.1, 1.3, 1., 6.6, 0., 0.], [0.087, 0.1, 0.22, 0.77, 0.94, 1.2, 1.1, 1.8, 2.6, 13.], [0.041, 0.042, 0.17, 0.6, 0.74, 0.52, 0.77, 2., 2.6, 9.8], [0., 0.029, 0.11, 0.43, 0.64, 0.9, 1.5, 2.1, 3.1, 6.2], [0., 0.022, 0.11, 0.32, 0.73, 0.82, 1.2, 1.4, 2.7, 4.4], [0., 0., 0.089, 0.19, 0.45, 0.7, 0.85, 1.1, 1.7, 3.2], [0.031, 0., 0.048, 0.12, 0.2, 0.2, 0.26, 0.76, 0.75, 1.7], [0., 0., 0., 0.045, 0.081, 0.19, 0.27, 0.39, 1.1, 0.], [0., 0.018, 0.064, 0.094, 0.098, 0.16, 0.39, 0.56, 1., 0.92], [0., 0.027, 0.073, 0.11, 0.053, 0.16, 0.3, 0.32, 0., 0.]])*1E-2
#            err_l = np.array([[0.25, 0.35, 1.3, 2.3, 5.1, 0., 0., 0., 0., 0.], [0.13, 0.17, 0.51, 1.3, 2.1, 5.7, 3.3, 0., 0., 0.], [0.21, 0.12, 0.35, 0.65, 1.9, 2.5, 4., 5.8, 0., 0.], [0.12, 0.1, 0.28, 0.44, 1.1, 1.3, 1., 6.6, 0., 0.], [0.087, 0.1, 0.22, 0.77, 0.94, 1.2, 1.1, 1.8, 2.6, 13.], [0.041, 0.042, 0.17, 0.6, 0.74, 0.52, 0.77, 2., 2.6, 9.8], [0., 0.029, 0.11, 0.43, 0.64, 0.9, 1.5, 2.1, 3.1, 6.2], [0., 0.022, 0.11, 0.32, 0.73, 0.82, 1.2, 1.4, 2.7, 4.4], [0., 0., 0.089, 0.19, 0.45, 0.7, 0.85, 1.1, 1.7, 3.2], [0.031, 0., 0.048, 0.12, 0.2, 0.2, 0.26, 0.76, 0.75, 1.7], [0., 0., 0., 0.045, 0.081, 0.19, 0.27, 0.39, 1.1, 0.], [0., 0.018, 0.064, 0.094, 0.098, 0.16, 0.39, 0.56, 1., 0.92], [0., 0.027, 0.073, 0.11, 0.053, 0.16, 0.3, 0.32, 0., 0.]])*1E-2
            rates_extrap = np.array([[0.003500, 0.003100, 0.026000, 0.077700, 0.083300, 0.127318, 0.205925, 0.306147, 0.418360, 0.523824], [0.002060, 0.002100, 0.010100, 0.037500, 0.050000, 0.082900, 0.052000, 0.161206, 0.219065, 0.272785], [0.002700, 0.001360, 0.002700, 0.024600, 0.037700, 0.048000, 0.060000, 0.112000, 0.140519, 0.174270], [0.001670, 0.002440, 0.009170, 0.019400, 0.032400, 0.025000, 0.011000, 0.022300, 0.100756, 0.124555], [0.000790, 0.000850, 0.006240, 0.013700, 0.021700, 0.028300, 0.021000, 0.021000, 0.032400, 0.180000], [0.000250, 0.000591, 0.002400, 0.007000, 0.011100, 0.005800, 0.012700, 0.013000, 0.022000, 0.170000], [0.000445, 0.000340, 0.002070, 0.010300, 0.033100, 0.055000, 0.053700, 0.064800, 0.056000, 0.054000], [0.000348, 0.000350, 0.000930, 0.010800, 0.025200, 0.037400, 0.061900, 0.037200, 0.077800, 0.075500], [0.000267, 0.000640, 0.000910, 0.003000, 0.007400, 0.022100, 0.022500, 0.036000, 0.029000, 0.051000], [0.000330, 0.000450, 0.000644, 0.001400, 0.003400, 0.003300, 0.003740, 0.015100, 0.009200, 0.026000], [0.000141, 0.000333, 0.000726, 0.000594, 0.000870, 0.002400, 0.006100, 0.005100, 0.011000, 0.015945], [0.000110, 0.000250, 0.000920, 0.001200, 0.001650, 0.002100, 0.003170, 0.009300, 0.018000, 0.020000], [0.000091, 0.000310, 0.001830, 0.001300, 0.000910, 0.002300, 0.005010, 0.004530, 0.007822, 0.009321]])
            err_u_extrap = np.array([[0.002500, 0.003500, 0.013000, 0.023000, 0.051000, 0.070336, 0.125336, 0.217702, 0.368586, 0.603218], [0.001300, 0.001700, 0.005100, 0.013000, 0.021000, 0.057000, 0.033000, 0.103845, 0.172654, 0.277563], [0.002100, 0.001200, 0.003500, 0.006500, 0.019000, 0.025000, 0.040000, 0.058000, 0.103413, 0.164077], [0.001200, 0.001000, 0.002800, 0.004400, 0.011000, 0.013000, 0.010000, 0.066000, 0.070802, 0.111169], [0.000870, 0.001000, 0.002200, 0.007700, 0.009400, 0.012000, 0.011000, 0.018000, 0.026000, 0.130000], [0.000410, 0.000420, 0.001700, 0.006000, 0.007400, 0.005200, 0.007700, 0.020000, 0.026000, 0.098000], [0.000391, 0.000290, 0.001100, 0.004300, 0.006400, 0.009000, 0.015000, 0.021000, 0.031000, 0.062000], [0.000317, 0.000220, 0.001100, 0.003200, 0.007300, 0.008200, 0.012000, 0.014000, 0.027000, 0.044000], [0.000257, 0.000474, 0.000890, 0.001900, 0.004500, 0.007000, 0.008500, 0.011000, 0.017000, 0.032000], [0.000310, 0.000362, 0.000480, 0.001200, 0.002000, 0.002000, 0.002600, 0.007600, 0.007500, 0.017000], [0.000167, 0.000298, 0.000517, 0.000450, 0.000810, 0.001900, 0.002700, 0.003900, 0.011000, 0.012052], [0.000150, 0.000180, 0.000640, 0.000940, 0.000980, 0.001600, 0.003900, 0.005600, 0.010000, 0.009200], [0.000145, 0.000270, 0.000730, 0.001100, 0.000530, 0.001600, 0.003000, 0.003200, 0.005425, 0.007560]])
            err_l_extrap = np.array([[0.002500, 0.003100, 0.013000, 0.023000, 0.051000, 0.070336, 0.125336, 0.217702, 0.368586, 0.523824], [0.001300, 0.001700, 0.005100, 0.013000, 0.021000, 0.057000, 0.033000, 0.103845, 0.172654, 0.272785], [0.002100, 0.001200, 0.002700, 0.006500, 0.019000, 0.025000, 0.040000, 0.058000, 0.103413, 0.164077], [0.001200, 0.001000, 0.002800, 0.004400, 0.011000, 0.013000, 0.010000, 0.022300, 0.070802, 0.111169], [0.000790, 0.000850, 0.002200, 0.007700, 0.009400, 0.012000, 0.011000, 0.018000, 0.026000, 0.130000], [0.000250, 0.000420, 0.001700, 0.006000, 0.007400, 0.005200, 0.007700, 0.013000, 0.022000, 0.098000], [0.000391, 0.000290, 0.001100, 0.004300, 0.006400, 0.009000, 0.015000, 0.021000, 0.031000, 0.054000], [0.000317, 0.000220, 0.000930, 0.003200, 0.007300, 0.008200, 0.012000, 0.014000, 0.027000, 0.044000], [0.000257, 0.000474, 0.000890, 0.001900, 0.004500, 0.007000, 0.008500, 0.011000, 0.017000, 0.032000], [0.000310, 0.000362, 0.000480, 0.001200, 0.002000, 0.002000, 0.002600, 0.007600, 0.007500, 0.017000], [0.000141, 0.000298, 0.000517, 0.000450, 0.000810, 0.001900, 0.002700, 0.003900, 0.011000, 0.012052], [0.000110, 0.000180, 0.000640, 0.000940, 0.000980, 0.001600, 0.003170, 0.005600, 0.010000, 0.009200], [0.000091, 0.000270, 0.000730, 0.001100, 0.000530, 0.001600, 0.003000, 0.003200, 0.005425, 0.007560]])
            
#            # Bivariate quadratic least squares fit to the data
#            x = (bins_Porb[:-1]+bins_Porb[1:])/2.
#            y = (bins_Rp[:-1]+bins_Rp[1:])/2.
#            xx, yy = np.meshgrid(x, y)
##            zz = np.log(rates)
##            zz = np.log(err_u)
#            zz = np.log(err_l)
#            
#            ww = np.isinf(zz) == False
#            XX = xx[ww]
#            YY = yy[ww]
#            ZZ = zz[ww]
#            
#            #
#            B = ZZ
#            
#            A_bqf = np.array([XX**2, YY**2, XX, YY, XX*YY, XX*0.+1.]).T
#            p_bqf = np.linalg.lstsq(A_bqf, B)
#            def bqf(p, xx, yy):
#                return p[0]*xx**2+p[1]*yy**2+p[2]*xx+p[3]*yy+p[4]*xx*yy+p[5]*(xx*0.+1.)
#            r_bqf = bqf(p_bqf[0], xx.flatten(), yy.flatten()).reshape(zz.shape)
#            
##            A_lin = np.array([XX, YY, XX*0.+1.]).T
##            p_lin = np.linalg.lstsq(A_lin, B)
##            def lin(p, xx, yy):
##                return p[0]*xx+p[1]*yy+p[2]*(xx*0.+1.)
##            r_lin = lin(p_lin[0], xx.flatten(), yy.flatten()).reshape(zz.shape)
#            
#            # Extrapolate empty bins
#            ww = np.isinf(zz) == True
##            rates_extrap = rates.copy()
##            rates_extrap[ww] = np.exp(bqf(p_bqf[0], xx[ww], yy[ww]))
##            err_u_extrap = err_u.copy()
##            err_u_extrap[ww] = np.exp(bqf(p_bqf[0], xx[ww], yy[ww]))
#            err_l_extrap = err_l.copy()
#            err_l_extrap[ww] = np.exp(bqf(p_bqf[0], xx[ww], yy[ww]))
#            err_l_extrap[rates_extrap < err_l_extrap] = rates_extrap[rates_extrap < err_l_extrap]
#            
#            #
#            fig = plt.figure()
#            ax = fig.add_subplot(111, projection='3d')
#            ax.plot_surface(xx, yy, zz)
#            ax.plot_wireframe(xx, yy, r_bqf, color='black')
#            ax.set_xlabel('Log(orbital period)')
#            ax.set_ylabel('Log(radius)')
##            ax.set_zlabel('Log($\eta$)')
#            ax.set_zlabel('Log($\Delta_\eta$)')
#            ax.view_init(elev=30., azim=130.)
#            plt.title('Hsu et al. 2019')
#            plt.show(block=True)
#            
#            np.set_printoptions(formatter={'float_kind':'{:.6f}'.format})
#            import pdb; pdb.set_trace()
#            
#            # Print planet occurrence rates for TeX
#            for i in range(rates.shape[0]):
#                string = ''
#                for j in range(rates.shape[1]):
#                    if (rates[i, j] > 1E-10):
#                        string += '$%.4f' % (rates_extrap[i, j]*100.)+'^{+%.4f' % (err_u_extrap[i, j]*100.)+'}_{-%.4f' % (err_l_extrap[i, j]*100.)+'}$'
#                    else:
#                        string += '\\textcolor{red}{$%.4f' % (rates_extrap[i, j]*100.)+'^{+%.4f' % (err_u_extrap[i, j]*100.)+'}_{-%.4f' % (err_l_extrap[i, j]*100.)+'}$}'
#                    if (j < rates.shape[1]-1):
#                        string += ' & '
#                    else:
#                        string += ' \\\\'
#                print(string)
#            
#            import pdb; pdb.set_trace()
            
            # Use baseline, min or max planet occurrence rates
            if (self.model == 0):
                F0 = rates_extrap
            elif (self.model == 1):
                F0 = rates_extrap-err_l_extrap
            elif (self.model == 2):
                F0 = rates_extrap+err_u_extrap
            else:
                raise UserWarning('Model must be 0, 1 or 2')
             
            # Go through all bins
            for i in range(len(bins_Rp)-1):
                for j in range(len(bins_Porb)-1):
                    
                    # Distribute number of planets per star Poissonian
                    for k in range(np.random.poisson(F0[i, j])):
                        
                        # Draw random planet radius and planet orbital period
                        dummy1 = np.exp(bins_Rp[i]+(bins_Rp[i+1]-bins_Rp[i])*np.random.rand())
                        dummy2 = np.exp(bins_Porb[j]+(bins_Porb[j+1]-bins_Porb[j])*np.random.rand())
                        
                        # Check if planet radius and planet orbital period lie inside range
                        if (float(Rp_range[0]) <= dummy1 <= float(Rp_range[1])):
                            if (float(Porb_range[0]) <= dummy2 <= float(Porb_range[1])):
                                
                                # Fill lists for planet radii and planet orbital periods
                                Rp_hsu2019 += [dummy1]
                                Porb_hsu2019 += [dummy2]
            
        # Return lists with planet radii and planet orbital periods
        return Rp_hsu2019, Porb_hsu2019
    
    def dressing2015(self,
                     Rp_range=[0.5, 4.0],
                     Porb_range=[0.5, 200.0]):
        """
        Draws planet radii and planet orbital periods according to Dressing &
        Charbonneau 2015 statistics
        (https://ui.adsabs.harvard.edu/abs/2015ApJ...807...45D/abstract)
        """
        
        # Create lists for planet radii and planet orbital periods
        Rp_dressing2015 = []
        Porb_dressing2015 = []
        
        # Check if host star type suits
        if (True):
            
            # Parameters from Dressing 2015
            bins_Rp = np.log(np.array([0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4.]))
            bins_Porb = np.log(np.array([0.5, 1.7, 5.5, 18.2, 60.3, 200.]))
#            rates = np.array([[1.38, 8.42, 20.59, 0., 0.], [1.95, 9.94, 0., 26.85, 28.85], [0.41, 4.15, 0., 24.59, 19.98], [0., 2.72, 18.73, 27.58, 18.08], [0., 1.59, 8.29, 14.51, 8.61], [0., 0.65, 3.25, 3.37, 1.97], [0., 0.38, 1.05, 0.56, 0.]])*1E-2
#            err_u = np.array([[0.93, 3.53, 8.7, 0., 0.], [0.93, 2.82, 0., 10.7, 28.66], [0.51, 1.94, 0., 9.08, 16.07], [0., 1.73, 5.39, 9.42, 13.21], [0., 1.39, 3.96, 7.71, 9.78], [0., 1., 2.72, 4.62, 5.87], [0., 0.77, 1.82, 2.32, 0.]])*1E-2
#            err_l = np.array([[0.53, 2.39, 5.57, 0., 0.], [0.61, 2.13, 0., 6.79, 10.34], [0.2, 1.28, 0., 6., 7.42], [0., 1.01, 3.95, 6.31, 6.57], [0., 0.69, 2.55, 4.62, 3.84], [0., 0.32, 1.37, 1.62, 0.85], [0., 0.19, 0.52, 0.19, 0.]])*1E-2
            rates_extrap = np.array([[0.0138, 0.0842, 0.2059, 0.2559, 0.1276], [0.0195, 0.0994, 0.3600, 0.2685, 0.2885], [0.0041, 0.0415, 0.2341, 0.2459, 0.1998], [0.0036, 0.0272, 0.1873, 0.2758, 0.1808], [0.0015, 0.0159, 0.0829, 0.1451, 0.0861], [0.0006, 0.0065, 0.0325, 0.0337, 0.0197], [0.0003, 0.0038, 0.0105, 0.0056, 0.0135]])
            err_u_extrap = np.array([[0.0093, 0.0353, 0.0870, 0.1839, 0.2879], [0.0093, 0.0282, 0.0769, 0.1070, 0.2866], [0.0051, 0.0194, 0.0583, 0.0908, 0.1607], [0.0050, 0.0173, 0.0539, 0.0942, 0.1321], [0.0039, 0.0139, 0.0396, 0.0771, 0.0978], [0.0030, 0.0100, 0.0272, 0.0462, 0.0587], [0.0023, 0.0077, 0.0182, 0.0232, 0.0489]])
            err_l_extrap = np.array([[0.0053, 0.0239, 0.0557, 0.0852, 0.0704], [0.0061, 0.0213, 0.0769, 0.0679, 0.1034], [0.0020, 0.0128, 0.0532, 0.0600, 0.0742], [0.0020, 0.0101, 0.0395, 0.0631, 0.0657], [0.0010, 0.0069, 0.0255, 0.0462, 0.0384], [0.0006, 0.0032, 0.0137, 0.0162, 0.0085], [0.0003, 0.0019, 0.0052, 0.0019, 0.0071]])
            
#            # Bivariate quadratic least squares fit to the data
#            x = (bins_Porb[:-1]+bins_Porb[1:])/2.
#            y = (bins_Rp[:-1]+bins_Rp[1:])/2.
#            xx, yy = np.meshgrid(x, y)
##            zz = np.log(rates)
##            zz = np.log(err_u)
#            zz = np.log(err_l)
#            
#            ww = np.isinf(zz) == False
#            XX = xx[ww]
#            YY = yy[ww]
#            ZZ = zz[ww]
#            
#            #
#            B = ZZ
#            
#            A_bqf = np.array([XX**2, YY**2, XX, YY, XX*YY, XX*0.+1.]).T
#            p_bqf = np.linalg.lstsq(A_bqf, B)
#            def bqf(p, xx, yy):
#                return p[0]*xx**2+p[1]*yy**2+p[2]*xx+p[3]*yy+p[4]*xx*yy+p[5]*(xx*0.+1.)
#            r_bqf = bqf(p_bqf[0], xx.flatten(), yy.flatten()).reshape(zz.shape)
#            
##            A_lin = np.array([XX, YY, XX*0.+1.]).T
##            p_lin = np.linalg.lstsq(A_lin, B)
##            def lin(p, xx, yy):
##                return p[0]*xx+p[1]*yy+p[2]*(xx*0.+1.)
##            r_lin = lin(p_lin[0], xx.flatten(), yy.flatten()).reshape(zz.shape)
#            
#            # Extrapolate empty bins
#            ww = np.isinf(zz) == True
##            rates_extrap = rates.copy()
##            rates_extrap[ww] = np.exp(bqf(p_bqf[0], xx[ww], yy[ww]))
##            err_u_extrap = err_u.copy()
##            err_u_extrap[ww] = np.exp(bqf(p_bqf[0], xx[ww], yy[ww]))
#            err_l_extrap = err_l.copy()
#            err_l_extrap[ww] = np.exp(bqf(p_bqf[0], xx[ww], yy[ww]))
#            err_l_extrap[rates_extrap < err_l_extrap] = rates_extrap[rates_extrap < err_l_extrap]
#            
#            #
#            fig = plt.figure()
#            ax = fig.add_subplot(111, projection='3d')
#            ax.plot_surface(xx, yy, zz)
#            ax.plot_wireframe(xx, yy, r_bqf, color='black')
#            ax.set_xlabel('Log(orbital period)')
#            ax.set_ylabel('Log(radius)')
##            ax.set_zlabel('Log($\eta$)')
##            ax.set_zlabel('Log($\Delta_+\eta$)')
#            ax.set_zlabel('Log($\Delta_-\eta$)')
#            ax.view_init(elev=30., azim=130.)
#            plt.title('Dressing & Charbonneau 2015')
#            plt.show(block=True)
#            
#            np.set_printoptions(formatter={'float_kind':'{:.4f}'.format})
#            import pdb; pdb.set_trace()
#            
#            # Print planet occurrence rates for TeX
#            for i in range(rates.shape[0]):
#                string = ''
#                for j in range(rates.shape[1]):
#                    if (rates[i, j] > 1E-10):
#                        string += '$%.2f' % (rates_extrap[i, j]*100.)+'^{+%.2f' % (err_u_extrap[i, j]*100.)+'}_{-%.2f' % (err_l_extrap[i, j]*100.)+'}$'
#                    else:
#                        string += '\\textcolor{red}{$%.2f' % (rates_extrap[i, j]*100.)+'^{+%.2f' % (err_u_extrap[i, j]*100.)+'}_{-%.2f' % (err_l_extrap[i, j]*100.)+'}$}'
#                    if (j < rates.shape[1]-1):
#                        string += ' & '
#                    else:
#                        string += ' \\\\'
#                print(string)
#            
#            import pdb; pdb.set_trace()
            
            # Use baseline, min or max planet occurrence rates
            if (self.model == 0):
                F0 = rates_extrap
            elif (self.model == 1):
                F0 = rates_extrap-err_l_extrap
            elif (self.model == 2):
                F0 = rates_extrap+err_u_extrap
            else:
                raise UserWarning('Model must be 0, 1 or 2')
            
            # Go through all bins
            for i in range(len(bins_Rp)-1):
                for j in range(len(bins_Porb)-1):
                    
                    # Distribute number of planets per star Poissonian
                    for k in range(np.random.poisson(F0[i, j])):
                        
                        # Draw random planet radius and planet orbital period
                        dummy1 = np.exp(bins_Rp[i]+(bins_Rp[i+1]-bins_Rp[i])*np.random.rand())
                        dummy2 = np.exp(bins_Porb[j]+(bins_Porb[j+1]-bins_Porb[j])*np.random.rand())
                        
                        # Check if planet radius and planet orbital period lie inside range
                        if (float(Rp_range[0]) <= dummy1 <= float(Rp_range[1])):
                            if (float(Porb_range[0]) <= dummy2 <= float(Porb_range[1])):
                                
                                # Fill lists for planet radii and planet orbital periods
                                Rp_dressing2015 += [dummy1]
                                Porb_dressing2015 += [dummy2]
        
        # Return lists with planet radii and planet orbital periods
        return Rp_dressing2015, Porb_dressing2015
    
    def fressin2013(self,
                    Rp_range=[0.8, 22.0],
                    Porb_range=[0.8, 418.0]):
        """
        Draws planet radii and planet orbital periods according to Fressin et
        al. 2013 statistics
        (https://ui.adsabs.harvard.edu/abs/2013ApJ...766...81F/abstract)
        """
        
        # Create lists for planet radii and planet orbital periods
        Rp_fressin2013 = []
        Porb_fressin2013 = []
        
        # Check if host star type suits
        if (True):
            
            # Parameters from Fressin 2013
            bins_Rp = np.log(np.array([0.8, 1.25, 2., 4., 6., 22.]))
            bins_Porb = np.log(np.array([0.8, 2., 3.4, 5.9, 10., 17., 29., 50., 85., 145., 245., 418.]))
#            rates = np.array([[0.18, 0.61, 1.72, 2.7, 2.7, 2.93, 4.08, 3.46, 0., 0., 0.], [0.17, 0.74, 1.49, 2.9, 4.3, 4.49, 5.29, 3.66, 6.54, 0., 0.], [0.035, 0.18, 0.73, 1.93, 3.67, 5.29, 6.45, 5.25, 4.31, 3.09, 0.], [0.004, 0.006, 0.11, 0.091, 0.29, 0.32, 0.49, 0.66, 0.43, 0.53, 0.24], [0.015, 0.067, 0.17, 0.18, 0.27, 0.23, 0.35, 0.71, 1.25, 0.94, 1.05]])*1E-2
#            err_u = np.array([[0.04, 0.15, 0.43, 0.6, 0.83, 1.05, 1.88, 2.81, 0., 0., 0.], [0.03, 0.13, 0.23, 0.56, 0.73, 1., 1.48, 1.21, 2.2, 0., 0.], [0.011, 0.03, 0.09, 0.19, 0.39, 0.64, 1.01, 1.05, 1.03, 0.9, 0.], [0.003, 0.006, 0.03, 0.03, 0.07, 0.08, 0.12, 0.16, 0.17, 0.21, 0.15], [0.007, 0.018, 0.03, 0.04, 0.06, 0.06, 0.1, 0.17, 0.29, 0.28, 0.3]])*1E-2
#            err_l = np.array([[0.04, 0.15, 0.43, 0.6, 0.83, 1.05, 1.88, 2.81, 0., 0., 0.], [0.03, 0.13, 0.23, 0.56, 0.73, 1., 1.48, 1.21, 2.2, 0., 0.], [0.011, 0.03, 0.09, 0.19, 0.39, 0.64, 1.01, 1.05, 1.03, 0.9, 0.], [0.003, 0.006, 0.03, 0.03, 0.07, 0.08, 0.12, 0.16, 0.17, 0.21, 0.15], [0.007, 0.018, 0.03, 0.04, 0.06, 0.06, 0.1, 0.17, 0.29, 0.28, 0.3]])*1E-2
            rates_extrap = np.array([[0.00180, 0.00610, 0.01720, 0.02700, 0.02700, 0.02930, 0.04080, 0.03460, 0.06672, 0.04862, 0.03049], [0.00170, 0.00740, 0.01490, 0.02900, 0.04300, 0.04490, 0.05290, 0.03660, 0.06540, 0.02737, 0.01775], [0.00035, 0.00180, 0.00730, 0.01930, 0.03670, 0.05290, 0.06450, 0.05250, 0.04310, 0.03090, 0.01017], [0.00004, 0.00006, 0.00110, 0.00091, 0.00290, 0.00320, 0.00490, 0.00660, 0.00430, 0.00530, 0.00240], [0.00015, 0.00067, 0.00170, 0.00180, 0.00270, 0.00230, 0.00350, 0.00710, 0.01250, 0.00940, 0.01050]])
            err_u_extrap = np.array([[0.00040, 0.00150, 0.00430, 0.00600, 0.00830, 0.01050, 0.01880, 0.02810, 0.04387, 0.04734, 0.04663], [0.00030, 0.00130, 0.00230, 0.00560, 0.00730, 0.01000, 0.01480, 0.01210, 0.02200, 0.01787, 0.01744], [0.00011, 0.00030, 0.00090, 0.00190, 0.00390, 0.00640, 0.01010, 0.01050, 0.01030, 0.00900, 0.00640], [0.00003, 0.00006, 0.00030, 0.00030, 0.00070, 0.00080, 0.00120, 0.00160, 0.00170, 0.00210, 0.00150], [0.00007, 0.00018, 0.00030, 0.00040, 0.00060, 0.00060, 0.00100, 0.00170, 0.00290, 0.00280, 0.00300]])
            err_l_extrap = np.array([[0.00040, 0.00150, 0.00430, 0.00600, 0.00830, 0.01050, 0.01880, 0.02810, 0.04387, 0.04734, 0.03049], [0.00030, 0.00130, 0.00230, 0.00560, 0.00730, 0.01000, 0.01480, 0.01210, 0.02200, 0.01787, 0.01744], [0.00011, 0.00030, 0.00090, 0.00190, 0.00390, 0.00640, 0.01010, 0.01050, 0.01030, 0.00900, 0.00640], [0.00003, 0.00006, 0.00030, 0.00030, 0.00070, 0.00080, 0.00120, 0.00160, 0.00170, 0.00210, 0.00150], [0.00007, 0.00018, 0.00030, 0.00040, 0.00060, 0.00060, 0.00100, 0.00170, 0.00290, 0.00280, 0.00300]])
            
#            # Bivariate quadratic least squares fit to the data
#            x = (bins_Porb[:-1]+bins_Porb[1:])/2.
#            y = (bins_Rp[:-1]+bins_Rp[1:])/2.
#            xx, yy = np.meshgrid(x, y)
##            zz = np.log(rates)
##            zz = np.log(err_u)
#            zz = np.log(err_l)
#            
#            ww = np.isinf(zz) == False
#            XX = xx[ww]
#            YY = yy[ww]
#            ZZ = zz[ww]
#            
#            #
#            B = ZZ
#            
#            A_bqf = np.array([XX**2, YY**2, XX, YY, XX*YY, XX*0.+1.]).T
#            p_bqf = np.linalg.lstsq(A_bqf, B)
#            def bqf(p, xx, yy):
#                return p[0]*xx**2+p[1]*yy**2+p[2]*xx+p[3]*yy+p[4]*xx*yy+p[5]*(xx*0.+1.)
#            r_bqf = bqf(p_bqf[0], xx.flatten(), yy.flatten()).reshape(zz.shape)
#            
##            A_lin = np.array([XX, YY, XX*0.+1.]).T
##            p_lin = np.linalg.lstsq(A_lin, B)
##            def lin(p, xx, yy):
##                return p[0]*xx+p[1]*yy+p[2]*(xx*0.+1.)
##            r_lin = lin(p_lin[0], xx.flatten(), yy.flatten()).reshape(zz.shape)
#            
#            # Extrapolate empty bins
#            ww = np.isinf(zz) == True
##            rates_extrap = rates.copy()
##            rates_extrap[ww] = np.exp(bqf(p_bqf[0], xx[ww], yy[ww]))
##            err_u_extrap = err_u.copy()
##            err_u_extrap[ww] = np.exp(bqf(p_bqf[0], xx[ww], yy[ww]))
#            err_l_extrap = err_l.copy()
#            err_l_extrap[ww] = np.exp(bqf(p_bqf[0], xx[ww], yy[ww]))
#            err_l_extrap[rates_extrap < err_l_extrap] = rates_extrap[rates_extrap < err_l_extrap]
#            
#            #
#            fig = plt.figure()
#            ax = fig.add_subplot(111, projection='3d')
#            ax.plot_surface(xx, yy, zz)
#            ax.plot_wireframe(xx, yy, r_bqf, color='black')
#            ax.set_xlabel('Log(orbital period)')
#            ax.set_ylabel('Log(radius)')
##            ax.set_zlabel('Log($\eta$)')
#            ax.set_zlabel('Log($\Delta_\eta$)')
#            ax.view_init(elev=30., azim=130.)
#            plt.title('Fressin et al. 2013')
#            plt.show(block=True)
#            
#            np.set_printoptions(formatter={'float_kind':'{:.5f}'.format})
#            import pdb; pdb.set_trace()
#            
#            # Print planet occurrence rates for TeX
#            for i in range(rates.shape[0]):
#                string = ''
#                for j in range(rates.shape[1]):
#                    if (rates[i, j] > 1E-10):
#                        string += '$%.3f' % (rates_extrap[i, j]*100.)+'^{+%.3f' % (err_u_extrap[i, j]*100.)+'}_{-%.3f' % (err_l_extrap[i, j]*100.)+'}$'
#                    else:
#                        string += '\\textcolor{red}{$%.3f' % (rates_extrap[i, j]*100.)+'^{+%.3f' % (err_u_extrap[i, j]*100.)+'}_{-%.3f' % (err_l_extrap[i, j]*100.)+'}$}'
#                    if (j < rates.shape[1]-1):
#                        string += ' & '
#                    else:
#                        string += ' \\\\'
#                print(string)
#            
#            import pdb; pdb.set_trace()
            
            # Use baseline, min or max planet occurrence rates
            if (self.model == 0):
                F0 = rates_extrap
            elif (self.model == 1):
                F0 = rates_extrap-err_l_extrap
            elif (self.model == 2):
                F0 = rates_extrap+err_u_extrap
            else:
                raise UserWarning('Model must be 0, 1 or 2')
            
            # Go through all bins
            for i in range(len(bins_Rp)-1):
                for j in range(len(bins_Porb)-1):
                    
                    # Distribute number of planets per star Poissonian
                    for k in range(np.random.poisson(F0[i, j])):
                        
                        # Draw random planet radius and planet orbital period
                        dummy1 = np.exp(bins_Rp[i]+(bins_Rp[i+1]-bins_Rp[i])*np.random.rand())
                        dummy2 = np.exp(bins_Porb[j]+(bins_Porb[j+1]-bins_Porb[j])*np.random.rand())
                        
                        # Check if planet radius and planet orbital period lie inside range
                        if (float(Rp_range[0]) <= dummy1 <= float(Rp_range[1])):
                            if (float(Porb_range[0]) <= dummy2 <= float(Porb_range[1])):
                                
                                # Fill lists for planet radii and planet orbital periods
                                Rp_fressin2013 += [dummy1]
                                Porb_fressin2013 += [dummy2]
        
        # Return lists with planet radii and planet orbital periods
        return Rp_fressin2013, Porb_fressin2013
    
    def sag13_dressing2015(self):
        """
        Draws planet radii and planet orbital periods according to
        - SAG13 for AFGK stars
        - Dressing & Charbonneau 2015 for M stars
        """
        
        if (self.stype in 'AFGK'):
            Rp, Porb = self.sag13()
        elif (self.stype in 'M'):
            Rp, Porb = self.dressing2015()
        else:
            Rp = []; Porb = []
        
        return Rp, Porb
    
    def hsu2019_dressing2015(self):
        """
        Draws planet radii and planet orbital periods according to
        - Hsu et al. 2019 for AFGK stars
        - Dressing & Charbonneau 2015 for M stars
        """
        
        if (self.stype in 'AFGK'):
            Rp, Porb = self.hsu2019()
        elif (self.stype in 'M'):
            Rp, Porb = self.dressing2015()
        else:
            Rp = []; Porb = []
        
        return Rp, Porb

class flux():
    
    def __init__(self,
                 Rp=1,
                 Tp=255,
                 rp=1,
                 AgeomMIR=0.367,
                 AgeomVIS=0.367,
                 f=1,
                 dist=10,
                 Rs=1,
                 Ts=5772,
                 nodes=np.linspace(9, 11, 20),
                 trans=np.ones((20)),
                 lam_eff=10,
                 W_eff=2,
                 mission='MIR',
                 mags=None,
                 F0=None):
        """
        Initializes flux instance
        """
        
        # Planet radius
        Rp = float(Rp)
        if (Rp <= 0):
            raise ValueError('Rp must be positive')
        self.Rp = Rp
        
        # Planet effective temperature
        Tp = float(Tp)
        if (Tp <= 0):
            raise ValueError('Tp must be positive')
        self.Tp = Tp
        
        # Planet host star separation
        rp = float(rp)
        if (rp <= 0):
            raise ValueError('rp must be positive')
        self.rp = rp
        
        # Planet geometric albedo (MIR)
        AgeomMIR = float(AgeomMIR)
        if (AgeomMIR < 0):
            raise ValueError('AgeomMIR must be positive')
        self.AgeomMIR = AgeomMIR
        
        # Planet geometric albedo (VIS)
        AgeomVIS = float(AgeomVIS)
        if (AgeomVIS < 0):
            raise ValueError('AgeomVIS must be positive')
        self.AgeomVIS = AgeomVIS
        
        # Lambertian reflectance
        f = float(f)
        if (f < 0):
            raise ValueError('f must be positive')
        self.f = f
        
        # Host star distance from Earth
        dist = float(dist)
        if (dist <= 0):
            raise ValueError('dist must be positive')
        self.dist = dist
        
        # Host star radius
        Rs = float(Rs)
        if (Rs <= 0):
            raise ValueError('Rs must be positive')
        self.Rs = Rs
        
        # Host star effective temperature
        Ts = float(Ts)
        if (Ts <= 0):
            raise ValueError('Ts must be positive')
        self.Ts = Ts
        
        # Filter nodes
        if (nodes is not None): # CAN BE NONE
            if (isinstance(nodes, np.ndarray) == False):
                raise TypeError('nodes must be of type numpy.ndarray')
        self.nodes = nodes
        
        # Filter transmission
        if (trans is not None): # CAN BE NONE
            if (isinstance(trans, np.ndarray) == False):
                raise TypeError('trans must be of type numpy.ndarray')
            if (len(nodes) != len(trans)):
                raise ValueError('trans must have the same length like nodes')
        self.trans = trans
        
        # Filter effective wavelength
        if (lam_eff is not None): # CAN BE NONE
            lam_eff = float(lam_eff)
            if (lam_eff <= 0):
                raise ValueError('lam_eff must be positive')
        self.lam_eff = lam_eff
        
        # Filter effective width
        if (W_eff is not None): # CAN BE NONE
            W_eff = float(W_eff)
            if (W_eff <= 0):
                raise ValueError('W_eff must be positive')
        self.W_eff = W_eff
        
        # Mission operating wavelength
        self.mission = mission
        
        # Host star magnitude
        if (mags is not None):  # CAN BE NONE
            mags = float(mags)
        self.mags = mags
        
        # Filter zero point (of the filter with which mags was measured)
        if (F0 is not None):  # CAN BE NONE
            F0 = float(F0)
            if (F0 <= 0):
                raise ValueError('F0 must be positive')
        self.F0 = F0
    
    def bb_therm_p(self,
                   lam = np.logspace(-1, 2, 100)):
        """
        Computes planet thermal blackbody flux (in Watts per square meter per
        micron)
        """
        
        # Check parameters
        if (isinstance(lam, np.ndarray) == False):
            raise TypeError('lam must be of type numpy.ndarray')
        
        # Define constants
        h = 6.62607004E-34
        c = 299792458.
        kB = 1.38064852E-23
        Rearth = 6371000.
        pc = 3.0856776E+16
        
        # Return planet thermal blackbody flux
        return 1E-6*(2.*np.pi*h*c**2./(lam*1E-6)**5./(np.exp(h*c/(lam*1E-6*kB*self.Tp))-1.)*((self.Rp*Rearth)/(self.dist*pc))**2.)
    
    def bb_therm_p_filter(self):
        """
        Computes planet thermal blackbody flux observed through filter (in
        Watts per square meter per micron)
        """
        
        # Check parameters
        if (self.nodes is None):
            raise UserWarning('This function is only available if a filter curve is provided')
        
        # Return planet thermal blackbody flux observed through filter
        return self.bb_therm_p(self.nodes)*self.trans
    
    def bb_therm_p_int(self):
        """
        Computes integrated planet thermal blackbody flux observed through
        filter (in Janskys)
        """
        
        # Check parameters
        if (self.nodes is None):
            raise UserWarning('This function is only available if a filter curve is provided')
        
        # Define constants
        c = 299792458.
        
        # Return integrated planet thermal blackbody flux observed through filter
        return 1E+6*si.simps(self.bb_therm_p(self.nodes)*self.trans, self.nodes)/self.W_eff*(self.lam_eff*1E-6)**2/c*1E+26
    
    def bb_therm_s(self,
                   lam = np.logspace(-1, 2, 100)):
        """
        Computes host star thermal blackbody flux (in Watts per square meter 
        per micron)
        """
        
        # Check parameters
        if (isinstance(lam, np.ndarray) == False):
            raise TypeError('lam must be of type numpy.ndarray')
        
        # Define constants
        h = 6.62607004E-34
        c = 299792458.
        kB = 1.38064852E-23
        Rsun = 695700000.
        pc = 3.0856776E+16
        
        # Return host star thermal blackbody flux
        return 1E-6*(2.*np.pi*h*c**2./(lam*1E-6)**5./(np.exp(h*c/(lam*1E-6*kB*self.Ts))-1.)*((self.Rs*Rsun)/(self.dist*pc))**2.)
    
    def bb_therm_s_filter(self):
        """
        Computes host star thermal blackbody flux observed through filter (in
        Watts per square meter per micron)
        """
        
        # Check parameters
        if (self.nodes is None):
            raise UserWarning('This function is only available if a filter curve is provided')
        
        # Return host star thermal blackbody flux observed through filter
        return self.bb_therm_s(self.nodes)*self.trans
    
    def bb_therm_s_int(self):
        """
        Computes integrated host star thermal blackbody flux observed through
        filter (in Janskys)
        """
        
        # Check parameters
        if (self.nodes is None):
            raise UserWarning('This function is only available if a filter curve is provided')
        
        # Define constants
        c = 299792458.
        
        # Return integrated host star thermal blackbody flux observed through filter
        return 1E+6*si.simps(self.bb_therm_s(self.nodes)*self.trans, self.nodes)/self.W_eff*(self.lam_eff*1E-6)**2/c*1E+26
    
    def refl_p(self,
               lam = np.logspace(-1, 2, 100)):
        """
        Computes planet reflected host star flux (in Watts per square meter per
        micron)
        """
        
        # Check parameters
        if (isinstance(lam, np.ndarray) == False):
            raise TypeError('lam must be of type numpy.ndarray')
        
        # Define constants
        Rearth = 6371000.
        pc = 3.0856776E+16        
        h = 6.62607004E-34
        c = 299792458.
        kB = 1.38064852E-23
        Rsun = 695700000.
        au = 149597870700.
        
        # Return planet reflected host star flux
        if (self.mission == 'MIR'):
            return self.AgeomMIR*self.f*((self.Rp*Rearth)/(self.dist*pc))**2.*(1E-6*(2.*np.pi*h*c**2./(lam*1E-6)**5./(np.exp(h*c/(lam*1E-6*kB*self.Ts))-1.)*((self.Rs*Rsun)/(self.rp*au))**2.))
        if (self.mission == 'VIS'):
            return self.AgeomVIS*self.f*((self.Rp*Rearth)/(self.dist*pc))**2.*(1E-6*(2.*np.pi*h*c**2./(lam*1E-6)**5./(np.exp(h*c/(lam*1E-6*kB*self.Ts))-1.)*((self.Rs*Rsun)/(self.rp*au))**2.))
        else:
            raise UserWarning('Mission operating wavelength must be either MIR or VIS')
    
    def refl_p_filter(self):
        """
        Computes planet reflected host star flux observed through filter (in
        Watts per square meter per micron)
        """
        
        # Check parameters
        if (self.nodes is None):
            raise UserWarning('This function is only available if a filter curve is provided')
        
        # Return planet reflected host star flux observed through filter
        return self.refl_p(self.nodes)*self.trans
    
    def refl_p_int(self):
        """
        Computes integrated planet reflected host star flux observed through
        filter (in Janskys)
        """
        
        # Check parameters
        if (self.nodes is None):
            raise UserWarning('This function is only available if a filter curve is provided')
        
        # Define constants
        c = 299792458.
        
        # Return integrated planet reflected host star flux observed through filter
        return 1E+6*si.simps(self.refl_p(self.nodes)*self.trans, self.nodes)/self.W_eff*(self.lam_eff*1E-6)**2/c*1E+26
    
    def s_uJy(self):
        """
        Computes observed host star flux (in micro-Janskys)
        """
        
        # Check parameters
        if (self.mags is None):
            raise UserWarning('This function is only available if a host star magnitude is provided')
        
        # Return observed host star flux
        return 1E+6*self.F0/10**(self.mags/2.5)
    
    def p_therm_uJy(self):
        """
        Computes observed planet thermal flux (in micro-Janskys)
        """
        
        # Check parameters
        if (self.mags is None):
            raise UserWarning('This function is only available if a host star magnitude is provided')
        
        # Return observed planet thermal flux
        return 0.

    
    def p_refl_uJy(self):
        """
        Computes observed planet reflected host star flux (in micro-Janskys)
        """
        
        # Check parameters
        if (self.mags is None):
            raise UserWarning('This function is only available if a host star magnitude is provided')
        
        # Define constants
        Rearth = 6371000.
        au = 149597870700.
        
        # Return observed planet reflected host star flux
        if (self.mission == 'MIR'):
            return self.AgeomMIR*self.f*(self.Rp*Rearth)**2*self.s_uJy()/(self.rp*au)**2
        if (self.mission == 'VIS'):
            return self.AgeomVIS*self.f*(self.Rp*Rearth)**2*self.s_uJy()/(self.rp*au)**2
        else:
            raise UserWarning('Mission operating wavelength must be either MIR or VIS')
