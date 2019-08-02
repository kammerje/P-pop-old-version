"""
P_pop, a Python package to simulate exoplanet populations based on Kepler
statistics. Information about the method can be found in this paper:
http://adsabs.harvard.edu/abs/2018A%26A...609A...4K. This library is maintained
on GitHub at https://github.com/kammerje/P_pop.

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

import inout as io
import p_pop as pp


# PARAMETERS - PLEASE CHANGE ACCORDING TO YOUR WISHES
#==============================================================================
# Name of the output table containing all the exoplanets
logfile_name = 'SAG13_exocat500.txt'

# Number of Monte-Carlo shoots
nMC = 500

# Star catalog
catalog = 'ExoCat'
    # 20pc_bright_sample - star catalog from Kammerer & Quanz 2018
    # ExoCat - Nearby Stellar Systems Catalog for Exoplanet Imaging Missions (Turnbull/NASA)

# Scenario for the planet occurrence rates
model = 0
    # 0 - baseline
    # 1 - minimum
    # 2 - maximum

# Planet occurrence rates
rates = 'sag13_dressing2015'
    # sag13 - NASA SAG13, Kopparapu et al. 2018
    # hsu2019 - Hsu et al. 2019
    # dressing2015 - Dressing & Charbonneau 2015
    # fressin2013 - Fressin et al. 2013
    #
    # sag13_dressing2015 - sag13 (AFGK) & dressing2015 (M)
    # hsu2019_dressing2015 - hsu2019 (AFGK) & dressing2015 (M)

# Albedo model
albedo_model = 'uniform'
    # uniform - draw uniformly distributed Bond & geometric albedos
    # cahoy2010 - draw geometric albedo from Table 4 of Cahoy et al. 2010

# Mission operating wavelength
mission = 'VIS'
    # MIR - mid-infrared (this will compute fluxes for LIFE)
    # VIS - visible (this will compute fluxes for LUVOIR)

# Flux computation model
fluxes_from = 'mags'
    # filters - computes fluxes using the filter curves provided below
    # mags - computes fluxes using the stellar magnitudes from the star catalog (so far, this only works for catalog=ExoCat & mission=VIS)


# READ FILTERS & STAR CATALOG
#==============================================================================
if (fluxes_from == 'filters'):
    # LIFE filters
    io.read_filter_MIRI(path='data/F560W.csv',
                        tag='F560W',
                        lam_eff=5.6,
                        W_eff=1.2)
    io.read_filter_MIRI(path='data/F1000W.csv',
                        tag='F1000W',
                        lam_eff=10.,
                        W_eff=2.)
    io.read_filter_MIRI(path='data/F1500W.csv',
                        tag='F1500W',
                        lam_eff=15.,
                        W_eff=3.)
    # LUVOIR filters
    io.read_filter_IRDIS(path='data/Vband.txt',
                         tag='Vband',
                         lam_eff=0.554,
                         W_eff=0.0806)
    io.read_filter_IRDIS(path='data/Jband.txt',
                         tag='Jband',
                         lam_eff=1.245,
                         W_eff=0.24)
    io.read_filter_IRDIS(path='data/Hband.txt',
                         tag='Hband',
                         lam_eff=1.625,
                         W_eff=0.29)

if (catalog == '20pc_bright_sample'):
    SC = io.read_20pc_bright_sample(path='data/20pc_bright_sample.tbl',
                                    s_type=['A', 'F', 'G', 'K', 'M'], # Spectral types which should be put into the star catalog
                                    max_dist=20, # Distance cut (pc) for the star catalog
                                    min_dec=-90, # Declination cut (deg) for the star catalog
                                    max_dec=90) # Declination cut (deg) for the star catalog
elif (catalog == 'ExoCat'):
    SC = io.read_exocat(path='data/ExoCat.csv',
                        s_type=['B', 'A', 'F', 'G', 'K', 'M', 'D'], # Spectral types which should be put into the star catalog
                        max_dist=30, # Distance cut (pc) for the star catalog
                        min_dec=-90, # Declination cut (deg) for the star catalog
                        max_dec=90) # Declination cut (deg) for the star catalog
else:
    raise UserWarning('Star catalog '+catalog+' is not known')


# GENERATE PLANET POPULATION
#==============================================================================
logfile = open(logfile_name, 'w')
logfile.write('nMC\tRp\tPorb\ta\trp\tang_sep\tang_sep_max\tinc\tOmega\tomega\ttheta\tecc\tFinc\tAbond\tAgeomMIR\tAgeomVIS\tf\tTp\tMp\tdist\tRs\tTs\tMs\tstype\tzodis\tra\tdec\tnstar\t\n')
logfile.close()

for i in range(len(SC)):
    print('Star '+str(int(i+1))+' of '+str(int(len(SC))))
    
    stype = SC[i]['stype']
    dist = SC[i]['dist']
    Rs = SC[i]['rad']
    Ts = SC[i]['teff']
    Ms = SC[i]['mass']
    ra = SC[i]['ra']
    dec = SC[i]['dec']
    
    for j in range(nMC):
        System = pp.system(stype=stype,
                           dist=dist,
                           Rs=Rs,
                           Ts=Ts,
                           Ms=Ms,
                           model=model,
                           rates=rates,
                           albedo_model=albedo_model)
        
        for k in range(System.nplanets):
            Exoplanet = System.planets[k]
            
            logfile = open(logfile_name, 'a')
            logfile.write('%.0f\t' % j+'%.5f\t' % Exoplanet.Rp+'%.5f\t' % Exoplanet.Porb+'%.5f\t' % Exoplanet.a()+'%.5f\t' % Exoplanet.rp()+'%.5f\t' % Exoplanet.ang_sep()+'%.5f\t' % Exoplanet.ang_sep_max()+'%.5f\t' % Exoplanet.inc+'%.5f\t' % Exoplanet.Omega+'%.5f\t' % Exoplanet.omega+'%.5f\t' % Exoplanet.theta+'%.5f\t' % Exoplanet.ecc+'%.5f\t' % Exoplanet.Finc()+'%.5f\t' % Exoplanet.Abond+'%.5f\t' % Exoplanet.AgeomMIR+'%.5f\t' % Exoplanet.AgeomVIS+'%.5f\t' % Exoplanet.f()+'%.5f\t' % Exoplanet.Tp()+'%.5f\t' % Exoplanet.Mp+'%.5f\t' % Exoplanet.dist+'%.5f\t' % Exoplanet.Rs+'%.5f\t' % Exoplanet.Ts+'%.5f\t' % Exoplanet.Ms+stype+'\t%.5f\t' % Exoplanet.zodis+'%.5f\t' % ra+'%.5f\t' % dec+'%.0f\t' % i+'\n')
            logfile.close()


# COMPUTE FLUXES
#==============================================================================
if (fluxes_from == 'filters'):
    if (mission == 'MIR'):
        # Read LIFE filters
        temp_name = logfile_name[:-4]+'_LIFE.txt'
        name1 = 'F560W'
        name2 = 'F1000W'
        name3 = 'F1500W'
    elif (mission == 'VIS'):
        # Read LUVOIR filters
        temp_name = logfile_name[:-4]+'_LUVOIR.txt'
        name1 = 'Vband'
        name2 = 'Jband'
        name3 = 'Hband'
    else:
        raise UserWarning('Mission operating wavelength '+mission+' is not known')
    filter1 = np.load('data/'+name1+'.npy')
    filter2 = np.load('data/'+name2+'.npy')
    filter3 = np.load('data/'+name3+'.npy')
    
    logfile = open(temp_name, 'w')
    logfile.write('Ftherm_star\tFtherm_planet\tFrefl_planet\t\n')
    logfile.close()
    planets = open(logfile_name, 'r')
    planets_lines = planets.readlines()
    
    i = 0
    for line in planets_lines:
        line_temp = line.split('\t')
        if (i == 0):
            ncol = len(line_temp)
            col_Rp = np.where(np.array(line_temp) == 'Rp')[0][0]
            col_Tp = np.where(np.array(line_temp) == 'Tp')[0][0]
            col_rp = np.where(np.array(line_temp) == 'rp')[0][0]
            col_AgeomMIR = np.where(np.array(line_temp) == 'AgeomMIR')[0][0]
            col_AgeomVIS = np.where(np.array(line_temp) == 'AgeomVIS')[0][0]
            col_f = np.where(np.array(line_temp) == 'f')[0][0]
            col_dist = np.where(np.array(line_temp) == 'dist')[0][0]
            col_Rs = np.where(np.array(line_temp) == 'Rs')[0][0]
            col_Ts = np.where(np.array(line_temp) == 'Ts')[0][0]
        elif (len(line_temp) == ncol):
            flux = pp.flux(Rp=float(line_temp[col_Rp]),
                           Tp=float(line_temp[col_Tp]),
                           rp=float(line_temp[col_rp]),
                           AgeomMIR=float(line_temp[col_AgeomMIR]),
                           AgeomVIS=float(line_temp[col_AgeomVIS]),
                           f=float(line_temp[col_f]),
                           dist=float(line_temp[col_dist]),
                           Rs=float(line_temp[col_Rs]),
                           Ts=float(line_temp[col_Ts]),
                           nodes=filter1[0],
                           trans=filter1[1],
                           lam_eff=filter1[2],
                           W_eff=filter1[3],
                           mission=mission)
            logfile = open(temp_name, 'a')
            logfile.write('%018.12f\t' % (flux.bb_therm_s_int()*1E6)+'%018.12f\t' % (flux.bb_therm_p_int()*1E6)+'%018.12f\t' % (flux.refl_p_int()*1E6))
            logfile.close()
            
            flux = pp.flux(Rp=float(line_temp[col_Rp]),
                           Tp=float(line_temp[col_Tp]),
                           rp=float(line_temp[col_rp]),
                           AgeomMIR=float(line_temp[col_AgeomMIR]),
                           AgeomVIS=float(line_temp[col_AgeomVIS]),
                           f=float(line_temp[col_f]),
                           dist=float(line_temp[col_dist]),
                           Rs=float(line_temp[col_Rs]),
                           Ts=float(line_temp[col_Ts]),
                           nodes=filter2[0],
                           trans=filter2[1],
                           lam_eff=filter2[2],
                           W_eff=filter2[3],
                           mission=mission)
            logfile = open(temp_name, 'a')
            logfile.write('%018.12f\t' % (flux.bb_therm_s_int()*1E6)+'%018.12f\t' % (flux.bb_therm_p_int()*1E6)+'%018.12f\t' % (flux.refl_p_int()*1E6))
            logfile.close()
            
            flux = pp.flux(Rp=float(line_temp[col_Rp]),
                           Tp=float(line_temp[col_Tp]),
                           rp=float(line_temp[col_rp]),
                           AgeomMIR=float(line_temp[col_AgeomMIR]),
                           AgeomVIS=float(line_temp[col_AgeomVIS]),
                           f=float(line_temp[col_f]),
                           dist=float(line_temp[col_dist]),
                           Rs=float(line_temp[col_Rs]),
                           Ts=float(line_temp[col_Ts]),
                           nodes=filter3[0],
                           trans=filter3[1],
                           lam_eff=filter3[2],
                           W_eff=filter3[3],
                           mission=mission)
            logfile = open(temp_name, 'a')
            logfile.write('%018.12f\t' % (flux.bb_therm_s_int()*1E6)+'%018.12f\t' % (flux.bb_therm_p_int()*1E6)+'%018.12f\t' % (flux.refl_p_int()*1E6)+'\n')
            logfile.close()
        else:
            print('Line %.0f: unappropriate data' % i)
            import pdb; pdb.set_trace()
        
        i += 1
        if (i % 1000 == 0):
            print('Processed '+str(i)+' planets')
    planets.close()

elif (fluxes_from == 'mags'):
    if (mission == 'VIS'):
        temp_name = logfile_name[:-4]+'_LUVOIR.txt'
        F0s = [3562.41, 1582.23, 1024.74]
        print('--> Using default filter zero points:\n'+str(F0s[0])+' Jy (ZIMPOL V)\n'+str(F0s[1])+' Jy (2MASS J)\n'+str(F0s[2])+' Jy (2MASS H)')
        
        logfile = open(temp_name, 'w')
        logfile.write('Ftherm_star\tFtherm_planet\tFrefl_planet\t\n')
        logfile.close()
        planets = open(logfile_name, 'r')
        planets_lines = planets.readlines()
        
        i = 0
        for line in planets_lines:
            line_temp = line.split('\t')
            if (i == 0):
                ncol = len(line_temp)
                col_Rp = np.where(np.array(line_temp) == 'Rp')[0][0]
                col_Tp = np.where(np.array(line_temp) == 'Tp')[0][0]
                col_rp = np.where(np.array(line_temp) == 'rp')[0][0]
                col_AgeomMIR = np.where(np.array(line_temp) == 'AgeomMIR')[0][0]
                col_AgeomVIS = np.where(np.array(line_temp) == 'AgeomVIS')[0][0]
                col_f = np.where(np.array(line_temp) == 'f')[0][0]
                col_dist = np.where(np.array(line_temp) == 'dist')[0][0]
                col_Rs = np.where(np.array(line_temp) == 'Rs')[0][0]
                col_Ts = np.where(np.array(line_temp) == 'Ts')[0][0]
                col_nstar = np.where(np.array(line_temp) == 'nstar')[0][0]
            elif (len(line_temp) == ncol):
                flux = pp.flux(Rp=float(line_temp[col_Rp]),
                               Tp=float(line_temp[col_Tp]),
                               rp=float(line_temp[col_rp]),
                               AgeomMIR=float(line_temp[col_AgeomMIR]),
                               AgeomVIS=float(line_temp[col_AgeomVIS]),
                               f=float(line_temp[col_f]),
                               dist=float(line_temp[col_dist]),
                               Rs=float(line_temp[col_Rs]),
                               Ts=float(line_temp[col_Ts]),
                               nodes=None,
                               trans=None,
                               lam_eff=None,
                               W_eff=None,
                               mission=mission,
                               mags=SC['vmag'][int(line_temp[col_nstar])],
                               F0=F0s[0])
                logfile = open(temp_name, 'a')
                logfile.write('%018.12f\t' % (flux.s_uJy())+'%018.12f\t' % (flux.p_therm_uJy())+'%018.12f\t' % (flux.p_refl_uJy()))
                logfile.close()
                
                flux = pp.flux(Rp=float(line_temp[col_Rp]),
                               Tp=float(line_temp[col_Tp]),
                               rp=float(line_temp[col_rp]),
                               AgeomMIR=float(line_temp[col_AgeomMIR]),
                               AgeomVIS=float(line_temp[col_AgeomVIS]),
                               f=float(line_temp[col_f]),
                               dist=float(line_temp[col_dist]),
                               Rs=float(line_temp[col_Rs]),
                               Ts=float(line_temp[col_Ts]),
                               nodes=None,
                               trans=None,
                               lam_eff=None,
                               W_eff=None,
                               mission=mission,
                               mags=SC['jmag'][int(line_temp[col_nstar])],
                               F0=F0s[1])
                logfile = open(temp_name, 'a')
                logfile.write('%018.12f\t' % (flux.s_uJy())+'%018.12f\t' % (flux.p_therm_uJy())+'%018.12f\t' % (flux.p_refl_uJy()))
                logfile.close()
                
                flux = pp.flux(Rp=float(line_temp[col_Rp]),
                               Tp=float(line_temp[col_Tp]),
                               rp=float(line_temp[col_rp]),
                               AgeomMIR=float(line_temp[col_AgeomMIR]),
                               AgeomVIS=float(line_temp[col_AgeomVIS]),
                               f=float(line_temp[col_f]),
                               dist=float(line_temp[col_dist]),
                               Rs=float(line_temp[col_Rs]),
                               Ts=float(line_temp[col_Ts]),
                               nodes=None,
                               trans=None,
                               lam_eff=None,
                               W_eff=None,
                               mission=mission,
                               mags=SC['hmag'][int(line_temp[col_nstar])],
                               F0=F0s[2])
                logfile = open(temp_name, 'a')
                logfile.write('%018.12f\t' % (flux.s_uJy())+'%018.12f\t' % (flux.p_therm_uJy())+'%018.12f\t' % (flux.p_refl_uJy())+'\n')
                logfile.close()
            else:
                print('Line %.0f: unappropriate data' % i)
                import pdb; pdb.set_trace()
            
            i += 1
            if (i % 1000 == 0):
                print('Processed '+str(i)+' planets')
        planets.close()
