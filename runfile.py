"""
p_pop, a Python package to simulate exoplanet populations based on Kepler
statistics. Information about the method can be found in this paper:
http://adsabs.harvard.edu/abs/2018A%26A...609A...4K. This library is maintained
on GitHub at https://github.com/kammerje/p_pop.

Author: Jens Kammerer
Version: 3.0.1
Last edited: 22.05.19
"""


# PREAMBLE
#==============================================================================
import matplotlib.pyplot as plt
import numpy as np

import inout as io
import p_pop as pp


# PARAMETERS
#==============================================================================
logfile_name = 'SAG13_planets5000.txt' # Name of the output table containing all the
                             # exoplanets
nMC = 5000 # Number of Monte-Carlo shoots
model = 0 # Scenario for the planet occurrence rates:
          # 0 - baseline
          # 1 - minimum
          # 2 - maximum
rates = 'sag13_dressing2015' # Planet occurrence rates:
                # sag13 - NASA SAG13, Kopparapu et al. 2018
                # hsu2019 - Hsu et al. 2019
                # dressing2015 - Dressing & Charbonneau 2015
                # fressin2013 - Fressin et al. 2013
                #
                # sag13_dressing2015 - sag13 (AFGK) & dressing2015 (M)
                # hsu2019_dressing2015 - hsu2019 (AFGK) & dressing2015 (M)
albedo = 'MIR' # Geometric albedo which is used to compute the reflected host star light
               # MIR - infrared
               # VIS - visible

# MAIN
#==============================================================================

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

SC = io.read_20pc_bright_sample(path='data/20pc_bright_sample.tbl',
                                s_type=['A', 'F', 'G', 'K', 'M'], # Spectral types which should be put into the star catalog
                                max_dist=20, # Distance cut (pc) for the star catalog
                                min_dec=-90, # Declination cut (deg) for the star catalog
                                max_dec=90) # Declination cut (deg) for the star catalog

#==============================================================================
# GENERATE PLANET POPULATION

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
                           rates=rates)
        
        for k in range(System.nplanets):
            Exoplanet = System.planets[k]
            
            logfile = open(logfile_name, 'a')
            logfile.write('%.0f\t' % j+'%.5f\t' % Exoplanet.Rp+'%.5f\t' % Exoplanet.Porb+'%.5f\t' % Exoplanet.a()+'%.5f\t' % Exoplanet.rp()+'%.5f\t' % Exoplanet.ang_sep()+'%.5f\t' % Exoplanet.ang_sep_max()+'%.5f\t' % Exoplanet.inc+'%.5f\t' % Exoplanet.Omega+'%.5f\t' % Exoplanet.omega+'%.5f\t' % Exoplanet.theta+'%.5f\t' % Exoplanet.ecc+'%.5f\t' % Exoplanet.Finc()+'%.5f\t' % Exoplanet.Abond+'%.5f\t' % Exoplanet.AgeomMIR+'%.5f\t' % Exoplanet.AgeomVIS+'%.5f\t' % Exoplanet.f()+'%.5f\t' % Exoplanet.Tp()+'%.5f\t' % Exoplanet.Mp+'%.5f\t' % Exoplanet.dist+'%.5f\t' % Exoplanet.Rs+'%.5f\t' % Exoplanet.Ts+'%.5f\t' % Exoplanet.Ms+stype+'\t%.5f\t' % Exoplanet.zodis+'%.5f\t' % ra+'%.5f\t' % dec+'%.0f\t' % i+'\n')
            logfile.close()

import pdb; pdb.set_trace()

#==============================================================================
# COMPUTE FLUXES

# Read LIFE filters
name1 = 'F560W'
name2 = 'F1000W'
name3 = 'F1500W'
filter1 = np.load('data/'+name1+'.npy')
filter2 = np.load('data/'+name2+'.npy')
filter3 = np.load('data/'+name3+'.npy')

logfile = open(logfile_name[:-4]+'_LIFE.txt', 'w')
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
                       albedo=albedo)
        logfile = open(logfile_name, 'a')
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
                       albedo=albedo)
        logfile = open(logfile_name, 'a')
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
                       albedo=albedo)
        logfile = open(logfile_name, 'a')
        logfile.write('%018.12f\t' % (flux.bb_therm_s_int()*1E6)+'%018.12f\t' % (flux.bb_therm_p_int()*1E6)+'%018.12f\t' % (flux.refl_p_int()*1E6)+'\n')
        logfile.close()
    else:
        print('Line %.0f: unappropriate data' % i)
        import pdb; pdb.set_trace()
    
    i += 1
    if (i % 1000 == 0):
        print('Processed '+str(i)+' planets')
planets.close()

import pdb; pdb.set_trace()
