"""
p_pop, a Python package to simulate exoplanet populations based on Kepler
statistics. Information about the method can be found in this paper:
http://adsabs.harvard.edu/abs/2018A%26A...609A...4K. This library is maintained
on GitHub at https://github.com/kammerje/p_pop.

Author: Jens Kammerer
Version: 3.0.0
Last edited: 05.04.19
"""


# PREAMBLE
#==============================================================================
import astropy.table as at
import astropy.io as io
import csv
import matplotlib.pyplot as plt
import numpy as np


# MAIN
#==============================================================================
def create_filter_rect(tag='F560W_rect',
                       lam_eff=5.6,
                       W_eff=1.2):
    """
    Parameters
    ----------
    tag : str
        Name of filter
    lam_eff : float
        Effective wavelength (microns)
    W_eff : float
        Effective width (microns)
    """
    
    print('Started create_filter_rect...')
    
    # Check parameters
    if (isinstance(tag, str) == False):
        raise TypeError('tag must be of type string')
    lam_eff = float(lam_eff)
    if (lam_eff <= 0):
        raise ValueError('lam_eff must be positive')
    W_eff = float(W_eff)
    if (W_eff <= 0):
        raise ValueError('W_eff must be positive')
    
    # Create lists for nodes and transmission
    x = np.linspace(lam_eff-0.5*W_eff, lam_eff+0.5*W_eff, 100)
    y = np.array([1.]*100)
    
    plt.figure(figsize=(12, 9))
    plt.plot(x, y, label=tag)
    plt.xlabel('Wavelength [microns]')
    plt.ylabel('Transmission')
    plt.legend()
    plt.title('Filter curve')
    plt.show(block=True)
    
    # Create and save output file
    out = [x, y, lam_eff, W_eff]
    np.save('data/'+tag, out)
    
    print('Finished create_filter_rect')
    
    # Return nothing
    return None

def read_20pc_bright_sample(path='data/20pc_bright_sample.tbl',
                            s_type=['A', 'F', 'G', 'K', 'M'],
                            max_dist=20,
                            min_dec=-90,
                            max_dec=90):
    """
    Parameters
    ----------
    path : str
        Path of "20pc_bright_sample.tbl"
    s_type : list
        List containing spectral types which should be considered
    max_dist : float
        Maximal distance which should be considered (parsec)
    min_dec : float
        Minimal declination which should be considered (degree)
    max_dec : float
        Maximal declination which should be considered (degree)
    
    Returns
    -------
    SC : table
        Star catalog in format of an astropy table
    """
    
    print('Started read_20pc_bright_sample...')
    
    # Check parameters
    if (isinstance(s_type, list) == False):
        raise TypeError('s_type must be of type list')
    for i in range(0, len(s_type)):
        if (isinstance(s_type[i], str) == False):
            raise TypeError('elements of s_type must be of type str')
        if (not s_type[i] in ['A', 'F', 'G', 'K', 'M']):
            raise ValueError('elements of s_type must be A, F, G, K or M')
    max_dist = float(max_dist)
    if (max_dist <= 0):
        raise ValueError('max_dist must be positive')
    min_dec = float(min_dec)
    if (min_dec < -90 or min_dec > 90):
        raise ValueError('min_dec must be between -90 and 90')
    max_dec = float(max_dec)
    if (max_dec < -90 or max_dec > 90):
        raise ValueError('max_dec must be between -90 and 90')
    if (min_dec >= max_dec):
        raise ValueError('min_dec must be smaller than max_dec')
    
    # Read input star catalog
    SC_in = at.Table.read(path, format = 'ipac')
    
    # Create output star catalog
    SC_out = at.Table(names = ('name', 'dist', 'stype', 'rad', 'teff', 'mass', 'ra', 'dec'), dtype = ('S32', 'd', 'c', 'd', 'd', 'd', 'd', 'd'))
    
    # Fill output star catalog
    for i in range(len(SC_in)):
        
        # Check whether stellar type fits
        if (str(SC_in[i]['Sptype'])[0] in s_type):
            
            # Check whether distance fits
            if (1./(1E-3*float(SC_in[i]['pi_trig'])) <= max_dist):
                
                # Convert declination from degrees-minutes-seconds to degrees
                dec = int(SC_in[i]['DEC_2000'][0:3])
                if (dec >= 0):
                    dec = float(dec)+int(SC_in[i]['DEC_2000'][4:6])/60.+float(SC_in[i]['DEC_2000'][7:14])/3600.
                else:
                    dec = float(dec)-int(SC_in[i]['DEC_2000'][4:6])/60.-float(SC_in[i]['DEC_2000'][7:14])/3600.
                
                # Check whether declination fits
                if (min_dec <= dec <= max_dec):
                    
                    # Convert right ascension from hours-minutes-seconds to degrees
                    ra = int(SC_in[i]['RA_2000'][0:2])*360./24.+int(SC_in[i]['RA_2000'][3:5])/60.+float(SC_in[i]['RA_2000'][6:13])/3600.
                    
                    SC_out.add_row([SC_in[i]['DiscoveryName'], 1./(1E-3*float(SC_in[i]['pi_trig'])), SC_in[i]['Sptype'], SC_in[i]['r'], SC_in[i]['teff'], SC_in[i]['m'], ra, dec])
    
    # Save output star catalog
    io.ascii.write(SC_out, 'data/SC_out.tbl', format = 'ipac')
    
    print('Finished read_20pc_bright_sample')
    
    # Return output star catalog
    return SC_out

def read_filter_MIRI(path='data/F560W.csv',
                     tag='F560W',
                     lam_eff=5.6,
                     W_eff=1.2):
    """
    Parameters
    ----------
    path : str
        Path of MIRI filter curve
    tag : str
        Name of filter
    lam_eff : float
        Effective wavelength (microns)
    W_eff : float
        Effective width (microns)
    """
    
    print('Started read_filter_MIRI...')
    
    # Check parameters
    if (isinstance(path, str) == False):
        raise TypeError('path must be of type string')
    if (isinstance(tag, str) == False):
        raise TypeError('tag must be of type string')
    lam_eff = float(lam_eff)
    if (lam_eff <= 0):
        raise ValueError('lam_eff must be positive')
    W_eff = float(W_eff)
    if (W_eff <= 0):
        raise ValueError('W_eff must be positive')
    
    # Create lists for nodes and transmission
    x = []
    y = []
    
    # Open file
    with open(path) as csvfile:
        data = csv.reader(csvfile, delimiter = ';')
        
        # Go through all rows
        for row in data:
            
            # Fill lists for nodes and transmission
            x += [float(row[0])]
            y += [float(row[2])/100.]
    
    # Convert lists for nodes and transmission to arrays
    x = np.array(x)
    y = np.array(y)
    
    # Define masks for arrays for nodes and transmission
    mask1 = (lam_eff-2.*W_eff < x)
    mask2 = (x < lam_eff+2.*W_eff)
    
    # Apply masks for arrays for nodes and transmission
    x = x[mask1 & mask2]
    y = y[mask1 & mask2]
    
#    plt.figure(figsize=(12, 9))
#    plt.plot(x, y, label=tag)
#    plt.xlabel('Wavelength [microns]')
#    plt.ylabel('Transmission')
#    plt.legend()
#    plt.title('Filter curve')
#    plt.show(block=True)
    
    # Create and save output file
    out = [x, y, lam_eff, W_eff]
    np.save('data/'+tag, out)
    
    print('Finished read_filter_MIRI')
    
    # Return nothing
    return None

def read_filter_IRDIS(path='data/Vband.txt',
                      tag='Vband',
                      lam_eff=0.7354,
                      W_eff=0.2905):
    """
    Parameters
    ----------
    path : str
        Path of MIRI filter curve
    tag : str
        Name of filter
    lam_eff : float
        Effective wavelength (microns)
    W_eff : float
        Effective width (microns)
    """
    
    print('Started read_filter_IRDIS...')
    
    # Check parameters
    if (isinstance(path, str) == False):
        raise TypeError('path must be of type string')
    if (isinstance(tag, str) == False):
        raise TypeError('tag must be of type string')
    lam_eff = float(lam_eff)
    if (lam_eff <= 0):
        raise ValueError('lam_eff must be positive')
    W_eff = float(W_eff)
    if (W_eff <= 0):
        raise ValueError('W_eff must be positive')
    
    # Create lists for nodes and transmission
    x = []
    y = []
    
    # Open file
    f = open(path, 'r')
    f = f.readlines()
    
    # Go through all lines
    for line in f[8:]:
        
        # Fill lists for nodes and transmission
        x += [float(line.split()[0])/float(1000)]
        y += [float(line.split()[1])]
    
    # Convert lists for nodes and transmission to arrays
    x = np.array(x)
    y = np.array(y)
    
#    plt.figure(figsize=(12, 9))
#    plt.plot(x, y, label=tag)
#    plt.xlabel('Wavelength [microns]')
#    plt.ylabel('Transmission')
#    plt.legend()
#    plt.title('Filter curve')
#    plt.show(block=True)
    
    # Create and save output file
    out = [x, y, lam_eff, W_eff]
    np.save('data/'+tag, out)
    
    print('Finished read_filter_IRDIS')
    
    # Return nothing
    return None
