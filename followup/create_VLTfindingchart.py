#!/usr/bin/env python

"""create_VLTfindingchart.py -- Input RA and DEC, 6' cutout, 2' cutout, 2' cutout's catalog (of stars with RA, DEC, g_mag, r_mag) and output finding chart suitable for VLT FORS2 LSS observations.

Usage: create_VLTfindingchart [-h] [-v] [--debug] [-o SAVELOC] [--pa FLOAT] [--norefstars INT] [--PI STRING] [--RunID STRING] [--OBname STRING] [--fitswaveband STRING] <RADEC> <6arcminfits> <2arcminfits> <catgaia>

Arguments:
    RADec (string)
        In decimals. E.g. 121.21,-78.31
    6arcminfits (string)
        Full path to 6arcmin fits (tested using DECam data)
    2arcminfits (string)
        Full path to 2arcmin fits (tested using DECam data). This should be cutout of 6arcminfits
    catgaia (string)
        Full path to catalog (of stars with RA, DEC, g_mag, r_mag) in the 2arcmin fits (cross matched with gaia)

Options:
    -h, --help                      Show this screen
    -v, --verbose                   Show extra information [default: False]     
    --debug                         Output more for debugging [default: False]
    -o SAVELOC, --out SAVELOC       Saved output as [default: ./finding_chart.jpg]
    --pa FLOAT                      Position angle of slit [default: 0.0]
    --norefstars INT                Number of reference stars to plot [default: 8]
    --PI STRING                     PI for VLT proposal (info to be printed) [default: 'Jeff Cooke']
    --RunID STRING                  VLT Observing Run ID (info to be printed) [default: '0103.D-0625(A)']
    --OBname STRING                 Observation Block or Target Name (info to be printed) [default: Transient1]
    --fitswaveband STRING           Finding chart background image wavelength band (info to be printed) [default: g-band]

Examples:
"""

import docopt
import os
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import aplpy
import astropy.visualization as astrovis
from astropy import wcs
import astropy
import pandas as pd
from pandas.plotting import table 

def draw_line(plot, theta, length, ra, dec, color='b', linewidth=1, alpha=0.7):
    theta = theta*np.pi/180.0
    length = length/2.0
    dx = np.sin(theta)*length/(np.cos(dec*np.pi/180.0)*60.0)
    dy = np.cos(theta)*length/60.0
    coords = np.array([[ra+dx, ra-dx], [dec+dy, dec-dy]])
    plot.show_lines([coords], color=color, linewidth=linewidth, alpha=alpha)
    return plot

def read_cat_gaia(f_cat,number_fstars,verbose=False):
    # Read in catalog
    cat = np.genfromtxt(f_cat,names=True,skip_header=0)
    catRA  = cat['RA']
    catDec = cat['DEC']
    catG   = cat['Gmag']
    catB   = cat['BPmag']
    catR   = cat['RPmag']

    # Sort catalog entries by g-band mag
    catG_sorted = np.sort(catG)
    catB_sorted = [x for _,x in sorted(zip(catG,catB))]
    catR_sorted = [x for _,x in sorted(zip(catG,catR))]
    catRA_sorted = [x for _,x in sorted(zip(catG,catRA))]
    catDec_sorted = [x for _,x in sorted(zip(catG,catDec))]

    # Pick X brightest stars
    catG_sorted_top    = catG_sorted[0:number_refstars]
    catB_sorted_top    = catB_sorted[0:number_refstars]
    catR_sorted_top    = catR_sorted[0:number_refstars]
    catRA_sorted_top   = catRA_sorted[0:number_refstars]
    catDec_sorted_top  = catDec_sorted[0:number_refstars]

    # Calculate other information required
    catGR_sorted_top   = catG_sorted_top-catR_sorted_top
    RAoffset_arcsec     = (np.array(catRA_sorted_top)-ra)*60.*60.
    DECoffset_arcsec    = (np.array(catDec_sorted_top)-dec)*60.*60.

    return catG_sorted_top,catGR_sorted_top,catRA_sorted_top,catDec_sorted_top,RAoffset_arcsec,DECoffset_arcsec

def create_refstardf(catRA, catDec, catg, catgr, RAoffset_arcsec, DECoffset_arcsec, catalog='gaia', verbose=False):

    # Define file
    number_refstars = len(catRA)

    gband_name = 'mag (gaia G)'
    color_name = 'color (gaia G-R)'

    if number_refstars > 0:
        labels = ['#'+str(x+1) for x in range(number_refstars)]
        d = {'Label': labels, 
             'RA': catRA, 
             'DEC': catDec, 
             gband_name:catg, 
             color_name:catgr,
             'RAoffset (arcsec)':RAoffset_arcsec,
             'DECoffset (arcsec)':DECoffset_arcsec
            }
    else:
        d = {'Label':['x'],
                'RA':['x'],
                'DEC':['x'],
                gband_name:['x'],
                color_name:['x'],
                'RAoffset (arcsec)':['x'],
                'DECoffset (arcsec':['x']}
        print('\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        print('ISSUE: There were no reference stars found in the catalog entry\n')
    df_stars = pd.DataFrame(data=d)

    return df_stars

def create_OBinfodf(PI='Jeff Cooke',RunID='?',OBname='Transient1',fitsband='g-band',pa=0.0,verbose=False):
    # Define info table
    # Table of Observing Block information
    d = {'Item': ['Observing Run ID', 'PI name', 'OB/ Target Name','slit PA (deg)','image wavelength band'], 
         'Entry':[RunID,PI,OBname,pa,fitsband]}
    df_info = pd.DataFrame(data=d)
    return df_info

def make_chart(ra,dec,f_6arcmin,f_2arcmin,df_stars,df_info,catRA,catDec,pa=0.0,saveloc='./finding_chart.jpg',verbose=False):
    color='red'
    fig = plt.figure(figsize=(12,24))

    # ====== 6 arcmin image ===

    # Load fits
    d, h = fits.getdata(f_6arcmin,header=True)

    # Visualize fits
    interval = astrovis.ZScaleInterval()
    vmin,vmax=interval.get_limits(d)
    img = aplpy.FITSFigure(f_6arcmin, figure=fig, subplot=[0.1,0.6,0.9,0.4])
    img.show_grayscale(stretch='linear', vmin=vmin, vmax=vmax)

    # Mark RA and DEC target
    img.show_markers([ra],[dec],coords_frame='world',s=10)

    # Add scalebar
    scalebar_length = 1./60. #arcmin
    img.add_scalebar(length=scalebar_length)
    img.scalebar.set_label('1 arcmin')
    img.scalebar.set_linewidth(3)
    img.scalebar.set_font_size(15)
    img.scalebar.set_font_weight('medium')

    # Draw box for 2'x2' zoom in
    boxsize_arcsec=2.*60
    w = wcs.WCS(h)
    [[x1,y1]]=w.wcs_world2pix([ [ra,dec] ], 1)
    pixelscale_arcsec = astropy.wcs.utils.proj_plane_pixel_scales(w)*60.*60.
    boxsize_pixel     = boxsize_arcsec/pixelscale_arcsec[0]
    img.show_rectangles([x1], [y1], boxsize_pixel, boxsize_pixel,
                        coords_frame='pixel',edgecolor='red',linewidth=2)

    # ====== 2 arc min image ===

    # Visualize image
    img2 = aplpy.FITSFigure(f_2arcmin, figure=fig, subplot=[0.1,0.2,0.9,0.4]) 
    img2.show_grayscale(stretch='linear', vmin=vmin, vmax=vmax)

    # Draw slit
    draw_line(img2,0.,4,ra,dec,color='red')

    # Draw scalebar
    scalebar_length = 10./60./60. #arcmin
    img2.add_scalebar(length=scalebar_length)
    img2.scalebar.set_label('10 arcsec')
    img2.scalebar.set_linewidth(3)
    img2.scalebar.set_font_size(15)
    img2.scalebar.set_font_weight('medium')

    # Mark RA and DEC target
    img2.show_markers([ra],[dec],coords_frame='world',s=100,edgecolor=color,linewidth=2)
    img2.add_label(ra+0.006,dec+0.0006,'TARGET',color=color,size='xx-large',weight='bold')

    # Show acquisition specific reference stars
    number_refstars = len(catRA)
    labels = ['#'+str(x+1) for x in range(number_refstars)]
    for refstar_ra,refstar_dec,label in zip(catRA,
                                                          catDec,
                                                          labels):
        img2.show_markers([refstar_ra],[refstar_dec],coords_frame='world',s=10,edgecolor=color,linewidth=3)
        img2.add_label(refstar_ra+0.002,refstar_dec+0.0002,label,color=color,size='xx-large',weight='bold')

    # ====== Tables to be printed ===

    # Table of acquisition reference stars
    ax2 = fig.add_axes([0.1,0.1,0.9,0.06]) 
    ax2.xaxis.set_visible(False)  # hide the x axis
    ax2.yaxis.set_visible(False)  # hide the y axis
    table(ax2, df_stars, rowLabels=['']*df_stars.shape[0], loc='center')  # where df is your data frame

    # Table of information needed to be printed
    ax3 = fig.add_axes([0.1,0.05,0.9,0.04]) 
    ax3.xaxis.set_visible(False)  # hide the x axis
    ax3.yaxis.set_visible(False)  # hide the y axis
    table(ax3, df_info, rowLabels=['']*df_info.shape[0], loc='center')  # where df is your data frame

    # Save figure
    fig.savefig(saveloc, bbox_inches='tight')
    full_path = os.path.abspath(saveloc)
    print('Finding chart saved: ',full_path)
    scpcommand='scp fstars@ozstar.swin.edu.au:'+full_path+' ./finding_chart.jpg'
    print(scpcommand)
    scpcommand='scp fstars@ozstar.swin.edu.au:'+full_path+' '+full_path.split('/')[-1]
    print(scpcommand)

    return scpcommand

def create_VLTfindingchart(ra,dec,f_6arcmin,f_2arcmin,f_catgaia,pa=0.0,number_refstars=8,
                            PI='Jeff Cooke',RunID='?',OBname='Transient1',fitsband='g-band',
                            saveloc='./finding_chart.jpg',
                            verbose=False):

    # Read in catalog
    catG,catGR,catRA,catDec,RAoffset_arcsec,DECoffset_arcsec = read_cat_gaia(f_catgaia,number_refstars,verbose=verbose)

    # Create pd table for acquisition specific reference stars
    df_stars = create_refstardf(catRA, catDec, catG, catGR, RAoffset_arcsec, DECoffset_arcsec, catalog='gaia', verbose=verbose)        
    
    # Create pd table of information
    df_info = create_OBinfodf(PI=PI,RunID=RunID,OBname=OBname,fitsband=fitsband,pa=pa,verbose=verbose) 

    # Create finding chart
    scpcommand = make_chart(ra,dec,f_6arcmin,f_2arcmin,df_stars,df_info,catRA,catDec,pa=pa,saveloc=saveloc,verbose=verbose)

    # Finish up
    print('#################################')
    print('Congrats, finding chart made! ')
    print('\nUse this to get your finding chart on whatever computer you want!\n')
    print(scpcommand)
    print('#################################')

    return None

####################### BODY OF PROGRAM STARTS HERE ########################

if __name__ == "__main__":

    # Read in arguments
    arguments       = docopt.docopt(__doc__)

    # Mandatory arguments
    RADEC           = arguments['<RADEC>']  
    ra              = float(RADEC.split(',')[0])
    dec             = float(RADEC.split(',')[1])
    f_6arcmin       = arguments['<6arcminfits>']
    f_2arcmin       = arguments['<2arcminfits>']
    f_catgaia       = arguments['<catgaia>']

    # Optional arguments
    verbose         = arguments['--verbose']
    debugmode       = arguments['--debug']
    saveloc         = arguments['--out']
    pa              = float(arguments['--pa'])
    number_refstars = int(arguments['--norefstars'])
    PI              = arguments['--PI']
    RunID           = arguments['--RunID']
    OBname          = arguments['--OBname']
    fitsband        = arguments['--fitswaveband']

    if debugmode:
        print(arguments)  

    create_VLTfindingchart(ra,dec,f_6arcmin,f_2arcmin,f_catgaia,pa=pa,number_refstars=number_refstars,
                            PI=PI,RunID=RunID,OBname=OBname,fitsband=fitsband,
                            saveloc=saveloc,
                            verbose=verbose)
