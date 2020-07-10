"""
create_findingChartInputs.py  -- Create a small cutout and a larger cutout of an input template image around an input RA and DEC. Also creates a catalog of stars overlapping with the cutout using the gaia catalog, downloaded within this script from Vizier. 

Usage: create_findingChartInputs [-h] [-v] [--debug] [--findingchartdir DIRECTORY] <RA> <DEC> <field> <source_name> <templateimage_path>

Arguments:
        RA (float)
            The RA in degrees.
        DEC (float)
            The DEC in degrees.
        field (string)
            The name of the DWF field. Only used to save output files in out directory/field/source_name
        source_name (string)
            The name you would like to call your source in the outputs. 
        templateimage_path (string)
            The template fits file you want to use for the finding chart.

Outputs:
        Saved in <out directory>/source_name/ (This overwrites)
            Small cutout (2' x 2' if DECam pixelscale) of template image centered around input RA and DEC
            Larger cutout (6' x 6' if DECam pixelscale) of template image centered around input RA and DEC
            Catalog of stars in the small cutout found in the gaia catalog

Options:
        -h, --help                        Show this screen
        -v, --verbose                     Show extra information [default: False]     
        --debug                           Output more for debugging [default: False]
        --findingchartdir DIRECTORY       Outputs are saved here/field/ [default: /fred/oz100/FINDING_CHARTS/]

Example:
       python create_VLTfindingchartInputs.py -v 246.50 -73.0 ngc6101 jielaitest /fred/oz100/pipes/DWF_PIPE/TEMPLATES/ngc6101/REQ16A_20160727_8c29d52-g-20160728T041310_si1_dither_seeing12_join.fits
"""
import docopt
import numpy as np
from astropy.table import Table 
from astropy import wcs
from astropy.wcs import WCS
from astropy.io import fits
import os
import datetime as dt
import os
from astropy.nddata.utils import Cutout2D
from astroquery.vizier import Vizier
from astropy.coordinates import Angle

def create_finding_chart_inputs(RA,DEC,
                                field,source_name,
                                templateimage_path,
                                findingchartdir='/fred/oz100/FINDING_CHARTS/',
                                verbose=False, debugmode=False):
    print('\n#############################################')
    print('Creating files needed to create finding chart image for RA: {} DEC: {}'.format(RA,DEC))
    print('#############################################\n')

    # Let's keep things organized, so always put the outputs of this script here:
    outputdir = findingchartdir+'/'+field+'/'+source_name+'/'

    if verbose:
        print('Template image used: {}\n'.format(templateimage_path))

    if not os.path.exists(outputdir):
        os.makedirs(outputdir, 0o755)
    else:
        pass

    with fits.open(templateimage_path) as hdu:
        size                    = 1370
        w                       = WCS(hdu[0].header)
        head                    = hdu[0].header
        date                    = dt.datetime.strptime(head['DATE'], '%Y-%m-%dT%H:%M:%S')
        xlim                    = head['NAXIS1']
        ylim                    = head['NAXIS2']
        pixcrd_im               = np.array([[xlim, ylim]], np.float_)
        world_im                = w.wcs_pix2world(pixcrd_im, 1)
        pixx_im, pixy_im        = world_im[0][0], world_im[0][1]

        pixcrd                  = np.array([[RA, DEC]], np.float_)
        worldpix                = w.wcs_world2pix(pixcrd, 1)
        pixx, pixy              = worldpix[0][0], worldpix[0][1]
    
        cutout                  = Cutout2D(hdu[0].data, (pixx, pixy), size, wcs= w)
        hdu[0].data             = cutout.data
        hdu[0].header['CRPIX1'] = cutout.wcs.wcs.crpix[0]
        hdu[0].header['CRPIX2'] = cutout.wcs.wcs.crpix[1]
        outfile_large           = outputdir +'/'+ source_name+'6arcmin.fits'
        hdu.writeto(outfile_large, overwrite = True)
        print('========SAVED===Large fits cutout file saved: \n     {}'.format(outfile_large))
  
    with fits.open(templateimage_path) as hdu:
        size                    = 457
        w                       = WCS(hdu[0].header)
        head                    = hdu[0].header
        date                    = dt.datetime.strptime(head['DATE'], '%Y-%m-%dT%H:%M:%S')
        xlim                    = head['NAXIS1']
        ylim                    = head['NAXIS2']
        pixcrd_im               = np.array([[xlim, ylim]], np.float_)
        world_im                = w.wcs_pix2world(pixcrd_im, 1)
        pixx_im, pixy_im        = world_im[0][0], world_im[0][1]

        pixcrd                  = np.array([[RA, DEC]], np.float_)
        worldpix                = w.wcs_world2pix(pixcrd, 1)
        pixx, pixy              = worldpix[0][0], worldpix[0][1]
        cutout                  = Cutout2D(hdu[0].data, (pixx, pixy), size, wcs= w)
        hdu[0].data             = cutout.data
        hdu[0].header['CRPIX1'] = cutout.wcs.wcs.crpix[0]
        hdu[0].header['CRPIX2'] = cutout.wcs.wcs.crpix[1]
        outfile_small           = outputdir + '/' + source_name+'2arcmin.fits'
        hdu.writeto(outfile_small, overwrite = True)
        print('========SAVED===Small fits cutout file saved: \n     {}'.format(outfile_small))
    
    # Create a catalog of sources within the field of view. 
    if debugmode:
        print('Starting to search for GAIA sources')
    GAIA_DR2 = "I/345"
    RA_DEC = str(f'{RA} {DEC}')

    result = Vizier.query_region(RA_DEC, radius=Angle('120"'), catalog=GAIA_DR2)
    GAIA_DR2          = Table()
    GAIA_DR2['RA']    = result[0]['RA_ICRS']
    GAIA_DR2['DEC']   = result[0]['DE_ICRS']
    GAIA_DR2['Gmag']  = result[0]['Gmag']
    GAIA_DR2['BPmag'] = result[0]['BPmag']
    GAIA_DR2['RPmag'] = result[0]['RPmag']
    outfile_cat        = outputdir + '/' + source_name + '_gaia_star_cat.ascii'
    GAIA_DR2.write(outfile_cat, format='ascii', overwrite=True)
    print('========SAVED===Small cutout star catalog saved: \n     {}'.format(outfile_cat))

    print('\nTo Create Finding Chart, run the following in terminal:')
    print('--pa 0.0 (position angle of slit)')
    print("--PI 'Jeff Cooke'")
    print("--RunID '0103.D-0625(A)'   (Take from proposal, VLT)")
    print("--OBname Transient1    (Observation Block or Target Name)")
    print("--fitswaveband g-band   (for background image)")
    print('RADEC={},{}'.format(RA,DEC))
    print('largefits={}'.format(outfile_large))
    print('smallfits={}'.format(outfile_small))
    print('catfile={}'.format(outfile_cat))
    print('outfile={}'.format(outputdir+'/'+source_name+'_finding_chart.jpg'))
    print('python ~/src/dwftools/FOLLOWUP/create_VLTfindingchart.py $RADEC $largefits $smallfits $catfile -o $outfile') 

    print('\n#############################################')
    print('# YOUR FINDING CHART INPUTS ARE DONE')
    print('# FIND THEM HERE: {}'.format(outputdir))
    print('#############################################\n')


if __name__ == "__main__":
    ## read in arguments 
    arguments = docopt.docopt(__doc__, options_first=True)

    # Mandatory arguments
    RA                  = arguments['<RA>']
    DEC                 = arguments['<DEC>']
    field               = arguments['<field>']
    source_name         = arguments['<source_name>']
    templateimage_path  = arguments['<templateimage_path>']
    
    # Optional arguments
    verbose         = arguments['--verbose']
    debugmode       = arguments['--debug']
    findingchartdir  = arguments['--findingchartdir']

    if debugmode:
            print(arguments)

    create_finding_chart_inputs(RA,DEC,
                                field,source_name,
                                templateimage_path,
                                findingchartdir=findingchartdir,
                                verbose=verbose,debugmode=debugmode)
    
