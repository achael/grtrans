# run_grtrans_3D.py
# runs grtrans on koral simulation images at very high resolution
# these will take a long time!

import grtrans_batch as gr

import numpy as np
import sys
import os
import subprocess
import time
import astropy.io.fits as fits
import scipy.ndimage.interpolation as interpolation

####################
#Constants
####################
pcG = 6.67259e-8
pcc2 = 8.98755179e20
msun = 1.99e33
cmperkpc=3.086e21

EP = 1.0e-10
C = 299792458.0
DEGREE = 3.141592653589/180.0
HOUR = 15.0*DEGREE
RADPERAS = DEGREE/3600.0
RADPERUAS = RADPERAS*1.e-6

####################
# Problem setup
####################
SOURCE = 'M87'
RA = 12.51373
DEC = 12.39112
MJD =  58211

hfile = './ana0990.h5' #jet, a=0.25
dfile = './ana0990.h5'
SPIN = 0.25
RESCALE = 2.46e-16
NPIX_IM = 128

#hfile = './ana1000.h5' #mks2, a=0.9375
#dfile = './ana1000.h5'
#SPIN=0.9375 
#RESCALE = 1.e-18 
#NPIX_IM = 128

MBH=6.5e9   # this was fixed in the simulation, but EHT measured 6.5
DTOBH = 16*1.e3*cmperkpc 
TGPERFILE = 10 # gravitational times per file


NGEO = 800
RAYTRACESIZE=50. # raytracing volume. For M87 at 17 degrees, 1000 rg = 1 milliarcsec


RERUN = True # rerun even if the output file exists
ROTATE = False # rotate the image before saving
ANGLE = 108 # rotation angle

# parameters to loop over (in serial, in this script)

ang = 20.            # inclination angle
sigma_cut = 10.
freq_ghz = 230.     # frequency
pixel_size_uas = 2  # microarcseconds per pixel

# derived quantities
lbh = pcG*msun*MBH / pcc2 
fac= (4*np.pi*lbh**2)
LumtoJy = 1.e23/(4*np.pi*DTOBH**2)
MuasperRg = lbh / DTOBH / RADPERUAS
secperg = MBH*4.927e-6 #seconds per tg 
hourperg = secperg / 3600.
weekperg = secperg / (604800.)
yrperg = secperg / (3.154e7)

size =  20#0.5*NPIX_IM*pixel_size_uas/MuasperRg


####################
# Functions
####################
def main():

    name = './test_hdf5'

    # skip over if output already exists, or delete and rerun

    # write input radiative transfer parameters
    mu = np.cos(ang*np.pi/180.)
    freq = freq_ghz*1.e9
    uout = 1./RAYTRACESIZE
    x=gr.grtrans()
    npol=1
    x.write_grtrans_inputs(name + '.in', oname=name+'.out',
                           fscalefac=RESCALE, sigcut=sigma_cut,
                           fname='KORALH5',phi0=0.,
                           nfreq=1,fmin=freq,fmax=freq,
                           ename='SYNCHTHAV',
                           nvals=npol,
                           gmin=-3, # confusingly, this is rhigh. #-2-->vladimir model, -3--> michael model
                           spin=SPIN,standard=1,
                           uout=uout,
                           mbh=MBH, 
                           mdotmin=1.,mdotmax=1.,nmdot=1,
                           #mdotmin=1.57e15,mdotmax=1.57e15,nmdot=1,
                           nmu=1,mumin=mu,mumax=mu,
                           gridvals=[-size,size,-size,size],
                           nn=[NPIX_IM,NPIX_IM,NGEO],
                           hhfile=hfile, hdfile=dfile,
                           hindf=1,hnt=1,
                           muval=1.)
    run=True
    if os.path.exists(name+'.out'):
        run = False
        if RERUN: 
            run=True
            os.remove(name+'.out') 
         
    # run grtrans 
    if run:
        x.run_grtrans()         

    # Read grtrans output
    try:    
        x.read_grtrans_output()
    except:# IOError:
        return None

    # pixel sizes
    da = x.ab[x.nx,0]-x.ab[0,0]
    db = x.ab[1,1]-x.ab[0,1]
    if (da!=db): raise Exception("pixel da!=db")
    psize = da*(lbh/DTOBH)

    #image values
    if npol==4:
        ivals = x.ivals[:,0,0]*fac*da*db*LumtoJy
        qvals = x.ivals[:,1,0]*fac*da*db*LumtoJy
        uvals = x.ivals[:,2,0]*fac*da*db*LumtoJy
        vvals = x.ivals[:,3,0]*fac*da*db*LumtoJy

        # mask nan failure points with zeros
        ivals = np.array(ivals)
        qvals = np.array(qvals)
        uvals = np.array(uvals)
        vvals = np.array(vvals)

        imask = np.isnan(ivals)
        qumask = ~(~imask * ~np.isnan(qvals) * ~np.isnan(uvals))
        vmask = ~(~imask * ~np.isnan(vvals))

        ivals[imask] = 0.
        qvals[qumask] = 0.
        uvals[qumask] = 0.
        vvals[vmask] = 0.

        ivals =  (np.flipud(np.transpose(ivals.reshape((NPIX_IM,NPIX_IM))))).flatten()
        qvals =  -(np.flipud(np.transpose(qvals.reshape((NPIX_IM,NPIX_IM))))).flatten()
        uvals =  -(np.flipud(np.transpose(uvals.reshape((NPIX_IM,NPIX_IM))))).flatten()
        vvals =  (np.flipud(np.transpose(vvals.reshape((NPIX_IM,NPIX_IM))))).flatten()
    else:
        ivals = x.ivals[:,0,0]*fac*da*db*LumtoJy
        ivals = np.array(ivals)
        imask = np.isnan(ivals)
        ivals[imask] = 0.
        ivals =  (np.flipud(np.transpose(ivals.reshape((NPIX_IM,NPIX_IM))))).flatten()
        qvals = uvals = vvals = 0*ivals
   
    print('total flux', np.sum(ivals))
    save_im_fits((ivals,qvals,uvals,vvals, psize))
    return

def save_im_fits(imdata,  mjd=MJD, source=SOURCE, ra=RA, dec=DEC, time=0):

    fname = 'hdf5test.fits'

    # unpack the image data
    (imvec, qvec, uvec, vvec, psize) = imdata

    # Create header and fill in some values
    header = fits.Header()
    header['OBJECT'] = source
    header['CTYPE1'] = 'RA---SIN'
    header['CTYPE2'] = 'DEC--SIN'
    header['CDELT1'] = -psize/DEGREE
    header['CDELT2'] =  psize/DEGREE
    header['OBSRA'] = ra * 180/12.
    header['OBSDEC'] = dec
    header['FREQ'] = freq_ghz * 1.e9

    #TODO these are the default values for centered images
    #TODO support for arbitrary CRPIX? 
    header['CRPIX1'] = NPIX_IM/2. + .5
    header['CRPIX2'] = NPIX_IM/2. + .5

    mjd += (time/24.)

    header['MJD'] = float(mjd)
    header['TELESCOP'] = 'VLBI'
    header['BUNIT'] = 'JY/PIXEL'
    header['STOKES'] = 'I'

    # Create the fits image
    image = np.reshape(imvec,(NPIX_IM,NPIX_IM))[::-1,:] #flip y axis!
    hdu = fits.PrimaryHDU(image, header=header)
    hdulist = [hdu]
    if len(qvec):
        qimage = np.reshape(qvec,(NPIX_IM,NPIX_IM))[::-1,:]
        uimage = np.reshape(uvec,(NPIX_IM,NPIX_IM))[::-1,:]
        header['STOKES'] = 'Q'
        hduq = fits.ImageHDU(qimage, name='Q', header=header)
        header['STOKES'] = 'U'
        hduu = fits.ImageHDU(uimage, name='U', header=header)
        hdulist = [hdu, hduq, hduu]
    if len(vvec):
        vimage = np.reshape(vvec,(NPIX_IM,NPIX_IM))[::-1,:]
        header['STOKES'] = 'V'
        hduv = fits.ImageHDU(vimage, name='V', header=header)
        hdulist.append(hduv)

    hdulist = fits.HDUList(hdulist)

    # Save fits
    try:
        hdulist.writeto(fname, overwrite=True)
    except:
        hdulist.writeto(fname, clobber=True)
    return



if __name__=='__main__':

    main()


