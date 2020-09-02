# Run Grtrans with rrjet model
# The rrjet model is defined in "fluid_model_rrjet.py"
# NOTE -- currently the power law emissivity is very slow because paralleization is off

# First make grtrans with 'make' 
# Then run this in python
from __future__ import division
from __future__ import print_function
import numpy as np
import argparse
import grtrans_batch as gr
import matplotlib.pyplot as plt
import scipy.ndimage.filters as filt
import sys, os, time
import astropy.io.fits as fits
import scipy.ndimage.interpolation as interpolation
from scipy.interpolate import interp1d

# Run parameters
FIND_BSCL = True
RUN_IMAGE = True      # run image
RUN_SPECTRUM = True   # run spectrum 
RERUN = True          # rerun 
SAVEOUT = True        # save output images
DISPLAYOUT = True     # display output image(s)

# RRJET parameters -- these can be changed in function args below
BETAECONST = 1.e-2      # constant bete0
FPOSITRON = 1           # 0 < npositron/nelectron < 1

PEGASRATIO = -1 #0.1    # ratio of electron to gas pressure (-1 to not use this model)
                        # THIS WILL OVERWRITE THE OTHER MODELS IF NOT EQUAL TO -1 
BETAECRIT = -1          # critical beta for exponential supression (-1 uses constant betae)
XIMAX = 10.             # maximum xi=s^2/z defining jet edge
GAMMAMIN = 10           # minimum gamma for power law distribution
GAMMAMAX = 5.e3         # maximum gamma for power law distribution
PNTH =  3.1             # nonthermal power law index
BSCL = 860.             # sets field at horizon
BSCLMIN= 1.             # bscl for search
BSCLMAX=100000. 

# Source parameters
SOURCE = 'M87'          # source for fits header
RA = 12.51373           # ra for fits header
DEC = 12.39112          # dec for fits header
MJD =  58211            # mjd for fits header

# Blackhole parameters
MBH = 6.5e9    # bh mass / Msun
DTOBH = 16.8e3 # bh distance / kpc
A = 0.5        # bh spin
ANG = 160.      # polar angle (degrees)
ROTANG = 117   # rotation angle in sky plane (degrees)

# Raytrace parameters - image
# TODO -- multiple frequencies? 
FLUX = 1.5          # desired flux in Jy
RFGHZ = 230.        # Frequency in Ghz
FOV = 120.          # FOV / Rg
NPIX = 128          # number of pixels
NGEO = 1000         # number of geodesic points

# Raytrace parameters - spectrum
NFREQ = 20         # number of frequencies
FOV_SPEC= FOV      # FOV / Rg
NPIX_SPEC = NPIX   # number of pixels in spectrum image
FMIN = 1.e9        # minimum freq in spectrum
FMAX = 1.e16       # maximum freq in spectrum 

DEPTH = 5 # raytracing outer volume is DEPTH*FOV/2 in Rg
          # might want to be large for nearly face on jet, but slower

# Output File Location
OUTDIR = '../rrjet_and_riaf'                            # output directory

# Constants
mu = np.cos(ANG*np.pi/180.)
cmperKpc = 3.086e21
degree = np.pi/180.
radperuas = degree/3600./1.e6
cmperrg = 147708.885 * MBH
bhdist = DTOBH * cmperKpc
psize_rg = FOV / float(NPIX)
psize_cm = psize_rg * cmperrg
psize_rad = psize_cm / bhdist
psize_uas = psize_rad / radperuas
LumFac = (4 * np.pi * cmperrg**2)
LumtoJy = 1.e23/(4*np.pi*bhdist**2)

# Plotting parameters
M87DIR =   './M87_data' # directory with Sgr A* data files
NEWDATA = False         # Use new m87 SED from 2017 (LOWER TOTAL FLUX)
cfun = 'afmhot'
cfun2 =  'Spectral'
xticks_maj = [1.e8,1.e10,1.e12,1.e14,1.e16]#1.e16,1.e18,1.e20,1.e22]
xticks_min = [1.e9,1.e11,1.e13,1.e15]#1.e17,1.e19,1.e21]
yticks_min = [5.e39,5.e40,5.e41]
yticks_maj = [1.e39,1.e40,1.e41,1.e42]
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def run_grtrans_image(fname, fpositron=FPOSITRON,pegasratio=PEGASRATIO,
                      betaeconst=BETAECONST,betaecrit=BETAECRIT,
                      bscl=BSCL,gmin=GAMMAMIN,gmax=GAMMAMAX,pnth=PNTH):
    """ run grtrans single image"""

    print("Run image, bscl=%.2f"%bscl)

    size  = 0.5*FOV         
    uout = 1./(10.*size)

    # pressure scale is fixed!
    pscl = (bscl**2)/(8*np.pi) 

    epcoefindx=[1,1,1,1,1,1,1]

    # TO TURN OFF FARADAY CONVERSTION
    #epcoefindx=[1,1,1,1,0,1,1]

    x=gr.grtrans()
    x.write_grtrans_inputs(fname+'_im.in', oname=fname+'_im.out',
                           fname='RRJET', phi0=0., pegasratio=pegasratio,
                           betaeconst=betaeconst, betaecrit=betaecrit, 
                           ximax=XIMAX, bscl=bscl, pscl=pscl,
                           nfreq=1,fmin=RFGHZ*1.e9,fmax=RFGHZ*1.e9,
                           gmin=gmin, gmax=gmax, p2=pnth, p1=pnth,
                           fpositron=fpositron,
                           epcoefindx=epcoefindx,
                           ename='POLSYNCHPL',
                           nvals=4,
                           spin=A, standard=1,
                           uout=uout,
                           mbh=MBH,
                           nmu=1,mumin=mu,mumax=mu,
                           gridvals=[-size,size,-size,size],
                           nn=[NPIX,NPIX,NGEO],
                           hindf=1,hnt=1,
                           muval=1.)

    # run grtrans
    if RERUN:
        x.run_grtrans()

    # load image data
    x.read_grtrans_output()

    # pixel sizes
    #da = x.ab[x.nx,0]-x.ab[0,0]
    da = x.ab[x.ny,0]-x.ab[0,0] ##TODO -- is this right ordering? 
    db = x.ab[1,1]-x.ab[0,1]
    if (da!=db): raise Exception("pixel da!=db")
    psize = da*(cmperrg/bhdist)

    #image values
    ivals = x.ivals[:,0,0]*LumFac*da*db*LumtoJy
    qvals = x.ivals[:,1,0]*LumFac*da*db*LumtoJy
    uvals = x.ivals[:,2,0]*LumFac*da*db*LumtoJy
    vvals = x.ivals[:,3,0]*LumFac*da*db*LumtoJy
    
    # mask pixels with  I < 0 
    imask = ivals < 0.
    ivals[imask] = 0.
    qvals[imask] = 0.
    uvals[imask] = 0.
    vvals[imask] = 0.

    # correct orientation for eht-imaging
    ivals =  (np.flipud(np.transpose(ivals.reshape((NPIX,NPIX))))).flatten()
    qvals =  -(np.flipud(np.transpose(qvals.reshape((NPIX,NPIX))))).flatten()
    uvals =  -(np.flipud(np.transpose(uvals.reshape((NPIX,NPIX))))).flatten()
    vvals =  (np.flipud(np.transpose(vvals.reshape((NPIX,NPIX))))).flatten()
    imdata = (ivals, qvals, uvals, vvals, psize)

    # Rotate the image
    if ROTANG!=0:
        imdata = rotateimdata(imdata, ROTANG, interp='cubic')

    # save image
    if SAVEOUT:
        save_im_fits(imdata, fname + ('%.0f.fits'%RFGHZ), freq_ghz=RFGHZ)

    # display images
    if DISPLAYOUT:
        #tmax=5.e10
        #pmax=2.e10
        tmax = np.max( ivals*3.254e13/((RFGHZ*1.e9)**2 * psize**2))
        pmax = np.max( np.sqrt(qvals**2 + uvals**2)*3.254e13/((RFGHZ*1.e9)**2 * psize**2))

        display_grtrans_image(imdata, tmax=tmax, pmax=pmax)


def findbscl(fname, flux, bsclmin, bsclmax, fpositron=FPOSITRON,pegasratio=PEGASRATIO,
             betaeconst=BETAECONST,betaecrit=BETAECRIT,
             gmin=GAMMAMIN,gmax=GAMMAMAX,pnth=PNTH):
    """ run grtrans single image to find the bscl that gives the correct flux with bisection, 
        for all other parameters fixed """

    # convergance parameters
    bedge_stop = 1
    fluxconvratio = .05
    itermax = 20

    # these image parameters are fixed for now
    fov_search = FOV/2.
    npix_search = int(NPIX/2)
    
    size  = 0.5*fov_search         
    uout = 1./(10.*size)

    bsclmin0 = bsclmin
    bsclmax0 = bsclmax
    for i in range(itermax): 

        # bscl by bisection
        bscl = (bsclmax+bsclmin)/2.

        if bsclmax0-bscl < bedge_stop:
            print("did not find solution -- to close to bsclmax!")
            break
        if bscl-bsclmin0 < bedge_stop:
            print("did not find solution -- to close to bsclmin!")
            break

        # pressure scale is fixed!
        pscl = (bscl**2)/(8*np.pi) 
        x=gr.grtrans()
        x.write_grtrans_inputs(fname+'_SEARCH.in', oname=fname+'_SEARCH.out',
                               fname='RRJET', phi0=0., pegasratio=pegasratio,
                               betaeconst=betaeconst, betaecrit=betaecrit, 
                               ximax=XIMAX, bscl=bscl, pscl=pscl,
                               nfreq=1,fmin=RFGHZ*1.e9,fmax=RFGHZ*1.e9,
                               gmin=gmin, gmax=gmax, p2=pnth, p1=pnth,
                               fpositron=fpositron,
                               ename='POLSYNCHPL',
                               nvals=1,
                               spin=A, standard=1,
                               uout=uout,
                               mbh=MBH,
                               nmu=1,mumin=mu,mumax=mu,
                               gridvals=[-size,size,-size,size],
                               nn=[npix_search,npix_search,NGEO],
                               hindf=1,hnt=1,
                               muval=1.)
        # run grtrans
        x.run_grtrans()

        # load image data
        x.read_grtrans_output()

        # pixel sizes
        #da = x.ab[x.nx,0]-x.ab[0,0]
        da = x.ab[x.ny,0]-x.ab[0,0] ##TODO -- is this right ordering? 
        db = x.ab[1,1]-x.ab[0,1]
        if (da!=db): raise Exception("pixel da!=db")
        psize = da*(cmperrg/bhdist)

        #image values
        ivals = x.ivals[:,0,0]*LumFac*da*db*LumtoJy
        imask = ivals < 0.
        ivals[imask] = 0.

        #total flux
        tflux = np.sum(ivals)
        tfluxdiff = tflux-flux
        tfluxdiff_rel = np.abs(tfluxdiff/flux)

        print("iter %i %.1f | %.3f/%.3f %.2f"%(i+1,bscl,tflux,flux,tfluxdiff_rel)) 

        # flux must monotonically increase with bscl,  if all other params fixed
        if tfluxdiff_rel < fluxconvratio:
            print("solution: bscl=%.2f, tflux=%.3f" %(bscl,tflux))
            break

        if (tflux<flux):
            bsclmin=bscl

        elif (tflux>flux):
            bsclmax=bscl

        if i==itermax-1:
            print("did not find solution -- reached itermax!")
            break

    return bscl

def run_grtrans_spectrum(fname, fpositron=FPOSITRON,pegasratio=PEGASRATIO,
                         betaeconst=BETAECONST,betaecrit=BETAECRIT,
                         bscl=BSCL,gmin=GAMMAMIN,gmax=GAMMAMAX,pnth=PNTH):
    """Run grtrans spectrum"""

    size_spec  = 0.5*FOV_SPEC      
    uout_spec = 1./(DEPTH*size_spec)

    npix_x = NPIX_SPEC
    npix_y = NPIX_SPEC

    size_x = size_spec
    size_y = size_spec

    # pressure scale is fixed!
    pscl = (bscl**2)/(8*np.pi) 

    x=gr.grtrans()
    x.write_grtrans_inputs(fname+'_spec.in', oname=fname+'_spec.out',
                           fname='RRJET', phi0=0., pegasratio=pegasratio,
                           betaeconst=betaeconst, betaecrit=betaecrit, 
                           ximax=XIMAX, bscl=bscl, pscl=pscl,
                           nfreq=NFREQ,fmin=FMIN,fmax=FMAX,
                           gmin=gmin, gmax=gmax, p2=pnth, p1=pnth,
                           fpositron=fpositron,
                           ename='POLSYNCHPL',
                           nvals=4,
                           spin=A, standard=1,
                           uout=uout_spec,
                           mbh=MBH,
                           nmu=1,mumin=mu,mumax=mu,
                           gridvals=[-size_x,size_x,-size_y,size_y],
                           nn=[npix_x,npix_y,NGEO],
                           hindf=1,hnt=1,
                           muval=1.)

    # run grtrans
    if RERUN:
        x.run_grtrans()

    # load spectrum
    x.read_grtrans_output()
    x.calc_freqs(NFREQ)
    x.convert_to_lum()
    spec = x.spec[0][0:NFREQ]
    qspec = x.spec[1][0:NFREQ]
    uspec = x.spec[2][0:NFREQ]
    vspec = x.spec[3][0:NFREQ]
    if npix_x==1 or npix_y==1:
        spec *= 0.5  # divide by 2 because we have +/- r in the strip
        qspec *= 0.5
        uspec *= 0.5
        vspec *= 0.5
    freqs = x.freqs

    # save spectrum
    if SAVEOUT:
        outdat = np.vstack([freqs,spec,qspec,uspec,vspec]).T
        np.savetxt(fname + '_spec.txt', outdat)

    # display spectrum
    if DISPLAYOUT:
        plot_grtrans_spectrum(freqs, spec, qspec, uspec, vspec)


def plot_grtrans_spectrum(freqs, spec, qspec, uspec, vspec):

    # plot Stokes I spectrum -- nu*Lnu
    f=plt.figure(111,figsize=(16,16))
    plt.clf()
    ax=f.add_subplot(111)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$\\nu$ (Hz)', size=26)
    plt.ylabel('$\\nu L_{\\nu}$ (erg s$^{-1}$)', size=26)
    plt.xlim([1.e9,1.e16])
    plt.ylim([1.e39,1.e42])
    plt.tick_params(axis='both',labelsize=22)
    plt.ion()
    ax = plot_m87_data(ax)

    #linestyles=['solid','dashdot','dashed']
    ls = 'solid'

    spec_interp = interp1d(np.log10(freqs), np.log10(spec), kind=2)
    logfreqs_plot = np.linspace(np.log10(FMIN), np.log10(FMAX), 500)
    logspec_plot = spec_interp(logfreqs_plot)

    plt.plot(10**logfreqs_plot, 10**(logfreqs_plot+logspec_plot), 'k-',
             linewidth=2, label=r'I, $f_p=%.1f$'%FPOSITRON, linestyle=ls)

    plt.legend()

    #########################
    # plot Stokes I spectrum -- Fnu
    f=plt.figure(222,figsize=(16,16))
    plt.clf()
    ax=f.add_subplot(111)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$\\nu$ (Hz)', size=26)
    plt.ylabel('$F_{\\nu}$ (Jy)', size=26)
    plt.xlim([1.e9,1.e16])
    plt.ylim([1.e-5,10.])
    plt.tick_params(axis='both',labelsize=22)
    plt.ion()
    ax = plot_m87_data_Jy(ax)

    #linestyles=['solid','dashdot','dashed']
    ls = 'solid'

    spec *= LumtoJy
    spec_interp = interp1d(np.log10(freqs), np.log10(spec), kind=2)
    logfreqs_plot = np.linspace(np.log10(FMIN), np.log10(FMAX), 500)
    logspec_plot = spec_interp(logfreqs_plot)

    #plt.plot(freqs, freqs*spec, 'k-', linewidth=2, label=r'I, $f_p=%.1f$'%FPOSITRON, linestyle=ls)
    plt.plot(10**logfreqs_plot, 10**(logspec_plot), 'k-',
             linewidth=2, label=r'I, $f_p=%.1f$'%FPOSITRON, linestyle=ls)

    plt.legend()

def display_grtrans_image(imdata,nvec=25,veccut=0.005,tmax=1.e10,pmax=1.e10,blur_kernel=0):

    (I_im, Q_im, U_im, V_im, psize) = imdata
    I_im =  I_im.reshape((NPIX,NPIX))
    Q_im =  Q_im.reshape((NPIX,NPIX))
    U_im =  U_im.reshape((NPIX,NPIX))
    V_im =  V_im.reshape((NPIX,NPIX))

    # convert to brightness temperature Tb
    totalflux = np.sum(I_im)
    factor = 3.254e13/((RFGHZ*1.e9)**2 * psize**2)
    I_im *= factor
    Q_im *= factor
    U_im *= factor
    V_im *= factor
    
    # Polarization Vectors
    P_im = np.abs(Q_im + 1j*U_im)
    m_im = P_im/I_im

    thin = NPIX//nvec
    mask = I_im > veccut * np.max(I_im)
    mask2 = mask[::thin, ::thin]

    m = m_im[::thin, ::thin][mask2]
    x = (np.array([[i for i in range(NPIX)] for j in range(NPIX)])[::thin, ::thin])
    x = x[mask2]
    y = (np.array([[j for i in range(NPIX)] for j in range(NPIX)])[::thin, ::thin])
    y = y[mask2]
    a = (-np.sin(np.angle(Q_im+1j*U_im)/2)[::thin, ::thin])
    a = a[mask2]
    b = ( np.cos(np.angle(Q_im+1j*U_im)/2)[::thin, ::thin])
    b = b[mask2]

    P_im[np.logical_not(mask)]=0.
    m_im[np.logical_not(mask)]=0.

    # ticks
    xticks = ticks(NPIX, FOV/float(NPIX))
    yticks = ticks(NPIX, FOV/float(NPIX))

    # display Stokes I 
    plt.figure(0)
    im = plt.imshow(I_im, cmap=plt.get_cmap(cfun), interpolation='gaussian',vmin=0,vmax=tmax)
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('Tb (K)', fontsize=14)
    plt.title(("Stokes I, %.2f GHz, %.2f Jy" % (RFGHZ, totalflux)), fontsize=16)
    plt.xticks(xticks[0], xticks[1])
    plt.yticks(yticks[0], yticks[1])
    plt.xlabel('x/rg')
    plt.ylabel('y/rg')

    # display Stokes Q 
    plt.figure(1)
    im = plt.imshow(Q_im, cmap=plt.get_cmap(cfun2), interpolation='gaussian',vmin=-pmax,vmax=pmax)
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('Tb (K)', fontsize=14)
    plt.title(("Stokes Q, %.2f GHz " % (RFGHZ)), fontsize=16)
    plt.xticks(xticks[0], xticks[1])
    plt.yticks(yticks[0], yticks[1])
    plt.xlabel('x/rg')
    plt.ylabel('y/rg')

    # display Stokes U
    plt.figure(2)
    im = plt.imshow(U_im, cmap=plt.get_cmap(cfun2), interpolation='gaussian',vmin=-pmax,vmax=pmax)
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('Tb (K)', fontsize=14)
    plt.title(("Stokes U, %.2f GHz " % (RFGHZ)), fontsize=16)
    plt.xticks(xticks[0], xticks[1])
    plt.yticks(yticks[0], yticks[1])
    plt.xlabel('x/rg')
    plt.ylabel('y/rg')

    # display Stokes V 
    plt.figure(3)
    im = plt.imshow(V_im, cmap=plt.get_cmap(cfun2), interpolation='gaussian',vmin=-pmax,vmax=pmax)
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('Tb (K)', fontsize=14)
    voi = np.sum(V_im)/np.sum(I_im)
    plt.title(("Stokes V, %.2f GHz , V/I=%.2f " % (RFGHZ, voi)), fontsize=16)
    plt.xticks(xticks[0], xticks[1])
    plt.yticks(yticks[0], yticks[1])
    plt.xlabel('x/rg')
    plt.ylabel('y/rg')

    # display P
    plt.figure(4)
    im = plt.imshow(P_im, cmap=plt.get_cmap(cfun), interpolation='gaussian',vmin=0,vmax=pmax)
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('Tb (K)', fontsize=14)
    poi = np.sum(P_im)/np.sum(I_im)
    plt.title(("P, %.2f GHz , P/I=%.2f" % (RFGHZ,poi)), fontsize=16)
    plt.xticks(xticks[0], xticks[1])
    plt.yticks(yticks[0], yticks[1])
    plt.xlabel('x/rg')
    plt.ylabel('y/rg')

    # display m
    plt.figure(5)
    im = plt.imshow(m_im, cmap=plt.get_cmap('viridis'), interpolation='gaussian')
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('P/I', fontsize=14)
    plt.title(("P/I, %.2f GHz , P/I=%.2f" % (RFGHZ,poi)), fontsize=16)
    plt.xticks(xticks[0], xticks[1])
    plt.yticks(yticks[0], yticks[1])
    plt.xlabel('x/rg')
    plt.ylabel('y/rg')

    # display I with pol ticks
    plt.figure(6)
    im = plt.imshow(I_im, cmap=plt.get_cmap(cfun), interpolation='gaussian',vmin=0,vmax=tmax)
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('Tb (K)', fontsize=14)
    plt.title(("I, %.2f GHz " % (RFGHZ)), fontsize=16)
    plt.xticks(xticks[0], xticks[1])
    plt.yticks(yticks[0], yticks[1])
    plt.xlabel('x/rg')
    plt.ylabel('y/rg')

    plt.quiver(x, y, a, b,
           headaxislength=20, headwidth=1, headlength=.01, minlength=0, minshaft=1,
           width=.01*NPIX, units='x', pivot='mid', color='k', angles='uv', 
           scale=1.0/thin)
    plt.quiver(x, y, a, b,
           headaxislength=20, headwidth=1, headlength=.01, minlength=0, minshaft=1,
           width=.005*NPIX, units='x', pivot='mid', color='w', angles='uv', 
           scale=1.1/thin)

def ticks(axisdim, psize, nticks=8):
    """Return a list of ticklocs and ticklabels
       psize should be in desired units
    """

    axisdim = int(axisdim)
    nticks = int(nticks)
    if not axisdim % 2: axisdim += 1
    if nticks % 2: nticks -= 1
    tickspacing = float((axisdim-1))/nticks
    ticklocs = np.arange(0, axisdim+1, tickspacing) - 0.5
    ticklabels= np.around(psize * np.arange((axisdim-1)/2.0, -(axisdim)/2.0, -tickspacing), decimals=1)

    return (ticklocs, ticklabels)

def rotateimdata(imdata, angle, interp='cubic'):

    (imvec, qvec, uvec, vvec, psize) = imdata

    order=3
    if interp=='linear': order=1
    elif interp=='cubic': order=3
    elif interp=='quintic': order=5
    
    # Define an interpolation function
    def rot_imvec(imvec):
        imarr_rot = interpolation.rotate(imvec.reshape((NPIX, NPIX)),
                                               angle, reshape=False, order=order,
                                               mode='constant', 
                                               cval=0.0, prefilter=True)
        return imarr_rot.flatten()

    # Make new image vectors
    imvec_rot = rot_imvec(imvec)
    vvec_rot  = rot_imvec(vvec)
    qvec_rot  = rot_imvec(qvec)
    uvec_rot  = rot_imvec(uvec)

    # correctly transform Q,U
    rlvec_rot = np.exp(1j*2*angle*degree) * (qvec_rot + 1j*uvec_rot)
    qvec_rot = np.real(rlvec_rot)
    uvec_rot = np.imag(rlvec_rot)

    imdata_out = (imvec_rot, qvec_rot, uvec_rot, vvec_rot, psize)

    return imdata_out

def save_im_fits(imdata, fname, freq_ghz=RFGHZ, 
                 mjd=MJD, source=SOURCE, ra=RA, dec=DEC):

    # unpack the image data
    (imvec, qvec, uvec, vvec, psize) = imdata

    # Create header and fill in some values
    header = fits.Header()
    header['OBJECT'] = source
    header['CTYPE1'] = 'RA---SIN'
    header['CTYPE2'] = 'DEC--SIN'
    header['CDELT1'] = -psize/degree
    header['CDELT2'] =  psize/degree
    header['OBSRA'] = ra * 180/12.
    header['OBSDEC'] = dec
    header['FREQ'] = freq_ghz * 1.e9

    #TODO these are the default values for centered images
    #TODO support for arbitrary CRPIX? 
    header['CRPIX1'] = NPIX/2. + .5
    header['CRPIX2'] = NPIX/2. + .5

    header['MJD'] = float(mjd)
    header['TELESCOP'] = 'VLBI'
    header['BUNIT'] = 'JY/PIXEL'
    header['STOKES'] = 'I'

    # Create the fits image
    image = np.reshape(imvec,(NPIX,NPIX))[::-1,:] #flip y axis!
    hdu = fits.PrimaryHDU(image, header=header)
    hdulist = [hdu]
    if len(qvec):
        qimage = np.reshape(qvec,(NPIX,NPIX))[::-1,:]
        uimage = np.reshape(uvec,(NPIX,NPIX))[::-1,:]
        header['STOKES'] = 'Q'
        hduq = fits.ImageHDU(qimage, name='Q', header=header)
        header['STOKES'] = 'U'
        hduu = fits.ImageHDU(uimage, name='U', header=header)
        hdulist = [hdu, hduq, hduu]
    if len(vvec):
        vimage = np.reshape(vvec,(NPIX,NPIX))[::-1,:]
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

def plot_m87_data(ax):

    CAPSIZE = .8
    if NEWDATA:
        # Data from 2017 eht campaign
        data1 =  np.loadtxt(M87DIR + '/m87_sed_2017.dat')
        data=data1
        freqsdat = data[:,1]
        lumval = freqsdat*data[:,3] / LumtoJy
        lumerr = freqsdat*data[:,4] / LumtoJy
        (_,caps,_) = plt.errorbar(freqsdat,lumval,yerr=lumerr,fmt='o',
                                  color='k',ecolor='k',markersize=4*CAPSIZE,
                                  markerfacecolor=None,capthick=CAPSIZE,capsize=CAPSIZE*4)

        for cap in caps:
            cap.set_color('black')
    else:

        # data table 1 -- core flux in quiescent --- lower resolution
        data1 =  np.loadtxt(M87DIR + '/m87_data_1.txt')
        data=data1
        freqsdat = data[:,0]
        lumval = freqsdat*data[:,1] / LumtoJy
        lumerr = freqsdat*data[:,2] / LumtoJy
        (_,caps,_) = plt.errorbar(freqsdat,lumval,yerr=lumerr,fmt='o',color='k',ecolor='k',
                                  markersize=4*CAPSIZE,
                                  markerfacecolor=None,capthick=CAPSIZE,capsize=CAPSIZE*4)
        for cap in caps:
            cap.set_color('black')

        # Michael total flux
        datam =  np.loadtxt(M87DIR + '/m87_data_michael.txt')
        data=datam
        freqsdat = data[:,0]
        lumval = freqsdat*data[:,1] / LumtoJy
        lumerr = freqsdat*data[:,2] / LumtoJy
        (_,caps,_) = plt.errorbar(freqsdat,lumval,yerr=lumerr,fmt="v",color='c',ecolor='c',
                                  markersize=6*CAPSIZE,capthick=2*CAPSIZE,capsize=CAPSIZE*4,zorder=10)
        for cap in caps:
            cap.set_color('c')
            #cap.set_markeredgewidth(2)

#        # Michael core flux
#        datam =  np.loadtxt(M87DIR + '/m87_data_michael_core.txt')
#        data=datam
#        freqsdat = data[:,0]
#        lumval = freqsdat*data[:,1] / LumtoJy
#        lumerr = freqsdat*data[:,2] / LumtoJy
#        (_,caps,_) = plt.errorbar(freqsdat,lumval,yerr=lumerr,fmt='^',color='m',ecolor='m',
#                                  markersize=6*CAPSIZE,capthick=2*CAPSIZE,
#                                  capsize=CAPSIZE*4, zorder=10)
#        for cap in caps:
#            cap.set_color('m')


    return ax

def plot_m87_data_Jy(ax):

    CAPSIZE = .8

    if NEWDATA:
        # data from 2017 eht campaign
        data1 =  np.loadtxt(M87DIR + '/m87_sed_2017.dat')
        data=data1
        freqsdat = data[:,1]
        lumval = data[:,3] 
        lumerr = data[:,4] 
        (_,caps,_) = plt.errorbar(freqsdat,lumval,yerr=lumerr,fmt='o',
                                  color='k',ecolor='k',markersize=4*CAPSIZE,
                                  markerfacecolor=None,capthick=CAPSIZE,capsize=CAPSIZE*4)
        for cap in caps:
            cap.set_color('black')
    else:

        # data table 1 -- core flux in quiescent --- lower resolution
        data1 =  np.loadtxt(M87DIR + '/m87_data_1.txt')
        data=data1
        freqsdat = data[:,0]
        lumval = data[:,1] 
        lumerr = data[:,2] 
        (_,caps,_) = plt.errorbar(freqsdat,lumval,yerr=lumerr,fmt='o',
                                   color='k',ecolor='k',markersize=4*CAPSIZE,
                                  markerfacecolor=None,capthick=CAPSIZE,capsize=CAPSIZE*4)
        for cap in caps:
            cap.set_color('black')


        # Michael total flux
        datam =  np.loadtxt(M87DIR + '/m87_data_michael.txt')
        data=datam
        freqsdat = data[:,0]
        lumval = data[:,1] 
        lumerr = data[:,2] 
        (_,caps,_) = plt.errorbar(freqsdat,lumval,yerr=lumerr,fmt="v",color='c',ecolor='c',
                                  markersize=6*CAPSIZE,capthick=2*CAPSIZE,capsize=CAPSIZE*4,zorder=10)
        for cap in caps:
            cap.set_color('c')
            #cap.set_markeredgewidth(2)

    return ax

def ticks(axisdim, psize, nticks=8):
    """Return a list of ticklocs and ticklabels
       psize should be in desired units
    """

    axisdim = int(axisdim)
    nticks = int(nticks)
    if not axisdim % 2: axisdim += 1
    if nticks % 2: nticks -= 1
    tickspacing = float((axisdim-1))/nticks
    ticklocs = np.arange(0, axisdim+1, tickspacing) - 0.5
    ticklabels= np.around(psize * np.arange((axisdim-1)/2.0, -(axisdim)/2.0, -tickspacing), decimals=1)

    return (ticklocs, ticklabels)


if __name__=='__main__':

    # parse parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('--betaeconst', type=float,default=BETAECONST)
    parser.add_argument('--fpositron',type=float,default=FPOSITRON)
    args = parser.parse_args()
    fpos = args.fpositron
    beconst = args.betaeconst

    print('betaeconst: %.2e fpositron: %.2f' % (beconst, fpos))

    # preliminaries for display/save output
    if DISPLAYOUT:
        plt.close('all')

    if SAVEOUT: 
        modeltag = 'betae%0.1e_fpos%.1f' % (beconst, fpos)
        fdir = OUTDIR + '/' + modeltag
        if not os.path.exists(fdir):
            os.mkdir(fdir)
        fname = fdir + '/' + modeltag
    else: 
        fname = OUTDIR + '/rrjet_tmp'

    # find the bscale factor
    if FIND_BSCL and RERUN:
        bscl = findbscl(fname, FLUX, BSCLMIN, BSCLMAX,
                        fpositron=fpos, betaeconst=beconst,
                        pegasratio=PEGASRATIO, betaecrit=BETAECRIT,
                        gmin=GAMMAMIN,gmax=GAMMAMAX,pnth=PNTH)
    else:
        bscl=BSCL

    if SAVEOUT: 
        np.savetxt(fname + '_bscl.txt', np.array([bscl]))
   
    # run the image
    # TODO -- multiple frequencies? 
    if RUN_IMAGE:
        run_grtrans_image(fname, fpositron=fpos, betaeconst=beconst,
                          pegasratio=PEGASRATIO, betaecrit=BETAECRIT,
                          bscl=bscl,gmin=GAMMAMIN,gmax=GAMMAMAX,pnth=PNTH)

    # run the spectrum    
    if RUN_SPECTRUM:
        run_grtrans_spectrum(fname, fpositron=fpos, betaeconst=beconst,
                             pegasratio=PEGASRATIO, betaecrit=BETAECRIT,
                             bscl=bscl,gmin=GAMMAMIN,gmax=GAMMAMAX,pnth=PNTH)

    print('BSCL', bscl)
    if DISPLAYOUT:
        plt.show()
