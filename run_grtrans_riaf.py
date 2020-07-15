# Run Grtrans with raif model
# NOTE -- currently the power law emissivity is slow because paralleization is off

# First make grtrans with 'make' 
# Then run this in python

import numpy as np
import grtrans_batch as gr
import matplotlib.pyplot as plt
import scipy.ndimage.filters as filt
import sys, os, time
import astropy.io.fits as fits
import scipy.ndimage.interpolation as interpolation
from scipy.interpolate import interp1d

# Run parameters
RUN_IMAGE = False    # run image
RUN_SPECTRUM = True # run spectrum 
RERUN = True      # rerun 
SAVEOUT = True    # save output images
DISPLAYOUT = True # display output image(s)

# Broderick&Loeb 06 parameters
NSCL = 1.0e7             # scaling factor for thermal electron number density
NNTHSCL = 1.e5          # scaling factor for nonthermal electron number density
TSCL = 1.5e11           # scaling factor for electron temperature
BETA = 10.              # plasma beta
BLO6 = 0                # Use Broderick 06 conventions (1) or Broderick 09 conventions (0)
NTH_RADIAL_PLAW = 2.02  # radial power law for nonthermal electrons (fixed to -2.02 in BL11)
GAMMAMIN = 100          # minimum gamma for power law distribution
PNTH = 3.5              # nonthermal power law index
FPOSITRON=0            # 0 < npositron/nelectron < 1

# Blackhole parameters
MBH = 4.1e6  # bh mass / Msun
DTOBH = 8.1  # bh distance / kpc
A = 0.1      # bh spin
ANG = 60.    # polar angle (degrees)
ROTANG = 156 # rotation angle in sky plane (degrees)

# Raytrace parameters - image
RFGHZ = 230.        # Frequency in Ghz
FOV = 30.           # FOV / Rg
NPIX = 128          # number of pixels
NGEO = 500          # number of geodesic points
size  = 0.5*FOV         
uout = 1./(2.*size)

# Raytrace parameters - spectrum
NFREQ = 20          # number of frequencies
FOV_SPEC = 250      # FOV / Rg
NPIX_SPEC = 64      # number of pixels in spectrum image
FMIN = 1.e9         # minimum freq in spectrum
FMAX = 1.e14        # maximum freq in spectrum 
size_spec  = 0.5*FOV_SPEC      
uout_spec = 1./(2.*size_spec)

# Output File names
OUTDIR = '../rrjet_and_riaf'                            # output directory
OUTNAME =  ('riaf_%0.1f_'%FPOSITRON) + str(ANG) + '_im' # output file name - fits
oname = 'riaf'                                          # output file name - grtrans internal
fname = OUTDIR + '/' + OUTNAME
SOURCE = 'SGRA'                                         # source for fits header
RA = 17.761122                                          # ra for fits header
DEC = -20.007                                           # dec for fits header
MJD =  58211                                            # mjd for fits header

# Constants
mu = np.cos(ANG*np.pi/180.)
cmperKpc = 3.086e21
degree = np.pi/180.
radperuas = degree/3600./1.e6
psize_rg = 2*size / NPIX
cmperrg = 147708.885 * MBH
bhdist = DTOBH * cmperKpc
psize_cm = psize_rg * cmperrg
psize_rad = psize_cm / bhdist
psize_uas = psize_rad / radperuas
LumFac = (4 * np.pi * cmperrg**2)
LumtoJy = 1.e23/(4*np.pi*bhdist**2)

# Plotting parameters
SGRADIR =   './SgrA_data' # directory with Sgr A* data files
sgra_radio_data =  np.loadtxt(SGRADIR + '/sgr_A_data_radio.txt') # Sgr A* radio data
sgra_ir_data = np.loadtxt(SGRADIR + '/sgr_A_data_ir.txt') # Sgr A* infrared data
cfun = 'afmhot'
cfun2 =  'Spectral'
xticks_maj = [1.e8,1.e10,1.e12,1.e14]#1.e16,1.e18,1.e20,1.e22]
xticks_min = [1.e9,1.e11,1.e13,1.e15]#1.e17,1.e19,1.e21]
yticks_min = [1.e31,1.e33,1.e35,1.e37]
yticks_maj = [1.e30,1.e32,1.e34,1.e36]
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def run_grtrans_image():
    """ run grtrans single image"""
    x=gr.grtrans()

    x.write_grtrans_inputs(oname+'_im.in', oname=oname+'_im.out',
                           fname='SARIAF', phi0=0.,
                           nfreq=1,fmin=RFGHZ*1.e9,fmax=RFGHZ*1.e9,
                           gmin=GAMMAMIN, p2=PNTH, p1=PNTH,
                           snscl=NSCL, ntscl=TSCL, snnthscl=NNTHSCL, 
                           snnthp=NTH_RADIAL_PLAW, sbeta=BETA, sbl06=BLO6,
                           fpositron=FPOSITRON,
                           ename='HYBRIDTHPL',
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
    da = x.ab[x.nx,0]-x.ab[0,0]
    db = x.ab[1,1]-x.ab[0,1]
    if (da!=db): raise Exception("pixel da!=db")
    psize = da*(cmperrg/bhdist)

    #image values
    ivals = x.ivals[:,0,0]*LumFac*da*db*LumtoJy
    qvals = x.ivals[:,1,0]*LumFac*da*db*LumtoJy
    uvals = x.ivals[:,2,0]*LumFac*da*db*LumtoJy
    vvals = x.ivals[:,3,0]*LumFac*da*db*LumtoJy
    
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
        save_im_fits(imdata, fname + '.fits', freq_ghz=RFGHZ)

    # display images
    if DISPLAYOUT:
        #tmax=5.e10
        #pmax=2.e10
        tmax = np.max( ivals*3.254e13/((RFGHZ*1.e9)**2 * psize**2))
        pmax = np.max( np.sqrt(qvals**2 + uvals**2)*3.254e13/((RFGHZ*1.e9)**2 * psize**2))

        display_grtrans_image(imdata, tmax=tmax, pmax=pmax)

def run_grtrans_spectrum():
    """Run grtrans spectrum"""

    #sizey =  size_spec / NPIX_SPEC

    x=gr.grtrans()
    x.write_grtrans_inputs(oname+'_spec.in', oname=oname+'_spec.out',
                           fname='SARIAF', phi0=0.,
                           nfreq=NFREQ,fmin=1.e8,fmax=1.e16,
                           gmin=GAMMAMIN, p2=PNTH, p1=PNTH,
                           snscl=NSCL, ntscl=TSCL, snnthscl=NNTHSCL, 
                           snnthp=NTH_RADIAL_PLAW, sbeta=BETA, sbl06=BLO6,
                           fpositron=FPOSITRON,
                           ename='HYBRIDTHPL',
                           nvals=1,
                           spin=A, standard=1,
                           uout=uout_spec,
                           mbh=MBH,
                           nmu=1,mumin=mu,mumax=mu,
                           #gridvals=[-size_spec,size_spec,-sizey,sizey],
                           #nn=[NPIX_SPEC,1,NGEO],
                           gridvals=[-size_spec,size_spec,-size_spec,size_spec],
                           nn=[NPIX_SPEC,NPIX_SPEC,NGEO],
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
    freqs = x.freqs

    # plot Stokes I spectrum
    f=plt.figure(111,figsize=(16,16))
    plt.clf()
    ax=f.add_subplot(111)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$\\nu$ (Hz)', size=26)
    plt.ylabel('$\\nu L_{\\nu}$ (erg s$^{-1}$)', size=26)
    plt.xlim([1.e9,1.e15])
    plt.ylim([1.e30,1.e36])
    plt.tick_params(axis='both',labelsize=22)
    plt.ion()
    ax = plot_sgra_data(ax)

    #linestyles=['solid','dashdot','dashed']
    ls = 'solid'

    spec_interp = interp1d(np.log10(freqs), np.log10(spec), kind=3)
    logfreqs_plot = np.linspace(np.log10(FMIN), np.log10(FMAX), 500)
    logspec_plot = spec_interp(logfreqs_plot)

    #plt.plot(freqs, freqs*spec, 'k-', linewidth=2, label=r'I, $f_p=%.1f$'%FPOSITRON, linestyle=ls)
    plt.plot(10**logfreqs_plot, 10**(logfreqs_plot+logspec_plot), 'k-',
             linewidth=2, label=r'I, $f_p=%.1f$'%FPOSITRON, linestyle=ls)

    plt.legend()

    # Plot linear polarization  and Stokes V spectra
#        specQ = x.spec[1][0:NFREQ]
#        specU = x.spec[2][0:NFREQ]
#        specP = np.sqrt(specQ**2 + specU**2)
#        specV = np.abs(x.spec[3][0:NFREQ])
#        plt.plot(freqs,freqs*specP,'b-',linewidth=2,label=r'P, $f_p=%.1f$'%FPOSITRON, linestyle=ls)
#        plt.plot(freqs,freqs*specV,'r-',linewidth=2,label=r'V, $f_p=%.1f$'%FPOSITRON, linestyle=ls)

#        f=plt.figure(2,figsize=(16,16))
#        plt.rc('text', usetex=True)
#        plt.rc('font', family='serif')
#        #plt.title('fpositron=%.1f'%FPOSITRON)
#        plt.xscale('log')
#        plt.yscale('log')
#        plt.xlim([1.e9,1.e15])
#        plt.ylim([1.e-3,1])

#        plt.figure(2)
#        plt.plot(freqs,specP/spec,'b-',linewidth=2,label=r'P/I, $f_p=%.1f$'%FPOSITRON, linestyle=ls)
#        plt.plot(freqs,specV/spec,'r-',linewidth=2,label=r'V/I, $f_p=%.1f$'%FPOSITRON, linestyle=ls)
#        plt.legend()

def display_grtrans_image(imdata,nvec=25,veccut=0.005,tmax=1.e10,pmax=1.e10,blur_kernel=0):

    (I_im, Q_im, U_im, V_im, psize) = imdata
    I_im =  I_im.reshape((NPIX,NPIX))
    Q_im =  Q_im.reshape((NPIX,NPIX))
    U_im =  U_im.reshape((NPIX,NPIX))
    V_im =  V_im.reshape((NPIX,NPIX))

    # convert to brightness temperature Tb
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
    xticks = ticks(NPIX, 2*size/NPIX)
    yticks = ticks(NPIX, 2*size/NPIX)

    # display Stokes I 
    plt.figure(0)
    im = plt.imshow(I_im, cmap=plt.get_cmap(cfun), interpolation='gaussian',vmin=0,vmax=tmax)
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('Tb (K)', fontsize=14)
    plt.title(("Stokes I, %.2f GHz " % (RFGHZ)), fontsize=16)
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

def plot_sgra_data(ax):
    CAPSIZE = .5

    # radio data
    data=sgra_radio_data
    freqsdat = data[:,0] * 1.e9
    lumval = freqsdat*data[:,1] / LumtoJy
    lumerr = freqsdat*data[:,2] / LumtoJy
    (_,caps,_) = plt.errorbar(freqsdat,lumval,yerr=lumerr,fmt='o',color='k',ecolor='k',
                              markersize=12*CAPSIZE,capthick=CAPSIZE,capsize=CAPSIZE*15)
    for cap in caps:
        cap.set_color('black')
        #cap.set_markeredgewidth(2)

    # infrared data
    data = sgra_ir_data
    freqsdat = 3.e8/ (data[:,0] * 1.e-6)
    lumval = 1.e-3*freqsdat*data[:,1] / LumtoJy
    lumerr = 1.e-3*freqsdat*data[:,2] / LumtoJy

    #not upper limits
    mask = data[:,3] == 0
    mask2 = data[mask][:,4]==1
    (_,caps,_) = plt.errorbar(freqsdat[mask][mask2],lumval[mask][mask2],yerr=lumerr[mask][mask2],fmt='o',color='r',
                              ecolor='r',markersize=12*CAPSIZE,capthick=CAPSIZE,capsize=CAPSIZE*15)
    for cap in caps:
        cap.set_color('red')
        #cap.set_markeredgewidth(2)
    (_,caps,_) = plt.errorbar(freqsdat[mask][~mask2],lumval[mask][~mask2],yerr=lumerr[mask][~mask2],
                              fmt='o',color='b',ecolor='b',markersize=12*CAPSIZE,capthick=CAPSIZE,capsize=CAPSIZE*15)
    for cap in caps:
        cap.set_color('blue')
        #cap.set_markeredgewidth(2)


    #upper limits
    mask = data[:,3] == 1
    (_,caps,_) = plt.errorbar(freqsdat[mask],lumval[mask],yerr=lumerr[mask],uplims=lumval[mask],
                              fmt='.',color='g',ecolor='g',capthick=CAPSIZE,capsize=CAPSIZE*15,markersize=12*CAPSIZE)
    for cap in caps:
        cap.set_color('green')
        #cap.set_markeredgewidth(10)


#    #NIR spectral slope (Witzel 2013)
#    nulnu_slope = 0.4
#    nirmin_freq = 1.e13
#    nirmax_freq = 1.e15
#    norm = 4.e34
#    plt.plot((nirmin_freq,nirmax_freq),(norm, norm*(nirmax_freq**(nulnu_slope))/(nirmin_freq**(nulnu_slope))),
#              'r--',linewidth=2)

#    # xray - 2-10 keV
#    #flare
#    #Neilsen 03
#    freqs = np.array([2,10])/4.14e-18
#    upper = np.array((2.0e35,2.0e35))
#    lower = np.array((1.0e34,1.0e34))
#    ax.add_patch(patches.Rectangle((2./4.14e-18,1.0e34),height=(2.0e35-1.0e34),width=8./4.14e-18,
#                 facecolor='k',alpha=0.25, linestyle='--',linewidth=1.3))
#    plt.fill_between(freqs, lower, upper, alpha=0.1, edgecolor='k', facecolor=None,)

#    #quiescent
#    #Baganoff 03
#    #approx 10% believed to come from inner region
#    upper = np.array((2.4e33,2.4e33))
#    lower = np.array((2.4e32,2.4e32))
#    ax.add_patch(patches.Rectangle((2./4.14e-18,2.4e32),height=(2.4e33-2.4e32),width=8./4.14e-18,
#                 facecolor='k',alpha=0.25, linestyle=':',linewidth=1.3))
#    plt.fill_between(freqs, lower, upper, alpha=0.1, edgecolor='k', facecolor=None, linestyle=':',linewidth=1.5)

    # ticks and return
    ax.set_xticks(xticks_min)
    ax.set_xticks(xticks_maj, minor=True)
    ax.set_xticklabels([], minor=True)

    ax.set_yticks(yticks_maj)
    ax.set_yticks(yticks_min, minor=True)
    ax.set_yticklabels([], minor=True)

    plt.tick_params(axis='both',which='minor',length=5)
    plt.tick_params(axis='both',which='major',length=8)
    return ax


if __name__=='__main__':
    plt.close('all')
    if RUN_IMAGE:
        run_grtrans_image()
    if RUN_SPECTRUM:
        run_grtrans_spectrum()
    if DISPLAYOUT:
        plt.show()
