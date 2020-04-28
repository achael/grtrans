# Run Grtrans with rrjet model
# The rrjet model is defined in "fluid_model_rrjet.py"
# NOTE -- currently the power law emissivity is very slow because paralleization is off

# First make grtrans with 'make' 
# Then run this in python

import numpy as np
import grtrans_batch as gr
import matplotlib.pyplot as plt
import scipy.ndimage.filters as filt

ang=20.
name = 'rrjet'+str(ang)
mu = np.cos(ang*np.pi/180.)
size  = 300.
uout = 1./(10*size)
npix = 100
ngeo = 1000

cmperMpc = 3.086e24
MBH = 6.7e9
DTOBH = 16.528*cmperMpc
RADPERUAS = np.pi/180./3600./1.e6

psize_rg = 2*size/npix
cmperrg = 147708.8 * MBH
psize_cm = psize_rg * cmperrg
psize_rad = psize_cm / DTOBH
psize_uas = psize_rad / RADPERUAS

pp= 2.001
RF = 230.e9
BETAECONST= -1 #1.e-10
FPOSITRON = 0
cfun = 'jet'
cfun2 = 'seismic'
RERUN = True
FNAME = 'grtrans_jet_compare_positrons.txt'

def main():

    # run grtrans
    x=gr.grtrans()

    x.write_grtrans_inputs(name+'.in', oname=name+'.out',
                           fname='RRJET',phi0=0.,
                           betaeconst=BETAECONST, ximax=10., 
                           nfreq=1,fmin=RF,fmax=RF,
                           gmin=10., gmax=1.e35, p2=pp, p1=pp,
                           #ename='SYNCHPL',
                           ename='POLSYNCHPL',
                           nvals=4, fpositron=FPOSITRON,
                           spin=0., standard=1,
                           uout=uout,
                           mbh=MBH,
                           epcoefindx=[1,1,1,1,1,1,1],
                           #epcoefindx=[1,1,1,1,0,0,0],
                           mdotmin=1.57e15,mdotmax=1.57e15,nmdot=1,
                           nmu=1,mumin=mu,mumax=mu,
                           gridvals=[-size,size,-size,size],
                           nn=[npix,npix,ngeo],
                           hindf=1,hnt=1,
                           muval=1.)


    if RERUN:
        x.run_grtrans()

    # load image
    x.read_grtrans_output()
    x.convert_to_Jy(DTOBH)

    #grt_obj=x
    save_grtrans_image(x)
    display_grtrans_image(x)

def save_grtrans_image(grt_obj):
    """quick save, not ehtim compatible"""
    I_im = grt_obj.ivals[:,0,0].reshape(npix,npix).flatten()
    Q_im = grt_obj.ivals[:,1,0].reshape(npix,npix).flatten()
    U_im = grt_obj.ivals[:,2,0].reshape(npix,npix).flatten()
    V_im = grt_obj.ivals[:,3,0].reshape(npix,npix).flatten()

    # convert to Tb
    factor = 3.254e13/(RF**2 * psize_rad**2)
    I_im *= factor
    Q_im *= factor
    U_im *= factor
    V_im *= factor

    x = np.array([[i for i in range(npix)] for j in range(npix)]).flatten().astype(float)
    y = np.array([[j for i in range(npix)] for j in range(npix)]).flatten().astype(float)

    x = x - npix/2
    y = y - npix/2
    x = x*psize_uas
    y = y*psize_uas

    outdat = np.vstack((x.T,y.T,I_im.T,Q_im.T,U_im.T,V_im.T)).T
    np.savetxt('../rrjet_and_riaf/'+FNAME,outdat)
    #np.savetxt('../rrjet_and_riaf/grtrans_jet_compare_positron_noconv.txt',outdat)
    return

def display_grtrans_image(grt_obj,nvec=20,veccut=0.005,blur_kernel=1.25):
    plt.close('all')

    I_im = grt_obj.ivals[:,0,0].reshape(npix,npix)
    Q_im = grt_obj.ivals[:,1,0].reshape(npix,npix)
    U_im = grt_obj.ivals[:,2,0].reshape(npix,npix)
    V_im = grt_obj.ivals[:,3,0].reshape(npix,npix)

    I_im = filt.gaussian_filter(I_im, (blur_kernel, blur_kernel))
    Q_im = filt.gaussian_filter(Q_im, (blur_kernel, blur_kernel))
    U_im = filt.gaussian_filter(U_im, (blur_kernel, blur_kernel))
    V_im = filt.gaussian_filter(V_im, (blur_kernel, blur_kernel))

    # convert to Tb
    factor = 3.254e13/(RF**2 * psize_rad**2)
    I_im *= factor
    Q_im *= factor
    U_im *= factor
    V_im *= factor
    
    # Polarization Vectors
    P_im = np.abs(Q_im + 1j*U_im)


    m_im = P_im/I_im
    voi_im = V_im/I_im

    thin = npix//nvec
    mask = I_im > veccut * np.max(I_im)
    mask2 = mask[::thin, ::thin]

    m = m_im[::thin, ::thin][mask2]
    x = (np.array([[i for i in range(npix)] for j in range(npix)])[::thin, ::thin])
    x = x[mask2]
    y = (np.array([[j for i in range(npix)] for j in range(npix)])[::thin, ::thin])
    y = y[mask2]
    a = (-np.sin(np.angle(Q_im+1j*U_im)/2)[::thin, ::thin])
    a = a[mask2]
    #a = m*a
    b = ( np.cos(np.angle(Q_im+1j*U_im)/2)[::thin, ::thin])
    b = b[mask2]
    #b = m*b

    P_im[np.logical_not(mask)]=0.
    m_im[np.logical_not(mask)]=0.
    voi_im[np.logical_not(mask)]=0.

    # ticks
    xticks = ticks(npix, 2*size/npix)
    yticks = ticks(npix, 2*size/npix)

    # display Stokes I 
#    plt.figure(0)
#    im = plt.imshow(I_im, cmap=plt.get_cmap(cfun), interpolation='gaussian')
#    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
#    cb.set_label('Tb (K)', fontsize=14)
#    plt.title(("Stokes I, %.2f GHz " % (RF/1e9)), fontsize=16)
#    plt.xticks(xticks[0], xticks[1])
#    plt.yticks(yticks[0], yticks[1])
#    plt.xlabel('x/rg')
#    plt.ylabel('y/rg')

    # display Stokes Q 
    plt.figure(1)
    im = plt.imshow(Q_im, cmap=plt.get_cmap(cfun2), interpolation='gaussian')
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('Tb (K)', fontsize=14)
    plt.title(("Stokes Q, %.2f GHz " % (RF/1e9)), fontsize=16)
    plt.xticks(xticks[0], xticks[1])
    plt.yticks(yticks[0], yticks[1])
    plt.xlabel('x/rg')
    plt.ylabel('y/rg')

    # display Stokes U
    plt.figure(2)
    im = plt.imshow(U_im, cmap=plt.get_cmap(cfun2), interpolation='gaussian')
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('Tb (K)', fontsize=14)
    plt.title(("Stokes U, %.2f GHz " % (RF/1e9)), fontsize=16)
    plt.xticks(xticks[0], xticks[1])
    plt.yticks(yticks[0], yticks[1])
    plt.xlabel('x/rg')
    plt.ylabel('y/rg')

    # display Stokes V 
    plt.figure(3)
    im = plt.imshow(V_im, cmap=plt.get_cmap(cfun2), interpolation='gaussian')
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('Tb (K)', fontsize=14)
    plt.title(("Stokes V, %.2f GHz " % (RF/1e9)), fontsize=16)
    plt.xticks(xticks[0], xticks[1])
    plt.yticks(yticks[0], yticks[1])
    plt.xlabel('x/rg')
    plt.ylabel('y/rg')

    # display P
#    plt.figure(4)
#    im = plt.imshow(P_im, cmap=plt.get_cmap(cfun), interpolation='gaussian')
#    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
#    cb.set_label('Tb (K)', fontsize=14)
#    plt.title(("P, %.2f GHz " % (RF/1e9)), fontsize=16)
#    plt.xticks(xticks[0], xticks[1])
#    plt.yticks(yticks[0], yticks[1])
#    plt.xlabel('x/rg')
#    plt.ylabel('y/rg')

#    # display m
    plt.figure(5)
    im = plt.imshow(m_im, cmap=plt.get_cmap('jet'), interpolation='gaussian')
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('P/I', fontsize=14)
    plt.title(("P/I, %.2f GHz " % (RF/1e9)), fontsize=16)
    plt.xticks(xticks[0], xticks[1])
    plt.yticks(yticks[0], yticks[1])
    plt.xlabel('x/rg')
    plt.ylabel('y/rg')

#    # display V/I
    plt.figure(6)
    im = plt.imshow(voi_im, cmap=plt.get_cmap('jet'), interpolation='gaussian')
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('V/I', fontsize=14)
    plt.title(("V/I, %.2f GHz " % (RF/1e9)), fontsize=16)
    plt.xticks(xticks[0], xticks[1])
    plt.yticks(yticks[0], yticks[1])
    plt.xlabel('x/rg')
    plt.ylabel('y/rg')


    # display I with pol ticks
    plt.figure(7)
    im = plt.imshow(I_im, cmap=plt.get_cmap(cfun), interpolation='gaussian')
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('Tb (K)', fontsize=14)
    plt.title(("I, %.2f GHz " % (RF/1e9)), fontsize=16)
    plt.xticks(xticks[0], xticks[1])
    plt.yticks(yticks[0], yticks[1])
    plt.xlabel('x/rg')
    plt.ylabel('y/rg')

    plt.quiver(x, y, a, b,
           headaxislength=20, headwidth=1, headlength=.01, minlength=0, minshaft=1,
           width=.01*npix, units='x', pivot='mid', color='k', angles='uv', 
           scale=1.0/thin)
    plt.quiver(x, y, a, b,
           headaxislength=20, headwidth=1, headlength=.01, minlength=0, minshaft=1,
           width=.005*npix, units='x', pivot='mid', color='w', angles='uv', 
           scale=1.1/thin)

    plt.show()

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
    main()


