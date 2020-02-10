# Run Grtrans with rrjet model
# The rrjet model is defined in "fluid_model_rrjet.py"
# NOTE -- currently the power law emissivity is very slow because paralleization is off

# First make grtrans with 'make' 
# Then run this in python

import numpy as np
import grtrans_batch as gr
import matplotlib.pyplot as plt
import scipy.ndimage.filters as filt

ang=45.
mu = np.cos(ang*np.pi/180.)
size  = 25.
uout = 1./(2*size)
npix = 100
ngeo = 800

cmperKpc = 3.086e21
MBH = 4.1e6
DTOBH = 8.1*cmperKpc
RADPERUAS = np.pi/180./3600./1.e6
psize_rg = 2*size/npix
cmperrg = 147708.8 * MBH
psize_cm = psize_rg * cmperrg
psize_rad = psize_cm / DTOBH
psize_uas = psize_rad / RADPERUAS

RF = 230.e9
cfun = 'afmhot'
cfun2 =  'Spectral'
RERUN = True

#Broderick&Loeb 06 parameters
NSCL=1.e7
TSCL=1.5e11
NNTHSCL=1.e5
BETA=10.
PNTH=2.8
FPOSITRON=0  #0 <npositron/nelectron <1
A = .998

name = ('riaf_%0.1f_'%FPOSITRON)+str(ang)+'_im'
FNAME = ('riaf_%0.1f_'%FPOSITRON)+str(ang)+'_im'

def main():
    # run grtrans

    x=gr.grtrans()

    x.write_grtrans_inputs(name+'.in', oname=name+'.out',
                           fname='SARIAF', phi0=0.,
                           #betaeconst=1.e-4, ximax=10., 
                           nfreq=1,fmin=RF,fmax=RF,
                           gmin=100., p2=4.5, p1=4.5,
                           snscl=NSCL, ntscl=TSCL, snnthscl=NNTHSCL, 
                           snnthp=PNTH, sbeta=BETA, sbl06=1,
                           fpositron=FPOSITRON,
                           ename='HYBRIDTHPL',
                           nvals=4,
                           spin=A,standard=1,
                           uout=uout,
                           mbh=MBH,
                           #mdotmin=1.57e15,mdotmax=1.57e15,nmdot=1,
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

    tmax = np.max( x.ivals[:,0,0]*3.254e13/(RF**2 * psize_rad**2))
    pmax = np.max( np.sqrt(x.ivals[:,1,0]**2 + x.ivals[:,2,0]**2)*3.254e13/(RF**2 * psize_rad**2))

    #tmax=5.e10
    #pmax=2.e10
    save_grtrans_image(x)
    display_grtrans_image(x,tmax=tmax,pmax=pmax)
    #x.disp_grtrans_image()
    #x.disp_pol_map()

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

    x = np.array([[i for i in range(npix)] for j in range(npix)]).flatten()
    y = np.array([[j for i in range(npix)] for j in range(npix)]).flatten()

    x -= npix/2
    y -= npix/2
    x = x*psize_uas
    y = y*psize_uas

    outdat = np.vstack((x.T,y.T,I_im.T,Q_im.T,U_im.T,V_im.T)).T
    np.savetxt('../rrjet_and_riaf/'+FNAME,outdat)
    return

def display_grtrans_image(x,nvec=25,veccut=0.005,tmax=1.e10,pmax=1.e10,blur_kernel=0):#1.25):
    plt.close('all')

    #I_im = np.transpose(x.ivals[:,0,0].reshape(npix,npix))
    #Q_im = np.transpose(x.ivals[:,1,0].reshape(npix,npix))
    #U_im = np.transpose(x.ivals[:,2,0].reshape(npix,npix))
    #V_im = np.transpose(x.ivals[:,3,0].reshape(npix,npix))

    #I_im = x.ivals[:,0,0].reshape(npix,npix)
    #Q_im = x.ivals[:,1,0].reshape(npix,npix)
    #U_im = x.ivals[:,2,0].reshape(npix,npix)
    #V_im = x.ivals[:,3,0].reshape(npix,npix)

    I_im =  (np.flipud(np.transpose(x.ivals[:,0,0].reshape((npix,npix)))))
    Q_im =  -(np.flipud(np.transpose(x.ivals[:,1,0].reshape((npix,npix)))))
    U_im =  -(np.flipud(np.transpose(x.ivals[:,2,0].reshape((npix,npix)))))
    V_im =  (np.flipud(np.transpose(x.ivals[:,3,0].reshape((npix,npix)))))

    #I_im = filt.gaussian_filter(I_im, (blur_kernel, blur_kernel))
    #Q_im = filt.gaussian_filter(Q_im, (blur_kernel, blur_kernel))
    #U_im = filt.gaussian_filter(U_im, (blur_kernel, blur_kernel))
    #V_im = filt.gaussian_filter(V_im, (blur_kernel, blur_kernel))

    # convert to Tb
    factor = 3.254e13/(RF**2 * psize_rad**2)
    I_im *= factor
    Q_im *= factor
    U_im *= factor
    V_im *= factor
    
    # Polarization Vectors
    P_im = np.abs(Q_im + 1j*U_im)
    m_im = P_im/I_im


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

    # ticks
    xticks = ticks(npix, 2*size/npix)
    yticks = ticks(npix, 2*size/npix)

    # display Stokes I 
    plt.figure(0)
    im = plt.imshow(I_im, cmap=plt.get_cmap(cfun), interpolation='gaussian',vmin=0,vmax=tmax)
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('Tb (K)', fontsize=14)
    plt.title(("Stokes I, %.2f GHz " % (RF/1e9)), fontsize=16)
    plt.xticks(xticks[0], xticks[1])
    plt.yticks(yticks[0], yticks[1])
    plt.xlabel('x/rg')
    plt.ylabel('y/rg')

    # display Stokes Q 
    plt.figure(1)
    im = plt.imshow(Q_im, cmap=plt.get_cmap(cfun2), interpolation='gaussian',vmin=-pmax,vmax=pmax)
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('Tb (K)', fontsize=14)
    plt.title(("Stokes Q, %.2f GHz " % (RF/1e9)), fontsize=16)
    plt.xticks(xticks[0], xticks[1])
    plt.yticks(yticks[0], yticks[1])
    plt.xlabel('x/rg')
    plt.ylabel('y/rg')

    # display Stokes U
    plt.figure(2)
    im = plt.imshow(U_im, cmap=plt.get_cmap(cfun2), interpolation='gaussian',vmin=-pmax,vmax=pmax)
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('Tb (K)', fontsize=14)
    plt.title(("Stokes U, %.2f GHz " % (RF/1e9)), fontsize=16)
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
    plt.title(("Stokes V, %.2f GHz , V/I=%.2f " % (RF/1e9, voi)), fontsize=16)
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
    plt.title(("P, %.2f GHz , P/I=%.2f" % (RF/1e9,poi)), fontsize=16)
    plt.xticks(xticks[0], xticks[1])
    plt.yticks(yticks[0], yticks[1])
    plt.xlabel('x/rg')
    plt.ylabel('y/rg')

    # display m
    plt.figure(5)
    im = plt.imshow(m_im, cmap=plt.get_cmap('viridis'), interpolation='gaussian')
    cb = plt.colorbar(im, fraction=0.046, pad=0.04, orientation="vertical")
    cb.set_label('P/I', fontsize=14)
    plt.title(("P/I, %.2f GHz , P/I=%.2f" % (RF/1e9,poi)), fontsize=16)
    plt.xticks(xticks[0], xticks[1])
    plt.yticks(yticks[0], yticks[1])
    plt.xlabel('x/rg')
    plt.ylabel('y/rg')

    # display I with pol ticks
    plt.figure(6)
    im = plt.imshow(I_im, cmap=plt.get_cmap(cfun), interpolation='gaussian',vmin=0,vmax=tmax)
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


