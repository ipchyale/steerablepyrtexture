import numpy as np
from numpy.fft import fft2, ifft2, fftshift
import matplotlib.pyplot as plt
import SPtools as SP
from scipy.stats import kurtosis
from scipy.stats import skew


def sp_expand(t,f):
    
    f = int(f)
    
    my = int(np.size(t,0))
    mx = int(np.size(t,1))
    
    mx *= f
    my *= f
    
    Te = np.zeros((my,mx), dtype=complex)
    T = np.power(f,2)*fftshift(fft2(t))
    
    y1=int(my/2+2-my/(2*f))-1
    y2=int(my/2+my/(2*f))-1
    x1=int(mx/2+2-mx/(2*f))-1
    x2=int(mx/2+mx/(2*f))-1
    Y = int(my/f)-1
    X = int(mx/f)-1
    
    Te[y1:y2+1,x1:x2+1] = T[1:Y+1,1:X+1]
    Te[y1-1,x1:x2+1] = T[0,1:X+1]/2
    Te[y2+1,x1:x2+1] = np.transpose(np.transpose(np.conjugate(T[0,list(range(X,0,-1))]/2)))
    Te[y1:y2+1,x1-1] = T[1:Y+1,0]/2;
    Te[y1:y2+1,x2+1] = np.transpose( np.transpose(np.conjugate(T[list(range(Y,0,-1)),0]/2)) )
    
    esq = T[0,0]/4
    Te[y1-1,x1-1] = esq
    Te[y1-1,x2+1] = esq
    Te[y2+1,x1-1] = esq
    Te[y2+1,x2+1] = esq
    
    Te=fftshift(Te)
    te=ifft2(Te)
    
    #tmp = np.imag(te) == 0
    #if np.sum(np.imag(te))==0:
        #te = np.real(te)
    
    return te


def texture_analyze(img,Nsc,Nor,Na, vector_output=False):
    '''Calculate texture parameters used in texture synthesis
    
    A steerable pyramid decomposition of img is created with the number of levels
    give by Nsc, the number of orientations, Nor. The number of neighbors kept
    around the central auto- and cross-correlations is determined by Na.

    The parameters are from Portilla and Simoncelli (2000). A Parametric 
    Texture Model Based on Joint Statistics of Complex Wavelet Coefficients
    
    This function is a python version of textureAnalysis.m from the texture 
    synthesis MATLAB package: 
        https://github.com/LabForComputationalVision/textureSynth


    Parameters
    ----------
    img : 2D array
        A grayscale texture image to analyze.
    Nsc : 'int'
        The number of levels in the wavelet decomposition.
    Nor : 'int'
        The number of orientation bands.
    Na : 'int'
        The number of neighbors kept in each direction around the central auto-
        or cross-correlation. MUST BE ODD. Na = 1 only incldes the central value.
    vector_output : 'bool'
        If False, the output parameters are in the form of a dictionary. 
        If True, the parameters are output as a feature vector. The default is
        False.



    Returns
    -------
    output : 
        The parameterization of the texture represented in img. See the 
        vector_output argument for the format of the output.

    '''


    if Na % 2 == 0:
       raise Exception('Na must be odd')
        
        
    la = int(np.floor((Na-1)/2))
    
    
    
    # Compute image statistics.
    mn0 = np.min(img)
    mx0 = np.max(img)
    mean0 = np.mean(img)
    var0 = np.var(img)
    skew0 = skew(img, axis=None)
    kurt0 = kurtosis(img, fisher=False, axis=None)
    
    # Add a little bit of noise to the original, in case it has been 
    # artificially generated, to avoid instability crated by symmetric
    # conditions at the synthesis stage.
    #imsz0 = np.size(img,0)
    #imsz1 = np.size(img,1)
    #im0 = img + 1e-4*(mx0 - mn0)*np.random.normal(size=(imsz0,imsz1))  
    # removed for comparisons and not intended for artificial texture analyses
    im0 = img
    
    #pyr = SP.pyramids.SteerablePyramidFreq(im0,is_complex=True)
    pyr0 = SP.pyramids.SteerablePyramidFreq(im0, height=Nsc, order=Nor-1, twidth=1, is_complex=True)
    
    
    
    for i in pyr0.pyr_size.values():
        if i[0]%2 == 1 or i[1]%2 == 1:
            raise Exception('Algorithm will fail: Some bands have odd dimensions!')
        
            
    #bandkeys = list()
    #for i in pyr0.pyr_coeffs.keys():
    #    bandkeys.append(i)
        
    bandkeys = []
    bandkeys.append('residual_highpass')
    for n in range(Nsc):
        for m in range(Nor):
            bandkeys.append((n,m))
            
    bandkeys.append('residual_lowpass')
    
    # Subtract mean of lowBand.
    pyr0.pyr_coeffs['residual_lowpass'] = np.real(pyr0.pyr_coeffs['residual_lowpass']) - np.mean(np.real(pyr0.pyr_coeffs['residual_lowpass']))
    
    # Subtract and output mean magnitude
    numbands = Nsc*(Nor)+2;
    magMeans0 = np.zeros(numbands)
    
    apyr = []
    rpyr = []
    for i in range(numbands):
        apyr.append(np.abs(pyr0.pyr_coeffs[bandkeys[i]]))
        rpyr.append(np.real(pyr0.pyr_coeffs[bandkeys[i]]))
        magMeans0[i] = np.mean(apyr[i])
        apyr[i] -= magMeans0[i]
        
    # Compute central autoCorr of lowband.
    mpyr = SP.pyramids.SteerablePyramidFreq(rpyr[numbands-1], height=0, order=0, is_complex=False)
    im = mpyr.pyr_coeffs['residual_lowpass']
    Nly = np.size(rpyr[numbands-1],0)
    Nlx = np.size(rpyr[numbands-1],1)
    Sch = np.min([Nly,Nlx])
    le = int(np.min([Sch/2-1, la]))
    cy = int(Nly/2)
    cx = int(Nlx/2)
    
    ac = fftshift(np.real(ifft2(np.power(np.abs(fft2(im)),2))))/np.size(rpyr[numbands-1])
    ac = ac[cy-le:cy+le+1, cx-le:cx+le+1]
    
    acr = np.nan * np.zeros((Na,Na,Nsc+1))
    acr[0:,0:,Nsc] = ac
    
    skew0p = np.zeros((Nsc+1,1))
    kurt0p = np.zeros((Nsc+1,1))
    
    vari = ac[le,le]
    if vari/var0 > 1e-6:
        skew0p[Nsc] = np.mean(np.power(im,3))/np.power(vari,1.5)
        kurt0p[Nsc] = np.mean(np.power(im,4))/np.power(vari,2)
    else:
        skew0p[Nsc] = 0
        kurt0p[Nsc] = 3
        
    # Compute  central autoCorr of each Mag band, and the autoCorr of the
    # combined (non-oriented) band.
    ace = np.zeros((Na,Na,Nsc,Nor))
    
    for n in range(Nsc,0,-1):
        for m in range(Nor):
            bandid = (n-1)*Nor+m+1
            ch = apyr[bandid]
            #print(bandid)
            #print(bandkeys[bandid])
            Nly = np.size(ch,0)
            Nlx = np.size(ch,1)
            Sch = np.min([Nly,Nlx])
            le = int(np.min([Sch/2-1, la]))
            cy = int(Nly/2)
            cx = int(Nlx/2)
            
            ac = fftshift(np.real(ifft2(np.power(np.abs(fft2(ch)),2))))/np.size(ch)
            ac = ac[cy-le:cy+le+1, cx-le:cx+le+1]
            ace[la-le:la+le+1,la-le:la+le+1,n-1,m] = ac
            
        bandnums = list(range((n-1)*Nor+1,n*Nor+1))
        fakepyr = SP.pyramids.SteerablePyramidFreq(np.zeros((Nly,Nlx)), height=1, order=3, twidth=1, is_complex=False)
        #print(fakepyr.pyr_size)
        for m in range(len(bandnums)):
            #print(np.size(rpyr[bandnums[m]],0))
            fakepyr.pyr_coeffs[(0,m)] = rpyr[bandnums[m]]
        
        ch = fakepyr.recon_pyr(levels=0)
        im = np.real(sp_expand(im,2))/4
        im += ch
        #plt.figure(n)
        #plt.imshow(im)
        
        ac = fftshift(np.real(ifft2(np.power(np.abs(fft2(im)),2))))/np.size(ch)
        ac = ac[cy-le:cy+le+1, cx-le:cx+le+1]
        #print(n-1)
        acr[0:,0:,n-1] = ac
        
        vari = ac[le,le]
        if vari/var0 > 1e-6:
            skew0p[n-1] = np.mean(np.power(im,3))/np.power(vari,1.5)
            kurt0p[n-1] = np.mean(np.power(im,4))/np.power(vari,2)
        else:
            skew0p[n-1] = 0
            kurt0p[n-1] = 3
            
    # Compute the cross-correlation matrices of the coefficient magnitudes
    # pyramid at the different levels and orientations
    
    C0 = np.zeros((Nor,Nor,Nsc+1))
    Cx0 = np.zeros((Nor,Nor,Nsc))
    Cr0 = np.zeros((2*Nor,2*Nor,Nsc+1))
    Crx0 = np.zeros((2*Nor,2*Nor,Nsc))
    
    for n in range(Nsc):
        firstBnum = (n)*Nor+1
        cousinSz = np.size(rpyr[firstBnum])
        
        if n<Nsc-1:
            parents = np.zeros((cousinSz,Nor))
            rparents = np.zeros((cousinSz,Nor*2))
            
            for m in range(Nor):
                bandid = (n+1)*Nor+m+1
                tmp = sp_expand(pyr0.pyr_coeffs[bandkeys[bandid]],2)/4
                rtmp = np.real(tmp)
                itmp = np.imag(tmp)
                # Double phase
                tmp = np.sqrt(np.power(rtmp,2) + np.power(itmp,2)) * np.exp(2 * 1j * np.arctan2(rtmp,itmp))
                
                rparents[:,m] = np.real(tmp.flatten())
                rparents[:,Nor+m] = np.imag(tmp.flatten())
                
                tmp = np.abs(tmp)
                tmp = tmp - np.mean(tmp)
                parents[:,m] = tmp.flatten()
                
        else:
            tmp = np.real(sp_expand(rpyr[numbands-1],2))/4
            rparents = np.zeros((np.size(tmp),5))
            rparents[:,0] = tmp.flatten()
            tmproll = np.roll(tmp, 1, axis=1)
            rparents[:,1] = tmproll.flatten()
            tmproll = np.roll(tmp, -1, axis=1)
            rparents[:,2] = tmproll.flatten()
            tmproll = np.roll(tmp, 1, axis=0)
            rparents[:,3] = tmproll.flatten()
            tmproll = np.roll(tmp, -1, axis=0)
            rparents[:,4] = tmproll.flatten()
            
            parents = []
            
        cousins = np.zeros((cousinSz,Nor))
        for m in range(Nor):
            bandid = (n)*Nor+m+1
            cousins[:,m] = apyr[bandid].flatten()
            
        nc = Nor
        if np.size(parents) > 0:
            nnp = np.size(parents,1)
        else:
            nnp = 0
                
        
        C0[0:nc,0:nc,n] = np.matmul(np.transpose(cousins),cousins)/cousinSz
        if nnp > 0:
            Cx0[0:nc,0:nnp,n] = np.matmul(np.transpose(cousins),parents)/cousinSz
            if n == Nsc-1:
                C0[0:nnp,0:nnp,Nsc] = np.matmul(np.transpse(parents),parents)/cousinSz*4
                
                
        cousins = np.zeros((cousinSz,Nor))
        for m in range(Nor):
            bandid = (n)*Nor+m+1
            cousins[:,m] = rpyr[bandid].flatten()
            
        nrc = Nor
        nrp = np.size(rparents,1)
        Cr0[0:nrc,0:nrc,n] = np.matmul(np.transpose(cousins),cousins)/cousinSz
        Crx0[0:nrc,0:nrp,n] = np.matmul(np.transpose(cousins),rparents)/cousinSz
        if n == Nsc-1:
            Cr0[0:nrp,0:nrp,Nsc] = np.matmul(np.transpose(rparents),rparents)/cousinSz*4
            
            
    ch = pyr0.pyr_coeffs[bandkeys[0]]
    vHPR0 = np.mean(np.power(ch,2))
            
    if not vector_output:
        output = dict()
        
        output['pixel_Stats'] = (mean0, var0, skew0, kurt0, mn0, mx0)
        output['pixel_LPstats'] = (skew0p, kurt0p)
        output['autocorr_Real'] = acr
        output['autocorr_Mag'] = ace
        output['mag_Means'] = magMeans0
        output['cousin_MagCorr'] = C0
        output['parent_MagCorr'] = Cx0
        output['cousin_RealCorr'] = Cr0
        output['parent_RealCorr'] = Crx0
        output['variance_HPR'] = vHPR0
        
    else:
        output = np.zeros((6,1))
        output[0] = mean0
        output[1] = var0
        output[2] = skew0
        output[3] = kurt0
        output[4] = mn0
        output[5] = mx0
        
        output = np.concatenate((output,skew0p), axis=None)
        output = np.concatenate((output,kurt0p), axis=None)
        output = np.concatenate((output,acr.transpose()), axis=None)
        output = np.concatenate((output,ace.transpose((2,3,1,0))), axis=None)
        output = np.concatenate((output,magMeans0), axis=None)
        output = np.concatenate((output,C0.transpose()), axis=None)
        output = np.concatenate((output,Cx0.transpose()), axis=None)
        output = np.concatenate((output,Cr0.transpose()), axis=None)
        output = np.concatenate((output,Crx0.transpose()), axis=None)
        output = np.concatenate((output,vHPR0), axis=None)        


    return output        
            
















    