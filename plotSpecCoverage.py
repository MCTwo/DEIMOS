import pylab

redshift = 0.092
lambda_central = 6300

fig = pylab.figure(figsize=(20,4.5))
pylab.xlim((lambda_central-2300, lambda_central+2300))
xl = pylab.xlim()
yl = (0,1)

lcl_low = lambda_central-1300
lcl_up = lambda_central+1300

#Plot the wavelength coverage
pylab.plot((lambda_central,lambda_central),yl,'-b',linewidth=3)
pylab.plot((lcl_low,lcl_low),yl,
           '-b',alpha=0.5,linewidth=3)
pylab.plot((lcl_up,lcl_up),yl,
           '-b',alpha=0.5,linewidth=3)
pylab.fill_between(pylab.arange(lcl_low,lcl_up+100,100),yl[0],yl[1],
                   facecolor='blue',alpha=0.25)

#Plot common spectral lines
x_Lyb = 1025.7*(1+redshift)
x_Lya = 1215.7*(1+redshift)
x_CIV = 1549.1*(1+redshift)
x_AlIII = 1858.7*(1+redshift)
x_FeII = 2600*(1+redshift)
x_MgII = 2799.8*(1+redshift)
x_MgI = 2852*(1+redshift)
x_OII = 3727.61*(1+redshift)
x_CalK = 3933.667*(1+redshift)
x_CalH = 3968.472*(1+redshift)
x_Hd = 4101.74*(1+redshift)
x_Gband = 4305*(1+redshift)
x_Hg = 4340.47*(1+redshift)
x_Hb = 4861.33*(1+redshift)
x_OIII_1 = 4960.3*(1+redshift)
x_OIII_2 = 5008.24*(1+redshift)
x_Mgb = 5176*(1+redshift)
x_FeI = 5269*(1+redshift)
x_NaD = 5893*(1+redshift)
x_NII_1 = 6548.06*(1+redshift)
x_Ha = 6562.799*(1+redshift)
x_NII_2 = 6585.2*(1+redshift)
x_SII = 6725.5*(1+redshift)

# Plot dashed lines at the respective redshifts
pylab.plot((x_Lyb,x_Lyb),yl,'--k')
pylab.plot((x_Lya,x_Lya),yl,'--k')
pylab.plot((x_CIV,x_CIV),yl,'--k')
pylab.plot((x_AlIII,x_AlIII),yl,'--k')
pylab.plot((x_FeII,x_FeII),yl,'--k')
pylab.plot((x_MgII,x_MgII),yl,'--k')
pylab.plot((x_MgI,x_MgI),yl,'--k')
pylab.plot((x_OII,x_OII),yl,'--k')
pylab.plot((x_CalK,x_CalK),yl,'--k')
pylab.plot((x_CalH,x_CalH),yl,'--k')
pylab.plot((x_Hd,x_Hd),yl,'--k')
pylab.plot((x_Gband,x_Gband),yl,'--k')
pylab.plot((x_Hg,x_Hg),yl,'--k')
pylab.plot((x_Hb,x_Hb),yl,'--k')
pylab.plot((x_OIII_1,x_OIII_1),yl,'--k')
pylab.plot((x_OIII_2,x_OIII_2),yl,'--k')
pylab.plot((x_Mgb,x_Mgb),yl,'--k')
pylab.plot((x_FeI,x_FeI),yl,'--k')
pylab.plot((x_NaD,x_NaD),yl,'--k')
pylab.plot((x_NII_1,x_NII_1),yl,'--k')
pylab.plot((x_Ha,x_Ha),yl,'--k')
pylab.plot((x_NII_2,x_NII_2),yl,'--k')
pylab.plot((x_SII,x_SII),yl,'--k')

labeloff = 0.5
pylab.text(lambda_central, labeloff*(yl[0]+yl[1]),
           '$\lambda_{central}$'+'={0}'.format(lambda_central), 
           horizontalalignment='right',verticalalignment='center', 
           rotation='vertical')

if x_Lyb > xl[0] and x_Lyb < xl[1]:    
    pylab.text(x_Lyb, labeloff*(yl[0]+yl[1]), 'Ly-beta', horizontalalignment='right',verticalalignment='center', rotation='vertical')
if x_Lya > xl[0] and x_Lya < xl[1]:    
    pylab.text(x_Lya, labeloff*(yl[0]+yl[1]), 'Ly-alpha', horizontalalignment='right',verticalalignment='center', rotation='vertical')
if x_CIV > xl[0] and x_CIV < xl[1]:    
    pylab.text(x_CIV, labeloff*(yl[0]+yl[1]), 'C IV', horizontalalignment='right',verticalalignment='center', rotation='vertical')
if x_AlIII > xl[0] and x_AlIII < xl[1]:    
    pylab.text(x_AlIII, labeloff*(yl[0]+yl[1]), 'Al III', horizontalalignment='right',verticalalignment='center', rotation='vertical')
if x_FeII > xl[0] and x_FeII < xl[1]:    
    pylab.text(x_FeII, labeloff*(yl[0]+yl[1]), 'Fe II', horizontalalignment='right',verticalalignment='center', rotation='vertical')
if x_MgII > xl[0] and x_MgII < xl[1]:    
    pylab.text(x_MgII, labeloff*(yl[0]+yl[1]), 'Mg II', horizontalalignment='right',verticalalignment='center', rotation='vertical')
if x_MgI > xl[0] and x_MgI < xl[1]:    
    pylab.text(x_MgI, labeloff*(yl[0]+yl[1]), 'Mg I', horizontalalignment='right',verticalalignment='center', rotation='vertical')
if x_OII > xl[0] and x_OII < xl[1]:    
    pylab.text(x_OII, labeloff*(yl[0]+yl[1]), '[O II]', horizontalalignment='right',verticalalignment='center', rotation='vertical')
if x_CalK > xl[0] and x_CalK < xl[1]:    
    pylab.text(x_CalK, labeloff*(yl[0]+yl[1]), 'Cal K', horizontalalignment='right',verticalalignment='center', rotation='vertical')
if x_CalH > xl[0] and x_CalH < xl[1]:    
    pylab.text(x_CalH, labeloff*(yl[0]+yl[1]), 'Cal H', horizontalalignment='right',verticalalignment='center', rotation='vertical')
if x_Hd > xl[0] and x_Hd < xl[1]:    
    pylab.text(x_Hd, labeloff*(yl[0]+yl[1]), 'Hd', horizontalalignment='right',verticalalignment='center', rotation='vertical')
if x_Gband > xl[0] and x_Gband < xl[1]:    
    pylab.text(x_Gband, labeloff*(yl[0]+yl[1]), 'G-band', horizontalalignment='right',verticalalignment='center', rotation='vertical')
if x_Hg > xl[0] and x_Hg < xl[1]:    
    pylab.text(x_Hg, labeloff*(yl[0]+yl[1]), 'Hg', horizontalalignment='right',verticalalignment='center', rotation='vertical')
if x_Hb > xl[0] and x_Hb < xl[1]:    
    pylab.text(x_Hb, labeloff*(yl[0]+yl[1]), 'Hb', horizontalalignment='right',verticalalignment='center', rotation='vertical')
if x_OIII_1 > xl[0] and x_OIII_1 < xl[1]:    
    pylab.text(x_OIII_1, labeloff*(yl[0]+yl[1]), '[OIII]', horizontalalignment='right',verticalalignment='center', rotation='vertical')
if x_OIII_2 > xl[0] and x_OIII_2 < xl[1]:    
    pylab.text(x_OIII_2, labeloff*(yl[0]+yl[1]), '[OIII]', horizontalalignment='right',verticalalignment='center', rotation='vertical')
if x_Mgb > xl[0] and x_Mgb < xl[1]:    
    pylab.text(x_Mgb, labeloff*(yl[0]+yl[1]), 'Mg I(b)', horizontalalignment='right',verticalalignment='center', rotation='vertical')
if x_FeI > xl[0] and x_FeI < xl[1]:    
    pylab.text(x_FeI, labeloff*(yl[0]+yl[1]), 'Fe I', horizontalalignment='right',verticalalignment='center', rotation='vertical')
if x_NaD > xl[0] and x_NaD < xl[1]:    
    pylab.text(x_NaD, labeloff*(yl[0]+yl[1]), 'Na I (D)', horizontalalignment='right',verticalalignment='center', rotation='vertical')
if x_NII_1 > xl[0] and x_NII_1 < xl[1]:    
    pylab.text(x_NII_1, (labeloff-0.1)*(yl[0]+yl[1]), '[NII]', horizontalalignment='right',verticalalignment='center', rotation='vertical')
if x_Ha > xl[0] and x_Ha < xl[1]:    
    pylab.text(x_Ha, labeloff*(yl[0]+yl[1]), 'Ha', horizontalalignment='right',verticalalignment='center', rotation='vertical')
if x_NII_2 > xl[0] and x_NII_2 < xl[1]:    
    pylab.text(x_NII_2, (labeloff+0.1)*(yl[0]+yl[1]), '[NII]', horizontalalignment='right',verticalalignment='center', rotation='vertical')
if x_SII > xl[0] and x_SII < xl[1]:    
    pylab.text(x_SII, labeloff*(yl[0]+yl[1]), '[SII]', horizontalalignment='right',verticalalignment='center', rotation='vertical')

pylab.xlim(xl)
frame1 = pylab.gca()
frame1.axes.get_yaxis().set_visible(False)
pylab.xlabel('$\lambda_{observed}$',fontsize=14)

pylab.show()