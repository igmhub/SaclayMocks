''' FGPA F = exp -a exp b g
    computes the pdf of F, <F> and <F^2>
'''
import numpy as np
import scipy as sp
from scipy import integrate
from iminuit import Minuit
#from iminuit.cost import LeastSquares


def F_pdf(F,ln_a,b):
    '''	
    pdf of F=exp[ -a exp(b g)] = exp[-exp(bg+ln(a))]  
	where g = gaus(0,1)  or bg+ln(a) = Gaus(ln(a),b)
	1/(F tau \sqrt{2PI} b) exp [-(ln(tau)-ln(a))^2/(2b^2)]    tau = -ln(F)
    '''
    if (b==0): print("erro b=0 in pdf\n") # exit(0)
    if (F<=0 or F>=1): print("error F=",F," in pdf\n") # exit(0)
    tau = - np.log(F)   # log is ln != log10
    xx = np.log(tau)-ln_a
    xx = np.exp(- xx*xx/2/b/b)
    xx /= F * tau * np.sqrt(2 * np.pi) * b
    return xx

def Ftimespdf(F, ln_a, b) :
    # F * pdf[F]
    return F*F_pdf(F,ln_a,b)

def F2timespdf(F, ln_a, b) :
    # F^2 * pdf[F]
    return F*F*F_pdf(F,ln_a,b)

def pdf_integrale(Fmax, ln_a, b) :
  # compute int_0^Fmax pdf(F)dF
  u = ln_a -np.log(-np.log(Fmax))
  return 0.5 * ( 1 + sp.special.erf( u/b/np.sqrt(2) ) )


def Fmean(ln_a,b):
    # <F> = \int_0^1 F pdf(F) dF
    # split integrale [0,1-eps] and [1-eps,1]
    # approximate F~1 in [1-eps,1] 
    epsilon = 1e-8
    F_mean = integrate.quad(Ftimespdf,0,1-epsilon,args=(ln_a,b))[0]  # int_0^1-eps
    xx= 1 - pdf_integrale(1-epsilon,ln_a,b) # int_{1-eps}^1
    #print (F_mean)
    #print (xx)
    #print (F_mean + xx)
    return F_mean + xx
  
def F2mean(ln_a,b):
    epsilon = 1e-8
    F2_mean = integrate.quad(F2timespdf,0,1-epsilon,args=(ln_a,b))[0]
    xx = 1 - pdf_integrale(1-epsilon,ln_a,b);
    return F2_mean + xx

def sig2F(ln_a,b):
    return F2mean(ln_a,b) - Fmean(ln_a,b)**2

def FF2sig2(ln_a,b):
    F_mean=Fmean(ln_a,b)    
    F2_mean=F2mean(ln_a,b)
    sig2=F2_mean-F_mean**2
    return F_mean,F2_mean,sig2

def chi2(ln_a,b,Fmean,sigF) :
    F_mean,_,sig2_F=FF2sig2(ln_a,b)
    xx = (Fmean-F_mean)**2 + (sigF - np.sqrt(sig2_F))**2
    return 1E6*xx

def fitab(F_mean,sigF,ln_a=-4,b=3):
    m = Minuit(chi2, ln_a=ln_a, b=b, Fmean=F_mean, sigF=sigF)
    m.errordef = 1 # Minuit.LEAST_SQUARES # i.e. 1
    m.fixed["Fmean"] = True
    m.fixed["sigF"] = True
    m.migrad()  # run optimiser
    result=m.values
    return result['ln_a'], result['b']

