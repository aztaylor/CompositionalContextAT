#!/bin/python3.11
import copy
import numpy as np 
import matplotlib.pyplot as plt

from scipy.integrate import solve_ivp

__author__ = "Aleczadner Taylor, Bahameh Banadaki, Encoh Yeung"
__data__ = "December 13th 2024"
__description__ = "Models Compositional Context of synthetic Genes."
__maintainer__ = "Aleczander Taylor"
__contact__="aztaylor76@fastmail.com"
__status__="Development"

"""This model is based off of the one developed in Biophysical 
Restraints Arising from Compositional Context in Synthetic Gene Networks
While we have written this code ourselves, it is based on similar 
code develped by Yeung et al. The original code is accesssible at 
the Yeung Lab Github Repository (). Or contributions include porting 
to an open source programming language, modification to an OOP approach 
(allowing for inheritence), and expolaration of the effect of 
DNA binding protiens and DNA binding protein-gyrase fusions on 
the simulation."""

# First we can define our parameters
hours = 3.1
atol=10**(-5) # Absolute error for ODE solver

#Create a dotdict class to make paramaeter and initial value handeling easier.
class dotdict(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

# and now the paparameters supplied to the system. See the system for the definition of each.

# And the intitial Conditions: Note that these need to allign with the states. 
ivs= {
      'stS0':-6,
      'stG0':-3,
      'spS0':-6,
      'spG0':-3,
      }

ivs = dotdict(ivs)

# These are the parameter supplied by the paper, not the github repo.
p  = {'h0': 10.5, # bps per right handed turn B-DNA.[bpturn]
      'TLs': 203, # mSpinach and T500 term len. [bp]
      'TLG': 68, # MG aptamer and T500 transcript len. [bp]
      'Ns': 150, # Intergenic Spacer Length. [bp]
      'PLs': 40, # Lac promotor Length. [bp]
      'PLG': 44, # Tet promoter length. [bp]
      'PLen': 2892, # Palsmid Length [bp]
      'kmsigma':50 , # Michailis-menten constant for supercoiling hill functions. [nM]
      'Rtot': 18.931, # Total RNAP concentration uM
      'Ribotot': 11.291, # Total Ribosome Concentration. [uM]
      'pTetTot': 11, # Total plamid consentration. [nM]
      'pLacTot': 11, # Total Plasmid Concentration. [nM]
      'G0': 12, # Gyrasae concentration. [nM]
      'T0': 2, # Topoisomerase Concentration. [nM]
      'LacItot': 10, # LacI Concentration. [nM]
      'TetRtot': 0, # TetR Concentration. [nM]
      'IPTGtot': 0, # IPTG Concentration. [nM]
      'aTcTot': 0 , # aTc Concentration. [nM]
      'kfmaxl': 7*10**(-2), # Leakynforward transcription initiation rate. [nM/s]
      'kr': 550, # Reverse Transcription Initiation Rate. [nM/s]
      'kcatG':85,# Per base-pair transcription rate. [nt/s]
      'kcatmax': lambda: p['kcatG']/105.5, # Averaged transcription rate for aptamers. []
      'kopen': 0.04,# Rate of open complex formation. [1/s]
      'k1': 0.02, # Fraction of terminator escaped transripts. [n/a]
      'rhoI': 0, # Constitutive production rate of lacI.[nM/s] (0 for cell free)
      'rhoT': 0, # Constitutive production rate of TetR. [nM/s] (0 for cell feee)
      'k0TL': 21, # Averaged Translation rate of CFP RFP.
      'TLC': 714, # Length of CFP ORF. [nt] )sourced from fpbase.org
      'TLR': 675, # Length of RFP ORF. [nt] sourced from fpbase.org
      'kfC': 1/(30*60), # Folding rate of CFP. [1/s]
      'kfR':1/(110*60), # Folding rate of RFP. [1/s]
      'kbgl':1.25, # Relative Ribosomal Affinity for BCGbgl [n/a]
      'kBCD16': 0.75,# Relative Ribosome Affinity for BCD16 (Mutalik) [n/a]
      'deltams': np.log(2)/(30*60), # mSpinach degredation rate. [1/s]
      'deltamg': np.log(2)/(60*60), # MG degredation rate. [1/s]
      'NS': 150, # Determined by the experimental design
      'tau': 0.5,  # Ngeative supercoils introduced per sec by TopoI. [turn/s]
      'gamma': 0.5, # Positive supercoils intridcued per sec by Gyrase. [turn/s]
      'sigma0': -0.065, # Natural Superhelical Density of DNA
      'kmGyr':200, # Hill coefficient for gyrase maintanence function. 
      'kal': 6*10**(3), # IPTG binding rate to LacI. [1/s]
      'kuaL': 1, # IPTG dissosociation rate with LacI. [1/s]
      'kseqL': 10, # Binding rate to RNAP to Lac promoter. [1/s]
      'kuL':.022, # Dissociation (w/o elongation) rate of RNAP to Lac promoter. [1/s]
      'kaT':lambda: p['kal'], # aTc binding rate to TetR.[1/s]
      'kuaT':lambda: p['kuaL'], # Bidning rate of TetR to operator. [1/s]
      'kseqT': lambda: p['kseqL'], # Dissociation w/o elongation of TetR to operator. [1/s]
      'kuT': lambda: p['kuL'],#*30, Fall of rate of TetR from DNA [1/s]
      'deltaP': 0, # General degredation rate of proteins in myTXTL [1/s]
      'kfolddimer': 1 # dimeriztion rate of RNA aptamers. [1/s]
    }

p = dotdict(p)
print(p.keys())


class Systems:
    def __init__(self, p:all, ivs:all):
        self.p = p
        self.ivs = ivs
    def initial_conditions(self):
        '''A similar function is used by the authors to define ininitial conditions, also ensures they are in the correct order '''
        x0 = np.zeros((19))
        x0[0] = 0 #mS0
        x0[1] = 0 #mMG0
        x0[2] = 0 #ECS0
        x0[3] = 0 #ECG0
        x0[4] = self.ivs['stS0'] #stS0
        x0[5] = self.ivs['stG0'] #stG0
        x0[6] = self.ivs['spS0'] #sps0
        x0[7] = self.ivs['spG0'] #spG0
        x0[8] = 0 #CCS0
        x0[9]  = 0 #CCG0
        x0[10] = self.p['LacItot'] #LacI0
        x0[11] = self.p['TetRtot'] #TetR0
        x0[12] = 1 #IPTG0
        x0[13] = 100 #aTc0
        x0[14] = 0 #PmS0
        x0[15] = 0 #PMG0
        x0[16] = 0 #fmSG0
        x0[17] = 0 # fMG0
        x0[18] = 0 # ffmG0
        return x0
    def  kf(self, t, sigma, promoter_len):
        """This equation gives the maturation time for the sighnal of interest,
        mSpinach, MG, CFP, RFP, etc... It is essenitally a rate from transcription to to
        experession"""
        p = self.p
        kf = (p.sigma0*p.kfmaxl)/(p.sigma0+(sigma-p.sigma0*promoter_len/p.PLen)**2)
        return kf

    def kcat(self, t, sigma, promoter_len):
        """This equation determines the rate of translation of each species."""
        p = self.p
        kcat = (p.sigma0*p.kcatmax())/(p.sigma0+(sigma - p.sigma0*promoter_len/p.PLen)**2)
        return kcat

    def m_gyrtopfunc(self, t, sigma):
        """Method which describes the activity of Gyrase, at hetero-quatramer enzyyme which relaxes positive suercoiling in 
        plasmids and chromosomes.
        Arguments: t[array-like]: Series of timepoints to evaluate.
                   sigma[float]: Supersoiling state. Continously updatated
                   pDNA[float]: Amount of plasmid DNA.
                   sigma[float]: Initital supercoiling state
                   t0[float]: Initial time.
                   g0[float]: Initial Gyrase Concentration.
                   tau:[float]:
                   gamma[float]:
                   
        Returns:
                gyrase activivity[float]:n*Hill_coeff*(t0+tau/pDNA)*(G0+gamma/pDNA)"""

        s0m = np.abs(p.sigma0)
        pX=0
        nX=0
        if sigma>0:
            pX = sigma# Positive supercoiling
        if sigma<0:    
            nX = sigma # Negative Supersoiling
        hill_coeff = np.abs(sigma-s0m)/p.kmGyr/(s0m+((sigma-s0m)/p.kmGyr)**2)
        m = nX*hill_coeff*(p.T0*p.tau/p.PLen)-pX*(hill_coeff)*(p.G0*p.gamma/p.PLen)
        return m
    def kTlfunc(self, t):
        nhours = 2
        p=self.p
        if t>=60*60*hours:
            half_life = 80*60
            alpha_decay = np.log(2)/half_life
            ksc = np.exp(-alpha_decay*(t-nhours*60*60))
            kTLeff = ksc+p.k0TL
        else:
            kTLeff = p.k0TL
        return kTLeff

    # Define the hill equations that defined the transcription rate.
    def convergent_TL_dynamics(self, t:all,y:all):
        '''The system of ordinary differential equations that satisfy the 
        convergent system expressing mSpinsach and MS.
        The states are y0=mS, y1=MG, y2=ECs, y3=ECG, y4=stS, y5=stG, y6=spS, y7=spG, y9=CCS,y9=CCG,y10=LacI, y11=TetR, y12=IPTG, y13=aTc, y14=PmS, y15=PMG, y16=fmS, y17=fmG, y18=ffMG.
        a[1] '''
        p = self.p
        # Conservation law equations for the active and deactive 
        # promoters and the repressors.
        pLacC = p.pLacTot-y[10]-p.IPTGtot+y[12] # Amount of bound and sequestered pLac
        pTetC = p.TetRtot-y[11]-p.aTcTot+y[13] # Similar for pTet
        pLac = p.pLacTot-y[2]-y[8]-pLacC # Free pLac
        pTet = p.TetRtot-y[3]-y[9]-pTetC # Same
        #R = p.Rtot-y[2]-y[3]-y[8]-y[9] # Ribosome present

        #Rates defined by the modifided hill functions (related to 
        #expressions/folding). Derived in our summary and the 
        #original paper.
        kf_mS = self.kf(t, y[6], p.PLs)
        kf_MG = self.kf(t, y[7], p.PLG)

        # Rates of transcription. Derived in our summary and in the paper.
        kcat_mS = self.kcat(t, y[4], p.TLs)
        kcat_MG = self.kcat(t, y[5], p.TLG)
        
        ydot = np.ndarray(self.initial_conditions().shape)

        # Now our ODEs.
        #dmSdt               
        ydot[0] = kcat_mS*y[2]-p.deltams*y[0]
        ydot[1] = kcat_MG*y[3]-p.deltamg*y[1]
        ydot[2] = p.kopen*y[8]-kcat_mS*y[2]
        ydot[3] = p.kopen*y[9]-kcat_MG*y[3]
        dKink = (y[4]+y[5])*p.h0
        nfs = p.PLs+p.TLs+p.NS/(2)*pTet/(pTet+pTetC)+(p.NS/2+p.TLG)*pTetC/(pTet+pTetC)-dKink
        nfG = p.PLG+p.TLG+p.NS/(2)*pLac/(pLac+pLacC)+(p.NS/2+p.TLs)*pLacC/(pLac+pLacC)-dKink
        ydot[8] = kf_mS*(p.Rtot-y[2]-y[3]-y[8]-y[9])*(pLac-y[8]-y[2]-pLacC)-(p.kr+p.kopen)*y[8]
        ydot[9] = kf_MG*(p.Rtot-y[2]-y[3]-y[8]-y[9])*(p.pTetTot-y[9]-y[3]-pTetC)-(p.kr+p.kopen)*y[9]
        ydot[10] = p.rhoI-p.kal*y[10]*y[12]+p.kuaL*(p.LacItot-y[10]-pLacC)+p.kuL*pLacC-p.kseqL*pLac*y[10]-p.deltaP*y[10]
        ydot[11] = p.rhoT-p.kal*y[11]*y[13]+p.kuaT()*(p.TetRtot-y[11]-pTetC)+p.kuT()*pTetC-p.kseqT()*pTet*y[11]-p.deltaP*y[11]
        ydot[12] = -p.kal*(y[10]+pLacC)*y[12]+p.kuaL*(p.LacItot-y[10]-pLacC) 
        ydot[13] = -p.kaT()*(y[11]+pTetC)*y[13]+p.kuaT()*(p.TetRtot-y[11]-pTetC)
        ydot[14] = p.kBCD16*self.kTlfunc(t)/(p.TLs/3)*p.Ribotot*y[0] - p.kfC*y[14]
        ydot[15] = p.kbgl*self.kTlfunc(t)/(p.TLG/3)*p.Ribotot*y[1]-p.kfR*y[15]
        ydot[16] = p.kfC*y[14]
        
        # Dimerization of MG. 
        ydot[17] = p.kfR*y[15] - p.kfolddimer*y[17]
        ydot[18] = p.kfolddimer*y[17]

        # Update the states of the DNA complexes as wekk as the rate 
        #of transcription of the aptamers.
        dmScrdt = ydot[0] + ydot[0]*y[0]
        dMGcrdt = ydot[1] + ydot[1]*y[1]


        # Updates t the super coiling states of the promoters, 
        #considering tge effect of gyrase.
        ydot[4] = 1/(2*p.h0)*((ydot[2]-ydot[8])*p.PLs/nfs+(dmScrdt-ydot[2])*p.TLs/(nfs)) + self.m_gyrtopfunc(t,y[4])
        ydot[5] = 1/(2*p.h0)*((ydot[3]-ydot[9])*p.PLG/nfG+(dMGcrdt-ydot[3])*p.TLG/nfG) + self.m_gyrtopfunc(t,y[5])

        ydot[6] = 1/(2*p.h0)*(-(ydot[2]-ydot[8])*p.PLs/nfs)+self.m_gyrtopfunc(t,y[6])
        ydot[7] = 1/(2*p.h0)*(-(ydot[2]-ydot[9])*p.PLG/nfG)+self.m_gyrtopfunc(t,y[7])

        return ydot

# Lets try it out with two concentrations of aTc:
IPTGconcs = [0,100]
tint = [0, hours*60*60]
p1 = copy.copy(p)
p2 = copy.copy(p)
p2['IPTGtot'] = 100
Sys1 = Systems(p1,ivs)
Sys2 = Systems(p2, ivs)
y0 = Sys1.initial_conditions()
sol1 = solve_ivp(Sys1.convergent_TL_dynamics, tint, y0, vectorized=True, method='LSODA', atol=atol)
sol2 = solve_ivp(Sys2.convergent_TL_dynamics, tint,
                 y0, vectorized=True, atol=atol)
fig, ax =plt.subplots(1,2,sharey=True)
ax[0].plot(sol1.t, sol1.y[0])
ax[0].plot(sol1.t, sol1.y[1])
ax[0].set_title('Convergent Orientation 100ng/ul aTx\n 0mM IPTG')
ax[0].set_ylabel('Fluorescence (RFU)')
ax[0].set_xlabel('Time[S]')
ax[0].legend(['mSpinach', 'MG'])
ax[1].set_title('Convergent Orientation 100ng/ul aTx\n 100mM IPTG')
ax[1].set_ylabel('Fluorescence (RFU)')
ax[1].set_xlabel('Time[S]')
ax[1].legend(['mSpinach', 'MG'])

ax[1].plot(sol2.t, sol2.y[0])
ax[1].plot(sol2.t, sol2.y[1])
plt.tight_layout()
plt.savefig('plot1.png')

class dCas12a_System(Systems):
    def __init__(self, p, ivs):
        super().__init__(p, ivs)
        self.p = p 
        self.ivs = ivs
    def convergent_TL_dynamics(self,t,y):
        # First lets define new state # Concentration of free guideRNA
        CasgRNA = y[19] # Concentration of bound gRNA
        print(self.p)
         ## %%
        kcas=30 
        dydt = super().convergent_TL_dynamics(t,y)
        #Ts = 150/2 #Spacer length
        pLen = self.p['pLen']
        h0 = self.p['h0']
        nfs = self.p['ns']
        nfg = self.p['ns']
        TLS = self.p['TLS']
        TLG = self.p['TLMG']
        mS = y[0]
        MG = y[1]
        stS = y[4]
        stG = y[5]
        dmSdt = dydt[0]
        dMGdt = dydt[1]
        dECSdt = dydt[2]
        dECGdt = dydt[3]
        dCCSdt = dydt[8]
        dCCGdt = dydt[9]
        pLacItot = self.p['pLactot']
        sigma0 = self.p['sigma0']
        T0 = self.p['T0']
        G0 = self.p['G0']
        tau = self.p['tau']
        gamma = self.p['gamma']
        kmGyr = self.p['kmGyr']
        dms = self.p['dms']
        dmg = self.p['dmg']
        PLS = self.p['PLmS']
        PLG = self.p['PLMG']
        pTetTot = self.p['pTettot']
        dCasgRNAdt = -self.p['dp']
        dsigmaSdt = -(dmSdt-dms*mS-dECSdt)*TLS/(2*h0*nfs)-(dECSdt-dCCSdt)*PLS/(2*h0*nfs)+kcas*dCasgRNAdt*CasgRNA*PLS/(2*h0*nfs)+self.m_gyrtopfunc(t,stS,pLacItot,sigma0,T0,G0,tau,gamma,TLS,pLen, kmGyr)
        dsigmaMGdt = -(dMGdt - dmg*MG-dECGdt)*TLG/(2*h0*nfg)-(dECGdt-dCCGdt)*PLG/(2*h0*nfg)+kcas*dCasgRNAdt*CasgRNA*PLG/(2*h0*nfg)+self.m_gyrtopfunc(t,stG,pTetTot,sigma0,T0,G0,tau,gamma,TLG,pLen, kmGyr)
        dydt[4]=dsigmaSdt
        dydt[5]=dsigmaMGdt
        dydt.append(dCasgRNAdt)
        return dydt
IPTGconcs = [0,100]
tint = [0, hours*60*60]
p1 = copy.copy(p)
p2 = copy.copy(p)
p2['IPTGtot'] = 100
#Sys1 =dCas12a_System(p1,ivs)
#Sys2 =dCas12a_System(p2, ivs)
#y0 = Sys1.initial_conditions()
#y0 = np.append(y0, [12])
#sol1 = solve_ivp(Sys1.convergent_TL_dynamics, tint, y0, vectorized=True, method='LSODA', atol=atol)
#sol2 = solve_ivp(Sys2.convergent_TL_dynamics, tint,
#                 y0, vectorized=True, method='LSODA',                 ato#l=atol)


