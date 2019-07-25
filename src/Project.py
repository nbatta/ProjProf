import numpy as np
from configparser import SafeConfigParser
from scipy.integrate import quad
from scipy import special
from scipy.special import j0
from scipy.interpolate import interp1d

from src.io import dict_from_section

class Profiles(object):
    def __init__(self,iniFile="input/params.ini",fwhm=1.4):

        Config = SafeConfigParser()
        Config.optionxform=str
        Config.read(iniFile)
        self.cc = dict_from_section(Config,'constants')
        self.pp = dict_from_section(Config,'params')

        self.TCMBmuK = self.cc['TCMB'] * 1e6

        self.XH = 0.76 ### FIX THIS INTO CONSTANTS ###
        self.NNR = 100
        self.disc_fac = np.sqrt(2)
        self.l0 = 30000. ### ell max
        self.fwhm = fwhm ### ARCMINS

        ### extrapolation variables
        self.inter_max = 5.
        self.inter_fine_fac = 40. #4.

        self.h = self.pp['H0']/100.
        self.om = (self.pp['ombh2']+self.pp['omch2'])/self.h**2

    #Hubble function
    def hub_func(self,z):
        Om = self.om
        Ol = 1 - Om
        #ans = np.sqrt(Om*(1.0 + z)**3 + Ol + (1 - O_tot)*(1 + z)**2)
        #FLAT LCDM
        ans = np.sqrt(Om*(1.0 + z)**3 + Ol)
        return ans

    #Comoving Distance Integrand
    def ComInt(self,z):
        ans = 1.0/self.hub_func(z)
        return ans

    #Comoving Distance
    def ComDist(self,z):
        Om = self.om 
        Ol = 1. - Om 
        O_tot = Om + Ol

        Dh = self.cc['C_OVER_HUBBLE']/self.h
        ans = Dh*quad(self.ComInt,0,z)[0]
        if (O_tot < 1.0): ans = Dh / np.sqrt(1.0-O_tot) *  np.sin(np.sqrt(1.0-O_tot) * quad(ComInt,0,z)[0])
        if (O_tot > 1.0): ans = Dh / np.sqrt(O_tot-1.0) *  np.sinh(np.sqrt(O_tot-1.0) * quad(ComInt,0,z)[0]) 
        return ans

    #Angular Distance
    def AngDist(self,z):
        ans = self.ComDist(z)/(1.0+z)
        return ans

    def project_prof_beam_interpol(self,tht,r200c,z,rho_sim,Pth_sim,test=False): #r200c in Mpc

        fwhm = self.fwhm

        drint = 1e-3 * (self.cc['MPC2CM'])
        AngDist = self.AngDist(z) * self.h
        disc_fac = self.disc_fac
        l0 = self.l0 
        NNR = self.NNR 
        NNR2 = 4*NNR
        
        fwhm *= np.pi / (180.*60.) #convert from arcmins to rads
        sigmaBeam = fwhm / np.sqrt(8.*np.log(2.))
        
        # JCH: convert r200c to Mpc/h units to match Nick's AngDist above
        r200c_Mpc_over_h = r200c * self.h

        sig = 0
        sig2 = 0
        sig_p = 0
        sig2_p = 0
        area_fac = 0
        
        r_ext = AngDist*np.arctan(np.radians(tht/60.))
        r_ext2 = AngDist*np.arctan(np.radians(tht*disc_fac/60.))
        
        rad = (np.arange(1e4) + 1.0)/1e3 #in MPC/h
        rad2 = (np.arange(1e4) + 1.0)/1e3 #in MPC/h
        
        radlim = r_ext
        radlim2 = r_ext2
        
        dtht = np.arctan(radlim/AngDist)/NNR # rads
        dtht2 = np.arctan(radlim2/AngDist)/NNR # rads
        
        thta = (np.arange(NNR) + 1.)*dtht
        thta2 = (np.arange(NNR) + 1.)*dtht2
        
        thta_smooth = (np.arange(NNR2) + 1.)*dtht
        thta2_smooth = (np.arange(NNR2) + 1.)*dtht2
        
        rho2D = thta_smooth * 0.0
        rho2D2 = thta_smooth * 0.0
        Pth2D = thta_smooth * 0.0
        Pth2D2 = thta_smooth * 0.0
        
        rho2D_beam = thta * 0.0
        rho2D2_beam = thta* 0.0
        Pth2D_beam = thta* 0.0
        Pth2D2_beam = thta* 0.0

        for kk in range(NNR2):
            rint  = np.sqrt((rad)**2  + thta_smooth[kk]**2 *AngDist**2)
            rint2 = np.sqrt((rad2)**2 + thta2_smooth[kk]**2*AngDist**2)

            if (test):
                rho2D[kk]  = np.sum(2.*self.rho_sim_test(theta_sim_rho,rint /r200c_Mpc_over_h)*drint)
                rho2D2[kk] = np.sum(2.*self.rho_sim_test(theta_sim_rho,rint2/r200c_Mpc_over_h)*drint)
                
                Pth2D[kk]  = np.sum(2.*self.Pth_sim_test(theta_sim_pth,rint /r200c_Mpc_over_h)*drint)
                Pth2D2[kk] = np.sum(2.*self.Pth_sim_test(theta_sim_pth,rint2/r200c_Mpc_over_h)*drint)
            else:
                rho2D[kk]  = np.sum(2.*rho_sim(rint /r200c_Mpc_over_h)*drint)
                rho2D2[kk] = np.sum(2.*rho_sim(rint2/r200c_Mpc_over_h)*drint)
                
                Pth2D[kk]  = np.sum(2.*Pth_sim(rint /r200c_Mpc_over_h)*drint)
                Pth2D2[kk] = np.sum(2.*Pth_sim(rint2/r200c_Mpc_over_h)*drint)
        
        for kk in range(NNR):

            special1 = special.iv(0, thta_smooth *thta[kk] / sigmaBeam**2)
            special2 = special.iv(0, thta2_smooth*thta2[kk]/ sigmaBeam**2)
            special1[np.where(special1 == np.inf)] = 1.e308
            special2[np.where(special2 == np.inf)] = 1.e308

            rho2D_beam[kk]  = np.sum(thta_smooth  * rho2D  * np.exp(-0.5*thta_smooth**2 /sigmaBeam**2)  
                                     * special1)*dtht
            rho2D2_beam[kk] = np.sum(thta2_smooth * rho2D2 * np.exp(-0.5*thta2_smooth**2/sigmaBeam**2)
                                     * special2)*dtht2
            Pth2D_beam[kk]  = np.sum(thta_smooth  * Pth2D  * np.exp(-0.5*thta_smooth**2 /sigmaBeam**2)  
                                     * special1)*dtht
            Pth2D2_beam[kk] = np.sum(thta2_smooth * Pth2D2 * np.exp(-0.5*thta2_smooth**2/sigmaBeam**2) 
                                     * special2)*dtht2

            area_fac += 2.0*np.pi*dtht*thta[kk]
        
            rho2D_beam[kk]  *= np.exp(-0.5*thta[kk]**2 /sigmaBeam**2) / sigmaBeam**2
            rho2D2_beam[kk] *= np.exp(-0.5*thta2[kk]**2/sigmaBeam**2) / sigmaBeam**2
            Pth2D_beam[kk]  *= np.exp(-0.5*thta[kk]**2 /sigmaBeam**2) / sigmaBeam**2
            Pth2D2_beam[kk] *= np.exp(-0.5*thta2[kk]**2/sigmaBeam**2) / sigmaBeam**2

        sig  = 2.0*np.pi*dtht *np.sum(thta *rho2D_beam) 
        sig2 = 2.0*np.pi*dtht2*np.sum(thta2*rho2D2_beam) 

        sig_all_beam = (2*sig - sig2) * 1e-3 * self.cc['SIGMA_T'] * self.TCMBmuK / self.cc['MP'] / (np.pi * np.radians(tht/60.)**2)  * ((2. + 2.*self.XH)/(3.+5.*self.XH)) 
        #sig_all_beam = (2*sig - sig2) * 1e-3 * self.cc.c['SIGMA_T'] / self.cc.c['ME'] / (np.pi * np.radians(tht/60.)**2) * ((2. + 2.*self.XH)/(3.+5.*self.XH)) 

        sig_p  = 2.0*np.pi*dtht*np.sum(thta*Pth2D_beam)
        sig2_p = 2.0*np.pi*dtht2*np.sum(thta2*Pth2D2_beam)

        # JCH: get rid of TCMB here so that results are in Compton-y units
        sig_all_p_beam = (2*sig_p - sig2_p) * self.cc['SIGMA_T']/(self.cc['ME']*self.cc['C']**2) / area_fac #* \
        #self.cc.c['TCMBmuK'] # muK #* ((2. + 2.*self.XH)/(3.+5.*self.XH))# muK

        sig_all_beam /= (self.h)
        sig_all_p_beam /= (self.h)

        return sig_all_beam, sig_all_p_beam #, sig_all_beam_cumul, sig_all_p_beam_cumul #JCH edit -- output cumulative profile as well

    def make_a_obs_profile_sim(self,thta_arc,r200c,z,rho_int,pres_int):
        rho = np.zeros(len(thta_arc))
        pth = np.zeros(len(thta_arc))
        for ii in range(len(thta_arc)):
            temp = self.project_prof_beam_interpol(thta_arc[ii],r200c,z,rho_int,pres_int)
            rho[ii] = temp[0]
            pth[ii] = temp[1]
        return rho,pth

    def interpol_sim_profile(self,x,prof):
        #Including extrapolation
        if (np.max(x) < self.inter_max):
            fine_fac = self.inter_fine_fac
            xtr_inds = np.ceil(self.inter_max - np.max(x))
            xtr_inds = np.floor(self.inter_max - np.max(x))*fine_fac
            str_xtr = np.ceil(np.max(x)) + 1./fine_fac
            xtr = np.arange(xtr_inds)/fine_fac + str_xtr
            extend = np.poly1d(np.polyfit(x, np.log(prof), 3))

            ans = interp1d(np.append(x,xtr),np.append(prof,np.exp(extend(xtr))),kind='slinear',bounds_error=False,fill_value=0)
        else: 
            ans = interp1d(x,prof,kind='slinear',bounds_error=False,fill_value=0)
        return ans

    #### JUST FOR TESTING ####
    def rho_sim_test(self, x):
        theta = np.array([3.6337402156859753, 1.0369351928324118, 3.3290812595973063])
        ###Battaglia 2016 profile parameters w/o M & z dependences
        fb = 0.0490086879038/0.314992203163
        rhoc = 2.77525e2
        Msol_cgs = 1.989e33
        kpc_cgs = 3.086e21
        a1,a2,a3 = theta
        gamma = 0.2
        ans = 10**a1 / ((x/0.5)**gamma * (1 + (x/0.5)**a2)**((a3 - gamma)/a2))
        ans *=  rhoc * Msol_cgs / kpc_cgs / kpc_cgs / kpc_cgs * fb 
        return ans
    
    def pth_sim_test(self,x):
        theta = np.array([18.1, 0.5, 4.35])
        ###Battaglia et al 2012b profile paramters w/o M & z dependences
        P0,xc,bt = theta
        al = 1.0
        gm = 0.3
        ans = P0 / ((x*xc)**gm * (1 + (x*xc)**al)**((bt-gm)/al))
        return ans
    #### JUST FOR TESTING ####
