from src.Project import Profiles
import numpy as np
import src.fisher_plots as fp

datahome = "../data/"

NGAL = np.array([38.e4,126.e4,333.e4,570.e4,442.e4,13.e4])
thta_arc = np.linspace(0.7, 5., 6.) *  1.5 / 1.4

P = Profiles()

xx = (np.arange(100) + 1) * 0.05

rho = P.rho_sim_test(xx)
pth = P.pth_sim_test(xx)

intr = P.interpol_sim_profile(xx,rho)
intp = P.interpol_sim_profile(xx,pth)

r200 = 0.5 ## Mpc/h                                                                                                           
z = 0.5

projr,projp = P.make_a_obs_profile_sim(theta,r200,z,intr,intp)

"""
This wont work, this is a template code.
You need to calculated the various derivates for density and pressure profiles below instead of what coded.
"""

projr_p1up = projr #param 1 density up
projr_p1dn = projr #param 1 density down
projp_p1up = projp #param 1 pressure up
projp_p1dn = projp #param 1 pressure down

projr_p2up = projr #param 2 density up
projr_p2dn = projr #param 2 density down
projp_p2up = projp #param 2 pressure up
projp_p2dn = projp #param 2 pressure down

projr_p3up = projr #param 3 density up
projr_p3dn = projr #param 3 density down
projp_p3up = projp #param 3 pressure up
projp_p3dn = projp #param 3 pressure down

p_intial = np.array([1,2,3])

p1diff = p1up - p1dn
p2diff = p2up - p2dn
p3diff = p3up - p3dn

#Now taking derivates
mu_rp1 = (projr_p1up - projr_p1up)/p1diff
mu_rp2 = (projr_p2up - projr_p2up)/p2diff
mu_rp3 = (projr_p3up - projr_p3up)/p3diff

mu_pp1 = (projp_p1up - projp_p1up)/p1diff
mu_pp2 = (projp_p2up - projp_p2up)/p2diff
mu_pp3 = (projp_p3up - projp_p3up)/p3diff

covrho = np.genfromtxt(datahome+"CovMat_CMB_S4-1.5arcmin_beam1.4_Lmax30000.0_2017-03-16.txt")
covpth = np.genfromtxt(datahome+"CovMat_y_S4-1.5arcmin_beam1.4_Lmax30000.0_2017-03-16.txt")

ii = 2 # index 2 in NGALs corresponds to redshift 

covrho /= NGAL[ii]
covpth /= NGAL[ii]

rho_mus = np.transpose([mu_rp1,mu_rp2,mu_rp3])
pth_mus = np.transpose([mu_pp1,mu_pp2,mu_pp3])

tot_fish_rho = np.dot(np.transpose(rho_mus),np.dot(np.linalg.inv(covrho),rho_mus))
tot_fish_pth = np.dot(np.transpose(pth_mus),np.dot(np.linalg.inv(covpth),pth_mus))

tot_fish = tot_fish_rho + tot_fish_pth

fish_inv = np.linalg.inv(tot_fish)

fp.marginal_pars(p_initial,fish_inv)
fp.tri_plot(p_initial,tot_fish,'test',saveFile=home+'test_fisher.png')
