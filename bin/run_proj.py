from src.Project import Profiles
import numpy as np

P = Profiles()
### test angular distance
print (P.AngDist(0.5))

### set up test profiles
xx = (np.arange(100) + 1) * 0.05

rho = P.rho_sim_test(xx)
pth = P.pth_sim_test(xx)

intr = P.interpol_sim_profile(xx,rho)
intp = P.interpol_sim_profile(xx,pth)

### set up for projection
theta = (np.arange(10) + 1) * 0.5
r200 = 0.5 ## Mpc/h
z = 0.5

projr,projp = P.make_a_obs_profile_sim(theta,r200,z,intr,intp)

print (projr,projp)
