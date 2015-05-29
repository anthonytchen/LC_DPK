import numpy as np
from scipy.optimize import minimize
from scipy.optimize import brute
from scipy import optimize
from runDPK import objFunCalib, calConcProfile, runDPKfunc
import matplotlib as mpl
import matplotlib.pyplot as plt

dat = np.loadtxt("Nicotine_Bannon89_Single.txt")
t_range = dat[:,0]
blood_conc_data = dat[:,1]
t_range *= 3600 # converting from hour to seconds

bnds = ((-4, 2), (-18, -8))
x0 = [-1, np.log10(1.2e-13)]
res = minimize(objFunCalib, x0, args=(t_range, blood_conc_data), method='TNC', bounds=bnds, options={'disp': True, 'maxiter': 50})
#res = minimize(objFunCalib, x0, args=(t_range, blood_conc_data), method='CG', options={'disp': True, 'maxiter': 100})
#SLSQP, L-BFGS-B, TNC
conc_simu = runDPKfunc( 10**res.x[0], 10**res.x[1], t_range )
blood_conc_simu = conc_simu[0]
print 1 - conc_simu[1][12] / conc_simu[1][0]

#rranges = (slice(-10, 2, 1), slice(-20, -5, 1))
#resbrute = brute(objFunCalib, rranges, args=(t_range, blood_conc_data), full_output=True, finish=optimize.fmin)
#blood_conc_simu = runDPKfunc( 10**resbrute[0,0], 10**resbrute[0,1], t_range )

plt.plot( t_range/3600, blood_conc_data, 'ro', t_range/3600, blood_conc_simu*1e6, 'b-')
plt.show()
#bnds = [(-10, 10), (-20, -5)]
#res = differential_evolution(objFunCalib, bnds, strategy='best1bin', maxiter=10)
