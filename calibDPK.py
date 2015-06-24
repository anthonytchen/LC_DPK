def test_bounds(f_new, x_new, f_old, x_old):
    ''' Function used in basin hopping global optimiser to define the bounds '''
    bnds = ((-0.5, 1.5), (-10.5, -8.5))
    # bnds = ((-0.3, 1), (-15, -11))
    if x_new[0]>bnds[0][0] and x_new[0]<bnds[0][1] and x_new[1]>bnds[1][0] and x_new[1]<bnds[1][1]:
        return True
    else:
        return False

def test_bounds_thickness(f_new, x_new, f_old, x_old):
    ''' Function used in basin hopping global optimiser to define the bounds '''
    bnds = ((20, 500),)
    if x_new[0]>bnds[0][0] and x_new[0]<bnds[0][1]:
        return True
    else:
        return False



def print_fun_basinhopping(x, f, accepted):
    print("    at minima %.4e accepted %d" % (f, int(accepted)))

import numpy as np
from scipy.optimize import minimize, basinhopping, brute
#from scipy import optimize
from runDPK import objFunCalib, calConcProfile, runDPKfunc
import matplotlib as mpl
import matplotlib.pyplot as plt

dat = np.loadtxt("Nicotine_Bannon89_Single.txt")
t_range = dat[:,0]
blood_conc_data = dat[:,2]
dose_factor = 2
t_range *= 3600 # converting from hour to seconds

#bnds = ((-0.3, 1), (-15, -11))
bnds = ((-0.5, 1), (-10.5, -8.5))
#x0 = [0, -12]
x0 = [0, -9.5]
#conc_simu = runDPKfunc( 10**x0[0]*0.5, 10**x0[1]*0.5, 100*1e-6, t_range,  dose_factor)

# optimising the depth of vehicle
#bnds = ((20, 500),)
#x0 = 100;
#res = basinhopping(objFunCalib, x0, niter=1000, stepsize=10, minimizer_kwargs={"method":"SLSQP", "options":{'maxiter':10}, "bounds":bnds, "args":(t_range, blood_conc_data, dose_factor)}, disp=True, niter_success=100, callback=print_fun_basinhopping, accept_test=test_bounds_thickness)

#conc_simu = runDPKfunc( 1, 9.727e-10, res.x[0]*1e-6, t_range,  dose_factor)

#print np.array(bnds, float).shape[0]

# bnds = ((-0.5, 1), (-15, -10))
# x0 = [-0.49974794, -12.8762165] # best results with the above bounds

#res = minimize(objFunCalib, x0, args=(t_range, blood_conc_data, dose_factor), method='L-BFGS-B', bounds=bnds, options={'disp': True, 'maxiter': 50})
res_brute = brute(objFunCalib, bnds, args=(t_range, blood_conc_data, dose_factor), Ns = 20, finish = None, full_output=True, disp=True)
res = minimize(objFunCalib, res_brute[0], args=(t_range, blood_conc_data, dose_factor), method='L-BFGS-B', bounds=bnds, options={'disp': True, 'maxiter': 100})

#res = basinhopping(objFunCalib, x0, niter=1000, stepsize=0.1, minimizer_kwargs={"method":"SLSQP", "options":{'maxiter':10}, "bounds":bnds, "args":(t_range, blood_conc_data, dose_factor)}, disp=True, niter_success=100, callback=print_fun_basinhopping, accept_test=test_bounds)

#SLSQP, L-BFGS-B, TNC
#conc_simu = runDPKfunc( 10**res[0][0], 10**res[0][1], 100*1e-6, t_range,  dose_factor)
conc_simu = runDPKfunc( 10**res.x[0], 10**res.x[1], 100*1e-6, t_range,  dose_factor)

blood_conc_simu = conc_simu[0]
print 1 - conc_simu[1][12] / conc_simu[1][0]

plt.plot( t_range/3600, blood_conc_data, 'ro', t_range/3600, blood_conc_simu*1e6, 'b-')
#matplotlib.rcParams.update({'font.size': 22})

plt.show()
#bnds = [(-10, 10), (-20, -5)]
#res = differential_evolution(objFunCalib, bnds, strategy='best1bin', maxiter=10)
