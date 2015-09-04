def test_bounds(f_new, x_new, f_old, x_old):
    ''' Function used in basin hopping global optimiser to define the bounds '''
    bnds = ((-0.5, 1.5), (-10.5, -8.5))
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

#bnds = ((-0.3, 0.3), (-9.313, -8.711))
x0 = [0, -9.012] # these are values as if the vehicle is water
conc_simu = runDPKfunc( 10**x0[0]*0.52, 10**x0[1], 100*1e-6, t_range,  dose_factor, clear_blood=23.3e-6)

# optimising the depth of vehicle
#bnds = ((-0.3, 0), (50, 200))
#x0 = [0, 100];
#res_brute = brute(objFunCalib, bnds, args=(t_range, blood_conc_data, dose_factor), Ns = 50, finish = None, full_output=True, disp=True)
#res = minimize(objFunCalib, res_brute[0], args=(t_range, blood_conc_data, dose_factor), method='L-BFGS-B', bounds=bnds, options={'disp': True, 'maxiter': 100})
# Options for the optimiser: SLSQP, L-BFGS-B, TNC

#conc_simu = runDPKfunc( 10**res.x[0], 9.727e-10, res.x[1]*1e-6, t_range,  dose_factor)
#conc_simu = runDPKfunc( 10**x0[0]*0.52, 9.727e-10, x0[1]*1e-6, t_range,  dose_factor)



blood_conc_simu = conc_simu[0]
print 1 - conc_simu[1][12] / conc_simu[1][0]

plt.plot( t_range/3600, blood_conc_data, 'ro', t_range/3600, blood_conc_simu*1e6, 'b-')
mpl.rcParams.update({'font.size': 22})
plt.show()

