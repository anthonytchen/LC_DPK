import numpy as np

x = np.array([0.0, -12.5])

#t_range = np.array(range(0,30))
t_range = np.linspace(0, 30, 61)
conc_simu_15 = runDPK.runDPKfunc( 10**x[0], 10**x[1], t_range*3600, 1 )
conc_simu_30 = runDPK.runDPKfunc( 10**x[0], 10**x[1], t_range*3600, 2 )
conc_simu_60 = runDPK.runDPKfunc( 10**x[0], 10**x[1], t_range*3600, 4 )

#conc_simu = runDPK.runDPKfunc( 10**-2.80, 10**-13.01202192, calibDPK.t_range )
#10**np.log10(5e-11), calibDPK.t_range )

dat = np.loadtxt("Nicotine_Bannon89_Single.txt")

plt.plot(dat[:,0], dat[:,1], 'ro', label='Data (15mg)')
plt.plot(t_range, conc_simu_15[0]*1e6, 'r-', label='Model (15mg)')
plt.plot(dat[:,0], dat[:,2], 'bs', label='Data (30mg)')
plt.plot(t_range, conc_simu_30[0]*1e6, 'b-', label='Model (30mg)')
plt.plot(dat[:,0], dat[:,3], 'c^', label='Data (60mg)')
plt.plot(t_range, conc_simu_60[0]*1e6, 'c-', label='Model (60mg)')

plt.annotate('Patch removed', xy=(24, 0), xytext=(24, 2), arrowprops=dict(facecolor='none', shrink=0.05),)
plt.legend(loc='upper right')
#plt.axis([0, 32, 0, 16])
plt.axis([0, 32, 0, 30])
plt.xlabel('Time (hr)')
plt.ylabel('Plasma concentration (ng/ml)')
matplotlib.rcParams.update({'font.size': 18})

plt.show()

#print 1 - conc_simu[1][12] / conc_simu[1][0]

dat = np.loadtxt("Nicotine_Bannon89_Multiple.txt")
t_range_multiple = dat[:,0]
blood_conc_data_multiple = dat[:,1]
t_range_multiple *= 3600 # converting from hour to seconds
#conc_simu_30_multiple = runDPK.runDPKfunc( 10**calibDPK.res.x[0], 10**calibDPK.res.x[1], 3600*range(169), 2 )
conc_simu_30_multiple = runDPK.runDPKfunc( 10**-3.20627986, 10**-12.99065003, 3600*np.array(range(169)), 2 )
plt.plot(t_range_multiple/3600, dat[:,1], 'bs', label='Data (30mg)')
plt.plot(range(169), conc_simu_30_multiple[0]*1e6, 'b-', label='Model (30mg)')
plt.legend(loc='upper right')
#plt.axis([0, 32, 0, 16])
plt.axis([0, 170, 0, 20])
plt.xlabel('Time (hr)')
plt.ylabel('Plasma concentration (ng/ml)')
matplotlib.rcParams.update({'font.size': 18})
plt.show()

