import numpy as np
import runDPK
import matplotlib as mpl
import matplotlib.pyplot as plt

x = np.array([-0.155, -9.012])  # 10**-0.155 = 0.7
#x = np.array([-0.097, -9.012]) # 10**-0.097 = 0.8
#x = np.array([-0.046, -9.012])  # 10**-0.046 = 0.9

t_range = np.linspace(0, 30, 181)
t_range_long = np.linspace(0, 60, 121)
#t_range = np.linspace(0, 30, 31)
conc_simu_15 = runDPK.runDPKfunc( 10**x[0], 10**x[1], 100*1e-6, t_range*3600, 1, clear_blood=20.8e-6 )
conc_simu_30 = runDPK.runDPKfunc( 10**x[0], 10**x[1], 100*1e-6, t_range*3600, 2, clear_blood=23.3e-6 )
conc_simu_60 = runDPK.runDPKfunc( 10**x[0], 10**x[1], 100*1e-6, t_range*3600, 4, clear_blood=25.7e-6 )

#conc_simu = runDPK.runDPKfunc( 10**-2.80, 10**-13.01202192, calibDPK.t_range )
#10**np.log10(5e-11), calibDPK.t_range )

dat = np.loadtxt("Nicotine_Bannon89_Single.txt")
dat_multi = np.loadtxt("Nicotine_Bannon89_Multiple.txt")

plt.figure(1, figsize=(18.5,10.5))
plt.plot(dat[:,0], dat[:,1], 'ro', label='Data (15mg)')
plt.plot(t_range, conc_simu_15[0]*1e6, 'r-', label='Model (15mg)')
plt.plot(dat[:,0], dat[:,2], 'bs', label='Data (30mg)')
plt.plot(dat_multi[0:12,0], dat_multi[0:12,1], 'bs')
plt.plot(t_range, conc_simu_30[0]*1e6, 'b-', label='Model (30mg)')
plt.plot(dat[:,0], dat[:,3], 'c^', label='Data (60mg)')
plt.plot(t_range, conc_simu_60[0]*1e6, 'c-', label='Model (60mg)')

plt.annotate('Patch removed', xy=(24, 0), xytext=(24, 2), arrowprops=dict(facecolor='none', shrink=0.05),)
plt.legend(loc='upper right')
#plt.axis([0, 32, 0, 16])
plt.axis([0, 32, 0, 30])
plt.xlabel('Time (hr)')
plt.ylabel('Plasma concentration (ng/ml)')
mpl.rcParams.update({'font.size': 20})

#plt.show()
plt.savefig('Nicotine_SingleDose.png', bbox_inches='tight')

#print 1 - conc_simu[1][12] / conc_simu[1][0]

conc_simu_30_long = runDPK.runDPKfunc( 10**x[0], 10**x[1], 100*1e-6, t_range_long*3600, 2, clear_blood=23.3e-6 )
total = conc_simu_30_long[2][1:]
#total[0] = 1 # to avoid division by zero
factor = (3.5*2*1e-4)/(40.075*0.01*1e-6)

fig = plt.figure(2, figsize=(18.5,10.5))
ax = fig.add_subplot(2,1,1)
ax.plot(t_range_long[1:], 100*conc_simu_30_long[3][1:]/total, 'r-', label='Stratum corneum', linewidth=3)
ax.plot(t_range_long[1:], 100*conc_simu_30_long[9][1:]/factor/total, 'm--', label='Cleared from blood', linewidth=3)
ax.annotate('Patch removed', xy=(24, 0), xytext=(24, 2), arrowprops=dict(facecolor='none', shrink=0.05),)
ax.axis([0.3, 50, 0, 105])
ax.legend(loc='best')
ax.set_xlabel('Time (hr)')
ax.set_ylabel('Percentage of applied dose (%)')

ax = fig.add_subplot(2,1,2)
ax.plot(t_range_long[1:], 100*conc_simu_30_long[6][1:]/total, 'b-', label='Viable epidermis', linewidth=3)
ax.plot(t_range_long[1:], 100*conc_simu_30_long[7][1:]/total, 'g-.', label='Dermis', linewidth=3)
ax.plot(t_range_long[1:], 100*conc_simu_30_long[8][1:]/factor/total, 'k--', label='Blood', linewidth=3)
ax.annotate('Patch removed', xy=(24, 0), xytext=(24, 2), arrowprops=dict(facecolor='none', shrink=0.05),)
#ax.set_yscale('log')
ax.legend(loc='upper right')
ax.axis([0.3, 50, 0, 1.5])
ax.set_xlabel('Time (hr)')
ax.set_ylabel('Percentage of applied dose (%)')
plt.savefig('Nicotine_Percent.png', bbox_inches='tight')
#plt.show()

fig = plt.figure(4, figsize=(18.5,10.5))
ax = fig.add_subplot(2,1,1)
ax.plot(t_range_long[1:], 1e9*conc_simu_30_long[4][1:], 'r-', label='Lipid', linewidth=3)
ax.plot(t_range_long[1:], 1e9*conc_simu_30_long[5][1:], 'b--', label='Corneocyte', linewidth=3)
ax.legend(loc='upper right')
ax.axis([0.3, 50, 0, 6])
ax.set_xlabel('Time (hr)')
ax.set_ylabel('Amount ($\mu$g)')

ax = fig.add_subplot(2,1,2)
ax.plot(t_range_long[1:], 100*conc_simu_30_long[4][1:]/conc_simu_30_long[3][1:], 'r-', label='Lipid', linewidth=3)
ax.plot(t_range_long[1:], 100*conc_simu_30_long[5][1:]/conc_simu_30_long[3][1:], 'b--', label='Corneocyte', linewidth=3)
ax.legend(loc='best')
ax.axis([0.3, 50, 0, 100])
ax.set_xlabel('Time (hr)')
#ax.set_yscale('log')
ax.set_ylabel('Percentage as in stratum corneum (%)')
mpl.rcParams.update({'font.size': 20})
plt.savefig('Nicotine_StraCorn.png', bbox_inches='tight')
#plt.show()

