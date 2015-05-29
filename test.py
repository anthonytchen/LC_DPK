#conc_simu = runDPK.runDPKfunc( 10**calibDPK.res.x[0], 10**calibDPK.res.x[1], calibDPK.t_range )

conc_simu = runDPK.runDPKfunc( 10**-1, 10**np.log10(1.2e-13), calibDPK.t_range )

blood_conc_simu = conc_simu[0]
plt.plot(calibDPK.t_range/3600, blood_conc_simu*1e6, 'b--', calibDPK.t_range/3600, calibDPK.dat[:,1], 'ro')
plt.show()

print 1 - conc_simu[1][12] / conc_simu[1][0]
