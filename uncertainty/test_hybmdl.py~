#from hybmdl import * 
#from qspr_Kcc import *
#from qspr_Klp import *
import hybmdl

dat_Plp = np.matrix( np.loadtxt("Kow_Plp_lg10.txt") )
dat_Ksc = np.matrix( np.loadtxt("Kow_Ksc_lg10.txt") )
dat_Kcc = np.matrix( np.loadtxt("Kow_Kcc_lg10.txt") )
paras0 = np.log( np.array([4.2, 0.31, 0.69]) )
bnds = ((-10, 10), (-10, 10), (-10, 10))

sig2_y = np.array([0.05])
sig2_z = np.array([0.05, 0.05])

Xy = dat_Ksc[:,0]
Y = dat_Ksc[:,1]
Xz = (dat_Plp[:,0], dat_Kcc[:,0])
Z = (dat_Plp[:,1], dat_Kcc[:,1])

paras = np.empty_like (paras0)
np.copyto(paras, paras0)

rlt_em = hybmdl.EMmain(hybmdl.testFunc_top, hybmdl.testFunc_low, Xy, Y, Xz, Z, paras0, sig2_y, sig2_z, 100, 10, bnds)
rlt_pred = hybmdl.pred(hybmdl.testFunc_top, hybmdl.testFunc_low, Xy, rlt_em[0], rlt_em[1], rlt_em[2], rlt_em[3])
plt.plot(Xy[:,0], np.squeeze(rlt_pred[2]), 'x', Xy[:,0], Y, 'o', Xy[:,0], np.squeeze(rlt_pred[2])+1.58*np.squeeze(np.sqrt(rlt_pred[3])), '.', Xy[:,0], np.squeeze(rlt_pred[2])-1.58*np.squeeze(np.sqrt(rlt_pred[3])), '.')


def qspr_lg10K_cc_lip(paras, lg10Kow):
    ''' Calculate log10 partition in CC-water and LP-water, and return a Nx2 matrix
    '''

    paras_Kcc = paras[:2]
    lg10Kcc = qspr_lg10Kcc(paras_Kcc, lg10Kow)

    paras_Klp = paras[2]
    lg10Klip = qspr_lg10Klp(paras_Klp, lg10Kow)

    return np.matrix( np.append(lg10Kcc, lg10Klip) )


sig2_y = np.array([0.05])
sig2_z = np.array([0.05, 0.05])

Samples = Estep(qspr_lg10Ksc, qspr_lg10K_cc_lip, paras, np.matrix(dat_Ksc[:,0]), dat_Ksc[:,1], sig2_y, sig2_z, N=10)

obj_val = Mstep_theta_obj(paras, qspr_lg10Ksc, qspr_lg10K_cc_lip, np.matrix(dat_Ksc[:,0]), dat_Ksc[:,1], 
                          Samples, (dat_Kcc[:,0], dat_Plp[:,0]), (dat_Kcc[:,1], dat_Plp[:,1]), sig2_y, sig2_z)

paras1 = Mstep_main(qspr_lg10Ksc, qspr_lg10K_cc_lip, paras, np.matrix(dat_Ksc[:,0]), dat_Ksc[:,1], 
                    Samples, (dat_Kcc[:,0], dat_Plp[:,0]), (dat_Kcc[:,1], dat_Plp[:,1]), sig2_y, sig2_z)
