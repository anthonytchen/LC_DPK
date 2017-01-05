
def plot_stracorn():

    import numpy as np
    import matplotlib.pyplot as plt

    nx_layers = 12

    nx_grids_lipid = 2; # of x-grids for lipid layer, 2
    nx_grids_cc = 4; # of x-grids for corneocyte layer, 4
    ny_grids_lipid = 2; # of y-grids for lipid layer, 2
    ny_grids_cc_dn = 2; # of y-grids for dn-part of the offset corneocyte layer, 2

    #plt.subplots(figsize=(8,8))

    # Plot the lines that mark the boundary of corneocytes
    for i in range(nx_layers): # starting from 0
        
        if np.remainder(i, 2) == 0: # if i is even number
            x_st = ny_grids_lipid - 1 - 0.5
            x_ed = ny_grids_lipid - 1 + ny_grids_cc_dn + ny_grids_lipid + ny_grids_cc_dn*8 - 0.5
            y_st = i*(nx_grids_lipid+nx_grids_cc) + nx_grids_lipid - 0.5
            y_ed = y_st + nx_grids_cc
            plt.plot([x_st, x_ed], [y_st, y_st], color='k', linestyle='-', linewidth=1)
            plt.plot([x_st, x_ed], [y_ed, y_ed], color='k', linestyle='-', linewidth=1)
            plt.plot([x_st, x_st], [y_st, y_ed], color='k', linestyle='-', linewidth=1)
            plt.plot([x_ed, x_ed], [y_st, y_ed], color='k', linestyle='-', linewidth=1)

        else:
            x_ed1 = ny_grids_lipid - 1 + ny_grids_cc_dn - 0.5
            x_st2 = ny_grids_lipid - 1 + ny_grids_cc_dn + ny_grids_lipid - 0.5
            x_ed2 = x_st2 + ny_grids_cc_dn*8 + ny_grids_lipid - 1
            y_st = i*(nx_grids_lipid+nx_grids_cc) + nx_grids_lipid - 0.5
            y_ed = y_st + nx_grids_cc
            plt.plot([-0.5, x_ed1], [y_st, y_st], color='k', linestyle='-', linewidth=1)
            plt.plot([x_st2, x_ed2], [y_st, y_st], color='k', linestyle='-', linewidth=1)
            plt.plot([-0.5, x_ed1], [y_ed, y_ed], color='k', linestyle='-', linewidth=1)
            plt.plot([x_st2, x_ed2], [y_ed, y_ed], color='k', linestyle='-', linewidth=1)
            
            plt.plot([x_ed1, x_ed1], [y_st, y_ed], color='k', linestyle='-', linewidth=1)
            plt.plot([x_st2, x_st2], [y_st, y_ed], color='k', linestyle='-', linewidth=1)

    #plt.imshow(dat, interpolation='none', aspect='auto')
    #cb = plt.colorbar(shrink=0.98, pad=0.05)
    #cb.set_clim(0, 2900)
    #plt.tight_layout()
    #plt.axis('off')

    #plt.show()

def calConcProfile(conc_obj, skin_obj, conc_raw, coord_x, coord_y):
    ''' The function to work out concentration at different depth    '''

    import numpy as np

    # calulcate concentration in vehicle
    idx = 0
    conc_obj.vehicle = np.matrix( [coord_x[ idx : idx+skin_obj.m_dim_vh ], conc_raw[ idx : idx+skin_obj.m_dim_vh ]] )
    idx += skin_obj.m_dim_vh


    # calculate concentration in stratum corneum
    n_layer_x_sc = skin_obj.getNLayerXSc(); # layer of cells in x direction in SC
    a1 = np.zeros(n_layer_x_sc);
    skin_obj.get1DCoordSC(a1);
    a2 = np.zeros(n_layer_x_sc);
    skin_obj.get1DConcSC(a2);
    conc_obj.stracorn1d = np.matrix( [a1, a2] );

    conc_obj.stracorn2d_coordx = np.array( coord_x[ idx : idx+skin_obj.m_dim_sc ] )
    conc_obj.stracorn2d_coordy = np.array( coord_y[ idx : idx+skin_obj.m_dim_sc ] )
    n_x_sc = skin_obj.getNGridsXSc()
    n_y_sc = skin_obj.getNGridsYSc()
    conc_obj.stracorn2d_data = np.zeros( (n_x_sc, n_y_sc) );
    for i in range(n_x_sc):
        st = idx + i*n_y_sc;
        ed = idx + (i+1)*n_y_sc;
        conc_obj.stracorn2d_data[i,:] = conc_raw[st:ed];

    idx += skin_obj.m_dim_sc;

    # calculate concentration in viable epidermis
    coord = coord_x[ idx : idx+skin_obj.m_dim_ve ] + coord_x[ idx+1 : idx+skin_obj.m_dim_ve+1 ]
    coord /= 2; # use the middle point as the coordinate
    conc_obj.viaepd = np.matrix( [coord, conc_raw[ idx : idx+skin_obj.m_dim_ve ]] )
    idx += skin_obj.m_dim_ve;

    # calculate concentration in dermis
    coord = coord_x[ idx : idx+skin_obj.m_dim_de ] + coord_x[ idx+1 : idx+skin_obj.m_dim_de+1 ]
    coord /= 2; # use the middle point as the coordinate
    conc_obj.dermis = np.matrix( [coord, conc_raw[ idx : idx+skin_obj.m_dim_de ]] )
    idx += skin_obj.m_dim_de;

    # calculate concentration in blood
    conc_obj.blood = np.matrix( [coord_x[ idx : idx+skin_obj.m_dim_bd ], conc_raw[ idx : idx+skin_obj.m_dim_bd ]] )

    # calculate amount of solute in differnet layers
    
    n_amount = 9
    dat = np.zeros(n_amount);
    skin_obj.getLayersAmount(dat);
    
    conc_obj.mass_applied = dat[0]-dat[1]
    conc_obj.mass_stracorn_total = dat[2]
    conc_obj.mass_stracorn_lipid = dat[3]
    conc_obj.mass_stracorn_corneocyte = dat[4]
    conc_obj.mass_viaepd = dat[5]
    conc_obj.mass_dermis = dat[6]
    conc_obj.mass_blood = dat[7]
    conc_obj.mass_cleared = dat[8]

def calcAUC(t, conc):
    '''Calculate AUC, t_max, C_max from time-profile of concentration values
    '''
    import numpy as np

    auc = .0
    for i in range(len(t)-1):
        auc += (conc[i] + conc[i+1])/2 * (t[i+1]-t[i])
    
    t_max = np.amax(t)
    t_interp = np.linspace(0, t_max, 1000)
    conc_interp = np.interp(t_interp, t, conc)

    t_max = t_interp[np.argmax(conc_interp)]
    C_max = np.amax(conc_interp)

    return (auc, t_max, C_max)

def runS3VDBfunc(t_range, t, fnchem):
    '''Read data and visualise
    '''
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    import matplotlib.patches as patches

    nt = len(t_range)
    #t = 50

    nx_sursb = 1
    ny_sursb = 20

    nx_sb_har = 10
    ny_sb_har = 1

    nx_layers = 12
    nx_sc = (4+2)*nx_layers+2 # 74
    ny_sc = 22
    ny_cc_cells = 1

    # load data
    fn_str = fnchem + "conc_sursb_chem0_comp0.txt"
    with open(fn_str) as f :
        dat_sursb = np.array(f.read().split(), dtype=float).reshape(nt, nx_sursb*ny_sursb+1)

    conc_max = np.max(dat_sursb)
    
    fn_str = fnchem + "conc_sursb_chem0_comp1.txt"
    dat_sursb_har = np.loadtxt(fn_str).reshape(nt, nx_sursb*ny_sb_har)

    conc1 = np.max(dat_sursb_har)
    if (conc1 > conc_max):
        conc_max = conc1[:]

    fn_str = fnchem + "conc_sc_chem0_comp0.txt"
    dat_sc = np.loadtxt(fn_str).reshape(nt, nx_sc*ny_sc*ny_cc_cells);

    conc1 = np.max(dat_sc)
    if (conc1 > conc_max):
        conc_max = conc1[:]
        
    fn_str = fnchem + "conc_sb_chem0_comp0.txt"
    dat_sb_har = np.loadtxt(fn_str).reshape(nt, nx_sb_har*ny_sb_har)

    conc1 = np.max(dat_sb_har)
    if (conc1 > conc_max):
        conc_max = conc1[:]
        
    # plot time course in dat_sursb_har
    fig = plt.figure(figsize=(8,6))
    plt.plot(np.array(t_range)/3600.0, dat_sursb_har*1e3)
    plt.xlabel('Time (hr)')
    plt.ylabel('Concentration ($\mu$g/mL)')
    plt.savefig("conc_hf.png", bbox_inches='tight', dpi=300)
    plt.close()
    
    # plot
    fig = plt.figure(figsize=(8,6))
    gs = gridspec.GridSpec(3, 2, width_ratios=[5, 0.4], height_ratios=[1, 1, 5]) 
    # plt.subplots(figsize=(8,8))

    solid_mass_ini = (0.5e-6*0.5e-6*0.01)*1.782*1e3
    #solid_mass_ini = (0.01e-6*10e-6*0.01)*1.782*1e3
    solid_mass_now = dat_sursb[t,0]
    ratio_mass = solid_mass_now / solid_mass_ini * 100
    ratio_len = np.sqrt(ratio_mass)
    tx = "Time = " + str(t_range[t]/3600.0) + \
         " hr; Concentration in $\mu$g/mL"#; \n Solid mass " + "{:.1f}".format(ratio_mass) + "% of the initial value"

    
    ax = plt.subplot(gs[0,0:1])
    if 0:
        ax.add_patch(
            patches.Rectangle(
                (0, 0),   # (x,y)
                0.1,          # width
                0.7,          # height
                facecolor="none"
            )
        )
        ax.add_patch(
            patches.Rectangle(
                (0, 0),   # (x,y)
                0.1*ratio_len,          # width
                0.7*ratio_len,          # height
                facecolor=None
            )
        )

    plt.text(0.2, 0.3, tx)
    plt.axis('off')

    conc_max *= 1e3
    
    ax = plt.subplot(gs[1,0])
    dat = dat_sursb[t,1:] *1e3
    dat = np.vstack( (dat, dat) )
    plt.imshow(dat, vmin=0, vmax=conc_max, interpolation='none', aspect='auto')
    #cb = plt.colorbar(shrink=0.98, pad=0.05)
    #cb.set_clim(0, 0.03)
    plt.tight_layout()
    plt.axis('off')

    ax = plt.subplot(gs[1,1])
    dat = dat_sursb_har[t,:].reshape(nx_sursb, ny_sb_har) *1e3
    #print dat
    dat = np.hstack( (dat, dat) )
    plt.imshow(dat,  vmin=0, vmax=conc_max, interpolation='none', aspect='auto')

    #cb.set_clim(0, 0.03)
    plt.tight_layout()
    plt.axis('off')


    ax = plt.subplot(gs[2,0])
    dat = dat_sc[t,:].reshape(nx_sc, ny_sc*ny_cc_cells) *1e3
    # plot_stracorn()
    plt.imshow(dat,  vmin=0, vmax=conc_max, interpolation='none', aspect='auto')
    #cb = plt.colorbar(shrink=0.98, pad=0.05)
    #cb.set_clim(0, 0.03)
    plt.tight_layout()
    plt.axis('off')

    ax = plt.subplot(gs[2,1])
    dat = dat_sb_har[t,:].reshape(nx_sb_har, ny_sb_har) *1e3
    dat = np.hstack( (dat, dat) )
    im = plt.imshow(dat,  vmin=0, vmax=conc_max, interpolation='none', aspect='auto')
    #cb = plt.colorbar(shrink=0.98, pad=0.05)
    plt.tight_layout()
    plt.axis('off')

    cax = fig.add_axes([1.05, 0.05, 0.03, 0.8])
    fig.colorbar(im, cax=cax)

#    fnsav = 
    plt.savefig("2D_" + str(t_range[t]/3600.0) + "hr.png", bbox_inches='tight', dpi=300)
    plt.close()

    return (dat_sursb, dat_sursb_har, dat_sc, dat_sb_har)



