class DPKInputParas:
    '''Class definition for the input parameters for DPK simulation'''

    # define input parameters and assign default values
    def __init__(self):       
        
        # skin geometric parameters
        self.n_layer_x_sc = 16
        self.n_layer_y_sc = 1
        self.x_len_viaepd = 70e-6 # depth of viable epidermis
        self.x_len_dermis = 1200e-6 # depth of dermis
        self.n_grids_x_ve = 10
        self.n_grids_x_de = 10
        self.offset_y_sc = 0
        self.b_inf_src = 0

        # solute parameters
        self.MW = []
        self.K_ow = []
        self.pKa = []
        self.AorB = [] # (A)cid or (B)ase
        self.frac_non_ion = [] # fraction of solute non-ionised at pH 7.4
        self.frac_unbound = [] # fraction of solute unbound in a 2.7% albumin solution at pH 7.4
        self.par_dermis2blood = [] # partition coefficient from dermis to blood
        self.k_clear_blood = [] # clearance rate in blood

        # vehicle parameters
        self.x_len_vehicle = [] # depth of vehicle
        self.diffu_vehicle = [] # diffusivity of solute in vehicle
        self.par_vehicle = [] # partition coefficient from vehicle to water
