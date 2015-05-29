def objFunCalib(paras, arg_list):
    ''' Function to calculate the objective function to be minimised for use with model calibration'''

    dx_vehicle = paras[0]
    diffu_vehicle = paras[1]

    # todo: the following data input should be independent of this function!
    t_range =         [0, 1, 2, 3, 4, 6, 8, 12, 16, 24, 27, 30]
    blood_conc_data = [0, 1, 2, 3, 5, 5, 5, 5, 5, 5, 2, 1] # these are just toy numbers for testing
    t_range *= 3600 # converting from hour to seconds

    blood_conc = runDKPfunc( dx_vehicle, diffu_vehicle, t_range )
    blood_conc *= 1e6 # converting from kg/m^3 to ng/ml
    err = np.mean( np.square( blood_conc_data - blood_conc ) )
    return err
