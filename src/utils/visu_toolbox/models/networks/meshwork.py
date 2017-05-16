import numpy as np




class Meshwork:
    '''
    Documentation
    '''






    def compute_signal(pop_params, indiv_params, indiv_number, Kxd, invKd):

        #delta = compute_interpolation(Kxd, invKd, pop_params["delta"])
        #nu = compute_interpolation(Kxd, invKd, pop_params["nu"])
        delta = pop_params["delta"]
        nu = pop_params["nu"]

        thickness = pop_params["thickness"][0]


        time_points = np.linspace(60,100,20)
        space_shift = [0 for i in range(0, len(delta))]
        if(indiv_number > -1):
            time_points = np.exp(indiv_params["ksi"][0]) * (time_points - indiv_params["tau"][0])
            space_shift = indiv_params["w"]

        signals = []
        for t in time_points:
            signal = np.divide(space_shift, thickness*np.exp(delta)) + delta + np.divide(nu, thickness * t)
            signal = thickness * np.exp(signal)
            signals.append(signal)

        return signals
