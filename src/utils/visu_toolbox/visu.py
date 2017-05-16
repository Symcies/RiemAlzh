import numpy as np
import matplotlib
import matplotlib.pyplot as plt


matplotlib.use('tkagg')

from models.scalar.univariate     import Univariate
from models.scalar.multivariate   import Multivariate

from utils.input_management import read_individual_parameters, read_population_parameters, extract_observations

class Visu:
    '''
    Documentation
    '''

    def __init__(self, type, population_params, individual_params, observations):
        self.fig = plt.figure()
        self.splot = self.fig.add_subplot(111)

        self.model_type = type.lower()
        self.pop_params = read_population_parameters(population_params)
        self.number_param, self.indiv_params = read_individual_parameters(individual_params)
        self.observations = extract_observations(observations)
        self.model = self.init_model()


    def init_model(self):
        if self.model_type == "univariate":
            return Univariate(self.pop_params, self.indiv_params, self.observations, self.splot)
        elif self.model_type == "multivariate":
            return Multivariate()



    def hist(self, variable_name, nb_beans):
        plt.hist(self.pop_params[variable_name], nb_beans)
        plt.show()

    def plot_mean(self):
        return self.model.plot_mean()


    def plot_patients(self, list_of_patients_to_plot, with_obs):
        return self.model.plot_patients(list_of_patients_to_plot, with_obs)

    def hold_plot(self):
        plt.show()


    ''' UTILS '''
    def print_population_params(self):
        return [key for key, value in self.pop_params.iteritems()]
