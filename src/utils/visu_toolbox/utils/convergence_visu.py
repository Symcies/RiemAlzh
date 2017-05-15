import matplotlib.pyplot as plt

class convergence_visu:
    '''
    Documentation
    '''

    def __init__(self, file_name):
        dict_convergence = self.read_convergence_output(file_name)
        self.iterations = dict_convergence["iter"]
        del dict_convergence["iter"]
        self.convergence_param = dict_convergence

    def name_variables(self):
        print [k for k, v in self.convergence_param.iteritems()]

    def plot(self, list_variables):

        f, subplots = plt.subplots(int(len(list_variables) + 1)/2,
                                   2,
                                   squeeze = False,
                                   figsize=(18, 4*len(list_variables)))

        for idx, var in enumerate(list_variables):
            self.display(var, subplots[int(idx/2)][idx%2])

        plt.show()

    def display(self, name, subplot):
        subplot.set_title(name, fontsize=18)

        # Todo : Assert if the name is not in the dict
        values = self.convergence_param[name]
        subplot.plot(self.iterations, values, linewidth=7)

    def read_convergence_output(self, file_name):
        f = open(file_name, 'r')
        labels = f.readline().split()

        values = [[] for _ in xrange(len(labels))]
        for line in f:
            line_vec = line.split()
            for idx, item in enumerate(line_vec):
                values[idx].append(item)
        f.close()

        dict_convergence = {}
        for idx, item in enumerate(labels):
            dict_convergence[item] = values[idx]

        return dict_convergence
