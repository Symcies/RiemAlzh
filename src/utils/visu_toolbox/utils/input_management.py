path_id = "group.csv"
path_ages = "X.csv"
path_scores = "Y.csv"

def extract_observations(path_to_folder):
    id_file = open(path_to_folder + path_id, 'r')
    age_file = open(path_to_folder + path_ages, 'r')
    score_file = open(path_to_folder + path_scores, 'r')

    #Check if all files have correct path
    #Check if all files have the same length

    dict_patients = {}
    for line in id_file:
        id = int(line.split()[0])
        age = age_file.readline()
        score = score_file.readline()
        if id not in dict_patients:
            dict_patients[id] = PatientObservation(id)
        dict_patients[id].add_observation(float(age), score.split(','))

    return dict_patients


class PatientObservation:
    def __init__(self, id):
        self.id = id
        self.obs = {}

    def add_observation(self, age, observation):
        self.obs[age] = observation


def read_population_parameters(path_to_file):
    params = {}
    f = open(path_to_file, 'r')

    for line in f:
        line = line.lower()

        # Extract parameters
        list_elements = line.split()
        param_name = list_elements[0]
        param_values = [float(list_elements[i]) for i in range(1, len(list_elements))]

        # Add to existing dict
        params[param_name] = param_values

    f.close()

    return params


def read_individual_parameters(filename):
    f_ind = open(filename, 'r')
    ordered_params = [] #TODO make better

    # Creation of dictionary matching the name of a parameter and the number of its realisations
    number_param = {}
    name_params = f_ind.readline().lower().split()
    num_params = f_ind.readline().split()
    for i in range(len(name_params)):
        number_param[name_params[i]] = int(num_params[i])
        for j in range(int(num_params[i])):
            ordered_params.append(name_params[i] + str(j))

    indiv_param = {}

    # Initialize dict
    for param in number_param.keys():
        indiv_param[param] = [[] for _ in xrange(number_param[param])]

    for whole_line in f_ind:
        line = whole_line.split()
        for param in number_param.keys():
            for i in range(number_param[param]):
                indiv_param[param][i].append(float(line[ordered_params.index(param + str(i))]))


    return number_param, indiv_param
