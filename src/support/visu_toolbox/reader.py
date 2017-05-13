def read_population_parameters(file_name):
    params = {}
    f = open(file_name, 'r')

    for line in f:
        line = line.lower()

        # Extract parameters
        list_elements = line.split()
        param_name = list_elements[0]
        param_values = []
        for i in range(1, len(list_elements)):
            try:
                param_values.append(float(list_elements[i]))
            except:
                pass

        # Add to existing dict
        params[param_name] = param_values

    return params



def add_individual_parameters(line, number_of_params):
    individual_parameters = {}

    list_elements = line.split(' ')
    idx_start = 0

    for tpl in number_of_params:
        idx_end = idx_start + tpl[1]
        individual_parameters[tpl[0].rstrip()] = [float(list_elements[i].rstrip()) for i in range(idx_start, idx_end)]
        idx_start = idx_end

    return individual_parameters

def read_individual_parameters(file_name):
    params = []
    labels = []
    number_per_label = []
    individuals = []

    f = open(file_name, 'r')

    line_count = -1

    while True:
        line_count += 1
        # Read line, continue if empty and break if last one
        line = f.readline()
        if not line:
            break
        if line == "":
            continue

        # Read line, continue if empty, break if last one
        if (line_count == 0):
            labels = line.split(' ')
            continue

        if (line_count == 1):
            number_per_label = line.split(' ')

            for idx, name in enumerate(labels):
                params.append((name, int(number_per_label[idx])))

            continue

        # Add individuals
        individuals.append(add_individual_parameters(line, params))

    return individuals
