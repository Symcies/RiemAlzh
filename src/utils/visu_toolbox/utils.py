from PatientObservation import PatientObservation

path_id = "group.csv"
path_ages = "X.csv"
path_scores = "Y.csv"
path = "/Users/clementine.fourrier/RiemAlzh/examples/scalar_models/univariate_model/data/CognitiveScores/MCI_Univar/"

def IsInt(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def readInput(path_to_folder):
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
