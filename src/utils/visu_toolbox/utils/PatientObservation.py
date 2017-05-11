class PatientObservation:
    def __init__(self, id):
        self.id = id
        self.obs = {}

    def add_observation(self, age, observation):
        self.obs[age] = observation
