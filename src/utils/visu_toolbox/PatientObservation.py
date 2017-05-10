class PatientObservation:
    def __init__(self, id):
        self.id = id
        self.map_obs = {}

    def add_observation(self, age, observation):
        self.map_obs[age] = observation
