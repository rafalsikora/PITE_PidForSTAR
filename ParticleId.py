# Class representing particle type (name, mass)

class ParticleId:
    PION = 0
    KAON = 1
    PROTON = 2

    mass = tuple([0.13956995, 0.493677, 0.93827231])
    name = tuple(["pion", "kaon", "proton"])

    def __init__(self):
        pass
