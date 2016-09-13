# Class describing main properties of the STAR central detector

from ParticleId import *
import numpy


class StarDetectorAcceptance:
    etaLimits = tuple([-1.0, 1.0])
    pTThreshold = numpy.empty(len(ParticleId.mass), dtype=float)
    pTThreshold[ParticleId.PION] = 0.15
    pTThreshold[ParticleId.KAON] = 0.2
    pTThreshold[ParticleId.PROTON] = 0.3
    zLimits = tuple([-215, 215])

    def __init__(self):
        pass
