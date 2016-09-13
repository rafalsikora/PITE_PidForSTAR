# Class representing 2-particle exclusive event

from ParticleId import *
import numpy
import ROOT

class ExclusiveEvent:
    def __init__(self, pid, vertex):
        assert(isinstance(pid, int))
        assert(isinstance(vertex, ROOT.TVector3))
        self.pairID = pid
        self.vrt = vertex
        self.charge = (1, -1)
        self.p = numpy.empty(2, dtype=float)
        self.dEdx = numpy.empty(2, dtype=float)
        self.trkLength = numpy.empty(2, dtype=float)
        self.R = numpy.empty(2, dtype=float)
        self.tofTime = numpy.empty(2, dtype=float)
        self.nSigma = numpy.empty([2, len(ParticleId.mass)], dtype=float)
        self.ToFhitPosition = []
        self.mSquared = numpy.empty(2, dtype=float)
