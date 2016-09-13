# Class performing simulation of a charged track propagation
# in the magnetic field of STAR central detector

import ROOT
import numpy
from StarDetectorAcceptance import *


class TrackingSimulation:
    # constructor
    def __init__(self):
        self.elCharge = 1.602176565e-19
        self.B = 0.5  # T
        self.BfieldVec = ROOT.TVector3(0, 0, self.B)
        self.tofBarrelRadius = 210  # cm
        self.cMperS = 299792458
        self.timeStep = 1.e-9  # s
        self.ToFhitPosition = [ROOT.TVector3(), ROOT.TVector3()]
        self.R = numpy.empty(2, dtype=float)

    # returns 3-vector of Lorentz force acting on a particle of given charge, mass
    # and momentum, passing through the longitudinal magnetic field self.B
    def F(self, pVec, charge, mass):
        return (charge * self.elCharge / mass) * pVec.Cross(self.BfieldVec)

    # computes shift of a track in the time interval self.timeStep
    def computeStep(self, rVec, pVec, charge, lorentzGamma, mass):
        speed = self.cMperS / ROOT.TMath.Sqrt(1 + pow(mass / pVec.Mag(), 2))
        rVec += pVec.Unit() * self.timeStep * speed
        pVec += self.F(pVec, charge, mass) * self.timeStep * (1. / lorentzGamma)

    # returns the 3-vector with (x,y,z) position of the TOF module, which was
    # hit by the particle of given charge and four-momentum
    def getTofHitPositionVector(self, fourVector, charge, vertexVector):
        rVec = ROOT.TVector3(vertexVector)
        while True:
            self.computeStep(rVec, ROOT.TVector3(fourVector.Px(), fourVector.Py(), fourVector.Pz()), charge,
                             fourVector.Gamma(), fourVector.M())
            if not (rVec.Perp() < self.tofBarrelRadius):
                break
        return rVec

    # returns True if both exclusive tracks have reached the TOF barrel,
    # otherwise returns False
    def bothTracksReachTOF(self, fourVectors, charges, vertex):
        skipEvent = False
        for it in range(2):
            self.R[it] = 100 * (ROOT.TVector3(fourVectors[it].Px(), fourVectors[it].Py(), fourVectors[it].Pz()).Mag()) / (0.3 * self.B)
            if self.R[it] < self.tofBarrelRadius / 2:
                skipEvent = True
                break
            self.ToFhitPosition[it] = self.getTofHitPositionVector(fourVectors[it], charges[it], vertex)
            if not StarDetectorAcceptance.zLimits[0] < self.ToFhitPosition[it].z() < StarDetectorAcceptance.zLimits[1]:
                skipEvent = True
                break
        return False if skipEvent else True

    # returns list with 3-vectors describing position of the TOF modules
    # thar have been hit by the tracks
    def getToFhitPosition(self):
        return self.ToFhitPosition

    # returns list with radii of track helices
    def getTrackRadius(self):
        return self.R
