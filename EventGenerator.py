# Class representing "data factory". Generates exclusive events
# according to GenEx Monte Carlo, accounting for detector acceptance

import ROOT
import numpy
from dEdxParametrisation import *
from StarDetectorAcceptance import *
from ExclusiveEvent import *
from TrackingSimulation import *


class EventGenerator:
    # constructor
    def __init__(self):
        self.dEdxEngine = dEdxParametrisation()
        self.tracking = TrackingSimulation()
        self.randomNumGen = ROOT.TRandom3(0)
        self.particleProbabilities = self.getProbabilities()
        self.particleLimits = self.getLimits(self.particleProbabilities)
        self.VertexParams = [0.015, 0.015, 50.]
        self.GenExFiles = []
        for pid in ParticleId.name:
            self.GenExFiles.append(open("GenEx/"+pid+".dat", "r"))

    # main method which generates event (ExclusiveEvent object)
    # and returns it
    def generateEvent(self):
        pid = self.generateParticleId()
        vertex = self.generateVertex()
        evt = ExclusiveEvent(pid, vertex)
        evt.p = self.generateMomentum(pid, evt.charge, vertex)
        evt.ToFhitPosition = self.tracking.getToFhitPosition()
        evt.R = self.tracking.getTrackRadius()
        for it in range(2):
            evt.dEdx[it] = self.generatedEdx(pid, evt.p[it])
        return evt

    # generates particles momenta according to GenEx output files,
    # checks for STAR detector acceptance - if particles are outside
    # the acceptance then new set of momenta is loaded (and so on)
    def generateMomentum(self, pid, charges, vertex):
        fourVectors = [ROOT.TLorentzVector(), ROOT.TLorentzVector()]
        while True:
            line = self.GenExFiles[pid].readline()
            data = map(float, line.split())
            fourVectors[0].SetXYZM(data[0], data[1], data[2], ParticleId.mass[pid])
            fourVectors[1].SetXYZM(data[3], data[4], data[5], ParticleId.mass[pid])
            if self.isInAcceptance(fourVectors, pid) and self.tracking.bothTracksReachTOF(fourVectors, charges, vertex):
                break
        return tuple([fourVectors[0].P(), fourVectors[1].P()])

    # generates particle momentum loss (dE/dx) based on particle
    # momentum, according to the Bichsel parametrisation tuned to
    # STAR detector response
    def generatedEdx(self, pid, p):
        return 1E-6*ROOT.TMath.Exp(self.randomNumGen.Gaus(
            self.dEdxEngine.GetMostProbableZ(ROOT.TMath.Log10(p/ParticleId.mass[pid]), 1.),
            self.dEdxEngine.GetRmsZ(ROOT.TMath.Log10(p/ParticleId.mass[pid]), 1.)/7))

    # generates particle ID (in fact, ID of particles in an exclusive pair) according to
    # the yields reconstructed in the data
    def generateParticleId(self):
        x = self.randomNumGen.Uniform()
        for upperBound in self.particleLimits:
            if x <= upperBound:
                return self.particleLimits.index(upperBound)
        return ParticleId.PION

    # extract relative yields of particles according to distributions from the data
    def getProbabilities(self):
        file = ROOT.TFile("dEdxData.root")
        eventCounts = []
        for pid in ParticleId.name:
            eventCounts.append(file.Get("DEdxVsMomentumPid_" + pid).GetEntries() / 2)
        return tuple(eventCounts / numpy.sum(eventCounts))

    # calculates numerical limits for the purpose of random PID generation
    def getLimits(self, pValues):
        limits = []
        for p in pValues:
            if p == pValues[0]:
                limits.append(p)
            else:
                limits.append(limits[len(limits) - 1] + p)
        return limits

    # checks whether the tracks are within acceptance of the central barrel detector
    def isInAcceptance(self, fourVec, pid):
        if StarDetectorAcceptance.etaLimits[0] < fourVec[0].Eta() < StarDetectorAcceptance.etaLimits[1] and \
           StarDetectorAcceptance.etaLimits[0] < fourVec[1].Eta() < StarDetectorAcceptance.etaLimits[1] and \
           fourVec[0].Pt() > StarDetectorAcceptance.pTThreshold[pid] and \
           fourVec[1].Pt() > StarDetectorAcceptance.pTThreshold[pid]:
            return True
        else:
            return False

    # generates vertex position
    def generateVertex(self):
        return ROOT.TVector3(self.randomNumGen.Gaus(self.VertexParams[0], self.VertexParams[1]),
                            self.randomNumGen.Gaus(self.VertexParams[0], self.VertexParams[1]),
                            self.randomNumGen.Gaus(self.VertexParams[0], self.VertexParams[2]))
