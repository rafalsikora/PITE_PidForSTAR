# Class reconstructing DST variables

import ROOT
from ParticleId import *
from ExclusiveEvent import *


class EventReconstruction:
    # constructor
    def __init__(self, dEdxEngine):
        self.c = 29.9792  # cm/ns
        self.dEdxEngine = dEdxEngine
        self.randomNumGen = ROOT.TRandom3(0)
        self.tofResolution = 0.1  # ns
        self.momentumResolution = 0.02  # %

    # method takes an ExclusiveEvent object, performs reconstruction and fills
    # object with reconstructed data
    def reconstructEvent(self, evt):
        assert isinstance(evt, ExclusiveEvent)
        self.getNSigma(evt)
        self.getTofPathLength(evt)
        self.getTofTime(evt)
        self.getSquaredMass(evt)

    # calculates nSigma variables accordind to Bichsel parametrisation
    # of dE/dx(p) for three PID assumptions: pi, K and p
    def getNSigma(self, evt):
        for i in range(2):
            for j in range(len(ParticleId.mass)):
                pOverM = evt.p[i]/ParticleId.mass[j]
                evt.nSigma[i][j] = ROOT.TMath.Log(evt.dEdx[i]/(1e-6*ROOT.TMath.Exp(self.dEdxEngine.GetMostProbableZ(
                                   ROOT.TMath.Log10(pOverM), 1.))))/(self.dEdxEngine.GetRmsZ(
                                   ROOT.TMath.Log10(pOverM), 1.)/7)

    # calculates length of a path that particle has followed from the vertex
    # to the TOF module
    def getTofPathLength(self, evt):
        for i in range(2):
            xdif = evt.ToFhitPosition[i].x() - evt.vrt.x()
            ydif = evt.ToFhitPosition[i].y() - evt.vrt.y()
            C = ROOT.TMath.Sqrt(xdif * xdif + ydif * ydif)
            s_perp = 2 * evt.R[i] * ROOT.TMath.ASin(C / (2 * evt.R[i]))
            s_z = abs(evt.ToFhitPosition[i].z() - evt.vrt.z())
            evt.trkLength[i] = ROOT.TMath.Sqrt(s_perp * s_perp + s_z * s_z)

    # calculates time of detection in the TOF module (accounts for resolution)
    def getTofTime(self, evt):
        for i in range(2):
            pOverM = evt.p[i] / ParticleId.mass[evt.pairID]
            evt.tofTime[i] = evt.trkLength[i] * ROOT.TMath.Sqrt(1./(pOverM * pOverM) + 1) / self.c +\
                             self.randomNumGen.Gaus(0, self.tofResolution)

    # calculates squared mass assuming equal masses of two particles,
    # making use of tracks momenta, lengths and times of detection in TOF
    # (accounts for resolution effects)
    def getSquaredMass(self, evt):
        LSq = []
        PSq = []
        for i in range(2):
            LSq.append(evt.trkLength[i]*evt.trkLength[i] / (self.c*self.c))
            momentumSmearing = self.randomNumGen.Gaus(1.0, self.momentumResolution)
            PSq.append(momentumSmearing * momentumSmearing * evt.p[i] * evt.p[i])
        dtHitSq = (evt.tofTime[0] - evt.tofTime[1]) * (evt.tofTime[0] - evt.tofTime[1])
        aEq = -2 * LSq[0] * LSq[1] / (PSq[0] * PSq[1]) + LSq[0] * LSq[0] / (PSq[0] * PSq[0]) +\
              LSq[1] * LSq[1] / (PSq[1] * PSq[1])
        bEq = -2 * LSq[0] * LSq[1] * (1 / PSq[0] + 1 / PSq[1]) + 2 * LSq[0] * LSq[0] / PSq[0] +\
              2 * LSq[1] * LSq[1] / PSq[1] - 2 * dtHitSq * (LSq[0] / PSq[0] + LSq[1] / PSq[1])
        cEq = dtHitSq * dtHitSq - 2 * dtHitSq * (LSq[0] + LSq[1]) +\
              LSq[0] * LSq[0] + LSq[1] * LSq[1] - 2 * LSq[0] * LSq[1]
        evt.mSquared = (-bEq + ROOT.TMath.Sqrt(bEq * bEq - 4 * aEq * cEq)) / (2 * aEq)
