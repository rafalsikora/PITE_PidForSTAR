# Main file of the simulation

import ROOT
from EventGenerator import *
from EventReconstruction import *
from PlottingHistograms import *
from numpy import *


# number of events to simulate
nEvents = 2e2

gen = EventGenerator()
reco = EventReconstruction(gen.dEdxEngine)
plotHist = PlottingHistograms(nEvents)

# main loop
for i in range(0, int(nEvents)):
    print i
    evt = gen.generateEvent()
    reco.reconstructEvent(evt)
    plotHist.h2dEdxVsMomentum.Fill(-evt.p[0], evt.dEdx[0])
    plotHist.h2dEdxVsMomentum.Fill(evt.p[1], evt.dEdx[1])
    squaredMass = evt.mSquared
    plotHist.mhSqMassTof.Fill(squaredMass)
    nParticles = 3
    nSigmaVal = zeros((2, 3), float)
    sqRootNSigma = []
    for w in range(0, nParticles):
        for j in range(0, 2):
            nSigmaVal[j][w] = evt.nSigma[j][w]
            plotHist.hNSigma[w].Fill(nSigmaVal[j][w])
        sqRootNSigma.append(ROOT.TMath.Sqrt(nSigmaVal[0][w]*nSigmaVal[0][w] + nSigmaVal[1][w]*nSigmaVal[1][w]))

    #determine the PID of pair
    pidReco = nParticles
    if sqRootNSigma[ParticleId.PION]>3 and sqRootNSigma[ParticleId.KAON]>3 and sqRootNSigma[ParticleId.PROTON]<3 and squaredMass > 0.7 and squaredMass < 1.1:
        pidReco = ParticleId.PROTON
    elif sqRootNSigma[ParticleId.PION]>3 and sqRootNSigma[ParticleId.KAON]<3 and sqRootNSigma[ParticleId.PROTON]>3 and squaredMass > 0.2 and squaredMass < 0.32:
        pidReco = ParticleId.KAON
    elif nSigmaVal[0][ParticleId.PION]<3 and nSigmaVal[1][ParticleId.PION]<3:
        pidReco = ParticleId.PION

    if sqRootNSigma[ParticleId.PION]>3 and sqRootNSigma[ParticleId.KAON]>3 and sqRootNSigma[ParticleId.PROTON]<3:
        plotHist.mhSqMassTofPid[ParticleId.PROTON].Fill( squaredMass )
    elif sqRootNSigma[ParticleId.PION]>3 and sqRootNSigma[ParticleId.KAON]<3 and sqRootNSigma[ParticleId.PROTON]>3:
        plotHist.mhSqMassTofPid[ParticleId.KAON].Fill(squaredMass)
    elif pidReco == ParticleId.PION:
        plotHist.mhSqMassTofPid[ParticleId.PION].Fill(squaredMass)

    for j in range(0, nParticles):
        plotHist.mhSqMassTofVsSqRootNSigma[j].Fill(sqRootNSigma[j], squaredMass)
    plotHist.mhSqRootNSigmaPionVsKaon.Fill(sqRootNSigma[ParticleId.KAON], sqRootNSigma[ParticleId.PION])
    plotHist.mhSqRootNSigmaPionVsProton.Fill(sqRootNSigma[ParticleId.PROTON], sqRootNSigma[ParticleId.PION])
    plotHist.mhSqRootNSigmaKaonVsProton.Fill(sqRootNSigma[ParticleId.PROTON], sqRootNSigma[ParticleId.KAON])
    plotHist.hPidRecoVsPidGenerated.Fill(evt.pairID, pidReco)

# print histograms in the Output directory
plotHist.plotHistograms()
