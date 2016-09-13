# Class describing Bichsel parametrisation of the
# energy loss dE/dx in the Time Projection Chamber
# in STAR main detector

import ROOT

class dEdxParametrisation:
    def __init__(self, Tag="p10", MostProbableZShift=0, AverageZShift=0, I70Shift=1, I60Shift=1):
        self.fTag = ROOT.TString(Tag)
        self.fP = 0
        self.fA = 0
        self.fI70 = 0
        self.fI60 = 0
        self.fD = 0
        self.fRms = 0
        self.fW = 0
        self.fPhi = 0
        self.fMostProbableZShift = MostProbableZShift
        self.fAverageZShift = AverageZShift
        self.fI70Shift = I70Shift
        self.fI60Shift = I60Shift
        self.fbgL10min = -1
        self.fbgL10max = 4
        self.fdxL2min = -0.3
        self.fdxL2max = 3
        self.fzmin = -4
        self.fzmax = 6

        dir = ROOT.gDirectory
        rootf = "P10T.root"
        if self.fTag.Contains("pai", ROOT.TString.kIgnoreCase):
            rootf = "PaiT.root"
        elif self.fTag.Contains("p10", ROOT.TString.kIgnoreCase):
            rootf = "P10T.root"
        elif self.fTag.Contains("bich", ROOT.TString.kIgnoreCase):
            rootf = "BichselT.root"

        path = "dEdxModel"
        file = ROOT.gSystem.Which(path, rootf, ROOT.kReadPermission)
        if not file:
            print("dEdxParameterization::GetFile: File " + rootf + " has not been found in path " + path)
        else:
            print("dEdxParameterization::GetFile: File " + rootf + " has been found as " + file)
            pFile = ROOT.TFile(file)

            self.fP = pFile.Get("bichP")
            self.fP.SetDirectory(0)
            self.fA = pFile.Get("bichA")
            self.fA.SetDirectory(0)
            self.fI70 = pFile.Get("bichI70")
            self.fI70.SetDirectory(0)
            self.fI60 = pFile.Get("bichI60")
            self.fI60.SetDirectory(0)
            self.fD = pFile.Get("bichD")
            self.fD.SetDirectory(0)
            self.fRms = pFile.Get("bichRms")
            self.fRms.SetDirectory(0)
            self.fW = pFile.Get("bichW")
            self.fW.SetDirectory(0)
            self.fPhi = pFile.Get("bichPhi")
            self.fPhi.SetDirectory(0)

            self.fbgL10min = self.fPhi.GetXaxis().GetBinCenter(1) + 1e-7
            self.fbgL10max = self.fPhi.GetXaxis().GetBinCenter(self.fPhi.GetXaxis().GetNbins()) - 1e-7
            self.fdxL2min = self.fPhi.GetYaxis().GetBinCenter(1) + 1e-7
            self.fdxL2max = self.fPhi.GetYaxis().GetBinCenter(self.fPhi.GetYaxis().GetNbins()) - 1e-7
            self.fzmin = self.fPhi.GetZaxis().GetBinCenter(1) + 1e-7
            self.fzmax = self.fPhi.GetZaxis().GetBinCenter(self.fPhi.GetZaxis().GetNbins()) - 1e-7

            self.fAXYZ = tuple([self.fPhi.GetXaxis(), self.fPhi.GetYaxis(), self.fPhi.GetZaxis()])
            self.fnBins = tuple([x.GetNbins() for x in self.fAXYZ])
            self.fbinW = tuple([x.GetBinWidth(1) for x in self.fAXYZ])

            # set normalization factor to 2.3976 keV/cm at beta*gamma = 4
            dEdxMIP = 2.39761562607903311 # [keV/cm]
            MIPBetaGamma10 = ROOT.TMath.Log10(4.)
            #  fMostProbableZShift = ROOT.TMath.Log(dEdxMIP) - Interpolation(self.fP, MIPBetaGamma10, 1, 0)
            #  fAverageZShift      = ROOT.TMath.Log(dEdxMIP) - Interpolation(self.fA, MIPBetaGamma10, 1, 0)
            self.fI70Shift *= dEdxMIP/self.GetI70(MIPBetaGamma10, 1)
            self.fI60Shift *= dEdxMIP/self.GetI60(MIPBetaGamma10, 1)
            self.fMostProbableZShift = ROOT.TMath.Log(self.fI70Shift)
            self.fAverageZShift = self.fMostProbableZShift

    def BetaGamma_Dx(self, log10bg, log2dx):
        log10bg = ROOT.TMath.Max(self.fbgL10min, ROOT.TMath.Min(self.fbgL10max, log10bg))
        log2dx = ROOT.TMath.Max(self.fdxL2min, ROOT.TMath.Min(self.fdxL2max, log2dx))
        return tuple([log10bg, log2dx])

    def GetMostProbableZ(self, log10bg, log2dx):
        log10bg = ROOT.TMath.Max(self.fbgL10min, ROOT.TMath.Min(self.fbgL10max, log10bg))
        log2dx = ROOT.TMath.Max(self.fdxL2min, ROOT.TMath.Min(self.fdxL2max, log2dx))
        return self.fMostProbableZShift + self.fP.Interpolate(log10bg, log2dx)

    def GetAverageZ(self, log10bg,  log2dx):
        (log10bg, log2dx) = self.BetaGamma_Dx(log10bg, log2dx)
        return self.fAverageZShift + self.MostProbableZCorrection(log10bg) + self.fA.Interpolate(log10bg, log2dx)
  
    def GetRmsZ(self, log10bg, log2dx):
        (log10bg, log2dx) = self.BetaGamma_Dx(log10bg, log2dx)
        return self.fRms.Interpolate(log10bg, log2dx)
  
    def GetI70(self, log10bg,  log2dx):
        (log10bg, log2dx) = self.BetaGamma_Dx(log10bg, log2dx)
        return self.fI70Shift*self.fI70.Interpolate(log10bg, log2dx)
  
    def GetI60(self, log10bg,  log2dx):
        (log10bg, log2dx) = self.BetaGamma_Dx(log10bg, log2dx)
        return self.fI60Shift*self.fI60.Interpolate(log10bg, log2dx)
  
    def GetMostProbabledEdx(self, log10bg,  log2dx):
        (log10bg, log2dx) = self.BetaGamma_Dx(log10bg, log2dx)
        return self.fD.Interpolate(log10bg,log2dx)
  
    def GetdEdxWidth(self, log10bg,  log2dx):
        (log10bg, log2dx) = self.BetaGamma_Dx(log10bg, log2dx)
        return self.fW.Interpolate(log10bg, log2dx)
  
    def GetMostProbableZM(self, log10bg,  log2dx):
        (log10bg, log2dx) = self.BetaGamma_Dx(log10bg, log2dx)
        return self.MostProbableZCorrection(log10bg) + self.GetMostProbableZ(log10bg, log2dx)
  
    def GetAverageZM(self, log10bg,  log2dx):
        (log10bg, log2dx) = self.BetaGamma_Dx(log10bg, log2dx)
        return self.MostProbableZCorrection(log10bg) + self.GetAverageZ(log10bg,log2dx)
  
    def GetI70M(self, log10bg, log2dx):
        (log10bg, log2dx) = self.BetaGamma_Dx(log10bg, log2dx)
        return self.I70Correction(log10bg)*self.GetI70(log10bg, log2dx)
  
    def GetProbability(self, log10bg,  log2dx,  z):
        (log10bg, log2dx) = self.BetaGamma_Dx(log10bg, log2dx)
        z = ROOT.TMath.Max(self.fzmin, ROOT.TMath.Min(self.fzmax, z))
        return self.fPhi.Interpolate(log10bg, log2dx, z)

    def MostProbableZCorrection(self, log10bg):
        pars = (-3.68846e-03, 4.72944e+00)
        return pars[0]*ROOT.TMath.Exp(-pars[1]*log10bg)

    def I70Correction(self, log10bg):
        pars = (-1.65714e-02, 3.27271e+00)
        return ROOT.TMath.Exp(pars[0]*ROOT.TMath.Exp(-pars[1]*log10bg))

    def Tag(self):
        return self.fTag.Data()

    def P(self):
        return self.fP

    def A(self):
        return self.fA

    def I70(self):
        return self.fI70

    def I60(self):
        return self.fI60

    def D(self):
        return self.fD

    def Rms(self):
        return self.fRms

    def W(self):
        return self.fW

    def Phi(self):
        return self.fPhi

    def bgL10min(self):
        return self.fbgL10min

    def bgL10max(self):
        return self.fbgL10max


if __name__ == '__main__':
    testObject = dEdxParametrisation()

    def bichsel70(x, par):
        pove = x[0]
        mass = par[0]
        if mass < 0:
            mass = -mass
        poverm = pove/mass
        return 1E-6*ROOT.TMath.Exp(testObject.GetMostProbableZ(ROOT.TMath.Log10(poverm), 1.))
        #return 1E-6*testObject.GetMostProbabledEdx(ROOT.TMath.Log10(poverm), 1.)
        #return 1E-6*testObject.GetI70(ROOT.TMath.Log10(poverm), 1.)

    file = ROOT.TFile("dEdxData.root")
    canv = ROOT.TCanvas("canv","",800,600)
    ROOT.gStyle.SetOptStat(0)
    h2dEdxVsMomentum = file.Get("DEdxVsMomentum")
    h2dEdxVsMomentum.Draw("colz")
    particleMass = tuple([0.13956995, 0.493677, 0.93827231])
    particleName = tuple(["Pion", "Kaon", "Proton"])
    particleColor = tuple([2, 1, 4])
    func = []
    for i in range(len(particleMass)):
        func.append(ROOT.TF1("dEdx_vs_p_"+particleName[i], bichsel70, 0.1, 3.5, 1))
        func[i].SetParameter(0, particleMass[i])
        func[i].SetLineColor(particleColor[i])
        func[i].Draw("same")
    canv.Print("TestOfDEdxParametrisation.pdf")