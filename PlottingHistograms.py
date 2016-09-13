from ParticleId import *
import ROOT
import os


class PlottingHistograms:

    def __init__(self, events):
        self.nEvents = events
        self.checkOutputDir()
        self.myfile = ROOT.TFile.Open("Output/analysisOutput.root", "RECREATE")
        self.h2dEdxVsMomentum = ROOT.TH2F("h2dEdxVsMomentum", "h2dEdxVsMomentum", 400, -4, 4, 200, 0, 2E-5)
        self.mhSqMassTofVsSqRootNSigma = []
        self.mhSqMassTofPid = []
        self.hNSigma = []
        for pid in ParticleId.name:
            self.mhSqMassTofVsSqRootNSigma.append(ROOT.TH2F(
                "m^{2}_{TOF} vs. #sqrt{n#sigma_{"+pid+",1}^{2} + n#sigma_{"+pid+",2}^{2}}",
                "SqMassTofVsSqRootNSigma_"+pid, 400, 0, 16, 400, -0.5, 1.5))
            self.mhSqMassTofPid.append(ROOT.TH1F("SqMassTofPid_"+pid,
                "m^{2}_{TOF} assuming same mase of two tracks [GeV^{2}/c^{4}], "+pid, 400, -0.5, 1.5))
            self.hNSigma.append(ROOT.TH1F("hNSigma"+pid, "hNSigma"+pid, 160, -20, 20))
        self.mhSqRootNSigmaPionVsKaon = ROOT.TH2F("SqRootNSigmaPionVsKaon",
                "#sqrt{n#sigma_{pion,1}^{2} + n#sigma_{pion,2}^{2}}  vs.  #sqrt{n#sigma_{kaon,1}^{2} + n#sigma_{kaon,2}^{2}}",
                100, 0, 50, 100, 0, 50)
        self.mhSqRootNSigmaPionVsProton = ROOT.TH2F("SqRootNSigmaPionVsProton",
                "#sqrt{n#sigma_{pion,1}^{2} + n#sigma_{pion,2}^{2}}  vs.  #sqrt{n#sigma_{proton,1}^{2} + n#sigma_{proton,2}^{2}}",
                100, 0, 50, 100, 0, 50)
        self.mhSqRootNSigmaKaonVsProton = ROOT.TH2F("SqRootNSigmaKaonVsProton",
                "#sqrt{n#sigma_{kaon,1}^{2} + n#sigma_{kaon,2}^{2}}  vs.  #sqrt{n#sigma_{proton,1}^{2} + n#sigma_{proton,2}^{2}}",
                100, 0, 50, 100, 0, 50)
        self.mhSqMassTof = ROOT.TH1F("SqMassTof", "m^{2}_{TOF} assuming same mase of two tracks [GeV^{2}/c^{4}]", 400, -0.5, 1.5)
        self.mhTofPathLengthVsP = ROOT.TH2F("mhTofPathLengthVsP", "L^{TOF} vs. p", 200, 0, 4, 200, 0, 500)
        self.hPidRecoVsPidGenerated = ROOT.TH2F("hPidRecoVsPidGenerated", "PID Reconstructed vs. PID True-level", 4, 0, 4, 4, 0, 4)

    def checkOutputDir(self):
        if not os.path.exists("Output"):
            os.makedirs("Output")

    def plotHistograms(self):
        file = ROOT.TFile("dEdxData.root")
        hdEdx = file.Get("DEdxVsMomentum")
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptTitle(0)
        c = ROOT.TCanvas("c", "c", 800, 600)
        c.SetLogy()
        c.SetLeftMargin(0.08)
        c.SetRightMargin(0.03)
        c.SetTopMargin(0.04)
        c.SetBottomMargin(0.1)

        colors = [ROOT.kBlack, ROOT.kRed, ROOT.kGreen+3, ROOT.kBlue]

        mhSqMassTofRebinned = self.mhSqMassTof.Rebin(4, "SqMassTofRebinned")
        mhSqMassTofRebinned.SetLineColor(colors[0])
        mhSqMassTofRebinned.SetLineWidth(2)
        nParticles = 3
        mhSqMassTofRebinnedPid = []
        for i in range(0, nParticles):
            mhSqMassTofRebinnedPid.append(self.mhSqMassTofPid[i].Rebin(4, "SqMassTofRebinnedPid_"+ParticleId.name[i]))
            mhSqMassTofRebinnedPid[i].SetLineColor( colors[i+1] )
            mhSqMassTofRebinnedPid[i].SetFillColorAlpha( colors[i+1], 0.3 )
            mhSqMassTofRebinnedPid[i].SetLineWidth( 1 )

        #--- scaling ---
        mhSqMassTofRebinned.Scale( (hdEdx.GetEntries()/2) / self.nEvents)
        for i in range(0, nParticles):
            mhSqMassTofRebinnedPid[i].Scale( (hdEdx.GetEntries()/2) / self.nEvents)
        #---------------

        mhSqMassTofRebinned.GetXaxis().SetTitle("m_{TOF}^{2} [GeV^{2}/c^{4}]")
        mhSqMassTofRebinned.GetXaxis().SetTitleOffset(1.1)
        mhSqMassTofRebinned.GetYaxis().SetRangeUser(0.5, 7.2e4)
        mhSqMassTofRebinned.GetYaxis().SetTitleOffset(1.03)
        binWidthSqMass = mhSqMassTofRebinned.GetBinWidth(1)
        mhSqMassTofRebinned.GetYaxis().SetTitle("Number of track pairs / "+format(binWidthSqMass, '.2f')+" GeV^{2}c^{-4}")
        mhSqMassTofRebinned.Draw("HIST")
        mhSqMassTofRebinnedPid[ParticleId.PION].Draw("HIST SAME")
        mhSqMassTofRebinnedPid[ParticleId.KAON].Draw("HIST SAME")
        mhSqMassTofRebinnedPid[ParticleId.PROTON].Draw("HIST SAME")

        legend = ROOT.TLegend(0.6, 0.65, 0.89, 0.89)
        legend.SetBorderSize(0)
        legend.SetTextSize(0.04)
        legend.SetHeader("Pair PID based on dE/dx")
        legend.AddEntry(mhSqMassTofRebinned, "All pairs (before PID)", "l")
        legend.AddEntry(mhSqMassTofRebinnedPid[ParticleId.PION], "#pi^{+}#pi^{-}", "fl")
        legend.AddEntry(mhSqMassTofRebinnedPid[ParticleId.KAON], "K^{+}K^{-}", "fl")
        legend.AddEntry(mhSqMassTofRebinnedPid[ParticleId.PROTON], "p#bar{p}", "fl")
        legend.Draw()

        line = []
        # kaon m^2 limits
        line.append(ROOT.TLine(0.2,  0.3, 0.2,  2*mhSqMassTofRebinnedPid[ParticleId.KAON].GetMaximum()))
        line.append(ROOT.TLine(0.32, 0.3, 0.32, 2*mhSqMassTofRebinnedPid[ParticleId.KAON].GetMaximum()))
        for i in range(0,2):
            line[i].SetLineColor( ROOT.kGreen+3 )
            line[i].SetLineWidth(2)
            line[i].SetLineStyle(7)
            line[i].Draw()

        #proton m^2 limits
        line[0] = ROOT.TLine(0.7, 0.3, 0.7, 2*mhSqMassTofRebinnedPid[ParticleId.PROTON].GetMaximum())
        line[1] = ROOT.TLine(1.1, 0.3, 1.1, 2*mhSqMassTofRebinnedPid[ParticleId.PROTON].GetMaximum())
        for i in range(0, 2):
            line[i].SetLineColor(ROOT.kBlue)
            line[i].SetLineWidth(2)
            line[i].SetLineStyle(7)
            line[i].Draw()

        c.Print("Output/sqMassFromTof.pdf", "pdf")



        for j in range(0, nParticles):
            c = ROOT.TCanvas("c", "c", 800, 800)
            c.SetLeftMargin(0.09)
            c.SetRightMargin(0.03)
            c.SetTopMargin(0.04)
            c.SetBottomMargin(0.09)
            c.SetLogz()
            self.mhSqMassTofVsSqRootNSigma[j].Rebin2D(2,2)
            self.mhSqMassTofVsSqRootNSigma[j].GetXaxis().SetTitle("n#sigma_{"+ParticleId.name[j]+"}^{pair}   ")
            self.mhSqMassTofVsSqRootNSigma[j].GetXaxis().SetTitleOffset(1.05)
            self.mhSqMassTofVsSqRootNSigma[j].GetYaxis().SetTitle("m_{TOF}^{2} [GeV^{2}/c^{4}]")
            self.mhSqMassTofVsSqRootNSigma[j].GetYaxis().SetTitleOffset(1.18)
            self.mhSqMassTofVsSqRootNSigma[j].Draw("col")
            line2 = ROOT.TLine(3, -0.6, 3, 1.55)
            line2.SetLineWidth(4)
            line2.SetLineStyle(9)
            line2.Draw()
            l = ROOT.TLatex(.7, .15, "(log z-scale)")
            l.SetNDC()
            l.SetTextFont(42)
            l.SetTextSize(.05)
            l.DrawLatex(.7, .15, "(log z-scale)")
            c.Print("Output/SqMassTofVsSqRootNSigma_"+ParticleId.name[j]+".pdf", "pdf")

        c = ROOT.TCanvas("c", "c", 800, 800)
        c.SetLeftMargin(0.09)
        c.SetRightMargin(0.03)
        c.SetTopMargin(0.04)
        c.SetBottomMargin(0.09)
        c.SetLogz()
        #   self.mhSqRootNSigmaPionVsKaon.Rebin2D(2,2);
        self.mhSqRootNSigmaPionVsKaon.GetXaxis().SetTitle("n#sigma_{kaon}^{pair}   ")
        self.mhSqRootNSigmaPionVsKaon.GetYaxis().SetTitle("n#sigma_{pion}^{pair}   ")
        self.mhSqRootNSigmaPionVsKaon.GetXaxis().SetTitleOffset(1.05)
        self.mhSqRootNSigmaPionVsKaon.GetYaxis().SetTitleOffset(1.18)
        self.mhSqRootNSigmaPionVsKaon.GetXaxis().SetRangeUser(0,35)
        self.mhSqRootNSigmaPionVsKaon.GetYaxis().SetRangeUser(0,35)
        self.mhSqRootNSigmaPionVsKaon.Draw("col")

        line[0] = ROOT.TLine(3, -2, 3, 37)
        line[0].SetLineWidth(4)
        line[0].SetLineStyle(9)
        line[0].Draw()
        line[1] = ROOT.TLine(-2, 3, 37, 3)
        line[1].SetLineWidth(4)
        line[1].SetLineStyle(9)
        line[1].Draw()

        l.DrawLatex(.52, .85, "Exclusive candidates")
        l.DrawLatex(.7, .8, "(log z-scale)")
        c.Print("Output/SqRootNSigmaPionVsKaon.pdf", "pdf")

        self.mhSqRootNSigmaKaonVsProton.GetXaxis().SetTitle("n#sigma_{proton}^{pair}   ")
        self.mhSqRootNSigmaKaonVsProton.GetYaxis().SetTitle("n#sigma_{kaon}^{pair}   ")
        self.mhSqRootNSigmaKaonVsProton.GetXaxis().SetTitleOffset(1.05)
        self.mhSqRootNSigmaKaonVsProton.GetYaxis().SetTitleOffset(1.18)
        self.mhSqRootNSigmaKaonVsProton.GetXaxis().SetRangeUser(0,35)
        self.mhSqRootNSigmaKaonVsProton.GetYaxis().SetRangeUser(0,35)
        self.mhSqRootNSigmaKaonVsProton.Draw("col")

        line[0] = ROOT.TLine(3, -2, 3, 37)
        line[0].SetLineWidth(4)
        line[0].SetLineStyle(9)
        line[0].Draw()
        line[1] = ROOT.TLine(-2, 3, 37, 3)
        line[1].SetLineWidth(4)
        line[1].SetLineStyle(9)
        line[1].Draw()

        l.DrawLatex(.52, .85, "Exclusive candidates")
        l.DrawLatex(.7, .8, "(log z-scale)")
        c.Print("Output/SqRootNSigmaKaonVsProton.pdf", "pdf")

        self.mhSqRootNSigmaPionVsProton.GetXaxis().SetTitle("n#sigma_{proton}^{pair}   ")
        self.mhSqRootNSigmaPionVsProton.GetYaxis().SetTitle("n#sigma_{pion}^{pair}   ")
        self.mhSqRootNSigmaPionVsProton.GetXaxis().SetTitleOffset(1.05)
        self.mhSqRootNSigmaPionVsProton.GetYaxis().SetTitleOffset(1.18)
        self.mhSqRootNSigmaPionVsProton.GetXaxis().SetRangeUser(0,35)
        self.mhSqRootNSigmaPionVsProton.GetYaxis().SetRangeUser(0,35)
        self.mhSqRootNSigmaPionVsProton.Draw("col")

        line[0] = ROOT.TLine(3, -2, 3, 37)
        line[0].SetLineWidth(4)
        line[0].SetLineStyle(9)
        line[0].Draw()
        line[1] = ROOT.TLine(-2, 3, 37, 3)
        line[1].SetLineWidth(4)
        line[1].SetLineStyle(9)
        line[1].Draw()

        l.DrawLatex(.52, .85, "Exclusive candidates")
        l.DrawLatex(.7, .8, "(log z-scale)")
        c.Print("Output/SqRootNSigmaPionVsProton.pdf", "pdf")

        c2 = ROOT.TCanvas("c", "c", 800, 800)
        ROOT.gStyle.SetOptTitle(1)
        c2.SetLogz()
        self.hPidRecoVsPidGenerated.GetXaxis().SetBinLabel(1, "Pion")
        self.hPidRecoVsPidGenerated.GetXaxis().SetBinLabel(2, "Kaon")
        self.hPidRecoVsPidGenerated.GetXaxis().SetBinLabel(3, "Proton")
        self.hPidRecoVsPidGenerated.GetXaxis().SetTitle("PID at true level")
        self.hPidRecoVsPidGenerated.GetYaxis().SetBinLabel(1, "Pion")
        self.hPidRecoVsPidGenerated.GetYaxis().SetBinLabel(2, "Kaon")
        self.hPidRecoVsPidGenerated.GetYaxis().SetBinLabel(3, "Proton")
        self.hPidRecoVsPidGenerated.GetYaxis().SetBinLabel(4, "PID failed")
        self.hPidRecoVsPidGenerated.Draw("colz text")
        c2.Print("Output/PidEfficiency.pdf")
        file.Close()
        self.myfile.Write()
        self.myfile.Close()
