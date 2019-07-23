import os,sys,string, time
import ROOT
from math import *
from ROOT import TTree, TObject, TFile, gDirectory, TH1D, TH2D, TH3D, TCanvas, gROOT, TGaxis, gStyle, TColor, TLegend, THStack, TChain, TLatex, TText, TLine
#from ROOT import *
from array import array
from glob import glob

ROOT.gStyle.SetOptStat(0);

# Importing rootlogon.C
ROOT.gROOT.SetMacroPath('~/');
ROOT.gROOT.Macro( os.path.expanduser( 'rootlogon.C' ) )

# set kBird palette
red   = [ 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764];
green = [ 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832];
blue = [ 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539];
stops = [ 0.0, 0.1, 0.2, 0.3, 0.4 ,0.5 ,0.6 ,0.7, 1.0 ];
stopsArray = array("d", stops)
redArray = array("d", red)
greenArray = array("d", green)
blueArray = array("d", blue)
ROOT.TColor.CreateGradientColorTable(9, stopsArray, redArray, greenArray, blueArray, 255)

# Opening root file
#f = TFile("files/NuMIFlux.root")
f = TFile("../NuMIFlux.root")
f.ls()

h_Kplus_long_trans = f.Get("Kplus_NA49")

c4 = TCanvas()

h_Kplus_long_trans.SetTitle("K^{+}")
h_Kplus_long_trans.GetXaxis().SetTitle("x_{F} ")
# h_Kplus_long_trans.GetXaxis().SetTitleFont(12)
h_Kplus_long_trans.GetYaxis().SetTitle("p_{T} [GeV/c]")
# h_Kplus_long_trans.GetYaxis().SetTitleFont(12)
h_Kplus_long_trans.GetXaxis().SetRangeUser(-0.04,1)
h_Kplus_long_trans.GetYaxis().SetRangeUser(0,3.5)
h_Kplus_long_trans.GetXaxis().SetLabelSize(0.05)
h_Kplus_long_trans.GetXaxis().SetTitleSize(0.05)
h_Kplus_long_trans.GetYaxis().SetLabelSize(0.05)
h_Kplus_long_trans.GetYaxis().SetTitleSize(0.05)
h_Kplus_long_trans.Draw("colz")
c4.SetLeftMargin(0.15)
c4.SetBottomMargin(0.12)
c4.SetRightMargin(0.14)
gStyle.SetTitleH(0.07)

# Gemma Tinti's (MINOS thesis) Kplus / Kminus measurements (NA49, thin target, @158 GeV,  so with tgptype particle) [FERMILAB-THESIS-2010-44]

#box 1
line1 = TLine( -0.0125, 0.05, 0.225, 0.05 )
line2 = TLine( 0.225, 0.05, 0.225, 0.75 )
line3 = TLine( 0.225, 0.75 , -0.0125, 0.75 )
line4 = TLine( -0.0125, 0.75 , -0.0125 ,0.05 )

#box 2
line5 = TLine( -0.0125, 0.8, 0.225, 0.8 )
line6 = TLine( 0.225, 0.8, 0.225, 1 )
line7 = TLine( 0.225, 1 , -0.0125, 1 )
line8 = TLine( -0.0125, 1 , -0.0125 ,0.8 )

lineVector= [ line1, line2, line3, line4, line5, line6, line7, line8]
for line in lineVector:
    line.SetLineColor(ROOT.kBlack)
    line.SetLineWidth(3)
    line.Draw()

leg = TLegend(.40, .75 , .44 , .79)
leg.SetFillStyle(0);
leg.SetLineWidth(4);
leg.Draw();


text = TLatex(.17, .755  , "NA49 Coverage");
text.SetTextColor(ROOT.kBlack);
text.SetNDC();
text.SetTextSize(1.4/30.);
text.SetTextAlign(11);
#text.DrawLatex(.48, .55, "#Square");
text.Draw();

c4.Print("plots/kaonplus_NA49.pdf");

# raw_input("Please press enter to exit.")

