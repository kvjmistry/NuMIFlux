#!/usr/bin/env python
#
# read a root anatree using Pyroot
#
# run with python pyrootMacroTPCActive.py 14 -1
# 14 is the nu flavour (nu-mu in this case)
# -1 means loop over all the entries

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

h_pionminus_long_trans = f.Get("pionminus_MIPP")

cpionminus_long_trans = TCanvas("c1","c1",1500,1200)

h_pionminus_long_trans.SetTitle("#pi^{-}")
h_pionminus_long_trans.GetXaxis().SetTitle("p_{Z} [GeV/c]")
h_pionminus_long_trans.GetXaxis().SetTitleFont(12)
h_pionminus_long_trans.GetYaxis().SetTitle("p_{T} [GeV/c]")
h_pionminus_long_trans.GetYaxis().SetTitleFont(12)
h_pionminus_long_trans.GetXaxis().SetRangeUser(0,130)
h_pionminus_long_trans.GetYaxis().SetRangeUser(0,5)
h_pionminus_long_trans.Draw("colz")

# Seun's K/pi measurements (MIPP, thick target, so with tptype particle) [FERMILAB-THESIS-2007-61]

line1 = TLine( 20, 0, 20, 1 )
line2 = TLine( 20, 1, 24, 1 )
line3 = TLine( 24, 1, 24, 1.2 )
line4 = TLine( 24, 1.2, 31, 1.2 )
line5 = TLine( 31, 1.2, 31, 1.55 )
line6 = TLine( 31, 1.55, 42, 1.55 )
line7 = TLine( 42, 1.55, 42, 2 )
line8 = TLine( 42, 2, 60, 2 )
line9 = TLine( 60, 2, 90, 2 )
line10 = TLine( 90, 2, 90, 0 )
line11 = TLine( 90, 0, 20, 0 )

lineVector= [ line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11 ]
for line in lineVector:
    line.SetLineColor(ROOT.kBlack)
    line.SetLineWidth(3)
    line.Draw()

leg = TLegend(.44, .55, .48, .59)
leg.SetFillStyle(0);
leg.SetLineWidth(4);
leg.Draw();

text = TLatex(.5, .55, "MIPP Coverage");
text.SetTextColor(ROOT.kBlack);
text.SetNDC();
text.SetTextSize(1.4/30.);
text.SetTextAlign(11);
#text.DrawLatex(.48, .55, "#Square");
text.Draw();

cpionminus_long_trans.Print("plots/pionminus_MIPP.pdf");

# raw_input("Please press enter to exit.")
