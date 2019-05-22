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
import csv

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

h_pionplus_long_trans = f.Get("pionplus_NA49")

c4 = TCanvas("c1","c1",1500,1100)

h_pionplus_long_trans.SetTitle(" All #nus from #pi^{+}")
h_pionplus_long_trans.GetXaxis().SetTitle("x_{F} ")
h_pionplus_long_trans.GetXaxis().SetTitleFont(12)
h_pionplus_long_trans.GetYaxis().SetTitle("p_{T} [GeV/c]")
h_pionplus_long_trans.GetYaxis().SetTitleFont(12)
h_pionplus_long_trans.GetXaxis().SetRangeUser(-0.15,1)
h_pionplus_long_trans.GetYaxis().SetRangeUser(0,3.5)
h_pionplus_long_trans.Draw("colz")

# FOR each row, for each column, fill TLine( 1,2,3,4)
#file_reader = csv.reader('NA49_pion_coverage.csv', 'rb', delimiter=' ,')
#    for row in file_reader:
#...         print ', '.join(row)
#xF    pT

#box 1
line1 = TLine( -0.1125 ,   0.35 , -0.0875  ,  0.35)
line2 = TLine( -0.0875 ,   0.35   ,  -0.0875  ,  0.25)
line3 = TLine( -0.0875 ,   0.25 ,   -0.0625 ,   0.25)
line4 = TLine( -0.0625 ,   0.25,    -0.0625  ,  0.75)
line5 = TLine( -0.0625  ,  0.75,    -0.055  ,  0.75)
line6 = TLine( -0.055  ,  0.75,  -0.055   ,  0.35)
line7 = TLine( -0.055   ,  0.35 ,    0.17  ,  0.35)
line8 = TLine( 0.17  ,  0.35,   0.17  ,  0.65)
line9 = TLine( 0.17  ,  0.65 ,   0.175  ,  0.65)
line10 = TLine( 0.175  ,  0.65   ,0.175  ,  0.35)
line11 = TLine( 0.175  ,  0.35, 0.225  ,  0.35)
line12 = TLine( 0.225  ,  0.35,  0.225  ,  0.325)
line13 = TLine( 0.225  ,  0.325, 0.175 ,   0.325)
line14 = TLine( 0.175 ,   0.325  ,   0.175  ,  0.025)
line15 = TLine( 0.175  ,  0.025, 0.225  ,  0.025)
line16 = TLine(0.225  ,  0.025  ,   0.225  ,  0.045)
line17 = TLine(0.225  ,  0.045,    0.325  ,  0.045)
line18 = TLine(0.325  ,  0.045, 0.325  ,  0.65)
line19 = TLine(0.325  ,  0.65,   0.225  ,  0.65)
line20 = TLine(0.225  ,  0.65  ,   0.225   , 0.7)
line21 = TLine(0.225   , 0.7  ,   0.325  ,  0.7)
line22 = TLine(0.325  ,  0.7   ,   0.325  ,  1.3)
line23 = TLine(0.325  ,  1.3,  0.35  ,  1.3)
line24= TLine(0.35  ,  1.3    ,0.35  ,  0.7)
line25 = TLine(0.35  ,  0.7 ,  0.55 ,   0.7)
line26 = TLine( 0.55 ,   0.7    ,   0.55  ,  1.5)
line27 = TLine( 0.55  ,  1.5    ,   0.45 ,   1.5)
line28 = TLine( 0.45 ,   1.5,  0.45 ,   1.7)
line29= TLine( 0.45 ,   1.7    ,0.25  ,  1.7)
line30 = TLine( 0.25  ,  1.7    ,   0.25 ,   1.9)
line31 = TLine( 0.25 ,   1.9    ,   0.15  ,  1.9)
line32 = TLine( 0.15  ,  1.9    ,0.15   , 1.3)
line33 = TLine( 0.15   , 1.3    ,0.125  ,  1.3)
line34 = TLine( 0.125  ,  1.3   ,   0.125  ,  1.9)
line35 = TLine( 0.125  ,  1.9   ,   -0.125  ,  1.9)
line36 = TLine( -0.125   , 1.9  ,   -0.125  ,  1.3)
line37 = TLine( -0.125  ,  1.3  ,   -0.1125 ,   1.3)
line38 = TLine( -0.1125 ,   1.3 ,   -0.1125    ,1.1)
line39 = TLine( -0.1125    ,1.1    ,0.1125   , 1.1)
line40 = TLine( 0.1125    ,1.1 ,   0.1125   , 1.3)
line41 = TLine( 0.1125   , 1.3 ,   0.125  ,  1.3)
line42 = TLine(0.125  ,  1.3   ,   0.125  ,  1.1)
line43 = TLine( 0.125  ,  1.1   ,   0.225   , 1.1)
line44 = TLine( 0.225   , 1.1   ,   0.225  ,  1.05)
line45 = TLine( 0.225  ,  1.05   ,   0.125  ,  1.05)
line46 = TLine( 0.125  ,  1.05  ,   0.125  ,  0.65)
line47 = TLine( 0.125  ,  0.65  ,   0.1125  ,  0.65)
line48 = TLine( 0.1125  ,  0.65 ,   0.1125  ,  1.05)
line49 = TLine( 0.1125  ,  1.05 ,   -0.1125 ,   1.05)
line50 = TLine(-0.1125 ,   1.05,    -0.1125, 0.35)

#box 2
line51 = TLine(-0.055 ,   0.325  ,  -0.055, 0.175)
line52 = TLine(-0.055, 0.175,    -0.045, 0.175)
line53 = TLine(-0.045, 0.175,  -0.045  ,  0.12)
line54 = TLine(-0.045  ,  0.12, -0.035  ,  0.125)
line55 = TLine( -0.035  ,  0.125,    -0.035 ,   0.075)
line56 = TLine( -0.035  ,  0.075,   -0.025  ,  0.075)
line57 = TLine( -0.025  ,  0.075,    -0.025  ,  0.025)
line58 = TLine( -0.025 ,   0.025,   0.1625  ,  0.025)
line59 = TLine( 0.1625  ,  0.025,   0.1625  ,  0.325)
line60 = TLine(  0.1625  ,  0.325, -0.055 ,   0.325 )

#box 3

line61 = TLine( 0.35, 0.65, 0.35, 0.045 )
line62 = TLine( 0.35, 0.045, 0.55, 0.045 )
line63 = TLine( 0.55, 0.045, 0.55, 0.65 )
line64 = TLine( 0.55, 0.65, 0.35, 0.65 )

lineVector= [ line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11, line12, line13, line14, line15, line16, line17, line18, line19, line20, line21, line22,line23,line24,line25,line26,line27,line28,line29,line30,line31,line32,line33,line34,line35,line36,line37,line38,line39,line40, line41,line42,line43,line44,line45,line46,line47,line48, line49, line50, line51, line52, line53,line54,line55,line56,line57,line58, line59,line60, line61, line62, line63, line64]

for line in lineVector:
    line.SetLineColor(ROOT.kBlack)
    line.SetLineWidth(3)
    line.Draw()


leg = TLegend(.54, .78 , .58, .82)
leg.SetFillStyle(0);
leg.SetLineWidth(4);
leg.Draw();


text = TLatex(.6, .78, "NA49 Coverage");
text.SetTextColor(ROOT.kBlack);
text.SetNDC();
text.SetTextSize(1.4/30.);
text.SetTextAlign(11);
#text.DrawLatex(.48, .55, "#Square");
text.Draw();

c4.Print("plots/pionplus_NA49.pdf");

raw_input("Please press enter to exit.")
