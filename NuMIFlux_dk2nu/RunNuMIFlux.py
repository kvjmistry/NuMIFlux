#!/usr/bin/env python
#

import os,sys,string, time
import ROOT
ROOT.gSystem.Load("NuMIFlux_cc.so")
from ROOT import NuMIFlux
from glob import glob

# fname = glob("/pnfs/uboone/persistent/uboonebeam/numi_dk2nu_zero_threshold/FHC/g4numiv6*.root")
fname = glob("/pnfs/uboone/persistent/uboonebeam/numi_dk2nu_zero_threshold/FHC/g4numiv6_minervame_me000z200i_11*.root")

f = NuMIFlux()
f.CalculateFlux()


raw_input("Please press enter to exit.")
