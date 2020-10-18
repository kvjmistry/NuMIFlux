
/*
This script will plot the flux at MicroBooNE for all the neutrino flavours on one plot

To run this script run the command root -l 'plot_total_flux.C("fhc")' 
where fhc/rhc are the available options.

This file depends on the functions.h script so make
sure this file is included in the same directory.

*/

#include "functions.h"



void plot_total_flux(const char* horn) {
    
    gStyle->SetOptStat(0); // say no to stats box

    TFile *file; 

    bool boolfile;

    // Get the flux file
    if (!strcmp(horn,"fhc")) {
        boolfile  = GetFile(file ,"/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold_v46/FHC/output_uboone_fhc_run0_set1.root");
        if (boolfile == false) gSystem->Exit(0);
    }
    else {
        boolfile  = GetFile(file ,"/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold_v46/RHC/output_uboone_rhc_run0_set1.root");
        if (boolfile == false) gSystem->Exit(0);
    }

    std::vector<std::string> flav = {"numu", "nue", "numubar", "nuebar" };


    // Resize the histogram vector
    std::vector<TH1D*> h_flux;
    h_flux.resize(flav.size());

    
    // Get the Histograms
    for (unsigned int f = 0; f < flav.size(); f++){
        h_flux.at(f) = (TH1D*) file->Get( Form("%s/Detsmear/%s_CV_AV_TPC_5MeV_bin", flav.at(f).c_str(), flav.at(f).c_str()) ); 
        
    }
    
    // define the energy threshold to integrate from 
    double energy_threshold = 0.125; // Energy threhsold
    double flux_int_tot = 0.0;
    std::vector<double> flux_int(flav.size(), 0.0);

    // Loop over the histograms, get the bin for the threshold and get the integral
    for (unsigned int f = 0; f < flav.size(); f++){
        
        double xbin_th = h_flux.at(f)->GetXaxis()->FindBin(energy_threshold);   // find the x bin to integrate from

        // Get the flux for the neutrino flavour
        flux_int.at(f) = h_flux.at(f)->Integral(xbin_th, h_flux.at(f)->GetNbinsX()+1);
        

        // Add to the total flux
        flux_int_tot += flux_int.at(f);

    }

    // Customise and normalise the histograms
    for (unsigned int f = 0; f < flav.size(); f++){
        
        // Rebin to slightly wider 10 MeV  bins
        h_flux.at(f) ->Rebin(2);
    
        // Normalise the histogram by bin width
        Normalise(h_flux.at(f));

        double POT_Scale = 1.0; // The POT to scale to // was 6.0e20
        
        // Get the POT from the file
        double fPOT = GetPOT(file);

        // Now scale, 1e-4 for m2->cm2
        h_flux.at(f)->Scale( (POT_Scale)/ (fPOT*1.0e4) );  

        // set generic histogram properties
        h_flux.at(f)->GetXaxis()->SetLabelSize(0.05);
        h_flux.at(f)->GetXaxis()->SetTitleSize(0.05);
        h_flux.at(f)->GetYaxis()->SetLabelSize(0.05);
        h_flux.at(f)->GetYaxis()->SetTitleSize(0.05);
        h_flux.at(f)->SetTitle(";Neutrino Energy [GeV];#nu / POT / GeV / cm^{2}");
        h_flux.at(f)->GetXaxis()->SetRangeUser(0,5);
        
    }
    
    TCanvas* c = new TCanvas();
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.12);
    
    TLegend* leg = new TLegend(0.74, 0.65, 0.89, 0.9);
    
    // NuMu
    h_flux.at(0)->SetLineColor(kRed+1);
    h_flux.at(0)->SetLineWidth(2);
    leg->AddEntry(h_flux.at(0), Form("#nu_{#mu} (%2.1f%%)",       100 * flux_int.at(0) / flux_int_tot),"l");
    h_flux.at(0)->Draw("hist,same");
    
    // Nue
    h_flux.at(1)->SetLineColor(kRed+1);
    h_flux.at(1)->SetLineWidth(2);
    h_flux.at(1)->SetLineStyle(2);
    leg->AddEntry(h_flux.at(1), Form("#nu_{e} (%2.1f%%)",         100 * flux_int.at(1) / flux_int_tot),"l");
    h_flux.at(1)->Draw("hist,same");
    
    // NuMu Bar
    h_flux.at(2)->SetLineColor(kBlue+1);
    h_flux.at(2)->SetLineWidth(2);
    leg->AddEntry(h_flux.at(2), Form("#bar{#nu}_{#mu} (%2.1f%%)", 100 * flux_int.at(2) / flux_int_tot),"l");
    h_flux.at(2)->Draw("hist,same");
    
    
    // Nue bar
    h_flux.at(3)->SetLineColor(kBlue+1);
    h_flux.at(3)->SetLineWidth(2);
    h_flux.at(3)->SetLineStyle(2);
    leg->AddEntry(h_flux.at(3), Form("#bar{#nu}_{e} (%2.1f%%)",   100 * flux_int.at(3) / flux_int_tot),"l");
    h_flux.at(3)->Draw("hist,same");

    c->SetLogy();

    leg->SetNColumns(1);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(62); 
    leg->Draw();

    Draw_Nu_Mode(c, horn); // Draw FHC Mode/RHC Mode Text

    c->Print(Form("NuMI_%s_Flux.pdf", horn));

}


