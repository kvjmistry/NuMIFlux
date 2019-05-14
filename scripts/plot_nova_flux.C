/*
This script will plot the nova flux for the intersection and detweights methods and compare them to novas flux.

To run this script execute the command root -l 'plot_nova_flux.C("nue")'
where nue, nuebar, numu and numubar are the available options.

*/

#include "plot_comp_functions.h"


void plot_nova_flux(const char* mode) { 


    TCanvas* c1 = new TCanvas();
    TLegend* lFlux = new TLegend(0.5, 0.65, 0.9, 0.9);
    bool boolfile, boolhist;

    TFile *f, *f_nova;

    // File in 
	// boolfile  = GetFile(f,"/uboone/data/users/kmistry/work/PPFX/nova/Norm_investigation/test/output.root"); if (boolfile == false) gSystem->Exit(0);
    boolfile  = GetFile(f,"/uboone/data/users/kmistry/work/PPFX/nova/Norm_investigation/test/output_v2.root"); if (boolfile == false) gSystem->Exit(0);
    boolfile  = GetFile(f_nova,"/uboone/app/users/kmistry/PPFX/numi-validation/nova_flux/FHC_Flux_NOvA_ND_2017.root"); if (boolfile == false) gSystem->Exit(0);

    double fPOT = GetPOT(f);
    // double normfactor = 1.0e6 / (fPOT*50);
    double normfactor =  1.0e6 / (fPOT);

    // Now get the flux for intersection
	TH1D* h_intersection_flux;
	boolhist = GetHist(f, h_intersection_flux, Form("%s/Window/%s_CV_Window", mode, mode)); if (boolhist == false) gSystem->Exit(0);
    Normalise(h_intersection_flux);
    h_intersection_flux->SetLineColor(kMagenta);
    h_intersection_flux->SetLineWidth(2);
    h_intersection_flux->Scale(normfactor);
    lFlux->AddEntry(h_intersection_flux, "Intersection/Window","l");
    h_intersection_flux->Draw("hist");
    h_intersection_flux->Draw("e,same");


    // Now get the flux for deteweights
	TH1D* h_detweights_flux;
	boolhist = GetHist(f, h_detweights_flux, Form("%s/Detsmear/%s_CV_AV_TPC", mode, mode)); if (boolhist == false) gSystem->Exit(0);
    Normalise(h_detweights_flux);
    h_detweights_flux->SetLineColor(kBlack);
    h_detweights_flux->SetLineWidth(2);
    // h_detweights_flux->Scale(normfactor* 1./(3.1415926)  );
    h_detweights_flux->Scale(normfactor );
    lFlux->AddEntry(h_detweights_flux, "Detector Smear","l");
    h_detweights_flux->Draw("hist, same");
    h_detweights_flux->Draw("e, same");

    // // Now get the flux for deteweights_notilt
	// TH1D* h_detweights_notilt_flux;
	// boolhist = GetHist(f, h_detweights_notilt_flux, Form("%s/%s_CV_AV_TPC_detweights_notilt", mode, mode)); if (boolhist == false) gSystem->Exit(0);
    // Normalise(h_detweights_notilt_flux);
    // h_detweights_notilt_flux->SetLineColor(kGreen+1);
    // h_detweights_notilt_flux->SetLineWidth(2);
    // // h_detweights_notilt_flux->Scale(normfactor* 1./(3.1415926) );
    // h_detweights_notilt_flux->Scale(normfactor );
    // lFlux->AddEntry(h_detweights_notilt_flux, "Detector Smear no Tilt","l");
    // h_detweights_notilt_flux->Draw("hist, same");
    // h_detweights_notilt_flux->Draw("e, same");

    // Now get Nova flux 
   TH1D* h_nova_flux;
   boolhist = GetHist(f_nova, h_nova_flux, Form("flux_%s", mode)); if (boolhist == false) gSystem->Exit(0);
   h_nova_flux->SetLineColor(kBlue+1);
   h_nova_flux->SetLineWidth(2);
   lFlux->AddEntry(h_nova_flux, "NOvA","l");
   h_nova_flux->Draw("hist, same");
   h_nova_flux->Draw("e, same");

//    h_intersection_flux->Scale(h_nova_flux->Integral(3,-1)/h_intersection_flux->Integral(3,-1));
//    h_detweights_flux->Scale(h_nova_flux->Integral(3,-1)/h_detweights_flux->Integral(3,-1));
//    h_detweights_notilt_flux->Scale(h_nova_flux->Integral(3,-1)/h_detweights_notilt_flux->Integral(3,-1));

   lFlux->Draw();
   gPad->SetLogy();
   h_intersection_flux->GetYaxis()->SetRangeUser(0.01, 1000);

   c1->Update();




}


