/*
This script will plot the flux at microboone for the flux made using flugg/gsimple files
and compare this with the dk2nu and ppfx predictions.

To run this script run the command root -l 'plot_gsimple_flux.C("fhc", "nue")' 
where fhc/rhc, nue/nuebar/numu/numubar are the available options.

This file depends on the functions.h script so make
sure this file is included in the same directory.

*/

#include "functions.h"

// ----------------------------------------------------------------------------
// Main
void plot_gsimple_flux(const char* horn, const char* mode) {
	gStyle->SetOptStat(0); // say no to stats box

	// Pre declare variables
	TString Gethist_TPC, g_simp_names, Gethist_TPC_dk2nu;
	TH1D *h_dk2nu_flux, *h_g_simp, *hnovafileflux, *hppfx, *hppfx_mod; 
	TFile *f_gsimple, *f_ppfx, *f_ppfx_mod;
	bool boolfile, boolhist;
	double rebin{5}; // number to rebin the histograms by

	// Select neutrino type to run with 
	if (!strcmp(mode, "numu")){
			std::cout << "\nUsing NuMu Mode!\n" << std::endl;
			Gethist_TPC = "numu/Detsmear/numu_CV_AV_TPC_5MeV_bin";			// AV in TPC flux prediction
			g_simp_names = "numuFluxHisto";									// G simple files
			Gethist_TPC_dk2nu = "numu/Detsmear/numu_UW_AV_TPC_5MeV_bin";	// uw AV in TPC flux prediction
	}
	else if (! strcmp(mode,"nue")){
			std::cout << "\nUsing Nue Mode!\n" << std::endl;
			Gethist_TPC = "nue/Detsmear/nue_CV_AV_TPC_5MeV_bin";
			g_simp_names = "nueFluxHisto";
			Gethist_TPC_dk2nu = "nue/Detsmear/nue_UW_AV_TPC_5MeV_bin";
	}
	else if (!strcmp(mode, "numubar")){
			std::cout << "\nUsing NuMubar Mode!\n" << std::endl;
			Gethist_TPC = "numubar/Detsmear/numubar_CV_AV_TPC_5MeV_bin";
			g_simp_names = "anumuFluxHisto";
			Gethist_TPC_dk2nu = "numubar/Detsmear/numubar_UW_AV_TPC_5MeV_bin";
	}
	else if (!strcmp(mode, "nuebar")){
			std::cout << "\nUsing Nuebar Mode!\n" << std::endl;
			Gethist_TPC = "nuebar/Detsmear/nuebar_CV_AV_TPC_5MeV_bin";
			g_simp_names = "anueFluxHisto";
			Gethist_TPC_dk2nu = "nuebar/Detsmear/nuebar_UW_AV_TPC_5MeV_bin";
	}
	else {
		std::cout << "Unknown nuetrino flavour type" << std::endl;
	}
	
	// ------------------------------------------------------------------------------------------------------------
	// Flux Comparisons
	// ------------------------------------------------------------------------------------------------------------
	TCanvas* c1 = new TCanvas();
	TLegend* lFlux = new TLegend(0.5, 0.65, 0.9, 0.9);

	// Dk2nu
	// boolfile  = GetFile(f_ppfx_mod ,"/uboone/data/users/kmistry/work/PPFX/uboone/DetectorWeights_withtilt/output.root"); if (boolfile == false) gSystem->Exit(0);
	// boolfile  = GetFile(f_ppfx_mod ,"/uboone/data/users/kmistry/work/PPFX/uboone/parent/v3/output_parent_all.root"); if (boolfile == false) gSystem->Exit(0); // file with all large weights kept and just the ones which are < 0 removed
	
	if (!strcmp(horn,"fhc")) {
		boolfile  = GetFile(f_ppfx_mod ,"/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold/output_uboone_run0.root");
		if (boolfile == false) gSystem->Exit(0);
		boolfile  = GetFile(f_gsimple , "/uboone/data/users/kmistry/work/PPFX/uboone/NuMIFlux_update_morebins.root"); 
		if (boolfile == false) gSystem->Exit(0);
	}
	else {
		boolfile  = GetFile(f_ppfx_mod ,"/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold/RHC/output_uboone_run0.root");
		if (boolfile == false) gSystem->Exit(0);
		boolfile  = GetFile(f_gsimple , "/uboone/app/users/kmistry/Flux/NuMIFlux/files/NuMIFlux_anti_morebins.root"); 
		if (boolfile == false) gSystem->Exit(0);
	}
	
	double fPOT = GetPOT(f_ppfx_mod);
	boolhist = GetHist(f_ppfx_mod, h_dk2nu_flux, Gethist_TPC_dk2nu); if (boolhist == false) gSystem->Exit(0);
	h_dk2nu_flux->Rebin(rebin);
	Normalise(h_dk2nu_flux);

	// Get Gsimple files	
	
	boolhist = GetHist(f_gsimple, h_g_simp, g_simp_names); if (boolhist == false) gSystem->Exit(0);
	h_g_simp->Rebin(rebin);
	Normalise(h_g_simp);

	// PPFX flux
	boolhist = GetHist(f_ppfx_mod, hppfx, Gethist_TPC); if (boolhist == false) gSystem->Exit(0);
	hppfx->Rebin(rebin);
	Normalise(hppfx);
	hppfx->SetDirectory(0);
	h_dk2nu_flux->SetDirectory(0);

	// Scalings
	h_dk2nu_flux->Scale((6.0e20)/ (fPOT*1.0e4) );  
	hppfx->Scale( (6.0e20)/ (fPOT*1.0e4) );  

	// Plottings
	h_dk2nu_flux->SetLineColor(kRed+1);
	h_dk2nu_flux->SetLineWidth(2);
	h_dk2nu_flux->SetTitle(";E_{#nu} [GeV];#nu / 6 #times 10^{20} POT / GeV / cm^{2}");
	// h_dk2nu_flux->SetTitle(";E_{#nu} [GeV];#nu / 6 #times 10^{20} POT / GeV / cm^{2}");
	IncreaseLabelSize(h_dk2nu_flux);
	h_dk2nu_flux->GetXaxis()->SetRangeUser(0,5);
	h_dk2nu_flux->Draw("hist");

	h_g_simp->SetLineWidth(2);
	h_g_simp->Draw("hist, same"); // Only Draw for FHC mode 

	hppfx->SetLineColor(kGreen+1);
	hppfx->SetLineWidth(2);
	hppfx->Draw("hist, same");

	gPad->SetLogy();
	gPad->Update();

	lFlux->SetNColumns(1);
	lFlux->SetBorderSize(0);
	lFlux->SetFillStyle(0);
	lFlux->SetTextFont(62); 
	
	lFlux->AddEntry(h_dk2nu_flux, "dk2nu","l");
	lFlux->AddEntry(h_g_simp, "flugg","l");
	lFlux->AddEntry(hppfx, "dk2nu PPFX Corrected","l");;
	lFlux->Draw();
	Draw_Nu_Mode(c1, horn); // Draw FHC Mode/RHC Mode Text

	if (!strcmp(mode, "numu"))		h_dk2nu_flux->SetTitle("#nu_{#mu}");
	if (!strcmp(mode, "nue"))		h_dk2nu_flux->SetTitle("#nu_{e}");
	if (!strcmp(mode, "numubar"))	h_dk2nu_flux->SetTitle("#bar{#nu_{#mu}}");
	if (!strcmp(mode, "nuebar"))	h_dk2nu_flux->SetTitle("#bar{#nu_{e}}");

	// h_dk2nu_flux->SetTitleSize(0.05);
	if (!strcmp(mode,"nue") || !strcmp(mode,"numu")) gStyle->SetTitleH(0.1);
	else gStyle->SetTitleH(0.07);

	// ------------------------------------------------------------------------------------------------------------
	// Draw all fluxes on one plot
	// ------------------------------------------------------------------------------------------------------------
	TCanvas* c_plotall = new TCanvas();
	TLegend* l_plotall = new TLegend(0.74, 0.65, 0.89, 0.9);

	c_plotall->cd();
	double tot_flux = GetTotalFlux(f_ppfx_mod); // Get the total integrated flux from all flavours
	PlotFluxSame(c_plotall, l_plotall, f_ppfx_mod, "numu", fPOT, "numu/Detsmear/numu_UW_AV_TPC_5MeV_bin", tot_flux );
	PlotFluxSame(c_plotall, l_plotall, f_ppfx_mod, "numubar", fPOT, "numubar/Detsmear/numubar_UW_AV_TPC_5MeV_bin", tot_flux );
	PlotFluxSame(c_plotall, l_plotall, f_ppfx_mod, "nue", fPOT, "nue/Detsmear/nue_UW_AV_TPC_5MeV_bin", tot_flux );
	PlotFluxSame(c_plotall, l_plotall, f_ppfx_mod, "nuebar", fPOT, "nuebar/Detsmear/nuebar_UW_AV_TPC_5MeV_bin", tot_flux );
	Draw_Nu_Mode(c_plotall, horn); // Draw FHC Mode/RHC Mode Text

	l_plotall->SetNColumns(1);
	l_plotall->SetBorderSize(0);
	l_plotall->SetFillStyle(0);
	l_plotall->SetTextFont(62); 
	l_plotall->Draw();

	// ------------------------------------------------------------------------------------------------------------
	// Save the plots as pdfs in the plots folder
	// ------------------------------------------------------------------------------------------------------------
	// create plots folder if it does not exist
	gSystem->Exec("if [ ! -d \"plots\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir plots; fi"); 
	gSystem->Exec("if [ ! -d \"plots/CV_Flux\" ]; then echo \"\n CV_Flux folder does not exist... creating\"; mkdir plots/CV_Flux; fi"); 
	
	c1->Print(Form("plots/CV_Flux/CV_Flux_Prediction_%s_%s.pdf", horn, mode));
	c_plotall->Print(Form("plots/CV_Flux/CV_Flux_%s.pdf", horn));

	// gSystem->Exit(0);

} // end of main




