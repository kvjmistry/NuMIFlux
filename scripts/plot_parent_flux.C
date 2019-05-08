/*
This script will plot the nuetrino flux broken down into the various parent types.

To run this script execute the command root -l 'plot_parent_flux.C("nue")'
where nue, nuebar, numu and numubar are the available options.

*/

#include "plot_comp_functions.h"

// ----------------------------------------------------------------------------
// Function for plotting
void SetColour(TH1D* hist, std::string parent, TLegend* leg){
	hist->SetLineWidth(2);
	
	if (parent == "PI_Plus"){
		hist->SetLineColor(42);
		leg->AddEntry(hist, "#pi^{+}","l");
		hist->Draw("hist,same");
	} 
	else if  (parent == "PI_Minus"){
		hist->SetLineColor(kMagenta+2);
		leg->AddEntry(hist, "#pi^{-}","l");
		hist->Draw("hist,same");
	}
	else if  (parent == "Mu_Plus"){ // only turn on if there is MIPP
		hist->SetLineColor(30);
		leg->AddEntry(hist, "#mu^{+}","l");
		hist->Draw("hist,same");
	}
	else if  ( parent == "Mu_Minus"){ // only turn on if there is MIPP
		hist->SetLineColor(38);
		leg->AddEntry(hist, "#mu^{-}","l");
		hist->Draw("hist,same");
	}
	else if  (parent == "Kaon_Plus"){
		hist->SetLineColor(28);
		leg->AddEntry(hist, "K^{+}","l");
		hist->Draw("hist,same");
	}
	else if  (parent == "Kaon_Minus"){
		hist->SetLineColor(36);
		leg->AddEntry(hist, "K^{-}","l");
		hist->Draw("hist,same");
	}
	else if  (parent == "K0L"){
		hist->SetLineColor(1001);
		leg->AddEntry(hist, "K^{0}_{L}","l");
		hist->Draw("hist,same");
	}
	else
		std::cout << "Unkown parent type"<< std::endl;
}
// ----------------------------------------------------------------------------
// Main
void plot_parent_flux(TString mode) { // (mippon/mippoff, input, Product/noThinKaon etc. numu/nue)
	gStyle->SetOptStat(0); // say no to stats box
	
	// Pre declare variables
	const char *Gethist_CV_AV;
	TFile *f;
	bool boolfile, boolhist;
	double rebin{5}; // number to rebin the histograms by
	const char * mode_char;

	std::vector<std::string> parent = {"PI_Plus", "PI_Minus", "Mu_Plus", "Mu_Minus", "Kaon_Plus", "Kaon_Minus" , "K0L"}; // parent names
	std::vector<TH1D*> h_parent;
	h_parent.resize(parent.size());

	// Select neutrino type to run with 
	switch (return_mode(mode)){
		case enumu:
			std::cout << "\nUsing Numu Mode!\n" << std::endl;
			Gethist_CV_AV = "numu/numu_CV_AV_TPC_rebin";      // CV Flux using the 5MeV bins for now 
			mode_char = "numu";
			break;

		case enue:
			std::cout << "\nUsing Nue Mode!\n" << std::endl;
			Gethist_CV_AV = "nue/nue_CV_AV_TPC_rebin";
			mode_char = "nue";
			break;

		case enumubar:
			std::cout << "\nUsing Numubar Mode!\n" << std::endl;
			Gethist_CV_AV = "numubar/numubar_CV_AV_TPC_rebin";
			mode_char = "numubar";
			break;

		case enuebar:
			std::cout << "\nUsing Nuebar Mode!\n" << std::endl;
			Gethist_CV_AV = "nuebar/nuebar_CV_AV_TPC_rebin";
			mode_char = "nuebar";
			break;

	}
	
	// ------------------------------------------------------------------------------------------------------------
	// Flux Comparisons
	// ------------------------------------------------------------------------------------------------------------
	TCanvas* c1 = new TCanvas();
	TLegend* lFlux = new TLegend(0.7, 0.45, 0.9, 0.9);

	// File in 
	boolfile  = GetFile(f,"/uboone/data/users/kmistry/work/PPFX/uboone/parent/output.root"); if (boolfile == false) gSystem->Exit(0);
	// boolfile  = GetFile(f,"/uboone/data/users/kmistry/work/g4numi/rm_KMinus_Capture/output_wDAR.root"); if (boolfile == false) gSystem->Exit(0); // turn off K- capture at rest
	
	// Get POT
	double fPOT = GetPOT(f);
	double normfactor = 6.0e20 / (fPOT*1e4);
	
	// Now get the cv flux
	TH1D* h_ppfx_flux;
	boolhist = GetHist(f, h_ppfx_flux, Gethist_CV_AV); if (boolhist == false) gSystem->Exit(0);
	
	// Rebin and normalise
	h_ppfx_flux->Rebin(rebin);
	Normalise(h_ppfx_flux);
	IncreaseLabelSize(h_ppfx_flux);
	h_ppfx_flux->Scale(normfactor);
	h_ppfx_flux->Draw("hist,same");
	h_ppfx_flux->SetLineWidth(2);
	h_ppfx_flux->SetLineColor(kBlack);

	if (mode == "numu")		h_ppfx_flux->SetTitle("#nu_{#mu}; Energy [GeV]; #nu / 6 #times 10^{20} POT / 25 MeV / cm^{2}");
	if (mode == "nue")		h_ppfx_flux->SetTitle("#nu_{e}; Energy [GeV]; #nu / 6 #times 10^{20} POT / 25 MeV / cm^{2}");
	if (mode == "numubar")	h_ppfx_flux->SetTitle("#bar{#nu_{#mu}}; Energy [GeV]; #nu / 6 #times 10^{20} POT / 25 MeV / cm^{2}");
	if (mode == "nuebar")	h_ppfx_flux->SetTitle("#bar{#nu_{e}}; Energy [GeV]; #nu / 6 #times 10^{20} POT / 25 MeV / cm^{2}");

	lFlux->AddEntry(h_ppfx_flux, "Total Flux","l");
	
	// Now do the same for each neutrino parent
	boolhist = GetHist(f, h_ppfx_flux, Gethist_CV_AV); if (boolhist == false) gSystem->Exit(0);

	for (unsigned int i = 0; i < parent.size(); i++){
		boolhist = GetHist(f, h_parent.at(i), Form("%s/%s/Enu_%s_%s_AV_TPC", mode_char, parent.at(i).c_str(),mode_char, parent.at(i).c_str() )); if (boolhist == false) gSystem->Exit(0);
		h_parent.at(i)->Rebin(rebin);
		Normalise(h_parent.at(i));
		h_parent.at(i)->Scale(normfactor);
		SetColour(h_parent.at(i), parent.at(i), lFlux);
		// lFlux->AddEntry(h_parent.at(i), parent.at(i).c_str(),"l");
	}
	h_ppfx_flux->Draw("hist,same"); // draw again so it is on top of all the other components

	lFlux->SetNColumns(1);
	lFlux->SetBorderSize(0);
	lFlux->SetFillStyle(0);
	// lFlux->SetTextFont(62); 
	lFlux->SetTextSize(0.05);
	lFlux->Draw();
	gPad->SetLogy();


	// ++++++++++++++++++++++++++++++++++
	// Save the plots as pdfs in the plots folder
	// ++++++++++++++++++++++++++++++++++
	// create plots folder if it does not exist
	gSystem->Exec("if [ ! -d \"plots\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir plots; fi"); 
	
		
	c1->Print(Form("plots/%s_parentflux.pdf", mode_char));
	std::cout << "\n"<< std::endl;
	
	
// gSystem->Exit(0);

} // end of main




