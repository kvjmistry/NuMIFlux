/*
This script will plot the nuetrino flux broken down into the various parent types.

To run this script execute the command root -l 'plot_parent_flux.C("fhc","nue")'
where fhc/rhc, nue/nuebar/numu/numubar are the available options.

*/

#include "functions.h"

// ----------------------------------------------------------------------------
// Function for plotting
void SetColour(TH1D* hist, std::string parent, TLegend* leg, char flux_pcent[15]){
	hist->SetLineWidth(2);
	
	if (parent == "PI_Plus"){
		hist->SetLineColor(42);
		leg->AddEntry(hist, Form("#pi^{+} (%s%%)", flux_pcent),"l");
		hist->Draw("hist,same");
	} 
	else if  (parent == "PI_Minus"){
		hist->SetLineColor(kMagenta+2);
		leg->AddEntry(hist, Form("#pi^{-} (%s%%)", flux_pcent),"l");
		hist->Draw("hist,same");
	}
	else if  (parent == "Mu_Plus"){ // only turn on if there is MIPP
		hist->SetLineColor(30);
		leg->AddEntry(hist, Form("#mu^{+} (%s%%)", flux_pcent),"l");
		hist->Draw("hist,same");
	}
	else if  ( parent == "Mu_Minus"){ // only turn on if there is MIPP
		hist->SetLineColor(38);
		leg->AddEntry(hist, Form("#mu^{-} (%s%%)", flux_pcent),"l");
		hist->Draw("hist,same");
	}
	else if  (parent == "Kaon_Plus"){
		hist->SetLineColor(28);
		leg->AddEntry(hist, Form("K^{+} (%s%%)", flux_pcent),"l");
		hist->Draw("hist,same");
	}
	else if  (parent == "Kaon_Minus"){
		hist->SetLineColor(36);
		leg->AddEntry(hist, Form("K^{-} (%s%%)", flux_pcent),"l");
		hist->Draw("hist,same");
	}
	else if  (parent == "K0L"){
		hist->SetLineColor(1001);
		leg->AddEntry(hist, Form("K^{0}_{L} (%s%%)", flux_pcent),"l");
		hist->Draw("hist,same");
	}
	else
		std::cout << "Unkown parent type"<< std::endl;
}
// ----------------------------------------------------------------------------
// Main
void plot_parent_flux(const char* horn, TString mode) { // (fhc/rhc, numu/nue)
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
			Gethist_CV_AV = "numu/Detsmear/numu_CV_AV_TPC_5MeV_bin";      // CV Flux using the 5MeV bins for now 
			mode_char = "numu";
			break;

		case enue:
			std::cout << "\nUsing Nue Mode!\n" << std::endl;
			Gethist_CV_AV = "nue/Detsmear/nue_CV_AV_TPC_5MeV_bin";
			mode_char = "nue";
			break;

		case enumubar:
			std::cout << "\nUsing Numubar Mode!\n" << std::endl;
			Gethist_CV_AV = "numubar/Detsmear/numubar_CV_AV_TPC_5MeV_bin";
			mode_char = "numubar";
			break;

		case enuebar:
			std::cout << "\nUsing Nuebar Mode!\n" << std::endl;
			Gethist_CV_AV = "nuebar/Detsmear/nuebar_CV_AV_TPC_5MeV_bin";
			mode_char = "nuebar";
			break;

	}
	
	// ------------------------------------------------------------------------------------------------------------
	// Flux Comparisons
	// ------------------------------------------------------------------------------------------------------------
	TCanvas* c1 = new TCanvas();
	TLegend* lFlux = new TLegend(0.69, 0.45, 0.89, 0.9);

	// FHC File in 
	if (!strcmp(horn, "fhc")) {
		boolfile  = GetFile(f,"/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold_v46/FHC/output_uboone_fhc_run0_set1.root");
		if (boolfile == false) gSystem->Exit(0);
	}
	else { // RHC file
		boolfile  = GetFile(f,"/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold_v46/RHC/output_uboone_rhc_run0_set1.root");
		if (boolfile == false) gSystem->Exit(0); // turn off K- capture at rest
	}

	// Get POT
	double fPOT = GetPOT(f);
	double normfactor = 1.0 / (fPOT*1e4);
	
	// Now get the cv flux
	TH1D* h_ppfx_flux;
	boolhist = GetHist(f, h_ppfx_flux, Gethist_CV_AV); if (boolhist == false) gSystem->Exit(0);
	

	// Define Integrals
  double energy_threshold = 0.06;
  std::cout <<"Using Energy Threshold of: " << energy_threshold << std::endl;
  double xbin_th = h_ppfx_flux->GetXaxis()->FindBin(energy_threshold);   // find the x bin to integrate from
	double flux_cv  = h_ppfx_flux->Integral(xbin_th,  h_ppfx_flux->GetNbinsX()+1) * normfactor;
	std::vector<double> flux_parent(parent.size());

	// Rebin and normalise
	h_ppfx_flux->Rebin(rebin);
	Normalise(h_ppfx_flux);
	IncreaseLabelSize(h_ppfx_flux);
	h_ppfx_flux->Scale(normfactor);
	h_ppfx_flux->Draw("hist,same");
	h_ppfx_flux->SetLineWidth(2);
	h_ppfx_flux->SetLineColor(kBlack);
	h_ppfx_flux->GetXaxis()->SetRangeUser(0,4.0);
	h_ppfx_flux->SetLineStyle(2);

	if (mode == "numu")		h_ppfx_flux->SetTitle("#nu_{#mu}; Energy [GeV]; #nu / POT / 25 MeV / cm^{2}");
	if (mode == "nue")		h_ppfx_flux->SetTitle("#nu_{e}; Energy [GeV]; #nu / POT / 25 MeV / cm^{2}");
	if (mode == "numubar")	h_ppfx_flux->SetTitle("#bar{#nu_{#mu}}; Energy [GeV]; #nu / POT / 25 MeV / cm^{2}");
	if (mode == "nuebar")	h_ppfx_flux->SetTitle("#bar{#nu_{e}}; Energy [GeV]; #nu / POT / 25 MeV / cm^{2}");
	
	if (mode == "nue" || mode == "numu") gStyle->SetTitleH(0.1);
	else gStyle->SetTitleH(0.07);

	lFlux->AddEntry(h_ppfx_flux, "Total Flux","l");
	
	// Now do the same for each neutrino parent
	boolhist = GetHist(f, h_ppfx_flux, Gethist_CV_AV); if (boolhist == false) gSystem->Exit(0);

	for (unsigned int i = 0; i < parent.size(); i++){
		boolhist = GetHist(f, h_parent.at(i), Form("%s/%s/Enu_%s_%s_AV_TPC", mode_char, parent.at(i).c_str(),mode_char, parent.at(i).c_str() )); if (boolhist == false) gSystem->Exit(0);
    xbin_th = h_parent.at(i)->GetXaxis()->FindBin(energy_threshold);   // find the x bin to integrate from
		flux_parent.at(i) =  h_parent.at(i)->Integral(xbin_th,  h_parent.at(i)->GetNbinsX()+1) * normfactor;
		h_parent.at(i)->Rebin(rebin);
		Normalise(h_parent.at(i));
		h_parent.at(i)->Scale(normfactor);

		// Convert the flux percentage to a char
		char flux_pcent[15];
		snprintf(flux_pcent, 15,"%2.1f" ,100 * flux_parent.at(i) / flux_cv);

		SetColour(h_parent.at(i), parent.at(i), lFlux, flux_pcent );
	}
	h_ppfx_flux->Draw("hist,same"); // draw again so it is on top of all the other components

	Draw_Nu_Mode(c1, horn); // Draw FHC Mode/RHC Mode Text

	lFlux->SetNColumns(1);
	lFlux->SetBorderSize(0);
	lFlux->SetFillStyle(0);
	// lFlux->SetTextFont(62); 
	lFlux->SetTextSize(0.05);
	lFlux->Draw();
	gPad->SetLogy();
	if (std::string(horn) == "numu" || std::string(horn) == "numubar")
		h_ppfx_flux->SetMinimum(1e-12);
	else 
		h_ppfx_flux->SetMinimum(1e-14);

	// Now print the percentages for the flux
	std::cout << std::setw(10);
	std::cout << "----------------------------------" << std::endl;
	std::cout << "\nFlux Percentages" << std::endl;
	std::cout << std::setw(10) <<"\nCV:\t" << flux_cv << "\n" << std::endl;
	for (unsigned int i = 0; i < parent.size(); i++){
		std::cout << std::setprecision(3) << std::setw(10) << parent.at(i)<< ":\t" << std::setw(10) << flux_parent.at(i) << "\tPercentage\t" << 100 * flux_parent.at(i) / flux_cv << "\%" << std::endl;
	}
	std::cout << "----------------------------------" << std::endl;


	// ++++++++++++++++++++++++++++++++++
	// Save the plots as pdfs in the plots folder
	// ++++++++++++++++++++++++++++++++++
	// create plots folder if it does not exist
	gSystem->Exec("if [ ! -d \"plots/parent\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir -p plots/parent; fi"); 
	
	c1->Print(Form("plots/parent/parentflux_%s_%s.pdf", horn, mode_char));
	std::cout << "\n"<< std::endl;
	
	
// gSystem->Exit(0);

} // end of main




