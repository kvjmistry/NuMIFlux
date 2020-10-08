/*
This script was intended to plot the beamline variations unwrapped 2d flux and 
also plot their ratio for inspection.

To run this script execute the command root -l plot_uboone_beamline.C(<nu flavour>)
The nuetrino flvaour to plot is currently hardcoded into the code,
just find the variable mode and change this to the desired neutrino flavour.

This file depends on the functions.h script so make
sure this file is included in the same  directory.

*/


// Function to plot the beamline uncertainties compared to the CV
#include "functions.h"
// ------------------------------------------------------------------------------------------------------------
// Function to unwrap the histogram
void UnwrapTH2D(TH2D* &hFlux2d,TH1D* &h_unwrap, double POT){

	// Get bins
	const int nBinsEnu = hFlux2d->GetXaxis()->GetNbins(); // Enu
	const int nBinsTh = hFlux2d->GetYaxis()->GetNbins(); // Theta

	// Normalise 2d hist by bin area / deg / GeV
	// Loop over rows
	for (int i=1; i<nBinsEnu+1; i++) {

		// Loop over columns
		for (int j=1; j<nBinsTh+1; j++) {
			hFlux2d->SetBinContent(i,j, hFlux2d->GetBinContent(i, j)/ ( hFlux2d->GetXaxis()->GetBinWidth(i) * hFlux2d->GetYaxis()->GetBinWidth(j) ));
		}
	} 

	hFlux2d->Scale((6.0e20)/ (POT*1.0e4)); // scale to right POT and m2
	
	h_unwrap = new TH1D("", "",nBinsEnu*nBinsTh, 0, nBinsEnu*nBinsTh );
	h_unwrap->SetDirectory(0);
	int counter{0};
	for (int i=1; i<nBinsEnu+1; i++) { // Loop over rows
		for (int j=1; j<nBinsTh+1; j++){// Loop over columns
			counter++;
			h_unwrap->SetBinContent(counter, hFlux2d->GetBinContent(i , j)  );
		}
	}

}
// ------------------------------------------------------------------------------------------------------------
// Legend and draw specifiers for beamine variations
void DrawSpecifiers(TH1D* &hist, TLegend* legend, std::string param, const char* specifier){
	// ----------------------
	//    Draw Specifiers
	// ----------------------

	hist->SetLineWidth(1);
	hist->GetXaxis()->SetLabelSize(0.05);
	hist->GetXaxis()->SetTitleSize(0.05);
	hist->GetYaxis()->SetLabelSize(0.05);
	hist->GetYaxis()->SetTitleSize(0.05);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.12);

	if (param == "CV"){
	hist->SetLineColor(kBlack);
	hist->SetLineStyle(1);
	legend->AddEntry(hist, "CV", "l");
	hist->Draw(specifier);
	} 
	else if  (param == "HP"){
		hist->SetLineColor(kMagenta+2);
		hist->SetLineStyle(1);
		legend->AddEntry(hist, "Hadron Prod.", "l");
	}
	else if  (param == "Horn_p2kA" ){ 
		hist->SetLineColor(30);
		legend->AddEntry(hist, "Horn +2kA", "l");
		hist->SetLineStyle(1);
	}
	else if  (param == "Horn_m2kA"){ 
		hist->SetLineColor(30);
		legend->AddEntry(hist,"Horn -2kA", "l");
		hist->SetLineStyle(2);
		hist->Draw(specifier);
	}
	else if  (param == "Horn1_x_p3mm"){
		hist->SetLineColor(28);
		legend->AddEntry(hist, "Horn1 x +3mm", "l");
		hist->SetLineStyle(1);
		
	}
	else if  (param == "Horm1_x_m3mm" ){
		hist->SetLineColor(28);
		legend->AddEntry(hist, "Horm1 x -3mm", "l");
		hist->SetLineStyle(2);
	}
	else if  (param == "Horn1_y_p3mm"){
		hist->SetLineColor(1001);
		legend->AddEntry(hist, "Horn1 y +3mm", "l");
		hist->SetLineStyle(1);
	}
	else if  (param == "Horn1_y_m3mm"){
		hist->SetLineColor(1001);
		legend->AddEntry(hist, "Horn1 y -3mm", "l");
		hist->SetLineStyle(2);
	}
	else if  (param == "Beam_spot_1_1mm"){
		hist->SetLineColor(kBlue+1);
		legend->AddEntry(hist, "Beam spot 1.1mm", "l");
		hist->SetLineStyle(1);
	}
	else if  (param == "Nominal"){
		hist->SetLineColor(kBlue+1);
		legend->AddEntry(hist, "Nominal", "l");
		hist->SetLineStyle(3);
	}
	else if  (param == "Beam_spot_1_5mm"){
		hist->SetLineColor(kBlue+1);
		legend->AddEntry(hist, "Beam spot 1.5mm", "l");
		hist->SetLineStyle(2);
	}
	else if  (param == "Horn2_x_p3mm"){
		hist->SetLineColor(kOrange+1);
		legend->AddEntry(hist, "Horn2 x +3mm", "l");
		hist->SetLineStyle(1);
	}
	else if  (param == "Horm2_x_m3mm"){
		hist->SetLineColor(kOrange+1);
		legend->AddEntry(hist, "Horm2 x -3mm", "l");
		hist->SetLineStyle(2);
	}
	else if  (param == "Horn2_y_p3mm"){
		hist->SetLineColor(kMagenta);
		legend->AddEntry(hist, "Horn2 y +3mm", "l");
		hist->SetLineStyle(1);
	}
	else if  (param == "Horn2_y_m3mm"){
		hist->SetLineColor(kMagenta);
		legend->AddEntry(hist, "Horn2 y -3mm", "l");
		hist->SetLineStyle(2);
	}
	else if  (param == "Horns_0mm_water"){
		hist->SetLineColor(kSpring-7);
		legend->AddEntry(hist, "Horns 0mm water", "l");
		hist->SetLineStyle(1);
	}
	else if  (param == "Horns_2mm_water"){
		hist->SetLineColor(kSpring-7);
		legend->AddEntry(hist, "Horns 2mm water", "l");
		hist->SetLineStyle(2);
	}
	else if  (param == "Old_Horn" ){ 
		hist->SetLineColor(kSpring+9);
		legend->AddEntry(hist, "Old Horn", "l");
		hist->SetLineStyle(1);
	}
	else if  (param == "Beam_shift_x_p1mm" ){
		hist->SetLineColor(36);
		legend->AddEntry(hist,"Beam shift x +1mm", "l");
		hist->SetLineStyle(1);
	}
	else if  (param == "Beam_shift_x_m1mm" ){
		hist->SetLineColor(36);
		legend->AddEntry(hist,"Beam shift x -1mm", "l");
		hist->SetLineStyle(2);
	}
	else if  (param == "Beam_shift_y_p1mm"){
		hist->SetLineColor(42);
		legend->AddEntry(hist, "Beam shift y +1mm", "l");
		hist->SetLineStyle(1);
	}
	else if  (param == "Beam_shift_y_m1mm"){
		hist->SetLineColor(42);
		legend->AddEntry(hist,"Beam shift y -1mm", "l");
		hist->SetLineStyle(2);
	}
	else if  (param == "Target_z_p7mm"){
		hist->SetLineColor(kOrange+10);
		legend->AddEntry(hist, "Target z +7mm", "l");
		hist->SetLineStyle(1);
	}
	else if  (param == "Target_z_m7mm"){
		hist->SetLineColor(kOrange+10);
		legend->AddEntry(hist, "Target z -7mm", "l");
		hist->SetLineStyle(2);
	}
	else if  (param == "Decay_pipe_Bfield"){
		hist->SetLineColor(50);
		legend->AddEntry(hist, "Decay pipe B field", "l");
		hist->SetLineStyle(1);
	}
	else if  (param == "Horn1_refined_descr"){
		hist->SetLineColor(kMagenta-10);
		legend->AddEntry(hist, "Horn1 refined descr.", "l");
		hist->SetLineStyle(1);
	}
	else if  (param == "Beam_divergence_54urad"){
		hist->SetLineColor(kTeal+6);
		legend->AddEntry(hist, "Beam divergence 54#murad", "l");
		hist->SetLineStyle(1);
	}
	else return;

	hist->Draw(specifier);
}
// ------------------------------------------------------------------------------------------------------------
void DivideHists(TH1D* hCV, TH1D* hUniv, TH1D* &h_1D){

	for (unsigned int i = 1; i < hCV->GetNbinsX()+1;i++){

		if (hUniv->GetBinContent(i) == 0) h_1D->SetBinContent(i, 0);
		else h_1D->SetBinContent(i, hCV->GetBinContent(i) / hUniv->GetBinContent(i));
		h_1D->SetBinError(i, 0.000000001);

	}
}

// ------------------------------------------------------------------------------------------------------------
void plot_beamline_flux(const char* mode, const char * horn){
	gStyle->SetOptStat(0); // say no to stats box

	std::cout << "Using horn mode: " << std::string(horn) << std::endl;

	std::vector<std::string> params = { // A vector with the variations NEW ONES with no threshold
		"CV",                
		"Horn_p2kA",         "Horn_m2kA",
		"Horn1_x_p3mm",      "Horm1_x_m3mm",
		"Horn1_y_p3mm",      "Horn1_y_m3mm",
		"Beam_spot_1_1mm",   "Beam_spot_1_5mm",
		"Horn2_x_p3mm",      "Horm2_x_m3mm",
		"Horn2_y_p3mm",      "Horn2_y_m3mm",
		"Horns_0mm_water",   "Horns_2mm_water",
		"Beam_shift_x_p1mm", "Beam_shift_x_m1mm",
		"Beam_shift_y_p1mm", "Beam_shift_y_m1mm",
		"Target_z_p7mm",     "Target_z_m7mm",
		"Horn1_refined_descr",
		"Decay_pipe_Bfield",
		"Old_Horn"
		};
	
	if (!strcmp(mode,"nue") || !strcmp(mode,"numu")) gStyle->SetTitleH(0.1);
	else gStyle->SetTitleH(0.07);



	bool boolfile, boolhist;
	TFile *f;
	
	// 1D Flux
	TH1D *h_1D_clone;

	std::vector<TH1D*> h_1D;
	h_1D.resize(params.size());

	//2D Flux
	TH2D *h_2D;
	TH1D *h_unwrap, *h_unwrap_clone, *h_2D_unwrap_divided;
	
	double POT;
	TCanvas *c_beamline = new TCanvas();
	TCanvas *c_beamline_ratio = new TCanvas();
	
	TCanvas *c_beamline_1D = new TCanvas();
	TCanvas *c_beamline_ratio_1D = new TCanvas();


	// c_beamline->SetWindowSize(1000, 1000);
	TLegend *lFlux = new TLegend(0.8, 0.10, 1.0, 0.91);
	lFlux->SetNColumns(1);
	lFlux->SetBorderSize(0);
	lFlux->SetFillStyle(0);
	lFlux->SetTextFont(62);

	TLegend *lFlux_ratio = new TLegend(0.8, 0.10, 1.0, 0.91);
	lFlux_ratio->SetNColumns(1);
	lFlux_ratio->SetBorderSize(0);
	lFlux_ratio->SetFillStyle(0);
	lFlux_ratio->SetTextFont(62); 

	TLegend *lFlux_1D = new TLegend(0.8, 0.10, 1.0, 0.91);
	lFlux_1D->SetNColumns(1);
	lFlux_1D->SetBorderSize(0);
	lFlux_1D->SetFillStyle(0);
	lFlux_1D->SetTextFont(62);
	
	TLegend *lFlux_ratio_1D = new TLegend(0.8, 0.10, 1.0, 0.91);
	lFlux_ratio_1D->SetNColumns(1);
	lFlux_ratio_1D->SetBorderSize(0);
	lFlux_ratio_1D->SetFillStyle(0);
	lFlux_ratio_1D->SetTextFont(62); 

	const char* mode_title;

	if (strncmp("numu", mode, 4) == 0)		mode_title = "#nu_{#mu}";
	if (strncmp("nue", mode, 3) == 0)		mode_title = "#nu_{e}";
	if (strncmp("numubar", mode, 7) == 0)	mode_title = "#bar{#nu_{#mu}}";
	if (strncmp("nuebar", mode, 6) == 0)	mode_title = "#bar{#nu_{e}}";

	// ------------------------------------------------------------------------------------------------------------
	// CV
	if (std::string(horn) == "fhc"){
		boolfile  = GetFile(f , "/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold_v46/FHC/output_uboone_fhc_run0_set1.root"); 
	}
	else { 
		boolfile  = GetFile(f , "/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold_v46/RHC/output_uboone_rhc_run0_set1.root"); 
	}
	
	
	if (boolfile == false) gSystem->Exit(0); // Most up to date version of CV
	
	boolhist  = GetHist(f, h_2D, Form("%s/Detsmear/%s_CV_AV_TPC_2D", mode, mode)); 
	if (boolhist == false) gSystem->Exit(0); // CV 2D Histogram
	
	boolhist  = GetHist(f, h_1D.at(0), Form("%s/Detsmear/%s_CV_AV_TPC", mode, mode));
	if (boolhist == false) gSystem->Exit(0); // CV 1D Histogram
	
	POT = GetPOT(f);

	double flux_cv  = h_1D.at(0)->Integral(0,  h_1D.at(0)->GetNbinsX()+1); // Get the CV Flux
	std::vector<double> beamline_flux(params.size(), 0);
	beamline_flux.at(0) = flux_cv * (6.0e20)/ (POT*1.0e4);

	// 2D
	// Get POT and Unwrap histogram
	UnwrapTH2D( h_2D, h_unwrap, POT );
	h_unwrap_clone = (TH1D*) h_unwrap->Clone("h_unwrap_clone");
	h_unwrap->SetTitle(Form("%s;Bin index; #nu / 6 #times 10^{20} POT / GeV / deg / cm^{2}", mode_title ));
	c_beamline->cd();
	DrawSpecifiers(h_unwrap, lFlux, "CV","hist,same");
	Draw_Nu_Mode(c_beamline, horn); // Draw FHC Mode/RHC Mode Text

	//1D
	c_beamline_1D->cd();
	Normalise(h_1D.at(0));
	h_1D.at(0)->Scale((6.0e20)/ (POT*1.0e4)); // scale to right POT and m2
	// h_1D.at(0)->Rebin(10);
	h_1D.at(0)->GetXaxis()->SetRangeUser(0, 6);
	h_1D.at(0)->SetTitle(Form("%s;Energy [GeV];#nu / 6 #times 10^{20} POT / GeV / cm^{2}", mode_title));
	DrawSpecifiers(h_1D.at(0), lFlux_1D, "CV","hist,same");
	Draw_Nu_Mode(c_beamline_1D, horn); // Draw FHC Mode/RHC Mode Text
	h_1D_clone = (TH1D*) h_1D.at(0)->Clone("h_1D_clone");


	// ------------------------------------------------------------------------------------------------------------
	// Loop over the beamline
	for (int i = 1; i < params.size(); i++){

		if (params.at(i) == "Horn1_refined_descr" || params.at(i) == "Old_Horn") continue; // Remove these variations
		
		if (std::string(horn) == "fhc"){
			boolfile  = GetFile(f , Form("/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold_v46/FHC/output_uboone_fhc_run%i.root",i)); 
			// std::cout << boolfile << std::endl;
			if (boolfile == false) continue; // Skip if the file does not exist
		}
		else {
			boolfile  = GetFile(f , Form("/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold_v46/RHC/output_uboone_rhc_run%i.root",i)); 
			if (boolfile == false) continue; // Skip if the file does not exist
		}
		
		if (boolfile == false) continue; // Skip if the file does not exist
		
		boolhist  = GetHist(f, h_2D, Form("%s/Detsmear/%s_CV_AV_TPC_2D", mode, mode)); if (boolhist == false) gSystem->Exit(0);
		boolhist  = GetHist(f, h_1D.at(i), Form("%s/Detsmear/%s_CV_AV_TPC", mode, mode)); if (boolhist == false) gSystem->Exit(0);

		beamline_flux.at(i) = h_1D.at(i)->Integral(0,  h_1D.at(i)->GetNbinsX()+1) * (6.0e20)/ (POT*1.0e4); // Get the beamline Flux

		// Get POT and Unwrap histogram
		std::cout << params[i] << " ";
		POT = GetPOT(f);
		UnwrapTH2D( h_2D, h_unwrap, POT );

		h_2D_unwrap_divided = (TH1D*) h_unwrap_clone->Clone("h_2D_unwrap_divided");
		DivideHists(h_unwrap_clone, h_unwrap, h_2D_unwrap_divided); // CV_clone, Beamline uni i, out

		// Normal FLux
		c_beamline->cd();
		gPad->SetRightMargin(0.2);
		DrawSpecifiers(h_unwrap, lFlux, params.at(i),"hist,same"); // Undivided
	
		// Ratio
		c_beamline_ratio->cd();
		h_2D_unwrap_divided->SetMarkerStyle(7);
		h_2D_unwrap_divided->GetYaxis()->SetRangeUser(0.0, 2);
		gPad->SetRightMargin(0.2);
		h_2D_unwrap_divided->SetTitle(Form("%s;Bin index; Ratio Beamline/CV", mode_title));

		for (unsigned int i=1; i<h_2D_unwrap_divided->GetNbinsX()+1; i++ ) {if (h_2D_unwrap_divided->GetBinContent(i) == 0) h_2D_unwrap_divided->SetBinContent(i, 1);} // Get zero bins to 1 
		// DrawSpecifiers(h_2D_unwrap_divided, lFlux_ratio, params.at(i),"E1,same");
		DrawSpecifiers(h_2D_unwrap_divided, lFlux_ratio, params.at(i),"hist,same"); // Divided flux

		c_beamline_1D->cd();
		gPad->SetRightMargin(0.2);
		h_1D.at(i)->SetDirectory(0);
		Normalise(h_1D.at(i));
		// h_1D.at(i)->Rebin(10);
		h_1D.at(i)->Scale((6.0e20)/ (POT*1.0e4)); // scale to right POT and m2
		DrawSpecifiers(h_1D.at(i), lFlux_1D, params.at(i), "hist,same");

		// Divide histograms and draw again
		TH1D* clone = (TH1D*) h_1D.at(i)->Clone("clone");
		c_beamline_ratio_1D->cd();
		gPad->SetRightMargin(0.2);
		clone->Divide(h_1D_clone);
		clone->GetXaxis()->SetRangeUser(0, 4);
		clone->GetYaxis()->SetRangeUser(0.6, 1.4);
		clone->SetTitle(Form("%s;Energy [GeV];Ratio with CV", mode_title));
		DrawSpecifiers(clone, lFlux_ratio_1D, params.at(i), "E1,same");
		
	}
	
	// ------------------------------------------------------------------------------------------------------------
	c_beamline->cd();
	lFlux->Draw();
	gPad->SetLogy();
	gPad->Update();
	
	c_beamline_ratio->cd();
	lFlux_ratio->Draw();
	TLine* flat = new TLine(0, 1, 42, 1);
	flat->SetLineStyle(7);
	flat->Draw();
	Draw_Nu_Mode(c_beamline_ratio, horn); // Draw FHC Mode/RHC Mode Text



	c_beamline_ratio_1D->cd();
	lFlux_ratio_1D->Draw();
	Draw_Nu_Mode(c_beamline_ratio_1D, horn); // Draw FHC Mode/RHC Mode Text


	c_beamline_1D->cd();
	lFlux_1D->Draw();
	gPad->SetLogy();
	gPad->Update();

	// Now print the percentages for the flux
	std::cout << std::setw(10);
	std::cout << "----------------------------------"  <<  std::endl;
	std::cout << "\nFlux Percentages" << std::endl;
	std::cout << std::setw(10) <<"\nCV:\t" << beamline_flux.at(0) << "\n" << std::endl;
	for (unsigned int i = 1; i < params.size(); i++){
		std::cout << std::setprecision(3) << std::setw(15) << params.at(i)<< ":\t" << std::setw(10) << beamline_flux.at(i) << "\tPercentage\t" << 100 * beamline_flux.at(i) / beamline_flux.at(0) - 100<< "\%" << std::endl;
	}
	std::cout << "----------------------------------" << std::endl;
	
	std::vector<std::string> params_tidy = { // A vector with the variations NEW ONES with no threshold tidied up for plot
		"CV",                
		"Horn +2kA",         "Horn -2kA",
		"Horn1 x +3mm",      "Horm1 x m3mm",
		"Horn1 y +3mm",      "Horn1 y m3mm",
		"Beam spot -2mm",    "Beam spot +2mm",
		"Horn2 x +3mm",      "Horm2 x -3mm",
		"Horn2 y +3mm",      "Horn2 y -3mm",
		"Horns 0mm water",   "Horns 2mm water",
		"Beam shift x +1mm", "Beam shift x -1mm",
		"Beam shift y +1mm", "Beam shift y -1mm",
		"Target z +7mm",     "Target z -7mm",
		"Horn1 refined descr.",
		"Decay pipe B field",
		"Old Horn"
		};
	// Now make a plot of the percentage changes from the CV for Integrated Flux
	TH1D *hBeamline_flux= new TH1D("Beamline",Form("%s Beamline Integrated Flux", mode_title), params.size(), 0, params.size());
	for (unsigned int i = 1; i < params.size(); i++){
		if (beamline_flux.at(i) == 0) {
			hBeamline_flux->Fill(params_tidy[i].c_str(),0);
			continue;
		}
		hBeamline_flux->Fill(params_tidy[i].c_str(), 100 * beamline_flux.at(i) / beamline_flux.at(0) - 100);

	}
	TCanvas* c_beamline_flux = new TCanvas();
	gStyle->SetTitleH(0.04);
	hBeamline_flux->SetLineColor(kViolet-6);
	hBeamline_flux->SetLineWidth(3);
	hBeamline_flux->GetYaxis()->SetTitle("Percentage Change From CV %");
	hBeamline_flux->LabelsOption("v");
	gPad->SetBottomMargin(0.33);

	hBeamline_flux->GetXaxis()->SetLabelSize(0.05);
	hBeamline_flux->GetXaxis()->SetTitleSize(0.05);
	hBeamline_flux->GetYaxis()->SetLabelSize(0.05);
	hBeamline_flux->GetYaxis()->SetTitleSize(0.05);
	hBeamline_flux->SetMarkerSize(1.8);
	gPad->SetLeftMargin(0.15);

	hBeamline_flux->Draw("hist");
	c_beamline_flux ->SetGrid();
	Draw_Nu_Mode(c_beamline_flux, horn); // Draw FHC Mode/RHC Mode Text


	// ++++++++++++++++++++++++++++++++++
	// Save the plots as pdfs in the plots folder
	// ++++++++++++++++++++++++++++++++++
	// create plots folder if it does not exist
	gSystem->Exec("if [ ! -d \"plots/beamline\" ]; then echo \"\nBeamline folder does not exist... creating\"; mkdir -p plots/beamline; fi"); 
	
	c_beamline->Print(Form("plots/beamline/%s_Beamline_2D_unwrapped_Flux_%s.pdf",mode, horn));
	c_beamline_1D->Print(Form("plots/beamline/%s_Beamline_1D_Flux_%s.pdf",mode, horn));

	c_beamline_ratio->Print(Form("plots/beamline/%s_Beamline_2D_unwrapped_Flux_ratio_%s.pdf",mode, horn));
	c_beamline_ratio_1D->Print(Form("plots/beamline/%s_Beamline_1D_Flux_ratio_%s.pdf",mode, horn));

	c_beamline_flux->Print(Form("plots/beamline/%s_Beamline_Integrated_Flux_Change_%s.pdf",mode, horn));

} // End
