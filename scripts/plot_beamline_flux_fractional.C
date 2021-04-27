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
// ------------------------------------------------------------------------------------------------------------
// Legend and draw specifiers for beamine variations
void DrawSpecifiers(TH1D* &hist, TLegend* legend, std::string param){
	// ----------------------
	//    Draw Specifiers
	// ----------------------

	hist->SetLineWidth(2);
	hist->GetXaxis()->SetLabelSize(0.05);
	hist->GetXaxis()->SetTitleSize(0.05);
	hist->GetYaxis()->SetLabelSize(0.05);
	hist->GetYaxis()->SetTitleSize(0.05);

	if (param == "Total"){
	hist->SetLineColor(kBlack);
	hist->SetLineStyle(1);
	legend->AddEntry(hist, "Total", "l");
	} 
	else if  (param == "HP"){
		hist->SetLineColor(kMagenta+2);
		hist->SetLineStyle(1);
		legend->AddEntry(hist, "Hadron Prod.", "l");
	}
	else if  (param == "Horn_curr" ){ 
		hist->SetLineColor(30);
		legend->AddEntry(hist, "Horn Current", "l");
		hist->SetLineStyle(1);
	}
	else if  (param == "Horn_m2kA"){ 
		hist->SetLineColor(30);
		legend->AddEntry(hist,"Horn -2kA", "l");
		hist->SetLineStyle(2);
	}
	else if  (param == "Horn1_x"){
		hist->SetLineColor(28);
		legend->AddEntry(hist, "Horn1 x", "l");
		hist->SetLineStyle(1);
		
	}
	else if  (param == "Horm1_x_m3mm" ){
		hist->SetLineColor(28);
		legend->AddEntry(hist, "Horn1 x -3mm", "l");
		hist->SetLineStyle(2);
	}
	else if  (param == "Horn1_y"){
		hist->SetLineColor(1001);
		legend->AddEntry(hist, "Horn1 y", "l");
		hist->SetLineStyle(1);
	}
	else if  (param == "Horn1_y_m3mm"){
		hist->SetLineColor(1001);
		legend->AddEntry(hist, "Horn1 y -3mm", "l");
		hist->SetLineStyle(2);
	}
	else if  (param == "Beam_spot"){
		hist->SetLineColor(kBlue+1);
		legend->AddEntry(hist, "Beam Spot Size", "l");
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
	else if  (param == "Horn2_x"){
		hist->SetLineColor(kOrange+1);
		legend->AddEntry(hist, "Horn2 x", "l");
		hist->SetLineStyle(1);
	}
	else if  (param == "Horm2_x_m3mm"){
		hist->SetLineColor(kOrange+1);
		legend->AddEntry(hist, "Horn2 x", "l");
		hist->SetLineStyle(2);
	}
	else if  (param == "Horn2_y"){
		hist->SetLineColor(kMagenta);
		legend->AddEntry(hist, "Horn2 y", "l");
		hist->SetLineStyle(1);
	}
	else if  (param == "Horn2_y_m3mm"){
		hist->SetLineColor(kMagenta);
		legend->AddEntry(hist, "Horn2 y -3mm", "l");
		hist->SetLineStyle(2);
	}
	else if  (param == "Horn_water"){
		hist->SetLineColor(kSpring-7);
		legend->AddEntry(hist, "Horn Water", "l");
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
	else if  (param == "Beam_shift_x" ){
		hist->SetLineColor(36);
		legend->AddEntry(hist,"Beam shift x", "l");
		hist->SetLineStyle(1);
	}
	else if  (param == "Beam_shift_x_m1mm" ){
		hist->SetLineColor(36);
		legend->AddEntry(hist,"Beam shift x -1mm", "l");
		hist->SetLineStyle(2);
	}
	else if  (param == "Beam_shift_y"){
		hist->SetLineColor(42);
		legend->AddEntry(hist, "Beam shift y", "l");
		hist->SetLineStyle(1);
	}
	else if  (param == "Beam_shift_y_m1mm"){
		hist->SetLineColor(42);
		legend->AddEntry(hist,"Beam shift y -1mm", "l");
		hist->SetLineStyle(2);
	}
	else if  (param == "Target_z"){
		hist->SetLineColor(kOrange+10);
		legend->AddEntry(hist, "Target z", "l");
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

}
// ------------------------------------------------------------------------------------------------------------
void DivideHists(TH1D* hCV, TH1D* hUniv, TH1D* &h_1D){

	for (unsigned int i = 1; i < hCV->GetNbinsX()+1;i++){

		if (hUniv->GetBinContent(i) == 0) h_1D->SetBinContent(i, 0);
		else h_1D->SetBinContent(i, hCV->GetBinContent(i) / hUniv->GetBinContent(i));
		h_1D->SetBinError(i, 0.000000001);

	}
}
void CalcCovariance(std::vector<TH1D*> h_universe, TH1D *h_CV, TH2D *h_cov){

    for (unsigned int uni = 0; uni < h_universe.size(); uni++){
        
        // Loop over the rows
        for (int row = 1; row < h_CV->GetNbinsX()+1; row++){
            
            double uni_row = h_universe.at(uni)->GetBinContent(row);
            double cv_row  = h_CV->GetBinContent(row);

            // Loop over the columns
            for (int col = 1; col < h_CV->GetNbinsX()+1; col++){

                double uni_col = h_universe.at(uni)->GetBinContent(col);
                double cv_col  = h_CV->GetBinContent(col);
                
                double c = (uni_row - cv_row) * (uni_col - cv_col);

                if (uni != h_universe.size()-1)    h_cov->SetBinContent(row, col, h_cov->GetBinContent(row, col) + c ); // Fill with variance 
                else h_cov->SetBinContent(row, col, (h_cov->GetBinContent(row, col) + c) / h_universe.size());          // Fill with variance and divide by nuni
            
            } // end loop over columns

        } // end loop over rows
    
    } // end loop over universes

}

// ------------------------------------------------------------------------------------------------------------
void plot_beamline_flux_fractional(const char* mode, const char * horn){
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
		"Target_z_p7mm",     "Target_z_m7mm"
		};
	
	if (!strcmp(mode,"nue") || !strcmp(mode,"numu")) gStyle->SetTitleH(0.1);
	else gStyle->SetTitleH(0.07);

	bool boolfile, boolhist;
	TFile *f;
	
	std::vector<TH1D*> h_1D;
	h_1D.resize(params.size());
	
	double POT;
	TLegend *lFlux = new TLegend(0.5, 0.65, 0.9, 0.9);
	lFlux->SetNColumns(3);
	lFlux->SetBorderSize(0);
	lFlux->SetFillStyle(0);
	lFlux->SetTextFont(62);

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
	
	TH1D* h_temp;
	boolhist  = GetHist(f, h_temp, Form("%s/Detsmear/%s_CV_AV_TPC", mode, mode));
	if (boolhist == false) gSystem->Exit(0); // CV 1D Histogram

	Double_t xbins[9] = {0.00 ,0.06, 0.125, 0.25, 0.5, 1.00, 1.50, 2.00, 5.00 };
	Double_t xbins_mu[12] = { 0.00, 0.025, 0.03, 0.235 ,0.24, 0.50, 0.75, 1.00, 2.00, 3.00, 6.00, 10.00};
	
	if (std::string(mode) == "nue" || std::string(mode) == "nuebar")
		h_1D.at(0) =  dynamic_cast<TH1D*>(h_temp->Rebin(8, "", xbins));
	else{
		
		h_1D.at(0) =  dynamic_cast<TH1D*>(h_temp->Rebin(11, "", xbins_mu));
	}

	POT = GetPOT(f);

	//1D CV
	Normalise(h_1D.at(0));
	h_1D.at(0)->Scale((1.0)/ (POT*1.0e4)); // scale to right POT and m2
	// h_1D.at(0)->Rebin(10);
	
	if (std::string(mode) == "nue" || std::string(mode) =="nuebar" )
		h_1D.at(0)->GetXaxis()->SetRangeUser(0, 5);
	else
		h_1D.at(0)->GetXaxis()->SetRangeUser(0, 6);


	h_1D.at(0)->SetTitle(Form("%s;Energy [GeV];#nu /  POT / GeV / cm^{2}", mode_title));


	// ------------------------------------------------------------------------------------------------------------
	// Loop over the beamline
	for (int i = 1; i < params.size(); i++){
		
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
		
		TH1D *htemp;
		boolhist  = GetHist(f, htemp, Form("%s/Detsmear/%s_CV_AV_TPC", mode, mode)); if (boolhist == false) gSystem->Exit(0);

		if (std::string(mode) == "nue" || std::string(mode) == "nuebar")
			h_1D.at(i) =  dynamic_cast<TH1D*>(htemp->Rebin(8, "", xbins));
		else
			h_1D.at(i) =  dynamic_cast<TH1D*>(htemp->Rebin(11, "", xbins_mu));

		// Get POT and Unwrap histogram
		std::cout << params[i] << " ";
		POT = GetPOT(f);

		Normalise(h_1D.at(i));
		// h_1D.at(i)->Rebin(10);
		
		h_1D.at(i)->Scale(1.0/ (POT*1.0e4)); // scale to right POT and m2

		// std::cout << h_1D.at(i)->GetNbinsX()<< std::endl;
	
	}



	//  Now we have the beamline variations, calculate the covariance matrices for each one
	std::vector<std::string> params2 = { 
		"Horn_curr",
		"Horn1_x",
		"Horn1_y",
		"Beam_spot",
		"Horn2_x",
		"Horn2_y",
		"Horn_water",
		"Beam_shift_x",
		"Beam_shift_y",
		"Target_z"
		};
	
	std::vector<TH2D*> h_cov_v(params2.size());
	
	int n_bins = h_1D.at(0)->GetNbinsX();

	std::vector<int> index = {1,3,5,7,9,11,13,15,17,19};
	std::vector<std::vector<TH1D*>> h_uni_beamline(index.size());

	// Create the Covariance Matrix
	for (int i = 0; i < index.size(); i++){
		h_cov_v.at(i) = new TH2D(Form("cov_%d", i), "Covariance Matrix ;Bin i; Bin j",  n_bins, 1, n_bins+1, n_bins, 1, n_bins+1);
		h_uni_beamline.at(i).resize(2);
	}

	// Store the beamline variations into pairs
	int counter = 0;
	for (int i = 1; i < params.size(); i++){
		h_uni_beamline.at(counter).at(0) = (TH1D*)h_1D.at(i)->Clone();
		h_uni_beamline.at(counter).at(1) = (TH1D*)h_1D.at(i+1)->Clone();
		i++;
		counter++;
	}

	// Now we can calculate the cov matrix
	for (int i = 0; i < index.size(); i++){
		CalcCovariance(h_uni_beamline.at(i), h_1D.at(0), h_cov_v.at(i));
	}

	// ------------------------------------------------------------------------------------------------------------

	// Now Get the fractional errors
	std::vector<TH1D*> h_err(index.size());
	TH1D* h_err_tot = (TH1D*)h_1D.at(0)->Clone();

	for (int bin = 1; bin < h_err_tot->GetNbinsX()+1; bin++){
		h_err_tot->SetBinContent(bin, 0);
	}

	for (int i = 0; i < index.size(); i++){
		h_err.at(i) = (TH1D*)h_1D.at(0)->Clone();

		for (int bin = 1; bin < h_err.at(i)->GetNbinsX()+1; bin++){
			h_err.at(i)->SetBinContent(bin, std::sqrt(h_cov_v.at(i)->GetBinContent(bin, bin)) /h_1D.at(0)->GetBinContent(bin));
			h_err_tot->SetBinContent(bin, h_err_tot->GetBinContent(bin) + h_cov_v.at(i)->GetBinContent(bin, bin)) ;
		}
		DrawSpecifiers(h_err.at(i), lFlux, params2.at(i));
		
		if (std::string(mode) == "nue" || std::string(mode) == "nuebar")
			h_err.at(i)->GetYaxis()->SetRangeUser(0, 0.15);
		else 
			h_err.at(i)->GetYaxis()->SetRangeUser(0, 0.15); // 0.15
	}

	for (int bin = 1; bin < h_err_tot->GetNbinsX()+1; bin++){
		h_err_tot->SetBinContent(bin,  std::sqrt(h_err_tot->GetBinContent(bin)) / h_1D.at(0)->GetBinContent(bin));
	}
	DrawSpecifiers(h_err_tot, lFlux, "Total");

	// ------------------------------------------------------------------------------------------------------------
	TCanvas *c = new TCanvas();
	IncreaseLabelSize(h_err.at(0));
	h_err.at(0)->GetYaxis()->SetTitle("Fractional Uncertainty");

	for (int i = 0; i < index.size(); i++){
		h_err.at(i)->Draw("hist,same");
	}
	lFlux->Draw();

	h_err_tot->Draw("hist,same");

	Draw_Nu_Mode(c, horn); // Draw FHC Mode/RHC Mode Text

	
	

	// ++++++++++++++++++++++++++++++++++
	// Save the plots as pdfs in the plots folder
	// ++++++++++++++++++++++++++++++++++
	// create plots folder if it does not exist
	gSystem->Exec("if [ ! -d \"plots/beamline\" ]; then echo \"\nBeamline folder does not exist... creating\"; mkdir -p plots/beamline; fi"); 
	
	c->Print(Form("plots/beamline/%s_Beamline_fractional_uncertainty_%s.pdf",mode, horn));

} // End
