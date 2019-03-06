/*
 * NuMI Validation Plotting
 *
 * Plots each individual weighting mode instead of a single one
 * does not plot the correlation matrix for each individual mode to speed up the time
 * 
 * A. Mastbaum <mastbaum@uchicago.edu> 2018/11
 * Modified by K. Mistry 12/18
*/

#include "plot_comp_functions.h"

// Main
void plot_comparison_all( TString mipp, TString inputfile, TString prodmode, TString wplot, TString mode) { // (mippon/mippoff, input, Product/noThinKaon etc. numu/nue)
	gStyle->SetOptStat(0); // say no to stats box
	bool overwrite_errors{true};

	std::vector<std::string> inputmode = loopdir(inputfile, mode); // Grab the names

	// Pre-decleaations
	TH1D* hCV_Flux;
	TH1D* hNOvA_CV_Flux;
	TString Getmode, Gethist_TPC, Getflux, Cov_names, Err_names;
	TFile* f1 = TFile::Open(inputfile);

	// Select neutrino type to run with 
	switch (return_mode(mode)){
		case enumu:
			std::cout << "\nUsing NuMu Mode!\n" << std::endl;
			Getmode = "numu"; 												// Folder name
			Gethist_TPC = "numu/numu_CV_AV_TPC";									// AV in TPC flux prediction
			Getflux = "flux_numu";											// CV flux from NOvA
			Cov_names = "numu/%s/Active_TPC_Volume/numu_%s_Uni_%i_AV_TPC";  // Covariance matrix names
			Err_names = "fractional_uncertainty_numu";						// NOvA fractional uncertainties plot
			break;

		case enue:
			std::cout << "\nUsing Nue Mode!\n" << std::endl;
			Getmode = "nue";
			Gethist_TPC = "nue/nue_CV_AV_TPC";
			Getflux = "flux_nue";
			Cov_names = "nue/%s/Active_TPC_Volume/nue_%s_Uni_%i_AV_TPC";
			Err_names = "fractional_uncertainty_nue";
			break;

		case enumubar:
			std::cout << "\nUsing NuMubar Mode!\n" << std::endl;
			Getmode = "numubar";
			Gethist_TPC = "numubar/numubar_CV_AV_TPC";
			Getflux = "flux_numubar";
			Cov_names = "numubar/%s/Active_TPC_Volume/numubar_%s_Uni_%i_AV_TPC";
			Err_names = "fractional_uncertainty_numubar";
			break;

		case enuebar:
			std::cout << "\nUsing Nuebar Mode!\n" << std::endl;
			Getmode = "nuebar";
			Gethist_TPC = "nuebar/nuebar_CV_AV_TPC";
			Getflux = "flux_nuebar";
			Cov_names = "nuebar/%s/Active_TPC_Volume/nuebar_%s_Uni_%i_AV_TPC";
			Err_names = "fractional_uncertainty_nuebar";
			break;

	}
	
	// ------------------------------------------------------------------------------------------------------------
	// Get the POT in the file
	// ------------------------------------------------------------------------------------------------------------
	double fPOT{0};
	TTree* TPOT;
	bool boolPOT = GetTree(f1, TPOT, "POT"); // Check if got ttree properly

	if (boolPOT == false)  fPOT = 2.5e8; // default to this POT
	else {
		TPOT->SetBranchAddress("POT", &fPOT); // Get the POT
		TPOT->GetEntry(0);
		std::cout << "TOTAL POT READ IN:\t" << fPOT << std::endl;
	}
	// ------------------------------------------------------------------------------------------------------------
	// CV Flux
	// ------------------------------------------------------------------------------------------------------------
	TCanvas* c1 = new TCanvas();

	// Check if sucessfully got histo
	bool boolhist = GetHist(f1, hCV_Flux, Gethist_TPC); if (boolhist == false) gSystem->Exit(0);

	hCV_Flux->SetDirectory(0);
	
	// Normalise flux by bin width (gives a flux/E [GeV])
	for (int i=1;i<hCV_Flux->GetNbinsX()+1;i++) {
		hCV_Flux->SetBinContent(i, hCV_Flux->GetBinContent(i)/hCV_Flux->GetBinWidth(i));
	}

	// Now get the NOvA CV Flux
	TFile* f2 = TFile::Open("/uboone/app/users/kmistry/PPFX/numi-validation/nova_flux/FHC_Flux_NOvA_ND_2017.root");
	hNOvA_CV_Flux = (TH1D*) f2->Get(Getflux);

	// Check if sucessfully got histo
	boolhist = GetHist(f2, hNOvA_CV_Flux, Getflux); if (boolhist == false) gSystem->Exit(0);

	TH1D* horig = (TH1D*) hCV_Flux->Clone("horig");   // Clone for multiverse calculations

	// Draw Specifiers
	hCV_Flux->Sumw2();
	hCV_Flux->SetLineColor(kRed);
	hCV_Flux->SetLineWidth(2);
	hCV_Flux->GetYaxis()->SetTitle(hNOvA_CV_Flux->GetYaxis()->GetTitle());
	if (mode == "numu")  hCV_Flux->SetTitle("#nu_{#mu}");
	if (mode == "nue")   hCV_Flux->SetTitle("#nu_{e}");
	if (mode == "numubar")  hCV_Flux->SetTitle("#bar{#nu_{#mu}}");
	if (mode == "nuebar")   hCV_Flux->SetTitle("#bar{#nu_{e}}");
	hCV_Flux->Draw("");
	hNOvA_CV_Flux->SetLineColor(kBlue);
	hNOvA_CV_Flux->SetLineWidth(2);
	hNOvA_CV_Flux->Draw("same");
	
	// Normalisation --  Need to divide by area front face is 12.39 m2, looks closer to 14.6, window area is 74.99992402
	hCV_Flux->Scale(1.0e6/(fPOT * 74.99992402)); 
	std::cout << "norm factor:\t" << hCV_Flux->Integral(2,-1)/hNOvA_CV_Flux->Integral(2,-1) << std::endl;
	hCV_Flux->Scale(hNOvA_CV_Flux->Integral(3,-1)/hCV_Flux->Integral(3,-1));
	
	// Legend
	TLegend* l = new TLegend(0.5, 0.6, 0.85, 0.8);
	l->AddEntry(hNOvA_CV_Flux, "NOvA","l");
	l->AddEntry(hCV_Flux, "Our Prediction (stat+sys)","l");
	l->Draw();
	l->SetBorderSize(0);
	l->SetFillStyle(0);
	gPad->SetLogy();
	gPad->Update();

	// ------------------------------------------------------------------------------------------------------------
	// Ratio of ours to NOvA
	// ------------------------------------------------------------------------------------------------------------
	TCanvas* c2 = new TCanvas();
	TH1D* hratio = (TH1D*) hCV_Flux->Clone("hratio"); // Clone for ratio plot
	hratio->SetDirectory(0);
	hratio->Divide(hNOvA_CV_Flux);
	hratio->SetLineColor(kBlack);
	hratio->SetMarkerStyle(7);
	hratio->Draw("e1");
	hratio->SetTitle(";E_{#nu} [GeV];Ratio to NOvA");
	hratio->GetYaxis()->SetRangeUser(0.85, 1.15);
	TLine* flat = new TLine(0, 1, 20, 1);
	flat->SetLineStyle(7);
	flat->Draw();

	// ------------------------------------------------------------------------------------------------------------
	// Correlations, Covariance & uncertainties
	// ------------------------------------------------------------------------------------------------------------

	f1->cd();

	// Varables and histograms
	const int nbins = hCV_Flux->GetNbinsX();
	double* edges = new double[nbins+1];
	TCanvas* c3 = new TCanvas();
	TCanvas* c4 = new TCanvas();
	TH1D* herr2; 				// Nova uncertainties
	std::vector<TH2D*> cov;		// Covariance
	std::vector<TH2D*> cor; 	// Correlation
	std::vector<TH1D*> herr ;	// Fractional Uncertenties 

	// Set bin widths to be the same as NOvA
	for (int i=1; i<nbins+1; i++) {
		edges[i-1] = hCV_Flux->GetBinLowEdge(i);
	}

	// Get bin Edges
	edges[nbins] = hCV_Flux->GetBinLowEdge(nbins-1) + 2 * (hCV_Flux->GetBinWidth(nbins-1));
	
	// Legened
	TLegend* lfrac = new TLegend(0.5, 0.65, 0.9, 0.9);
	lfrac->SetNColumns(3);
	lfrac->SetBorderSize(0);
	lfrac->SetFillStyle(0);
	lfrac->SetTextFont(62); 

	// Resize
	cor.resize(inputmode.size());
	cov.resize(inputmode.size());
	herr.resize(inputmode.size());
	
	// Create histograms
	for (unsigned int l = 0; l < inputmode.size(); l++){
		cor[l]  = new TH2D(Form("%s_cor",inputmode[l].c_str()), ";E_{#nu} (GeV);E_{#nu} (GeV)", nbins, edges, nbins, edges);
		cov[l]  = new TH2D(Form("%s_cov",inputmode[l].c_str()), ";E_{#nu} (GeV);E_{#nu} (GeV)", nbins, edges, nbins, edges);
		herr[l] = new TH1D(Form("%s_herr",inputmode[l].c_str()),";E_{#nu} (GeV);Fractional Uncertainty", nbins, edges);
	}

	// Loop over all input modes, get cov matrix and then get fractional uncertainties
	for (unsigned int l = 0; l < inputmode.size(); l++){
		
		// +++++++++++++++++
		// Covariance matrix
		// +++++++++++++++++
		CalcCovariance(inputmode[l], Cov_names, f1, cov[l], horig, nbins);

		double cii{0}, cjj{0}, n{1}, cor_bini{0}, horig_cont{0};
		
		// ++++++++++++++++++++++++++++++++++++++
		// Correlation matrix & fractional errors
		// ++++++++++++++++++++++++++++++++++++++
		cor[l] = (TH2D*) cov[l]->Clone("cor");
		
		// loop over rows
		for (int i=1; i<nbins+1; i++) {
			cii = cov[l]->GetBinContent(i, i);

			// Catch zeros , set to arbitary 1.0
			if (horig->GetBinContent(i) <= 0) horig_cont = 1.0;
			else horig_cont = horig->GetBinContent(i);

			herr[l]->SetBinContent(i, sqrt(cii) / horig_cont);

			// Loop over columns
			for (int j=1; j<nbins+1; j++) {
				cjj = cov[l]->GetBinContent(j, j);
				n = sqrt(cii * cjj);

				// Catch Zeros, set to arbitary 1.0
				if (n == 0) cor_bini=1.0;
				else cor_bini = cov[l]->GetBinContent(i, j) / n;

				cor[l]->SetBinContent(i, j, cor_bini );
			}
		}

		// +++++++++++++++++
		// Plot correlations
		// +++++++++++++++++
		c3->cd();
		
		// Draw MasterWeight Cor plot only
		if (inputmode.size() == 2) { // checks if long or short list
			cor[0]->SetTitle("Correlation Master Weight");
			cor[0]->Draw("colz");
		} 
		else  {
			cor[11]->SetTitle("Correlation Master Weight");
			cor[11]->Draw("colz");
		}

		gStyle->SetPalette(55); // kRainbow

		// Plot fractional errors overlaid with official NOvA plot
		c4->cd();

		// Make the plot
		if (overwrite_errors == false) legDraw(lfrac, herr[l], prodmode, mipp, inputmode[l], mode);
		
		// Set the axes
		herr[l]->GetYaxis()->SetRangeUser(0,0.35);
		// herr[l]->GetYaxis()->SetRangeUser(0,1.75);
		
	}

	f2->cd();

	// Get the nova  histogram
	boolhist = GetHist(f2, herr2, Err_names); if (boolhist == false) gSystem->Exit(0);

	// Set the axes
	if (mode == "numu"){
		herr2->GetYaxis()->SetRangeUser(0,0.35);
	}
	if (mode == "nue"){
		herr2->GetYaxis()->SetRangeUser(0,0.40);
	}
	if (mode == "numubar"){
		herr2->GetYaxis()->SetRangeUser(0,0.65);
	}
	if (mode == "nuebar"){
		herr2->GetYaxis()->SetRangeUser(0,0.65);
	}

	// Draw Nova errors
	herr2->SetLineColor(kBlue);
	herr2->SetLineWidth(2);
	herr2->Draw("same");
	lfrac->AddEntry(herr2, "NOvA","l");

	// Draw the legend
	lfrac->Draw();

	// ------------------------------------------------------------------------------------------------------------
	// Override the errors to use Leos method
	// ------------------------------------------------------------------------------------------------------------
	// Decide if the errors need overwriting
	if (overwrite_errors == true ){
		c4->cd();
		std::cout << "Overwriting the errors" << std::endl;
		for (unsigned int l = 0; l < inputmode.size(); l++){
			// herr[l]->Reset();
			HPUncertainties_Leo(f1, herr[l], inputmode[l], mode);
			if (l == 11)legDraw(lfrac, herr[l], prodmode, mipp, inputmode[l], mode);
		}
		c4->Update();
	}

	// ------------------------------------------------------------------------------------------------------------
	// Create a tfile with the nova HP uncertainties
	// ------------------------------------------------------------------------------------------------------------
	std::string mode_str;
	if (mode == "numu") mode_str = "numu";
	else if (mode == "numubar") mode_str = "numubar";
	else if (mode == "nue") mode_str = "nue";
	else mode_str = "numubar";

	TFile* output = new TFile("nova_numu_hp_uncertainties.root", "RECREATE");
	
	for (int i=0; i<inputmode.size();i++){
		herr[i]->SetOption("hist"); // Overwrite the histo option
		herr[i]->SetName(Form("%s_%s", mode_str.c_str(), inputmode[i].c_str())); // Rename as names got grblled.. who knows why!
		herr[i]->Write();
	}
	herr2->Write();
	output->Close();

	// ------------------------------------------------------------------------------------------------------------
	// Make a plot with the uncertainties with errorbands
	// ------------------------------------------------------------------------------------------------------------
	TCanvas* cband = new TCanvas();
	cband->cd();
	TLegend* leg = new TLegend(0.60,0.70,0.90,0.90);
	DrawErrorBand(f1, mode, leg, "PPFXMaster"); // Plot for masterweight
	leg->Draw();

	// ------------------------------------------------------------------------------------------------------------
	// Make the weight histogram
	// ------------------------------------------------------------------------------------------------------------

	TCanvas* c5 = new TCanvas();
	if (wplot == "wplot"){
	
		TLegend* lwght = new TLegend(0.5, 0.65, 0.9, 0.9);
		lwght->SetNColumns(3);
		lwght->SetBorderSize(0);
		lwght->SetFillStyle(0);
		lwght->SetTextFont(62);

		weight_plots(mode, inputmode, f1, prodmode, mipp, c5, lwght);
	}
	else c5->Close(); // option is off so close the histogram

	double sigma{0}, stat{0}, sys{0};
	// ------------------------------------------------------------------------------------------------------------
	// Update the error to add in the beamline uncertainties
	// ------------------------------------------------------------------------------------------------------------
	// Choose one of these
	BeamlineUncertainties(herr, hCV_Flux, "file",mode); // likely to not be implemented anymore
	// BeamlineUncertainties(herr, hCV_Flux, "stdev",mode);
	// BeamlineUncertainties(herr, hCV_Flux, "quad",mode);

	// ------------------------------------------------------------------------------------------------------------
	// Update the CV flux prediction to include stat+sys errors
	// ------------------------------------------------------------------------------------------------------------
	c1->cd();
	

	std::cout << "bins:\t" << hCV_Flux->GetNbinsX() << std::endl;

	// Loop over the bins
	for (int bin=1; bin<hCV_Flux->GetNbinsX()+1; bin++){
		
		// Get the bin error (stat)
		stat = hCV_Flux->GetBinError(bin);

		// Get the bin error (sys)
		// sys = herr[11]->GetBinContent(bin);
		sys = herr[0]->GetBinContent(bin);
		sys = hCV_Flux->GetBinContent(bin) * sys; 

		// add in quadrature
		sigma =std::sqrt( stat*stat + sys*sys );

		// std::cout << "stat:\t" << stat << "\t" << "sys:\t" <<sys << "\tsigma:\t"<< sigma<< "\tbin content:\t"<< hCV_Flux->GetBinContent(bin) <<  std::endl;

		// Update the error on the plot
		hCV_Flux->SetBinError(bin, sigma);
	}

	// ------------------------------------------------------------------------------------------------------------
	// Take the ratio of novas uncertainties to ours
	// ------------------------------------------------------------------------------------------------------------
	TCanvas* c6;
	TH1D* hratio_sig;
	
	c6 = new TCanvas();
	// Clone for ratio plot
	hratio_sig = (TH1D*) herr[12]->Clone("hratio_sig");
	
	hratio_sig->Divide(herr2);
	hratio_sig->SetLineColor(kBlack);
	hratio_sig->SetMarkerStyle(7);
	hratio_sig->Draw("e1");
	hratio_sig->SetTitle("Ratio of uncertainties to NOvA;E_{#nu} (GeV);Ratio");
	hratio_sig->GetYaxis()->SetRangeUser(0, 3);
	hratio_sig->Draw("hist");
	flat->Draw();

	// Redraw the plot with the new errors
	c1->Update();
	c2->Update();
	c4->Update();

	// ------------------------------------------------------------------------------------------------------------
	// Save the plots as pdfs in the plots folder
	// ------------------------------------------------------------------------------------------------------------
	// create plots folder if it does not exist
	
	gSystem->Exec("if [ ! -d \"plots\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir plots; fi"); 
	
	// using mipp
	if (mipp == "mippon"){
		if (mode == "numu"){ 	
			c1->Print("plots/CV_Flux_Prediction_NuMu_MIPPOn.pdf");
			c2->Print("plots/Ratio_FLux_Prediction_NuMu_MIPPOn.pdf");
			c3->Print("plots/Correlation_Matrix_NuMu_MIPPOn.pdf");
			c4->Print("plots/Fractional_Uncertainties_NuMu_MIPPOn.pdf");
			if (wplot == "wplot") c5->Print("plots/Weightplot_NuMu_MIPPOn.pdf");
			std::cout << "\n"<< std::endl;
		}
		else if (mode == "nue") {
			c1->Print("plots/CV_Flux_Prediction_Nue_MIPPOn.pdf");
			c2->Print("plots/Ratio_FLux_Prediction_Nue_MIPPOn.pdf");
			c3->Print("plots/Correlation_Matrix_Nue_MIPPOn.pdf");
			c4->Print("plots/Fractional_Uncertainties_Nue_MIPPOn.pdf");
			if (wplot == "wplot") c5->Print("plots/Weightplot_Nue_MIPPOn.pdf");
			std::cout << "\n"<< std::endl;

		}
		else if (mode == "numubar"){
			c1->Print("plots/CV_Flux_Prediction_Numubar_MIPPOn.pdf");
			c2->Print("plots/Ratio_FLux_Prediction_Numubar_MIPPOn.pdf");
			c3->Print("plots/Correlation_Matrix_Numubar_MIPPOn.pdf");
			c4->Print("plots/Fractional_Uncertainties_Numubar_MIPPOn.pdf");
			if (wplot == "wplot") c5->Print("plots/Weightplot_NuMubar_MIPPOn.pdf");
			std::cout << "\n"<< std::endl;

		}
		else {
			c1->Print("plots/CV_Flux_Prediction_Nuebar_MIPPOn.pdf");
			c2->Print("plots/Ratio_FLux_Prediction_Nuebar_MIPPOn.pdf");
			c3->Print("plots/Correlation_Matrix_Nuebar_MIPPOn.pdf");
			c4->Print("plots/Fractional_Uncertainties_Nuebar_MIPPOn.pdf");
			if (wplot == "wplot") c5->Print("plots/Weightplot_Nuebar_MIPPOn.pdf");
			std::cout << "\n"<< std::endl;

		}
	}
	// no mipp
	else {
		if (mode == "numu"){ 	
			c1->Print("plots/CV_Flux_Prediction_NuMu_MIPPOff.pdf");
			c2->Print("plots/Ratio_FLux_Prediction_NuMu_MIPPOff.pdf");
			c3->Print("plots/Correlation_Matrix_NuMu_MIPPOff.pdf");
			c4->Print("plots/Fractional_Uncertainties_NuMu_MIPPOff.pdf");
			if (wplot == "wplot") c5->Print("plots/Weightplot_NuMu_MIPPOff.pdf");
			std::cout << "\n"<< std::endl;
		}
		else if (mode == "nue"){
			c1->Print("plots/CV_Flux_Prediction_Nue_MIPPOff.pdf");
			c2->Print("plots/Ratio_FLux_Prediction_Nue_MIPPOff.pdf");
			c3->Print("plots/Correlation_Matrix_Nue_MIPPOff.pdf");
			c4->Print("plots/Fractional_Uncertainties_Nue_MIPPOff.pdf");
			if (wplot == "wplot") c5->Print("plots/Weightplot_Nue_MIPPOff.pdf");
			std::cout << "\n"<< std::endl;

		}
		else if (mode == "numubar"){
			c1->Print("plots/CV_Flux_Prediction_Numubar_MIPPOff.pdf");
			c2->Print("plots/Ratio_FLux_Prediction_Numubar_MIPPOff.pdf");
			c3->Print("plots/Correlation_Matrix_Numubar_MIPPOff.pdf");
			c4->Print("plots/Fractional_Uncertainties_Numubar_MIPPOff.pdf");
			if (wplot == "wplot") c5->Print("plots/Weightplot_NuMubar_MIPPOff.pdf");
			std::cout << "\n"<< std::endl;

		}
		else {
			c1->Print("plots/CV_Flux_Prediction_Nuebar_MIPPOff.pdf");
			c2->Print("plots/Ratio_FLux_Prediction_Nuebar_MIPPOff.pdf");
			c3->Print("plots/Correlation_Matrix_Nuebar_MIPPOff.pdf");
			c4->Print("plots/Fractional_Uncertainties_Nuebar_MIPPOff.pdf");
			if (wplot == "wplot") c5->Print("plots/Weightplot_Nuebar_MIPPOff.pdf");
			std::cout << "\n"<< std::endl;

		}
	}

} // end of main