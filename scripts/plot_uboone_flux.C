/**
 * NuMI Flux at uboone Plotting
 *
 * Plots each individual weighting mode instead of a single one
 * does not plot the correlation matrix for each individual mode to speed up the time
 * 
 * A. Mastbaum <mastbaum@uchicago.edu> 2018/11
 * Modified by K. Mistry 12/18
 */

#include "plot_comp_functions.h"

// ----------------------------------------------------------------------------
// Main
void plot_uboone_flux( TString mipp, TString inputfile, TString prodmode, TString wplot, TString mode) { // (mippon/mippoff, input, Product/noThinKaon etc. numu/nue)
	gStyle->SetOptStat(0); // say no to stats box
	bool overwrite_errors{false};
	// bool overwrite_errors{true};

	std::vector<std::string> inputmode = loopdir(inputfile, mode); // Grab the names of the input reweighters

	// Pre declare variables
	TString Getmode, Gethist_TPC, Getflux, Cov_names, g_simp_names;
	TH1D* hCV_Flux;
	TH1D*  h_g_simp; 
	TFile* f_gsimple;
	TFile* f1 = TFile::Open(inputfile);

	// Select neutrino type to run with 
	switch (return_mode(mode)){
		case enumu:
			std::cout << "\nUsing NuMu Mode!\n" << std::endl;
			Getmode = "numu"; 												// Folder name
			Gethist_TPC = "numu/numu_CV_AV_TPC";							// AV in TPC flux prediction
			Getflux = "flux_numu";											// CV flux from NOvA
			Cov_names = "numu/%s/Active_TPC_Volume/numu_%s_Uni_%i_AV_TPC";  // Covariance matrix names
			g_simp_names = "numuFluxHisto";									// G simple files
			break;

		case enue:
			std::cout << "\nUsing Nue Mode!\n" << std::endl;
			Getmode = "nue";
			Gethist_TPC = "nue/nue_CV_AV_TPC";
			Getflux = "flux_nue";
			Cov_names = "nue/%s/Active_TPC_Volume/nue_%s_Uni_%i_AV_TPC";
			g_simp_names = "nueFluxHisto";
			break;

		case enumubar:
			std::cout << "\nUsing NuMubar Mode!\n" << std::endl;
			Getmode = "numubar";
			Gethist_TPC = "numubar/numubar_CV_AV_TPC";
			Getflux = "flux_numubar";
			Cov_names = "numubar/%s/Active_TPC_Volume/numubar_%s_Uni_%i_AV_TPC";
			g_simp_names = "anumuFluxHisto";
			break;

		case enuebar:
			std::cout << "\nUsing Nuebar Mode!\n" << std::endl;
			Getmode = "nuebar";
			Gethist_TPC = "nuebar/nuebar_CV_AV_TPC";
			Getflux = "flux_nuebar";
			Cov_names = "nuebar/%s/Active_TPC_Volume/nuebar_%s_Uni_%i_AV_TPC";
			g_simp_names = "anueFluxHisto";
			break;

	}
	
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Get the POT in the file
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	TTree* TPOT = (TTree*) f1->Get("POT");
	if (TPOT == NULL) std::cout << "Error cant get POT info" << std::endl;

	double fPOT{0};
	TPOT->SetBranchAddress("POT", &fPOT); // Get the POT
	TPOT->GetEntry(0);
	std::cout << "TOTAL POT READ IN:\t" << fPOT << std::endl;
	// ------------------------------------------------------------------------------------------------------------
	// CV Flux vs gsimple flux
	// ------------------------------------------------------------------------------------------------------------
	TCanvas* c1 = new TCanvas();

	// Check if sucessfully got histo
	bool boolhist = GetHist(f1, hCV_Flux, Gethist_TPC); if (boolhist == false) gSystem->Exit(0);
	hCV_Flux->SetDirectory(0);
	
	// Normalise flux by bin width (gives a flux/E [GeV])
	for (int i=1;i<hCV_Flux->GetNbinsX()+1;i++) {
		hCV_Flux->SetBinContent(i, hCV_Flux->GetBinContent(i)/hCV_Flux->GetBinWidth(i));		
	}


	TH1D* horig = (TH1D*) hCV_Flux->Clone("horig"); // Clone for plotting so dont need to norm the ms histograms
		
	bool boolfile  = GetFile(f_gsimple , "/uboone/data/users/kmistry/work/PPFX/uboone/NuMIFlux_update_morebins.root"); if (boolfile == false) gSystem->Exit(0);
	boolhist = GetHist(f_gsimple, h_g_simp, g_simp_names); if (boolhist == false) gSystem->Exit(0);

	// Norm
	// 20 is to get the bins in 50 MeV from 1GeV, POT counting done wrong becuase of >1 file per job
	hCV_Flux->Scale( (3* 6.0e20)/ (fPOT*1.0e4) * (50./1000.) );  
	// hCV_Flux->Scale( 3.14159* (6.0e20)/ (100000*950*1.0e4*20) );  // 671.36 is the window area, 20 is to get the bins in 50 MeV from 1GeV pi is fudge factor

	hCV_Flux->Sumw2();
	hCV_Flux->SetLineColor(kRed);
	hCV_Flux->SetLineWidth(2);
	hCV_Flux->SetTitle(";E_{#nu} (GeV);#nu / 6 #times 10^{20} POT / 5 MeV / cm^{2}");
	h_g_simp->SetLineWidth(2);
	gPad->SetLogy();
	gPad->Update();
	hCV_Flux->Draw("hist,same");
	h_g_simp->Draw("hist, same");

	TLegend* lFlux = new TLegend(0.5, 0.65, 0.9, 0.9);
	lFlux->SetNColumns(3);
	lFlux->SetBorderSize(0);
	lFlux->SetFillStyle(0);
	lFlux->SetTextFont(62); 
	lFlux->AddEntry(hCV_Flux, "DK2NU Flux","l");
	lFlux->AddEntry(h_g_simp, "G Simple Flux","l");
	lFlux->Draw();
	
	if (mode == "numu")		hCV_Flux->SetTitle("#nu_{#mu}");
	if (mode == "nue")		hCV_Flux->SetTitle("#nu_{e}");
	if (mode == "numubar")	hCV_Flux->SetTitle("#bar{#nu_{#mu}}");
	if (mode == "nuebar")	hCV_Flux->SetTitle("#bar{#nu_{e}}");

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Correlations, Covariance & uncertainties
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
		if (inputmode.size() == 1) { // checks if long or short list
		std::cout << "Drawing correlation plot" << std::endl;
			cor[0]->SetTitle("Correlation Master Weight");
			cor[0]->Draw("colz");
		} 
		else  {
			// cor[11]->SetTitle("Correlation Master Weight");
			// cor[11]->Draw("colz");
		}
		gStyle->SetPalette(55); // kRainbow

		// Plot fractional errors overlaid with official NOvA plot
		c4->cd();

		// Make the plot
		if (overwrite_errors == false) legDraw(lfrac, herr[l], prodmode, mipp, inputmode[l], mode);
		
		herr[l]->GetYaxis()->SetRangeUser(0,0.5);
		if (mode == "numu")		herr[l]->SetTitle("#nu_{#mu}; Energy [GeV];Fractional Uncertainty");
		if (mode == "nue")		herr[l]->SetTitle("#nu_{e}; Energy [GeV];Fractional Uncertainty");
		if (mode == "numubar")	herr[l]->SetTitle("#bar{#nu_{#mu}}; Energy [GeV];Fractional Uncertainty");
		if (mode == "nuebar")	herr[l]->SetTitle("#bar{#nu_{e}}; Energy [GeV];Fractional Uncertainty");
		// herr[l]->GetYaxis()->SetRangeUser(0,1.75);
		
	}


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
			herr[l]->Reset();
			HPUncertainties_Leo(f1, herr[l], inputmode[l], mode);
			// legDraw(lfrac, herr[l], prodmode, mipp, inputmode[l], mode);

			herr[l]->SetLineColor(kBlack);
			herr[l]->SetLineWidth(2);
			lfrac->AddEntry(herr[l], "PPFXMaster", "l");
			herr[l]->Draw("hist");
			lfrac->Draw();
		}
		c4->Update();
	}

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
	// BeamlineUncertainties(herr, hCV_Flux, "file"); // likely to not be implemented anymore
	// BeamlineUncertainties(herr, hCV_Flux, "stdev");
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

	// Redraw the plot with the new errors
	c1->Update();
	c4->Update();
	// ++++++++++++++++++++++++++++++++++
	// Save the plots as pdfs in the plots folder
	// ++++++++++++++++++++++++++++++++++
	// create plots folder if it does not exist
	gSystem->Exec("if [ ! -d \"plots\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir plots; fi"); 
	// using mipp
	if (mipp == "mippon"){
		if (mode == "numu"){ 	
			c1->Print("plots/CV_Flux_Prediction_NuMu_MIPPOn.pdf");
			c3->Print("plots/Correlation_Matrix_NuMu_MIPPOn.pdf");
			c4->Print("plots/Fractional_Uncertainties_NuMu_MIPPOn.pdf");
			if (wplot == "wplot") c5->Print("plots/Weightplot_NuMu_MIPPOn.pdf");
			std::cout << "\n"<< std::endl;
		}
		else if (mode == "nue") {
			c1->Print("plots/CV_Flux_Prediction_Nue_MIPPOn.pdf");
			c3->Print("plots/Correlation_Matrix_Nue_MIPPOn.pdf");
			c4->Print("plots/Fractional_Uncertainties_Nue_MIPPOn.pdf");
			if (wplot == "wplot") c5->Print("plots/Weightplot_Nue_MIPPOn.pdf");
			std::cout << "\n"<< std::endl;

		}
		else if (mode == "numubar"){
			c1->Print("plots/CV_Flux_Prediction_Numubar_MIPPOn.pdf");
			c3->Print("plots/Correlation_Matrix_Numubar_MIPPOn.pdf");
			c4->Print("plots/Fractional_Uncertainties_Numubar_MIPPOn.pdf");
			if (wplot == "wplot") c5->Print("plots/Weightplot_NuMubar_MIPPOn.pdf");
			std::cout << "\n"<< std::endl;

		}
		else {
			c1->Print("plots/CV_Flux_Prediction_Nuebar_MIPPOn.pdf");
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
			c3->Print("plots/Correlation_Matrix_NuMu_MIPPOff.pdf");
			c4->Print("plots/Fractional_Uncertainties_NuMu_MIPPOff.pdf");
			if (wplot == "wplot") c5->Print("plots/Weightplot_NuMu_MIPPOff.pdf");
			std::cout << "\n"<< std::endl;
		}
		else if (mode == "nue"){
			c1->Print("plots/CV_Flux_Prediction_Nue_MIPPOff.pdf");
			c3->Print("plots/Correlation_Matrix_Nue_MIPPOff.pdf");
			c4->Print("plots/Fractional_Uncertainties_Nue_MIPPOff.pdf");
			if (wplot == "wplot") c5->Print("plots/Weightplot_Nue_MIPPOff.pdf");
			std::cout << "\n"<< std::endl;

		}
		else if (mode == "numubar"){
			c1->Print("plots/CV_Flux_Prediction_Numubar_MIPPOff.pdf");
			c3->Print("plots/Correlation_Matrix_Numubar_MIPPOff.pdf");
			c4->Print("plots/Fractional_Uncertainties_Numubar_MIPPOff.pdf");
			if (wplot == "wplot") c5->Print("plots/Weightplot_NuMubar_MIPPOff.pdf");
			std::cout << "\n"<< std::endl;

		}
		else {
			c1->Print("plots/CV_Flux_Prediction_Nuebar_MIPPOff.pdf");
			c3->Print("plots/Correlation_Matrix_Nuebar_MIPPOff.pdf");
			c4->Print("plots/Fractional_Uncertainties_Nuebar_MIPPOff.pdf");
			if (wplot == "wplot") c5->Print("plots/Weightplot_Nuebar_MIPPOff.pdf");
			std::cout << "\n"<< std::endl;

		}
	}

} // end of main




