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
	// bool novafiles{false};
	bool novafiles{false};
	bool unweighted{false};
	// bool novafiles{true};

	std::vector<std::string> inputmode = loopdir(inputfile, mode); // Grab the names of the input reweighters

	// Pre declare variables
	TString Getmode, Gethist_TPC, Getflux, Cov_names, g_simp_names, Gethist_TPC_uw, Gethist_TPC_Th;
	TH1D *hCV_Flux, *hUW_Flux, *h_g_simp, *hnovafileflux;
	TFile *f_gsimple, *f_novafiles;
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
			Gethist_TPC_uw = "numu/numu_unweighted_AV_TPC";					// AV in TPC flux prediction unweighted
			Gethist_TPC_Th = "numu/Th_numu_CV_AV_TPC";						// AV in TPC flux prediction in theta
			break;

		case enue:
			std::cout << "\nUsing Nue Mode!\n" << std::endl;
			Getmode = "nue";
			Gethist_TPC = "nue/nue_CV_AV_TPC";
			Getflux = "flux_nue";
			Cov_names = "nue/%s/Active_TPC_Volume/nue_%s_Uni_%i_AV_TPC";
			g_simp_names = "nueFluxHisto";
			Gethist_TPC_uw = "nue/nue_unweighted_AV_TPC";
			Gethist_TPC_Th = "nue/Th_nue_CV_AV_TPC";
			break;

		case enumubar:
			std::cout << "\nUsing NuMubar Mode!\n" << std::endl;
			Getmode = "numubar";
			Gethist_TPC = "numubar/numubar_CV_AV_TPC";
			Getflux = "flux_numubar";
			Cov_names = "numubar/%s/Active_TPC_Volume/numubar_%s_Uni_%i_AV_TPC";
			g_simp_names = "anumuFluxHisto";
			Gethist_TPC_uw = "numubar/numubar_unweighted_AV_TPC";
			Gethist_TPC_Th = "numubar/Th_numubar_CV_AV_TPC";
			break;

		case enuebar:
			std::cout << "\nUsing Nuebar Mode!\n" << std::endl;
			Getmode = "nuebar";
			Gethist_TPC = "nuebar/nuebar_CV_AV_TPC";
			Getflux = "flux_nuebar";
			Cov_names = "nuebar/%s/Active_TPC_Volume/nuebar_%s_Uni_%i_AV_TPC";
			g_simp_names = "anueFluxHisto";
			Gethist_TPC_uw = "nuebar/nuebar_unweighted_AV_TPC";
			Gethist_TPC_Th = "nuebar/Th_nuebar_CV_AV_TPC";
			break;

	}
	
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Get the POT in the file
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	double fPOT = GetPOT(TFile* f);
	// ------------------------------------------------------------------------------------------------------------
	// CV Flux vs gsimple flux vs nova files flux
	// ------------------------------------------------------------------------------------------------------------
	TCanvas* c1 = new TCanvas();
	TLegend* lFlux = new TLegend(0.5, 0.65, 0.9, 0.9);

	// Check if sucessfully got histo
	bool boolhist = GetHist(f1, hCV_Flux, Gethist_TPC); if (boolhist == false) gSystem->Exit(0);
	hCV_Flux->SetDirectory(0);
	
	// Normalise flux by bin width (gives a flux/E [GeV])
	for (int i=1;i<hCV_Flux->GetNbinsX()+1;i++) {
		hCV_Flux->SetBinContent(i, hCV_Flux->GetBinContent(i)/hCV_Flux->GetBinWidth(i));		
	}
	Normalise(hCV_Flux);


	TH1D* horig = (TH1D*) hCV_Flux->Clone("horig"); // Clone for plotting so dont need to norm the ms histograms
		
	bool boolfile  = GetFile(f_gsimple , "/uboone/data/users/kmistry/work/PPFX/uboone/NuMIFlux_update_morebins.root"); if (boolfile == false) gSystem->Exit(0);
	boolhist = GetHist(f_gsimple, h_g_simp, g_simp_names); if (boolhist == false) gSystem->Exit(0);

	// Norm
	// 20 is to get the bins in 50 MeV from 1GeV, POT counting done wrong becuase of >1 file per job
	hCV_Flux->Scale( (6.0e20)/ (fPOT*1.0e4) );  
	// hCV_Flux->Scale( 3.14159* (6.0e20)/ (100000*950*1.0e4*20) );  // 671.36 is the window area, 20 is to get the bins in 50 MeV from 1GeV pi is fudge factor

	hCV_Flux->Sumw2();
	hCV_Flux->SetLineColor(kRed+1);
	hCV_Flux->SetLineWidth(2);
	hCV_Flux->SetTitle(";E_{#nu} [GeV];#nu / 6 #times 10^{20} POT / 50 MeV / cm^{2}");
	h_g_simp->SetLineWidth(2);
	gPad->SetLogy();
	gPad->Update();
	hCV_Flux->Draw("");
	hCV_Flux->Draw("hist,same");
	//h_g_simp->Draw("hist, same");

	// choose whether to draw the flux using nova files which have a threshold
	if (novafiles){
		bool boolfile  = GetFile(f_novafiles, "/uboone/data/users/kmistry/work/PPFX/uboone/bugfix_release_novafiles/output.root"); if (boolfile == false) gSystem->Exit(0);
		boolhist = GetHist(f_novafiles, hnovafileflux, Gethist_TPC); if (boolhist == false) gSystem->Exit(0);
		hnovafileflux->SetLineColor(kGreen+3);
		hnovafileflux->SetLineWidth(2);

		// Normalise flux by bin width (gives a flux/E [GeV])
		for (int i=1;i<	hnovafileflux->GetNbinsX()+1;i++) {
			hnovafileflux->SetBinContent(i, hnovafileflux->GetBinContent(i)/hnovafileflux->GetBinWidth(i));		
		}
		
		// Norm
		hnovafileflux->Scale( (6.0e20)/ (2.5e8*1.0e4) * (50./1000.) );  
		
		lFlux->AddEntry(hnovafileflux, "PPFX Flux with NOvA files","l");
		hnovafileflux->Draw("hist,same");
		hnovafileflux->Draw("same");

	}

	lFlux->SetNColumns(1);
	lFlux->SetBorderSize(0);
	lFlux->SetFillStyle(0);
	lFlux->SetTextFont(62); 
	lFlux->AddEntry(hCV_Flux, "PPFX Flux","l");
	//lFlux->AddEntry(h_g_simp, "G Simple Flux (no tilt wght)","l");
	lFlux->Draw();
	
	if (mode == "numu")		hCV_Flux->SetTitle("#nu_{#mu}");
	if (mode == "nue")		hCV_Flux->SetTitle("#nu_{e}");
	if (mode == "numubar")	hCV_Flux->SetTitle("#bar{#nu_{#mu}}");
	if (mode == "nuebar")	hCV_Flux->SetTitle("#bar{#nu_{e}}");

	// ------------------------------------------------------------------------------------------------------------
	// Draw all fluxes on one plot
	// ------------------------------------------------------------------------------------------------------------
	TCanvas* c_plotall = new TCanvas();
	TLegend* l_plotall = new TLegend(0.8, 0.65, 0.95, 0.9);

	c_plotall->cd();
	PlotFluxSame(c_plotall, l_plotall, f1, "numu", fPOT, "numu/numu_CV_AV_TPC" );
	PlotFluxSame(c_plotall, l_plotall, f1, "numubar", fPOT, "numubar/numubar_CV_AV_TPC" );
	PlotFluxSame(c_plotall, l_plotall, f1, "nue", fPOT, "nue/nue_CV_AV_TPC" );
	PlotFluxSame(c_plotall, l_plotall, f1, "nuebar", fPOT, "nuebar/nuebar_CV_AV_TPC" );

	l_plotall->SetNColumns(1);
	l_plotall->SetBorderSize(0);
	l_plotall->SetFillStyle(0);
	l_plotall->SetTextFont(62); 
	l_plotall->Draw();

	// ------------------------------------------------------------------------------------------------------------
	// Draw weighted flux vs unweighted flux
	// ------------------------------------------------------------------------------------------------------------
	TCanvas* c_uw_v_w = new TCanvas();
	TLegend* l_uw_v_w = new TLegend(0.6, 0.65, 0.75, 0.9);

	c_uw_v_w->cd();
	// Check if sucessfully got histo
	boolhist = GetHist(f1, hUW_Flux, Gethist_TPC_uw); if (boolhist == false) gSystem->Exit(0);
	
	// Normalise flux by bin width (gives a flux/E [GeV])
	for (int i=1;i<hUW_Flux->GetNbinsX()+1;i++) {
		hUW_Flux->SetBinContent(i, hUW_Flux->GetBinContent(i)/hUW_Flux->GetBinWidth(i));		
	}

	// 6e20 POT, 50/1000 for 50 MeV from 1GeV, 1e-4 for m2->cm2
	hUW_Flux->Scale( (6.0e20)/ (fPOT*1.0e4) * (50./1000.) );  


	hUW_Flux->SetLineColor(kBlue+1);
	hUW_Flux->SetLineWidth(2);
	hUW_Flux->SetTitle(";Energy [GeV];#nu / 6 #times 10^{20} POT / 50 MeV / cm^{2}");
	gPad->SetLogy();
	gPad->Update();
	hCV_Flux->GetXaxis()->SetRangeUser(0,5);
	
	hCV_Flux->Draw("hist");
	hUW_Flux->Draw("hist,same");
	
	l_uw_v_w->AddEntry(hUW_Flux, "Unweighted","l");
	l_uw_v_w->AddEntry(hCV_Flux, "PPFX","l");
	l_uw_v_w->SetNColumns(1);
	l_uw_v_w->SetBorderSize(0);
	l_uw_v_w->SetFillStyle(0);
	l_uw_v_w->SetTextFont(62); 
	l_uw_v_w->Draw();

	// Lets integrate the Gsimple flux and the dk2nu flux and see how differnent they really are!
	
	
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

			// // Catch zeros , set to arbitary 1.0
			// if (horig->GetBinContent(i) <= 0) horig_cont = 1.0;
			// else horig_cont = horig->GetBinContent(i);

			if (horig->GetBinContent(i) == 0) herr[l]->SetBinContent(i, 0);
			else herr[l]->SetBinContent(i, sqrt(cii) / horig->GetBinContent(i));
			
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
			// herr[l]->Reset();
			HPUncertainties_Leo(f1, herr[l], inputmode[l], mode);
			// legDraw(lfrac, herr[l], prodmode, mipp, inputmode[l], mode);

			herr[l]->SetLineColor(kBlack);
			herr[l]->SetLineWidth(2);
			lfrac->AddEntry(herr[l], "PPFXMaster", "l");
			if (mode == "numu")		herr[l]->SetTitle("#nu_{#mu}; Energy [GeV];Fractional Uncertainty");
			if (mode == "nue")		herr[l]->SetTitle("#nu_{e}; Energy [GeV];Fractional Uncertainty");
			if (mode == "numubar")	herr[l]->SetTitle("#bar{#nu_{#mu}}; Energy [GeV];Fractional Uncertainty");
			if (mode == "nuebar")	herr[l]->SetTitle("#bar{#nu_{e}}; Energy [GeV];Fractional Uncertainty");
			herr[l]->Draw("hist");
			// lfrac->Draw();
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
	// Make a 4D Covariance matrix for re-weighing
	// ------------------------------------------------------------------------------------------------------------
	// Call caclulate covzriance matrix function. For now write function here.
	TH1D *hCV2d, *hu, *hCV_unwrap, *hu_unwrap; 	// Flux hist for each universe
	TH2D *cov4d;
	TFile* f2d;
	std::string mode_str;
	if (mode == "numu") mode_str = "numu";
	else if (mode == "numubar") mode_str = "numubar";
	else if (mode == "nue") mode_str = "nue";
	else mode_str = "nuebar";

	boolfile  = GetFile(f2d , "/uboone/data/users/kmistry/work/PPFX/uboone/with_tilt_2Dhists/output.root"); if (boolfile == false) gSystem->Exit(0);
	// boolfile  = GetFile(f2d , "/uboone/data/users/kmistry/work/PPFX/uboone/beamline/run0008/output.root"); if (boolfile == false) gSystem->Exit(0);


	// Get the CV 2D matrix
	boolhist = GetHist(f2d, hCV2d, Gethist_TPC);   if (boolhist == false) gSystem->Exit(0);
	const int nBinsEnu = hCV2d->GetXaxis()->GetNbins(); // Enu
	const int nBinsTh = hCV2d->GetYaxis()->GetNbins(); // Theta
	int nuni{100};
	
	//------------------------------
	// Normalise 2d hist by bin area / deg / GeV
	// Loop over rows
	for (int i=1; i<nBinsEnu+1; i++) {

		// Loop over columns
		for (int j=1; j<nBinsTh+1; j++) {
			hCV2d->SetBinContent(i,j, hCV2d->GetBinContent(i, j)/ ( hCV2d->GetXaxis()->GetBinWidth(i) * hCV2d->GetYaxis()->GetBinWidth(j) ));
		}

	} 
	hCV2d->Scale((6.0e20)/ (2.5e8*1.0e4)); // scale to right POT and m2
	//------------------------------
	// Unwrap the histogram to binindex
	hCV_unwrap = new TH1D("", "",nBinsEnu*nBinsTh, 0, nBinsEnu*nBinsTh );
	int counter{0};
	for (int i=1; i<nBinsEnu+1; i++) { // Loop over rows
		for (int j=1; j<nBinsTh+1; j++){// Loop over columns
			counter++;
			hCV_unwrap->SetBinContent(counter, hCV2d->GetBinContent(i , j)  );
		}
	}

	// Draw the CV Unwrapped
	TCanvas* c_CV_unwrap= new TCanvas();
	c_CV_unwrap->cd();
	hCV_unwrap->SetTitle("CV Unwrapped; Bin index; Flux ");
	hCV_unwrap->SetLineWidth(2);
	hCV_unwrap->SetLineColor(kBlack);
	hCV_unwrap->Draw("hist");

	// Draw the CV in 2D
	TCanvas* c_CV2d= new TCanvas();
	c_CV2d->cd();
	hCV2d->SetTitle("CV 2D; Energy [GeV]; Theta [deg] ");
	gPad->SetLogz();
	gPad->Update();
	hCV2d->Draw("colz");

	//------------------------------
	// 4D Covariance matrix	
	cov4d  = new TH2D(Form("PPFXMaster_cov_4d"), ";Bin i; Bin j", nBinsEnu*nBinsTh, 0, nBinsEnu*nBinsTh , nBinsEnu*nBinsTh, 0, nBinsEnu*nBinsTh); // Bin i  Bin j

	// Loop over universes
	for (int k=0; k < nuni; k++) {
		char name[500];
		snprintf(name, 500,"%s/PPFXMaster/Active_TPC_Volume/%s_PPFXMaster_Uni_%i_AV_TPC" ,mode_str.c_str(), mode_str.c_str(), k); // Get uni i
		
		
		// Check if sucessfully got histo
		boolhist = GetHist(f2d, hu, name); if (boolhist == false) gSystem->Exit(0);
		//------------------------------
		// Normalise 2d hist by bin area / deg / GeV
		// Loop over rows
		for (int p=1; p<nBinsEnu+1; p++) {
			// Loop over columns
			for (int q=1; q<nBinsTh+1; q++) {
				hu->SetBinContent(p,q, hu->GetBinContent(p, q)/ ( hu->GetXaxis()->GetBinWidth(p) * hu->GetYaxis()->GetBinWidth(q) ));
			}
		} 
		hu->Scale((6.0e20)/ (2.5e8*1.0e4)); // scale to right POT and m2
		//------------------------------
		// Unwrap the histogram to binindex
		hu_unwrap = new TH1D("", "",nBinsEnu*nBinsTh, 0, nBinsEnu*nBinsTh );

		counter = 0;
		// Loop over rows
		for (int i=1; i<nBinsEnu+1; i++) { 
			// Loop over columns
			for (int j=1; j<nBinsTh+1; j++){ 
				counter ++;
				hu_unwrap->SetBinContent( counter, hu->GetBinContent(i , j)  );
			}
		}
		//------------------------------
		// Now calculate the covariance matrix
		CalcCovariance_4D(cov4d, hCV_unwrap, hu_unwrap, nBinsEnu*nBinsTh, nBinsEnu*nBinsTh, k  );

	} // End cov calc for universe i
	
	TCanvas* c_cov= new TCanvas();
	c_cov->cd();
	cov4d->SetTitle("4D Covariance Matrix ; Bin i; Bin j");
	gStyle->SetPalette(kDeepSea);
	cov4d->Draw("colz");
	gPad->SetLogz();
	gPad->Update();

	//------------------------------
	// Calculate the correlation matrix
	TCanvas *c_corr4d = new TCanvas();
	c_corr4d->cd();
	TH2D *hcorr4d = (TH2D*) cov4d->Clone("hCorr4d");
	CalcCorrelation(hcorr4d, cov4d, nBinsEnu*nBinsTh );
	hcorr4d->SetTitle("4D Correlation Matrix ; Bin i; Bin j");
	hcorr4d->Draw("colz");
	// gPad->SetLogz();
	// gPad->Update();
	//------------------------------
	// Calculate the fractional covariance matrix
	TCanvas *c_fraccov4d = new TCanvas();
	c_fraccov4d->cd();
	TH2D *hfraccov4d = (TH2D*) cov4d->Clone("hfraccov4d");
	CalcFracCovariance(hCV_unwrap, hfraccov4d, nBinsEnu*nBinsTh );
	hfraccov4d->SetTitle("4D Fractional Covariance Matrix ; Bin i; Bin j");
	hfraccov4d->Draw("colz");

	//------------------------------
	// Create a histogram with the bin indexes
	TCanvas *c_binidx = new TCanvas();
	c_binidx->cd();
	TH2D *hbinidx = (TH2D*) hCV2d->Clone("hCVClone");
	counter = 0;
	// Loop over rows
	for (int i=1; i<nBinsEnu+1; i++) {  // Enu
		// Loop over columns
		for (int j=1; j<nBinsTh+1; j++){ // Theta
			counter ++;
			hbinidx->SetBinContent( i, j, counter );
		}
	}
	hbinidx->SetTitle("Bin Indexes; Energy [GeV]; Theta [deg]");
	hbinidx->Draw("text00");
	gPad->SetLogx();
	gPad->Update();


	// ------------------------------------------------------------------------------------------------------------
	// Compare the percentage difference between the mean and the CV in each unwrapped histogram bin right now we output the pull
	// ------------------------------------------------------------------------------------------------------------
	TCanvas *cMeanCV = new TCanvas();
	TH1D* hRatioCVMean;
	TH1D *hMean_unwrap = (TH1D*) hCV_unwrap->Clone("hMean_unwrap");
	CalcMeanHist(f2d, hMean_unwrap, nBinsEnu, nBinsTh,mode);

	// Now call a function that calcuates the ratio between them
	CalcRatioMeanCV(hCV_unwrap, hMean_unwrap, hRatioCVMean);

	// Fill a histogram with the Ratio values to better visualise
	TH1D* hRatioMeanCVHist = new TH1D("hRatioMeanCVHist", "Ratio CV to Mean; Ratio; Entries", 30, 0.6, 1.4);
	TH1D* hPull = new TH1D("hPull", "CV - Mean / (stdev/#sqrt{n}); Pull; Entries", 25,-10, 10); // histogram with the pull 

	// Get rid of zeros
	for (unsigned int i =1; i <  hRatioCVMean->GetNbinsX()+1; i++ ){
		if (hRatioCVMean->GetBinContent(i) == 0) continue;
		hRatioMeanCVHist->Fill(hRatioCVMean->GetBinContent(i));
	}

	CalcPull(hCV_unwrap, hMean_unwrap, hPull);
	hPull->Draw("hist");

	// hCV_unwrap->Draw("hist");
	// hMean_unwrap->Draw("histsame");
	// hRatioCVMean->Draw("E2");
	// hRatioMeanCVHist->SetLineWidth(2);
	// hRatioMeanCVHist->SetLineColor(kBlack);
	// hRatioMeanCVHist->Draw("hist");
	// TLine* flat = new TLine(0, 1, 400, 1);
	// flat->SetLineStyle(7);
	// flat->Draw();

	// ------------------------------------------------------------------------------------------------------------
	// Make a plot of the fractional uncertainties from the 4d covariance matrix
	// ------------------------------------------------------------------------------------------------------------
	TCanvas *c_FracError4d = new TCanvas();
	TH1D *hFracError4d = (TH1D*) hCV_unwrap->Clone("hFracError4d");

	// Function to caluclate the fractional uncertainties
	CalcFractionalError(cov4d, hCV_unwrap, hFracError4d );
	c_FracError4d->cd();
	hFracError4d->Draw("hist");



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

	c4->Update();
	cband->Update();
	

	// Redraw the plot with the new errors
	c1->Update();
	c_uw_v_w->Update();
	

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
			c_uw_v_w->Print("plots/Unweighted_vs_ppfx_NuMu_MIPPOff.pdf");
			if (wplot == "wplot") c5->Print("plots/Weightplot_NuMu_MIPPOff.pdf");
			c_cov->Print("plots/CovarianceMarix4D_NuMu_MIPPOff.pdf");			
			c_corr4d->Print("plots/CorrelationMarix4D_NuMu_MIPPOff.pdf");
			c_FracError4d->Print("plots/FracCovarianceMarix4D_NuMu_MIPPOff.pdf");
			c_binidx->Print("plots/BinIndex_NuMu.pdf");
			c_fraccov4d->Print("plots/FracCovariance4d_NuMu_MIPPOff.pdf");
			c_CV2d->Print("plots/CV2d_NuMu_MIPPOff.pdf");
			std::cout << "\n"<< std::endl;
		}
		else if (mode == "nue"){
			c1->Print("plots/CV_Flux_Prediction_Nue_MIPPOff.pdf");
			c3->Print("plots/Correlation_Matrix_Nue_MIPPOff.pdf");
			c4->Print("plots/Fractional_Uncertainties_Nue_MIPPOff.pdf");
			c_uw_v_w->Print("plots/Unweighted_vs_ppfx_Nue_MIPPOff.pdf");
			if (wplot == "wplot") c5->Print("plots/Weightplot_Nue_MIPPOff.pdf");
			c_cov->Print("plots/CovarianceMarix4D_Nue_MIPPOff.pdf");
			c_corr4d->Print("plots/CorrelationMarix4D_Nue_MIPPOff.pdf");
			c_FracError4d->Print("plots/FracCovarianceMarix4D_Nue_MIPPOff.pdf");
			c_binidx->Print("plots/BinIndex_Nue.pdf");
			c_fraccov4d->Print("plots/FracCovariance4d_Nue_MIPPOff.pdf");
			c_CV2d->Print("plots/CV2d_Nue_MIPPOff.pdf");
			std::cout << "\n"<< std::endl;

		}
		else if (mode == "numubar"){
			c1->Print("plots/CV_Flux_Prediction_Numubar_MIPPOff.pdf");
			c3->Print("plots/Correlation_Matrix_Numubar_MIPPOff.pdf");
			c4->Print("plots/Fractional_Uncertainties_Numubar_MIPPOff.pdf");
			c_uw_v_w->Print("plots/Unweighted_vs_ppfx_Numubar_MIPPOff.pdf");
			if (wplot == "wplot") c5->Print("plots/Weightplot_NuMubar_MIPPOff.pdf");
			c_cov->Print("plots/CovarianceMarix4D_Numubar_MIPPOff.pdf");
			c_corr4d->Print("plots/CorrelationMarix4D_Numubar_MIPPOff.pdf");
			c_FracError4d->Print("plots/FracCovarianceMarix4D_Numubar_MIPPOff.pdf");
			c_binidx->Print("plots/BinIndex_Numubar.pdf");
			c_fraccov4d->Print("plots/FracCovariance4d_Numubar_MIPPOff.pdf");
			c_CV2d->Print("plots/CV2d_Numubar_MIPPOff.pdf");
			std::cout << "\n"<< std::endl;

		}
		else {
			c1->Print("plots/CV_Flux_Prediction_Nuebar_MIPPOff.pdf");
			c3->Print("plots/Correlation_Matrix_Nuebar_MIPPOff.pdf");
			c4->Print("plots/Fractional_Uncertainties_Nuebar_MIPPOff.pdf");
			c_uw_v_w->Print("plots/Unweighted_vs_ppfx_Nuebar_MIPPOff.pdf");
			if (wplot == "wplot") c5->Print("plots/Weightplot_Nuebar_MIPPOff.pdf");
			c_cov->Print("plots/CovarianceMarix4D_Nuebar_MIPPOff.pdf");
			c_corr4d->Print("plots/CorrelationMarix4D_Nuebar_MIPPOff.pdf");
			c_FracError4d->Print("plots/FracCovarianceMarix4D_Nuebar_MIPPOff.pdf");
			c_binidx->Print("plots/BinIndex_Nuebar.pdf");
			c_fraccov4d->Print("plots/FracCovariance4d_Nuebar_MIPPOff.pdf");
			c_CV2d->Print("plots/CV2d_Nuebar_MIPPOff.pdf");
			std::cout << "\n"<< std::endl;

		}
		c_plotall->Print("plots/flux_all_flavours.pdf");
	}

} // end of main




