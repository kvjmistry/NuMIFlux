/**
 * NuMI Flux at uboone Plotting
 *
 * Plots each individual weighting mode instead of a single one
 * does not plot the correlation matrix for each individual mode to speed up the time
 * 
 * A. Mastbaum <mastbaum@uchicago.edu> 2018/11
 * Modified by K. Mistry 12/18
 */

#include "/uboone/app/users/kmistry/PPFX/numi-validation/scripts/plot_comp_functions.h"

// ----------------------------------------------------------------------------
// Main
void plot_gsimple_flux(TString mode) { // (mippon/mippoff, input, Product/noThinKaon etc. numu/nue)
	gStyle->SetOptStat(0); // say no to stats box
	bool overwrite_errors{false};
	// bool overwrite_errors{true};
	bool novafiles{false};
	//bool novafiles{true};
	bool unweighted{false};
	// bool novafiles{true};

	// Pre declare variables
	TString Getmode, Gethist_TPC, Getflux, Cov_names, g_simp_names, Gethist_TPC_dk2nu;
	TH1D *h_dk2nu_flux;
	TH1D *h_g_simp, *hnovafileflux, *hppfx, *hppfx_mod; 
	TFile *f_gsimple, *f_novafiles, *f_ppfx, *f_ppfx_mod;
	bool boolfile, boolhist;
	double rebin{10}; // number to rebin the histograms by

	// Select neutrino type to run with 
	switch (return_mode(mode)){
		case enumu:
			std::cout << "\nUsing NuMu Mode!\n" << std::endl;
			Getmode = "numu"; 												// Folder name
			Gethist_TPC = "numu/numu_CV_AV_TPC";							// AV in TPC flux prediction
			Getflux = "flux_numu";											// CV flux from NOvA
			Cov_names = "numu/%s/Active_TPC_Volume/numu_%s_Uni_%i_AV_TPC";  // Covariance matrix names
			g_simp_names = "numuFluxHisto";									// G simple files
			Gethist_TPC_dk2nu = "numu/numu_unweighted_AV_TPC";				// uw AV in TPC flux prediction
			break;

		case enue:
			std::cout << "\nUsing Nue Mode!\n" << std::endl;
			Getmode = "nue";
			Gethist_TPC = "nue/nue_CV_AV_TPC";
			Getflux = "flux_nue";
			Cov_names = "nue/%s/Active_TPC_Volume/nue_%s_Uni_%i_AV_TPC";
			g_simp_names = "nueFluxHisto";
			Gethist_TPC_dk2nu = "nue/nue_unweighted_AV_TPC";
			break;

		case enumubar:
			std::cout << "\nUsing NuMubar Mode!\n" << std::endl;
			Getmode = "numubar";
			Gethist_TPC = "numubar/numubar_CV_AV_TPC";
			Getflux = "flux_numubar";
			Cov_names = "numubar/%s/Active_TPC_Volume/numubar_%s_Uni_%i_AV_TPC";
			g_simp_names = "anumuFluxHisto";
			Gethist_TPC_dk2nu = "numubar/numubar_unweighted_AV_TPC";
			break;

		case enuebar:
			std::cout << "\nUsing Nuebar Mode!\n" << std::endl;
			Getmode = "nuebar";
			Gethist_TPC = "nuebar/nuebar_CV_AV_TPC";
			Getflux = "flux_nuebar";
			Cov_names = "nuebar/%s/Active_TPC_Volume/nuebar_%s_Uni_%i_AV_TPC";
			g_simp_names = "anueFluxHisto";
			Gethist_TPC_dk2nu = "nuebar/nuebar_unweighted_AV_TPC";
			break;

	}
	
	// ------------------------------------------------------------------------------------------------------------
	// Flux Comparisons
	// ------------------------------------------------------------------------------------------------------------
	TCanvas* c1 = new TCanvas();
	TLegend* lFlux = new TLegend(0.5, 0.65, 0.9, 0.9);

	// Dk2nu
	boolfile  = GetFile(f_ppfx_mod ,"/uboone/data/users/kmistry/work/PPFX/uboone/DetectorWeights_withtilt/output.root"); if (boolfile == false) gSystem->Exit(0);
	double fPOT = GetPOT(f_ppfx_mod);
	boolhist = GetHist(f_ppfx_mod, h_dk2nu_flux, Gethist_TPC_dk2nu); if (boolhist == false) gSystem->Exit(0);
	h_dk2nu_flux->Rebin(rebin);
	Normalise(h_dk2nu_flux);

	// Get Gsimple files	
	boolfile  = GetFile(f_gsimple , "/uboone/data/users/kmistry/work/PPFX/uboone/NuMIFlux_update_morebins.root"); if (boolfile == false) gSystem->Exit(0);
	boolhist = GetHist(f_gsimple, h_g_simp, g_simp_names); if (boolhist == false) gSystem->Exit(0);
	h_g_simp->Rebin(rebin);
	Normalise(h_g_simp);

	// PPFX flux
	boolhist = GetHist(f_ppfx_mod, hppfx, Gethist_TPC); if (boolhist == false) gSystem->Exit(0);
	hppfx->Rebin(rebin);
	Normalise(hppfx);
	hppfx->SetDirectory(0);

	// Scalings
	h_dk2nu_flux->Scale( (6.0e20)/ (fPOT*1.0e4) );  
	hppfx->Scale( (6.0e20)/ (fPOT*1.0e4) );  

	// Plottings
	h_dk2nu_flux->SetLineColor(kRed+1);
	h_dk2nu_flux->SetLineWidth(2);
	h_dk2nu_flux->SetTitle(";E_{#nu} [GeV];#nu / 6 #times 10^{20} POT / 50 MeV / cm^{2}");
	// h_dk2nu_flux->SetTitle(";E_{#nu} [GeV];#nu / 6 #times 10^{20} POT / GeV / cm^{2}");
	IncreaseLabelSize(h_dk2nu_flux);
	h_dk2nu_flux->GetXaxis()->SetRangeUser(0,10);
	h_dk2nu_flux->Draw("hist");

	h_g_simp->SetLineWidth(2);
	h_g_simp->Draw("hist, same");

	hppfx->SetLineColor(kGreen+1);
	hppfx->SetLineWidth(2);
	hppfx->Draw("hist, same");

	gPad->SetLogy();
	gPad->Update();

	


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
		hnovafileflux->Scale( (3* 6.0e20)/ (2.5e8*1.0e4) * (50./1000.) );  
		
		lFlux->AddEntry(hnovafileflux, "PPFX Flux with NOvA files","l");
		//hnovafileflux->Draw("hist,same");
		//hnovafileflux->Draw("same");

	}

	lFlux->SetNColumns(1);
	lFlux->SetBorderSize(0);
	lFlux->SetFillStyle(0);
	lFlux->SetTextFont(62); 
	
	lFlux->AddEntry(h_dk2nu_flux, "Dk2Nu Flux (no ppfx)","l");
	lFlux->AddEntry(h_g_simp, "G Simple Flux","l");
	lFlux->AddEntry(hppfx, "PPFX Flux","l");;
	lFlux->Draw();
	
	if (mode == "numu")		h_dk2nu_flux->SetTitle("#nu_{#mu}");
	if (mode == "nue")		h_dk2nu_flux->SetTitle("#nu_{e}");
	if (mode == "numubar")	h_dk2nu_flux->SetTitle("#bar{#nu_{#mu}}");
	if (mode == "nuebar")	h_dk2nu_flux->SetTitle("#bar{#nu_{e}}");

	// ------------------------------------------------------------------------------------------------------------
	// Draw all fluxes on one plot
	// ------------------------------------------------------------------------------------------------------------
	TCanvas* c_plotall = new TCanvas();
	TLegend* l_plotall = new TLegend(0.8, 0.65, 0.95, 0.9);

	c_plotall->cd();
	PlotFluxSame(c_plotall, l_plotall, f_ppfx_mod, "numu", fPOT, "numu/numu_CV_AV_TPC" );
	PlotFluxSame(c_plotall, l_plotall, f_ppfx_mod, "numubar", fPOT, "numubar/numubar_CV_AV_TPC" );
	PlotFluxSame(c_plotall, l_plotall, f_ppfx_mod, "nue", fPOT, "nue/nue_CV_AV_TPC" );
	PlotFluxSame(c_plotall, l_plotall, f_ppfx_mod, "nuebar", fPOT, "nuebar/nuebar_CV_AV_TPC" );

	l_plotall->SetNColumns(1);
	l_plotall->SetBorderSize(0);
	l_plotall->SetFillStyle(0);
	l_plotall->SetTextFont(62); 
	l_plotall->Draw();

	
	c1->Update();
	c_plotall->Update();

	// ++++++++++++++++++++++++++++++++++
	// Save the plots as pdfs in the plots folder
	// ++++++++++++++++++++++++++++++++++
	// create plots folder if it does not exist
	gSystem->Exec("if [ ! -d \"plots\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir plots; fi"); 
	
	if (mode == "numu"){ 	
		c1->Print("plots/CV_Flux_Prediction_NuMu_MIPPOff.pdf");
		std::cout << "\n"<< std::endl;
	}
	else if (mode == "nue"){
		c1->Print("plots/CV_Flux_Prediction_Nue_MIPPOff.pdf");
		std::cout << "\n"<< std::endl;

	}
	else if (mode == "numubar"){
		c1->Print("plots/CV_Flux_Prediction_Numubar_MIPPOff.pdf");
		std::cout << "\n"<< std::endl;

	}
	else {
		c1->Print("plots/CV_Flux_Prediction_Nuebar_MIPPOff.pdf");
		std::cout << "\n"<< std::endl;

	}
	
	c_plotall->Print("plots/CVFlux_All_Flavours.pdf");

// gSystem->Exit(0);

} // end of main




