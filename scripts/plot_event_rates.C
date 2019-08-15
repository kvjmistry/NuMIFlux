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
void DrawSpecifiers(TCanvas* c, TH1D* &h, const char* mode){
	c->cd();

	// Plottings
	if (!strcmp(mode, "numu")){
		h->SetLineColor(kRed+1);
		h->SetLineWidth(2);
	}
	if (!strcmp(mode, "nue")){
		h->SetLineColor(kRed+1);
		h->SetLineWidth(2);
		h->SetLineStyle(2);
	}
	if (!strcmp(mode, "numubar")){
		h->SetLineColor(kBlue+1);
		h->SetLineWidth(2);
	}
	if (!strcmp(mode, "nuebar")){
		h->SetLineColor(kBlue+1);
		h->SetLineWidth(2);
		h->SetLineStyle(2);
	}
	
	IncreaseLabelSize(h);
	h->GetXaxis()->SetRangeUser(0,5);
}
// ----------------------------------------------------------------------------
constexpr double RELATIVE_TOLERANCE = 1e-5;

double bin_average_total_xsec(TH1D* flux_hist, TGraph *xsec_spline, int bin_number) {
	double E_min = flux_hist->GetBinLowEdge( bin_number );
	double E_max = flux_hist->GetBinLowEdge( bin_number + 1 );

	// Make sure that the TF1 that we define is valid up to the maximum energy
	// used in either the spline or the flux histogram
	double spline_max_energy = xsec_spline->GetX()[ xsec_spline->GetN() - 1 ];
	double flux_max_energy   = flux_hist->GetBinLowEdge( flux_hist->GetNbinsX() + 1 );
	double max_energy        = std::max( spline_max_energy, flux_max_energy );

	// Function to use for numerical integration. Element x[0] is the
	// neutrino energy
	TF1 temp_func("temp_func", [&](double* x, double*)
	{ double total_xsec = xsec_spline->Eval( x[0] ); return total_xsec; },
	0., max_energy, 0);

	double integ = temp_func.Integral(E_min, E_max, RELATIVE_TOLERANCE);
	return integ / (E_max - E_min);
}
// ----------------------------------------------------------------------------
// Main
void plot_event_rates(const char* horn) {
	gStyle->SetOptStat(0); // say no to stats box

	// Pre declare variables
	TString Gethist_TPC, Gethist_TPC_dk2nu;
	TFile *f, *f_gsimp, *f_genie, *f_genie_nue, *f_gevgen_numu, *f_gevgen_numubar, *f_gevgen_nue, *f_gevgen_nuebar;
	bool boolfile, boolhist;
	double rebin{10}; // number to rebin the histograms by

	// double Ntarget = 4.76e31/56.41e6* 256.35*233*1036.8; //TPC active!!! Marco
	// double Ntarget = (1.3836*6.022e23*40*256.35*233*1036.8) / 39.95; //TPC active!!! Colton
	double Ntarget = (1.3954*6.022e23*40*256.35*233*1036.8) / 39.95; //TPC active with MC lar density -- Krishan
	std::cout << "N_Targ:\t" << Ntarget << std::endl;


	// std::cout << ((1.3836*6.022e23*40*256.35*233*1036.8) / 39.95)/(4.76e31/56.41e6* 256.35*233*1036.8)<< std::endl;

	

	double histMin = 0;
	double histMax = 20;
	int histNbins = 4000;
	
	TGraph *genieXsecNumuCC;
	TGraph *genieXsecNumubarCC;
	TGraph *genieXsecNueCC;
	TGraph *genieXsecNuebarCC;
	
	TH1D* numuCCHisto  = new TH1D("numuCCHisto", "#nu_{#mu} CC; #nu_{#mu} Energy [GeV]; #nu_{#mu} CC / 79 t / 6 #times 10^{20} POT / 50 MeV",histNbins,histMin,histMax);
	TH1D* anumuCCHisto = new TH1D("anumuCCHisto", "#bar{#nu}_{#mu} CC; #bar{#nu}_{#mu} Energy [GeV]; #bar{#nu}_{#mu} CC / 79 t / 6 #times 10^{20} POT / 50 MeV",histNbins,histMin,histMax);
	TH1D* nueCCHisto   = new TH1D("nueCCHisto", "#nu_{e} CC; #nu_{e} Energy [GeV]; #nu_{e} CC / 79 t / 6 #times 10^{20} POT / 50 MeV",histNbins,histMin,histMax);
	TH1D* anueCCHisto  = new TH1D("anueCCHisto", "#bar{#nu}_{e} CC; #bar{#nu}_{e} Energy [GeV]; #bar{#nu}_{e} CC / 79 t / 6 #times 10^{20} POT / 50 MeV",histNbins,histMin,histMax);

	// Gsimple 
	TH1D* numuCCHisto_gsimp  = new TH1D("numuCCHisto_gsimp", "#nu_{#mu} CC; #nu_{#mu} Energy [GeV]; #nu_{#mu} CC / 79 t / 6 #times 10^{20} POT / 50 MeV",histNbins,histMin,histMax);
	TH1D* anumuCCHisto_gsimp = new TH1D("anumuCCHisto_gsimp", "#bar{#nu}_{#mu} CC; #bar{#nu}_{#mu} Energy [GeV]; #bar{#nu}_{#mu} CC / 79 t / 6 #times 10^{20} POT / 50 MeV",histNbins,histMin,histMax);
	TH1D* nueCCHisto_gsimp   = new TH1D("nueCCHisto_gsimp", "#nu_{e} CC; #nu_{e} Energy [GeV]; #nu_{e} CC / 79 t / 6 #times 10^{20} POT / 50 MeV",histNbins,histMin,histMax);
	TH1D* anueCCHisto_gsimp  = new TH1D("anueCCHisto_gsimp", "#bar{#nu}_{e} CC; #bar{#nu}_{e} Energy [GeV]; #bar{#nu}_{#mu} CC / 79 t / 6 #times 10^{20} POT / 50 MeV",histNbins,histMin,histMax);
	
	// Gevgen
	TH1D* numuCCHisto_gevgen  = new TH1D("numuCCHisto_gevgen",  "#nu_{#mu} CC; #nu_{#mu} Energy [GeV]; #nu_{#mu} CC", 4000,histMin,histMax);
	TH1D* anumuCCHisto_gevgen = new TH1D("anumuCCHisto_gevgen", "#bar{#nu}_{#mu} CC; #bar{#nu}_{#mu} Energy [GeV]; #bar{#nu}_{#mu} CC", 4000, histMin,histMax);
	TH1D* nueCCHisto_gevgen   = new TH1D("nueCCHisto_gevgen",  "#nu_{e} CC; #nu_{e} Energy [GeV]; #nu_{e} CC", 4000,histMin,histMax);
	TH1D* anueCCHisto_gevgen  = new TH1D("anueCCHisto_gevgen", "#bar{#nu}_{e} CC; #bar{#nu}_{e} Energy [GeV]; #bar{#nu}_{e} CC", 4000, histMin,histMax);

	TH1D *h_nue_genie, *h_nuebar_genie, *h_numu_genie, *h_numubar_genie;


	// ------------------------------------------------------------------------------------------------------------
	// Make a plot of flux x genie spline
	// ------------------------------------------------------------------------------------------------------------
	if (!strcmp(horn,"fhc")) {
		boolfile  = GetFile(f ,"/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold/output_uboone_run0.root");
		// boolfile  = GetFile(f ,"/uboone/data/users/kmistry/work/PPFX/uboone/test/output_run0_lowstats.root");
		if (boolfile == false) gSystem->Exit(0);
		boolfile  = GetFile(f_gsimp , "/uboone/data/users/kmistry/work/PPFX/uboone/NuMIFlux_update_morebins.root");
		if (boolfile == false) gSystem->Exit(0);
	}
	else {
		boolfile  = GetFile(f ,"/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold/RHC/output_uboone_run0.root");
		if (boolfile == false) gSystem->Exit(0);
		boolfile  = GetFile(f_gsimp , "/uboone/app/users/kmistry/Flux/NuMIFlux/files/NuMIFlux_anti_morebins.root"); 
		if (boolfile == false) gSystem->Exit(0);
		
	}

	// Get the event rate distribution generated through GENIE
	boolfile  = GetFile(f_genie ,"../files/NuMI_EventRate.root");
	if (boolfile == false) gSystem->Exit(0);
	
	// Get the event rate distribution generated through GENIE -- sample with nue enhanced for stats
	boolfile  = GetFile(f_genie_nue ,"../files/NuMI_EventRate_nue.root");
	if (boolfile == false) gSystem->Exit(0);
	
	// Get the event rate distribution generated through gevgen for numu
	boolfile  = GetFile(f_gevgen_numu ,"/uboone/data/users/kmistry/work/PPFX/uboone/genie/my_events.gst_numu.root");
	if (boolfile == false) gSystem->Exit(0);

	// Get the event rate distribution generated through gevgen for numubar
	boolfile  = GetFile(f_gevgen_numubar ,"/uboone/data/users/kmistry/work/PPFX/uboone/genie/my_events.gst_numubar.root");
	if (boolfile == false) gSystem->Exit(0);

	// Get the event rate distribution generated through gevgen for nue
	boolfile  = GetFile(f_gevgen_nue ,"/uboone/data/users/kmistry/work/PPFX/uboone/genie/my_events.gst_nue.root");
	if (boolfile == false) gSystem->Exit(0);

	// Get the event rate distribution generated through gevgen for nuebar
	boolfile  = GetFile(f_gevgen_nuebar ,"/uboone/data/users/kmistry/work/PPFX/uboone/genie/my_events.gst_nuebar.root");
	if (boolfile == false) gSystem->Exit(0);

	// Now get the associated TTree
	bool booltree;
	TTree *gevgen_tree_numu, *gevgen_tree_numubar, *gevgen_tree_nue, *gevgen_tree_nuebar;
	booltree = GetTree(f_gevgen_numu, gevgen_tree_numu, "gst");
	booltree = GetTree(f_gevgen_numubar, gevgen_tree_numubar, "gst");
	booltree = GetTree(f_gevgen_nue, gevgen_tree_nue, "gst");
	booltree = GetTree(f_gevgen_nuebar, gevgen_tree_nuebar, "gst");
	if (booltree == false) gSystem->Exit(0);

	
	double E_numu, E_numubar, E_nue, E_nuebar;
	bool CC_numu, CC_numubar, CC_nue, CC_nuebar;
	gevgen_tree_numu   ->SetBranchAddress("Ev", &E_numu);
	gevgen_tree_numubar->SetBranchAddress("Ev", &E_numubar);
	gevgen_tree_nue    ->SetBranchAddress("Ev", &E_nue);
	gevgen_tree_nuebar ->SetBranchAddress("Ev", &E_nuebar);
	
	gevgen_tree_numu   ->SetBranchAddress("cc", &CC_numu);
	gevgen_tree_numubar->SetBranchAddress("cc", &CC_numubar);
	gevgen_tree_nue    ->SetBranchAddress("cc", &CC_nue);
	gevgen_tree_nuebar ->SetBranchAddress("cc", &CC_nuebar);

	for ( int l = 0; l < gevgen_tree_numu->GetEntries(); l++) {
		gevgen_tree_numu->GetEntry(l);
		if (CC_numu == 1) numuCCHisto_gevgen->Fill(E_numu);
	}
	
	for ( int l = 0; l < gevgen_tree_numubar->GetEntries(); l++) {
		gevgen_tree_numubar->GetEntry(l);
		if (CC_numubar == 1) anumuCCHisto_gevgen->Fill(E_numubar);
	}

	for ( int l = 0; l < gevgen_tree_nue->GetEntries(); l++) {
		gevgen_tree_nue->GetEntry(l);
		if (CC_nue == 1) nueCCHisto_gevgen->Fill(E_nue);
	}
	
	for ( int l = 0; l < gevgen_tree_nuebar->GetEntries(); l++) {
		gevgen_tree_nuebar->GetEntry(l);
		if (CC_nuebar == 1) anueCCHisto_gevgen->Fill(E_nuebar);
	}


	// Get the POT
	double fPOT       = GetPOT(f);
	double fPOT_genie = GetPOT(f_genie, "NuMIEventRates/pottree", "pot");
	double fPOT_genie_nue = GetPOT(f_genie_nue, "NuMIEventRates/pottree", "pot");

	// Get Histograms
	TH1D *h_nue, *h_nuebar, *h_numu, *h_numubar;
	boolhist = GetHist(f, h_nue,     "nue/Detsmear/nue_UW_AV_TPC_5MeV_bin");         if (boolhist == false) gSystem->Exit(0);
	boolhist = GetHist(f, h_nuebar,  "nuebar/Detsmear/nuebar_UW_AV_TPC_5MeV_bin");   if (boolhist == false) gSystem->Exit(0);
	boolhist = GetHist(f, h_numu,    "numu/Detsmear/numu_UW_AV_TPC_5MeV_bin");       if (boolhist == false) gSystem->Exit(0);
	boolhist = GetHist(f, h_numubar, "numubar/Detsmear/numubar_UW_AV_TPC_5MeV_bin"); if (boolhist == false) gSystem->Exit(0);

	// Get Gsimple histograms	
	TH1D *h_nue_gsimp, *h_nuebar_gsimp, *h_numu_gsimp, *h_numubar_gsimp;
	boolhist = GetHist(f_gsimp, h_nue_gsimp,     "nueFluxHisto");   if (boolhist == false) gSystem->Exit(0);
	boolhist = GetHist(f_gsimp, h_nuebar_gsimp,  "anueFluxHisto");  if (boolhist == false) gSystem->Exit(0);
	boolhist = GetHist(f_gsimp, h_numu_gsimp,    "numuFluxHisto");  if (boolhist == false) gSystem->Exit(0);
	boolhist = GetHist(f_gsimp, h_numubar_gsimp, "anumuFluxHisto"); if (boolhist == false) gSystem->Exit(0);
	
	// Get the GENIE histograms
	boolhist = GetHist(f_genie, h_nue_genie,     "NuMIEventRates/Nue_dir/Nue_Energy");           if (boolhist == false) gSystem->Exit(0);
	boolhist = GetHist(f_genie_nue, h_nuebar_genie,  "NuMIEventRates/Nue_bar_dir/Nue_bar_Energy");   if (boolhist == false) gSystem->Exit(0);
	boolhist = GetHist(f_genie, h_numu_genie,    "NuMIEventRates/NuMu_dir/NuMu_Energy");         if (boolhist == false) gSystem->Exit(0);
	boolhist = GetHist(f_genie, h_numubar_genie, "NuMIEventRates/NuMu_bar_dir/NuMu_bar_Energy"); if (boolhist == false) gSystem->Exit(0);


	// Normalise
	// Normalise(hist);

	// Scale
	h_nue->Scale((6.0e20)/ (fPOT*1.0e4) );
	h_nuebar->Scale((6.0e20)/ (fPOT*1.0e4) );
	h_numu->Scale((6.0e20)/ (fPOT*1.0e4) );
	h_numubar->Scale((6.0e20)/ (fPOT*1.0e4) );

	// Scale the genie histograms
	h_nue_genie->Scale(  (6.0e20)/ (fPOT_genie) );
	h_nuebar_genie->Scale(  (6.0e20)/ (fPOT_genie_nue) );
	h_numu_genie->Scale(  (6.0e20)/ (fPOT_genie) );
	h_numubar_genie->Scale(  (6.0e20)/ (fPOT_genie) );
	
	const char* genieXsecPath = gSystem->ExpandPathName("$(GENIEXSECPATH)");
	if ( !genieXsecPath ) {
		std::cout << "$(GENIEXSECPATH) not defined." << std::endl;
		std::cout << "Please setup *genie_xsec*." << std::endl; 
	}

	if ( genieXsecPath ) {
		TString genieXsecFileName = genieXsecPath;
		genieXsecFileName += "/xsec_graphs.root";
		TFile *genieXsecFile = new TFile(genieXsecFileName,"READ");
		genieXsecNumuCC    = (TGraph *) genieXsecFile->Get("nu_mu_Ar40/tot_cc");
		genieXsecNumubarCC = (TGraph *) genieXsecFile->Get("nu_mu_bar_Ar40/tot_cc");
		genieXsecNueCC     = (TGraph *) genieXsecFile->Get("nu_e_Ar40/tot_cc");
		genieXsecNuebarCC  = (TGraph *) genieXsecFile->Get("nu_e_bar_Ar40/tot_cc");
		genieXsecFile->Close();

		TSpline3* genieXsecSplineNumuCC = new TSpline3("genieXsecSplineNumuCC", genieXsecNumuCC, "", 0,20);

		double value;
		for(int i=1; i<histNbins+1; i++) {

			// std::cout << i<< "  "<< h_numu->GetBinContent(i) << "  " << h_numu->GetBinCenter(i) << "   " << h_numu->GetBinContent(i+10)  << "  "<< h_numu->GetBinCenter(i+10)<< std::endl;
			
			// Nue
			value = h_nue->GetBinContent(i);
			value *= genieXsecNueCC->Eval(h_nue->GetBinCenter(i)); // Eval implies linear interpolation
			value *= (1e-38 * Ntarget/40.); // 1/40 is due to I'm considering nu_mu_Ar40.
			nueCCHisto->SetBinContent(i, value);

			// Nuebar
			value = h_nuebar->GetBinContent(i);
			value *= genieXsecNuebarCC->Eval(h_nuebar->GetBinCenter(i)); // Eval implies linear interpolation
			value *= (1e-38 * Ntarget/40.); // 1/40 is due to I'm considering nu_mu_Ar40.
			anueCCHisto->SetBinContent(i, value);

			// Numu
			value = h_numu->GetBinContent(i);
			value *= genieXsecNumuCC->Eval(h_numu->GetBinCenter(i)); // Eval implies linear interpolation
			// value *= bin_average_total_xsec(h_numu, genieXsecNumuCC, i);
			// value *= genieXsecSplineNumuCC->Eval(h_numu->GetBinCenter(i));
			// value *= genieXsecSplineNumuCC->Eval(h_numu->GetBinLowEdge(i+10));
			value *= (1e-38 * Ntarget/40.); // 1/40 is due to I'm considering nu_mu_Ar40.
			numuCCHisto->SetBinContent(i, value);

			// Numubar
			value = h_numubar->GetBinContent(i);
			value *= genieXsecNumubarCC->Eval(h_numubar->GetBinCenter(i)); // Eval implies linear interpolation
			value *= (1e-38 * Ntarget/40.); // 1/40 is due to I'm considering nu_mu_Ar40.
			anumuCCHisto->SetBinContent(i, value);

			// ------------------------------  GSimple  -------------------------------------------------
			// Nue
			value = h_nue_gsimp->GetBinContent(i);
			value *= genieXsecNueCC->Eval(h_nue_gsimp->GetBinCenter(i)); // Eval implies linear interpolation
			value *= (1e-38 * Ntarget/40.); // 1/40 is due to I'm considering nu_mu_Ar40.
			nueCCHisto_gsimp->SetBinContent(i, value);

			// Nuebar
			value = h_nuebar_gsimp->GetBinContent(i);
			value *= genieXsecNuebarCC->Eval(h_nuebar_gsimp->GetBinCenter(i)); // Eval implies linear interpolation
			value *= (1e-38 * Ntarget/40.); // 1/40 is due to I'm considering nu_mu_Ar40.
			anueCCHisto_gsimp->SetBinContent(i, value);

			// Numu
			value = h_numu_gsimp->GetBinContent(i);
			value *= genieXsecNumuCC->Eval(h_numu_gsimp->GetBinCenter(i)); // Eval implies linear interpolation
			value *= (1e-38 * Ntarget/40.); // 1/40 is due to I'm considering nu_mu_Ar40.
			numuCCHisto_gsimp->SetBinContent(i, value);

			// Numubar
			value = h_numubar_gsimp->GetBinContent(i);
			value *= genieXsecNumubarCC->Eval(h_numubar_gsimp->GetBinCenter(i)); // Eval implies linear interpolation
			value *= (1e-38 * Ntarget/40.); // 1/40 is due to I'm considering nu_mu_Ar40.
			anumuCCHisto_gsimp->SetBinContent(i, value);
			// ---------------------------------------------------------------------------------------------

		}
	} // end if ( genieXsecPath )
	
	TLegend* leg = new TLegend(0.75, 0.65, 0.95, 0.9);
	leg->SetNColumns(1);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetTextFont(62);
	leg->AddEntry(numuCCHisto, "#nu_{#mu}","l");
	leg->AddEntry(anumuCCHisto, "#bar{#nu_{#mu}}","l");
	leg->AddEntry(nueCCHisto, "#nu_{e} ","l");
	leg->AddEntry(anueCCHisto, "#bar{#nu_{e}} ","l");

	TLegend* leg_gsimp = new TLegend(0.75, 0.65, 0.95, 0.9);
	leg_gsimp->SetNColumns(1);
	leg_gsimp->SetBorderSize(0);
	leg_gsimp->SetFillStyle(0);
	leg_gsimp->SetTextFont(62);
	leg_gsimp->SetTextSize(0.04);

	// Rebin
	nueCCHisto        ->Rebin(rebin);
	anueCCHisto       ->Rebin(rebin);
	numuCCHisto       ->Rebin(rebin);
	anumuCCHisto      ->Rebin(rebin);
	nueCCHisto_gsimp  ->Rebin(rebin);
	anueCCHisto_gsimp ->Rebin(rebin);
	numuCCHisto_gsimp ->Rebin(rebin);
	anumuCCHisto_gsimp->Rebin(rebin);
	nueCCHisto_gevgen  ->Rebin(rebin);
	anueCCHisto_gevgen ->Rebin(rebin);
	numuCCHisto_gevgen ->Rebin(rebin);
	anumuCCHisto_gevgen->Rebin(rebin);

	std::cout <<nueCCHisto->Integral(0, -1) / h_nue_genie->Integral(0, -1)<< std::endl;

	// Area normalise to check the shape 
	h_nue_genie    ->Scale( nueCCHisto  ->Integral(0, -1) / h_nue_genie    ->Integral(0, -1) );
	h_nuebar_genie ->Scale( anueCCHisto ->Integral(0, 30) / h_nuebar_genie ->Integral(0, 30) );
	h_numu_genie   ->Scale( numuCCHisto ->Integral(0, -1) / h_numu_genie   ->Integral(0, h_numu_genie->GetNbinsX()-5) );
	h_numubar_genie->Scale( anumuCCHisto->Integral(0, -1) / h_numubar_genie->Integral(0, -1) );

	numuCCHisto_gevgen ->Scale( numuCCHisto ->Integral(0, -1)  / numuCCHisto_gevgen ->Integral(0,-1) );
	anumuCCHisto_gevgen->Scale( anumuCCHisto->Integral(0, -1)  / anumuCCHisto_gevgen->Integral(0,-1) );
	nueCCHisto_gevgen  ->Scale( nueCCHisto  ->Integral(0, -1)  / nueCCHisto_gevgen  ->Integral(0,-1) );
	anueCCHisto_gevgen ->Scale( anueCCHisto ->Integral(0, -1)  / anueCCHisto_gevgen ->Integral(0,-1) );

	

	// create plots folder if it does not exist
	gSystem->Exec("if [ ! -d \"plots/Event_Rates\" ]; then echo \"\n Event_Rates folder does not exist... creating\"; mkdir -p plots/Event_Rates; fi"); 

	// -----------------------------  Nue --------------------------------------
	TCanvas* c_nue = new TCanvas();
	DrawSpecifiers(c_nue, nueCCHisto, "nue");
	DrawSpecifiers(c_nue, nueCCHisto_gsimp, "nue");
	DrawSpecifiers(c_nue, h_nue_genie, "nue");
	DrawSpecifiers(c_nue, nueCCHisto_gevgen, "nue"); // Gevgen
	
	nueCCHisto_gevgen->SetLineColor(41);
	h_nue_genie->SetLineColor(40);
	nueCCHisto_gsimp->SetLineColor(kGreen+1);
	nueCCHisto->GetYaxis()->SetRangeUser(0,200);
	
	if (!strcmp(horn,"fhc")) nueCCHisto_gevgen->Draw("his,same");
	nueCCHisto->Draw("his, same");
	nueCCHisto_gsimp->Draw("his,same");
	if (!strcmp(horn,"fhc")) h_nue_genie->Draw("his,same");
	

	Draw_Nu_Mode(c_nue, horn); // Draw FHC Mode/RHC Mode Text
	leg_gsimp->AddEntry(nueCCHisto, "dk2nu","l");
	leg_gsimp->AddEntry(nueCCHisto_gsimp, "flugg","l");
	if (!strcmp(horn,"fhc")) leg_gsimp->AddEntry(h_nue_genie, "genie","l");
	if (!strcmp(horn,"fhc")) leg_gsimp->AddEntry(nueCCHisto_gevgen, "gevgen","l");
	leg_gsimp->Draw();
	gStyle->SetTitleH(0.07);
	c_nue->Print(Form("plots/Event_Rates/Event_Rate_Prediction_%s_%s.pdf", horn, "nue"));
	leg_gsimp->Clear();
	
	
	// --------------------------- Nuebar --------------------------------------
	TCanvas* c_nuebar = new TCanvas();
	DrawSpecifiers(c_nuebar, anueCCHisto, "nuebar");
	DrawSpecifiers(c_nuebar, anueCCHisto_gsimp, "nuebar");
	DrawSpecifiers(c_nuebar, h_nuebar_genie, "nuebar");
	DrawSpecifiers(c_nuebar, anueCCHisto_gevgen, "nuebar"); // Gevgen
	
	anueCCHisto_gevgen->SetLineColor(41);
	h_nuebar_genie->SetLineColor(40);
	anueCCHisto_gsimp->SetLineColor(kGreen+1);
	anueCCHisto->GetYaxis()->SetRangeUser(0,30);
	
	if (!strcmp(horn,"fhc")) anueCCHisto_gevgen->Draw("his,same");
	anueCCHisto->Draw("his, same");
	anueCCHisto_gsimp->Draw("his,same");
	if (!strcmp(horn,"fhc")) h_nuebar_genie->Draw("his,same");
	
	Draw_Nu_Mode(c_nuebar, horn); // Draw FHC Mode/RHC Mode Text
	leg_gsimp->AddEntry(anueCCHisto, "dk2nu","l");
	leg_gsimp->AddEntry(anueCCHisto_gsimp, "flugg","l");
	if (!strcmp(horn,"fhc")) leg_gsimp->AddEntry(h_nuebar_genie, "genie","l");
	if (!strcmp(horn,"fhc")) leg_gsimp->AddEntry(anueCCHisto_gevgen, "gevgen","l");
	leg_gsimp->Draw();
	gStyle->SetTitleH(0.07);
	c_nuebar->Print(Form("plots/Event_Rates/Event_Rate_Prediction_%s_%s.pdf", horn, "nuebar"));
	leg_gsimp->Clear();
	
	// ----------------------------  Numu --------------------------------------
	TCanvas* c_numu = new TCanvas();
	DrawSpecifiers(c_numu, numuCCHisto, "numu");
	DrawSpecifiers(c_numu, numuCCHisto_gsimp, "numu");
	DrawSpecifiers(c_numu, h_numu_genie, "numu");
	DrawSpecifiers(c_numu, numuCCHisto_gevgen, "numu"); // Gevgen
	
	numuCCHisto_gevgen->SetLineColor(41);
	h_numu_genie->SetLineColor(40);
	numuCCHisto_gsimp->SetLineColor(kGreen+1);
	numuCCHisto->GetYaxis()->SetRangeUser(0,2000);
	
	if (!strcmp(horn,"fhc")) numuCCHisto_gevgen->Draw("his,same");
	numuCCHisto->Draw("his, same");
	numuCCHisto_gsimp->Draw("his,same");
	if (!strcmp(horn,"fhc")) h_numu_genie->Draw("his,same");
	
	
	Draw_Nu_Mode(c_numu, horn); // Draw FHC Mode/RHC Mode Text
	leg_gsimp->AddEntry(numuCCHisto, "dk2nu","l");
	leg_gsimp->AddEntry(numuCCHisto_gsimp, "flugg","l");
	if (!strcmp(horn,"fhc")) leg_gsimp->AddEntry(h_numu_genie, "genie","l");
	if (!strcmp(horn,"fhc")) leg_gsimp->AddEntry(numuCCHisto_gevgen, "gevgen","l");
	leg_gsimp->Draw();
	gStyle->SetTitleH(0.07);
	c_numu->Print(Form("plots/Event_Rates/Event_Rate_Prediction_%s_%s.pdf", horn, "numu"));
	leg_gsimp->Clear();
	
	// ---------------------------  Numubar ------------------------------------
	TCanvas* c_numubar = new TCanvas();
	DrawSpecifiers(c_numubar, anumuCCHisto, "numubar");
	DrawSpecifiers(c_numubar, anumuCCHisto_gsimp, "numubar");
	DrawSpecifiers(c_numubar, h_numubar_genie, "numubar");
	DrawSpecifiers(c_numubar, anumuCCHisto_gevgen, "numubar"); // Gevgen

	anumuCCHisto_gevgen->SetLineColor(41);
	h_numubar_genie->SetLineColor(40);
	anumuCCHisto_gsimp->SetLineColor(kGreen+1);
	anumuCCHisto->GetYaxis()->SetRangeUser(0,400);
	
	if (!strcmp(horn,"fhc")) anumuCCHisto_gevgen->Draw("his,same");
	anumuCCHisto->Draw("his, same");
	anumuCCHisto_gsimp->Draw("his,same");
	if (!strcmp(horn,"fhc")) h_numubar_genie->Draw("his,same");
	
	Draw_Nu_Mode(c_numubar, horn); // Draw FHC Mode/RHC Mode Text
	leg_gsimp->AddEntry(anumuCCHisto, "dk2nu","l");
	leg_gsimp->AddEntry(anumuCCHisto_gsimp, "flugg","l");
	if (!strcmp(horn,"fhc")) leg_gsimp->AddEntry(h_numubar_genie, "genie","l");
	if (!strcmp(horn,"fhc")) leg_gsimp->AddEntry(anumuCCHisto_gevgen, "gevgen","l");
	leg_gsimp->Draw();
	gStyle->SetTitleH(0.07);
	c_numubar->Print(Form("plots/Event_Rates/Event_Rate_Prediction_%s_%s.pdf", horn, "numubar"));
	leg_gsimp->Clear();
	
	// -----------------------------  ALL --------------------------------------
	TCanvas* c_all = new TCanvas();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.12);
	numuCCHisto->SetTitle(";#nu Energy [GeV]; #nu CC / 79 t / 6 #times 10^{20} POT / 25 MeV");
	numuCCHisto->Draw("hist,same");
	anumuCCHisto->Draw("hist,same");
	nueCCHisto->Draw("hist,same");
	anueCCHisto->Draw("hist,same");
	Draw_Nu_Mode(c_all, horn); // Draw FHC Mode/RHC Mode Text
	leg->Draw();
	c_all->Print(Form("plots/Event_Rates/Event_Rate_Prediction_%s_all.pdf", horn));
	gPad->SetLogy();
	c_all->Print(Form("plots/Event_Rates/Event_Rate_Prediction_%s_all_logy.pdf", horn));

	// gSystem->Exit(0);

} // end of main




