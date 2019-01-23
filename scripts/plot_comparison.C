/**
 * NuMI Validation Plotting
 *
 * A. Mastbaum <mastbaum@uchicago.edu> 2018/11
 * Updated to plot multiple nu modes in a single plot
 * by Krishan Mistry
 */

// Fucnction that grabs the reweighted histogram names for plotting
std::vector<std::string> loopdir(TString  inputfile, TString mode) {
	std::vector<std::string> inputmode;

	TFile *f1 = TFile::Open(inputfile);
	f1->cd(mode);
	
	TKey *key;
	TIter nextkey(gDirectory->GetListOfKeys());

	std::cout << "\n=================================================" << std::endl;	
	std::cout << "Using input modes:" << std::endl;	
  	while ( ( key =  (TKey*)nextkey()) ) { // extra brackets to omit a warning 
    	if (key->IsFolder()) {
			std::cout << key->GetName() << std::endl; // print hte input modes
			inputmode.push_back(key->GetName());
		}
	}
	std::cout << "=================================================\n" << std::endl;

	return (inputmode);
}

// function that makes a legend for multiple histograms and draws them to the canvas
void legDraw(TLegend *legend, TH1D *hist, TString prodmode, std::string inputmode, TString mode){
	
	if (inputmode == "PPFXMaster"){

		if (mode == "numu"){
			hist->SetLineColor(kOrange+8);
			hist->SetLineWidth(2);
			hist->SetLineStyle(2);
			legend->AddEntry(hist, "#nu_{#mu}", "l");
			hist->Draw("hist,same");
		}
		else if (mode == "nue") {
			hist->SetLineColor(kOrange+1);
			hist->SetLineWidth(2);
			hist->SetLineStyle(2);
			legend->AddEntry(hist, "#nu_{e}", "l");
			hist->Draw("hist,same");
		}
		else if (mode == "numubar") {
			hist->SetLineColor(kMagenta-3);
			hist->SetLineWidth(2);
			hist->SetLineStyle(2);
			legend->AddEntry(hist, "#bar{#nu_{#mu}}", "l");
			hist->Draw("hist,same");

		}
		else if(mode == "nuebar"){
			hist->SetLineColor(kMagenta+3);
			hist->SetLineWidth(2);
			hist->SetLineStyle(2);
			legend->AddEntry(hist, "#bar{#nu_{e}}", "l");
			hist->Draw("hist,same");
		}
		else return;
	} 
	else return;

}

// Enumbers for the input mode 
enum e_mode{ enumu, enue, enumubar, enuebar};

// Function to retun enum from mode label
e_mode return_mode(TString mode){
		if (mode == "numu")    return enumu;
		if (mode == "nue")     return enue;
		if (mode == "numubar") return enumubar;
		if (mode == "nuebar")  return enuebar;
		else return enumu;

}

// Function to get the fractional errors for each nu mode
void GetFracErrors(TString mipp, TString inputfile, TString prodmode, TString mode, TCanvas* c, TLegend* leg ){
  	
 	std::vector<std::string> inputmode = loopdir(inputfile, mode); // Grab the names

	TString Getmode;
	TString Gethist_TPC;
	TString Getflux;
	TString Cov_names;
	TString Err_names;

	// Select neutrino type to run with 
	switch (return_mode(mode)){
		case enumu:
			std::cout << "\nUsing NuMu Mode!\n" << std::endl;
			Getmode = "numu"; 												// Folder name
			Gethist_TPC = "numu_CV_AV_TPC";									// AV in TPC flux prediction
			Getflux = "flux_numu";											// CV flux from NOvA
			Cov_names = "numu/%s/Active_TPC_Volume/numu_%s_Uni_%i_AV_TPC";  // Covariance matrix names
			Err_names = "fractional_uncertainty_numu";						// NOvA fractional uncertainties plot
			break;

		case enue:
			std::cout << "\nUsing Nue Mode!\n" << std::endl;
			Getmode = "nue";
			Gethist_TPC = "nue_CV_AV_TPC";
			Getflux = "flux_nue";
			Cov_names = "nue/%s/Active_TPC_Volume/nue_%s_Uni_%i_AV_TPC";
			Err_names = "fractional_uncertainty_nue";
			break;

		case enumubar:
			std::cout << "\nUsing NuMubar Mode!\n" << std::endl;
			Getmode = "numubar";
			Gethist_TPC = "numubar_CV_AV_TPC";
			Getflux = "flux_numubar";
			Cov_names = "numubar/%s/Active_TPC_Volume/numubar_%s_Uni_%i_AV_TPC";
			Err_names = "fractional_uncertainty_numubar";
			break;

		case enuebar:
			std::cout << "\nUsing Nuebar Mode!\n" << std::endl;
			Getmode = "nuebar";
			Gethist_TPC = "nuebar_CV_AV_TPC";
			Getflux = "flux_nuebar";
			Cov_names = "nuebar/%s/Active_TPC_Volume/nuebar_%s_Uni_%i_AV_TPC";
			Err_names = "fractional_uncertainty_nuebar";
			break;

	} 
	

	// Root is dumb and so need to pre-decalre  some stuff here
	TDirectory* d;
	TH1D* hCV_Flux;

	// ++++++++++++++++++++++++++++++++++
	// Get the CV flux
	// ++++++++++++++++++++++++++++++++++
 
	TFile* f1 = TFile::Open(inputfile);

	d = (TDirectory*)f1->Get(Getmode);

	d->cd();
	hCV_Flux = (TH1D*) (gDirectory->Get(Gethist_TPC)->Clone("fx"));

	hCV_Flux->SetDirectory(0);
	
	// Normalise flux by bin width (gives a flux/E [GeV])
	for (int i=0;i<hCV_Flux->GetNbinsX()+1;i++) {
		hCV_Flux->SetBinContent(i, hCV_Flux->GetBinContent(i)/hCV_Flux->GetBinWidth(i));
	}
	
	hCV_Flux->Sumw2();
	TH1D* horig = (TH1D*) hCV_Flux->Clone("horig");
	
  	// ++++++++++++++++++++++++++++++++++
	// Correlations, Covariance & uncertainties
	// ++++++++++++++++++++++++++++++++++

	f1->cd();
	const int nuni = 100; // num universes
	const int nbins = hCV_Flux->GetNbinsX();
	double* edges = new double[nbins+1];

	// Set bin widths to be the same as NOvA
	for (int i=1; i<nbins+1; i++) {
		edges[i-1] = hCV_Flux->GetBinLowEdge(i);
	}
	edges[nbins] = hCV_Flux->GetBinLowEdge(nbins-1) + 2 * (hCV_Flux->GetBinWidth(nbins-1));
	
	// More pre-declarations go here
	TH1D* hu; // Flux hist for each universe
	
	TH2D* cov;	// Covariance
	TH1D* herr ; // Fractional Uncertenties 

	cov  = new TH2D("cov" , ";E_{#nu} (GeV);E_{#nu} (GeV)", nbins, edges, nbins, edges);
	herr = new TH1D("herr", ";E_{#nu};Fractional Uncertainty", nbins, edges);
	
		
  // Covariance matrix only for Masterweight mode
  for (int k=0; k<nuni; k++) {
    char name[500];
    snprintf(name, 500, Cov_names ,"PPFXMaster","PPFXMaster", k); 

    hu = (TH1D*) f1->Get(name);
    
    for (int m=0; m<nbins; m++) {
      hu->SetBinContent(m, hu->GetBinContent(m) / hu->GetBinWidth(m));
    }

    for (int i=1; i<nbins+1; i++) {

		double cvi = horig->GetBinContent(i); // CV bin i
		double uvi = hu->GetBinContent(i); // Univ bin i 

		for (int j=1; j<nbins+1; j++) {
			
			double cvj = horig->GetBinContent(j); // CV bin j
			double uvj = hu->GetBinContent(j);    // Univ bin j 

			double c = (uvi - cvi) * (uvj - cvj);

			cov->SetBinContent(i, j, cov->GetBinContent(i, j) + c / nuni); // Fill with variance
		}
    }
    hu->Reset();
  }

  double cii{0}, cjj{0}, n{1}, horig_cont{0};

  //  Fractional errors
  // loop over rows
  for (int i=1; i<nbins+1; i++) {
    cii = cov->GetBinContent(i, i);

    // Catch zeros , set to arbitary 1.0
    if (horig->GetBinContent(i) <= 0) horig_cont = 1.0;
    else horig_cont = horig->GetBinContent(i);

    herr->SetBinContent(i, sqrt(cii) / horig_cont);

  }

	// Plot fractional errors overlaid with official NOvA plot
	c->cd();

	// Make the plot
	legDraw(leg, herr, prodmode, "PPFXMaster", mode);
	
	herr->GetYaxis()->SetRangeUser(0,0.35);
	// herr->GetYaxis()->SetRangeUser(0,1.75);
		
}

TString func_return_mode(TString mode){
	
	TString Err_names;
	// Select neutrino type to run with 
	switch (return_mode(mode)){
		case enumu:
			Err_names = "fractional_uncertainty_numu";						// NOvA fractional uncertainties plot
			break;

		case enue:
			Err_names = "fractional_uncertainty_nue";
			break;

		case enumubar:
			Err_names = "fractional_uncertainty_numubar";
			break;

		case enuebar:
			Err_names = "fractional_uncertainty_nuebar";
			break;
	} 

	return Err_names;
}

// Main
void plot_comparison( TString mipp, TString inputfile, TString prodmode, TString wplot, TString mode) { // (mippon/mippoff, input, Product/noThinKaon etc. numu/nue)
	gStyle->SetOptStat(0); // say no to stats box

	TCanvas* c1 = new TCanvas();
	TCanvas* c2 = new TCanvas();
	TH1D* herr2;

	TString Err_names;

	TLegend* lfrac1 = new TLegend(0.2, 0.65, 0.6, 0.9);
	lfrac1->SetNColumns(1);
	lfrac1->SetBorderSize(0);
	lfrac1->SetFillStyle(0);
	lfrac1->SetTextFont(62); 

	TLegend* lfrac2 = new TLegend(0.2, 0.65, 0.6, 0.9);
	lfrac2->SetNColumns(1);
	lfrac2->SetBorderSize(0);
	lfrac2->SetFillStyle(0);
	lfrac2->SetTextFont(62); 

	// Now get the NOvA uncertainties
	TFile* f2 = TFile::Open("/uboone/app/users/kmistry/PPFX/numi-validation/nova_flux/FHC_Flux_NOvA_ND_2017.root");
	f2->cd();

	mode = "numu";
	Err_names = func_return_mode(mode);
	GetFracErrors(mipp, inputfile, prodmode, mode, c1, lfrac1 );
	herr2 = (TH1D*) f2->Get(Err_names);
	herr2->SetLineColor(kOrange+8);
	herr2->SetLineWidth(2);
	herr2->Draw("same");
	lfrac1->AddEntry(herr2, "#nu_{#mu} NOvA","l");
	
	mode = "numubar";
	Err_names = func_return_mode(mode);
	GetFracErrors(mipp, inputfile, prodmode, mode, c1, lfrac1 );
	herr2 = (TH1D*) f2->Get(Err_names);
	herr2 = (TH1D*) f2->Get(Err_names);
	herr2->SetLineColor(kMagenta-3);
	herr2->SetLineWidth(2);
	herr2->Draw("same");
	lfrac1->AddEntry(herr2, "#bar{#nu_{#mu}} NOvA","l");


	mode = "nue";
	Err_names = func_return_mode(mode);
	GetFracErrors(mipp, inputfile, prodmode, mode, c2, lfrac2 );
	herr2 = (TH1D*) f2->Get(Err_names);
	herr2->SetLineColor(kOrange+1);
	herr2->SetLineWidth(2);
	herr2->Draw("same");
	lfrac2->AddEntry(herr2, "#nu_{e} NOvA","l");
	
	mode = "nuebar";
	Err_names = func_return_mode(mode);
	GetFracErrors(mipp, inputfile, prodmode, mode, c2, lfrac2 );
	herr2 = (TH1D*) f2->Get(Err_names);
	herr2 = (TH1D*) f2->Get(Err_names);
	herr2->SetLineColor(kMagenta+3);
	herr2->SetLineWidth(2);
	herr2->Draw("same");
	lfrac2->AddEntry(herr2, "#bar{#nu_{e}} NOvA","l");

	// Draw the legend
	c1->cd();
	lfrac1->Draw();

	c2->cd();
	lfrac2->Draw();

	// create plots folder if it does not exist
	gSystem->Exec("if [ ! -d \"plots\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir plots; fi"); 


	c1->Print("plots/Fractional_Uncertainties_Master_NuMu_NuMubar_MIPPOff.pdf");
	c2->Print("plots/Fractional_Uncertainties_Master_Nue_Nuebar_MIPPOff.pdf");

} // end of main

