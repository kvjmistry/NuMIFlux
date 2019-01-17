/**
 * NuMI Validation Plotting
 *
 * Plots each individual weighting mode instead of a single one
 * does not plot the correlation matrix for each individual mode to speed up the time
 * 
 * A. Mastbaum <mastbaum@uchicago.edu> 2018/11
 * Modified by K. Mistry 12/18
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

enum e_mode{ enumu, enue, enumubar, enuebar};


e_mode return_mode(TString mode){
		if (mode == "numu")    return enumu;
		if (mode == "nue")     return enue;
		if (mode == "numubar") return enumubar;
		if (mode == "nuebar")  return enuebar;
		else return enumu;

}


// Main
void plot_comparison_all( TString mipp, TString inputfile, TString prodmode, TString mode) { // (mippon/mippoff, input, Product/noThinKaon etc. numu/nue)
	gStyle->SetOptStat(0); // say no to stats box

	std::vector<std::string> inputmode = loopdir(inputfile, mode); // Grab the names

	TString Getmode;
	TString Gethist_TPC;
	TString Getflux;
	TString Cov_names;
	TString Err_names;

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
	TH1D* hNOvA_CV_Flux;

	// ++++++++++++++++++++++++++++++++++
	// Overlay of output plot with official NOvA FHC numu flux
	// ++++++++++++++++++++++++++++++++++
 
	TFile* f1 = TFile::Open(inputfile);
	TCanvas* c1 = new TCanvas();

	d = (TDirectory*)f1->Get(Getmode);

	d->cd();
	hCV_Flux = (TH1D*) (gDirectory->Get(Gethist_TPC)->Clone("fx"));

	hCV_Flux->SetDirectory(0);
	
	// Normalise flux by bin width (gives a flux/E [GeV])
	for (int i=0;i<hCV_Flux->GetNbinsX()+1;i++) {
		hCV_Flux->SetBinContent(i, hCV_Flux->GetBinContent(i)/hCV_Flux->GetBinWidth(i));
	}
	
	hCV_Flux->Sumw2();
	hCV_Flux->SetLineColor(kRed);
	hCV_Flux->SetLineWidth(2);
	hCV_Flux->SetTitle(";E_{#nu} (GeV);Fraction/GeV");
	hCV_Flux->Draw("");
	TH1D* horig = (TH1D*) hCV_Flux->Clone("horig");

	// Now get the NOvA CV Flux
	TFile* f2 = TFile::Open("/uboone/app/users/kmistry/PPFX/numi-validation/nova_flux/FHC_Flux_NOvA_ND_2017.root");
	// f2->ls();
	
	hNOvA_CV_Flux = (TH1D*) f2->Get(Getflux);

	hNOvA_CV_Flux->SetLineColor(kBlue);
	hNOvA_CV_Flux->SetLineWidth(2);

	hNOvA_CV_Flux->Draw("same");
	hCV_Flux->GetYaxis()->SetTitle(hNOvA_CV_Flux->GetYaxis()->GetTitle());
	
	// Need to divide by area front face is 12.39 m2, looks closer to 14.6
	hCV_Flux->Scale(1.0e6/(12.39*5.0e5*497)); 
	
	std::cout << "norm factor:\t" << hCV_Flux->Integral(2,-1)/hNOvA_CV_Flux->Integral(2,-1) << std::endl;
	
	hCV_Flux->Scale(hNOvA_CV_Flux->Integral(3,-1)/hCV_Flux->Integral(3,-1));
	

	TH1D* hratio = (TH1D*) hCV_Flux->Clone("hratio");
	hratio->SetDirectory(0);

	TLegend* l = new TLegend(0.5, 0.6, 0.7, 0.8);
	l->AddEntry(hNOvA_CV_Flux, "NOvA","l");
	l->AddEntry(hCV_Flux, "Our Prediction","l");
	l->Draw();
	l->SetBorderSize(0);
	l->SetFillStyle(0);

	gPad->SetLogy();
	gPad->Update();

	// ++++++++++++++++++++++++++++++++++
	// Ratio of ours to NOvA
	// ++++++++++++++++++++++++++++++++++
	TCanvas* c2 = new TCanvas();
	hratio->Divide(hNOvA_CV_Flux);
	hratio->SetLineColor(kBlack);
	hratio->SetMarkerStyle(7);
	hratio->Draw("e1");
	hratio->SetTitle(";E_{#nu} (GeV);Ratio to NOvA");
	hratio->GetYaxis()->SetRangeUser(0.85, 1.15);
	TLine* flat = new TLine(0, 1, 20, 1);
	flat->SetLineStyle(7);
	flat->Draw();

	// ++++++++++++++++++++++++++++++++++
	// Correlations & uncertainties
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
	TCanvas* c3 = new TCanvas();
	TCanvas* c4 = new TCanvas();
	TH1D* hu; // Flux hist for each universe
	TH1D* herr2;
	std::vector<TH2D*> cov;
	std::vector<TH2D*> cor;
	std::vector<TH1D*> herr ; // Fractional Uncertenties 
	TLegend* lfrac = new TLegend(0.5, 0.65, 0.9, 0.9);
	lfrac->SetNColumns(3);
	lfrac->SetBorderSize(0);
	lfrac->SetFillStyle(0);
	lfrac->SetTextFont(62); 

	cor.resize(inputmode.size());
	cov.resize(inputmode.size());
	herr.resize(inputmode.size());
	
	for (unsigned int l = 0; l < inputmode.size(); l++){
		cor[l]  = new TH2D(Form("%s_cor",inputmode[l].c_str()), ";E_{#nu} (GeV);E_{#nu} (GeV)", nbins, edges, nbins, edges);
		cov[l]  = new TH2D(Form("%s_cov",inputmode[l].c_str()), ";E_{#nu} (GeV);E_{#nu} (GeV)", nbins, edges, nbins, edges);
		herr[l] = new TH1D(Form("%s_herr",inputmode[l].c_str()),";E_{#nu};Fractional Uncertainty", nbins, edges);
	}

	// Loop over all input modes, get cov matrix and then get fractional uncertainties
	for (unsigned int l = 0; l < inputmode.size(); l++){
		
		// Covariance matrix
		for (int k=0; k<nuni; k++) {
			char name[500];
			snprintf(name, 500, Cov_names ,inputmode[l].c_str(),inputmode[l].c_str(), k); 
  

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

					cov[l]->SetBinContent(i, j, cov[l]->GetBinContent(i, j) + c / nuni); // Fill with variance
				}
			}
			hu->Reset();
		}

		double cii{0}, cjj{0}, n{1}, cor_bini{0}, horig_cont{0};

		// Correlation matrix & fractional errors
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

		// Plot correlations
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

		if (inputmode[l] == "PPFXMaster"){
			herr[l]->SetLineColor(kBlack);
			herr[l]->SetLineWidth(2);
			inputmode[l].erase(0,4); // Trim PPFX off the name
			lfrac->AddEntry(herr[l], inputmode[l].c_str(), "l");
			herr[l]->Draw("hist,same");
		} 
		else if  (inputmode[l] == "ms_PPFX"){
			herr[l]->SetLineColor(kMagenta+2);
			herr[l]->SetLineWidth(2);
			lfrac->AddEntry(herr[l], prodmode, "l"); // prodmode is overridden in the terminal input
			herr[l]->Draw("hist,same");
		}
		else if  (inputmode[l] == "PPFXMIPPKaon" && mipp =="mippon"){ // only turn on if there is MIPP
			herr[l]->SetLineColor(30);
			herr[l]->SetLineWidth(2);
			inputmode[l].erase(0,4);
			lfrac->AddEntry(herr[l], inputmode[l].c_str(), "l");
			herr[l]->SetLineStyle(2);
			herr[l]->Draw("hist,same");
		}
		else if  ( inputmode[l] == "PPFXMIPPPion" && mipp =="mippon"){ // only turn on if there is MIPP
			herr[l]->SetLineColor(38);
			herr[l]->SetLineWidth(2);
			inputmode[l].erase(0,4);
			lfrac->AddEntry(herr[l], inputmode[l].c_str(), "l");
			herr[l]->SetLineStyle(2);
			herr[l]->Draw("hist,same");
		}
		else if  (inputmode[l] == "PPFXOther"){
			herr[l]->SetLineColor(28);
			herr[l]->SetLineWidth(2);
			inputmode[l].erase(0,4);
			lfrac->AddEntry(herr[l], inputmode[l].c_str(), "l");
			herr[l]->SetLineStyle(2);
			herr[l]->Draw("hist,same");
			
		}
		else if  (inputmode[l] == "PPFXTargAtten" ){
			herr[l]->SetLineColor(36);
			herr[l]->SetLineWidth(2);
			inputmode[l].erase(0,4);
			lfrac->AddEntry(herr[l], inputmode[l].c_str(), "l");
			herr[l]->SetLineStyle(2);
			herr[l]->Draw("hist,same");
		}
		else if  (inputmode[l] == "PPFXThinKaon"){
			herr[l]->SetLineColor(1001);
			herr[l]->SetLineWidth(2);
			inputmode[l].erase(0,4);
			lfrac->AddEntry(herr[l], "pC #rightarrow KX", "l");
			herr[l]->SetLineStyle(2);
			herr[l]->Draw("hist,same");
		}
		else if  (inputmode[l] == "PPFXThinMeson"){
			herr[l]->SetLineColor(46);
			herr[l]->SetLineWidth(2);
			inputmode[l].erase(0,4);
			lfrac->AddEntry(herr[l], "Meson Incident.", "l");
			herr[l]->SetLineStyle(2);
			herr[l]->Draw("hist,same");
		}
		else if  (inputmode[l] == "PPFXThinNeutron"){
			herr[l]->SetLineColor(42);
			herr[l]->SetLineWidth(2);
			inputmode[l].erase(0,4);
			lfrac->AddEntry(herr[l], "nC #rightarrow #piX", "l");
			herr[l]->SetLineStyle(2);
			herr[l]->Draw("hist,same");
		}
		else if  (inputmode[l] == "PPFXThinNucA"){
			herr[l]->SetLineColor(50);
			herr[l]->SetLineWidth(2);
			inputmode[l].erase(0,4);
			lfrac->AddEntry(herr[l], "Nucleon-A", "l");
			herr[l]->SetLineStyle(2);
			herr[l]->Draw("hist,same");
		}
		else if  (inputmode[l] == "PPFXThinNuc"){
			herr[l]->SetLineColor(41);
			herr[l]->SetLineWidth(2);
			inputmode[l].erase(0,4);
			lfrac->AddEntry(herr[l], "pC #rightarrow NucleonX", "l");
			herr[l]->SetLineStyle(2);
			herr[l]->Draw("hist,same");
		}
		else if  (inputmode[l] == "PPFXThinPion"){
			herr[l]->SetLineColor(8);
			herr[l]->SetLineWidth(2);
			inputmode[l].erase(0,4);
			lfrac->AddEntry(herr[l], "pC #rightarrow #piX", "l");
			herr[l]->SetLineStyle(2);
			herr[l]->Draw("hist,same");
		}
		else if  (inputmode[l] == "PPFXTotAbsorp"){
			herr[l]->SetLineColor(kMagenta-7);
			herr[l]->SetLineWidth(2);
			inputmode[l].erase(0,4);
			lfrac->AddEntry(herr[l], inputmode[l].c_str(), "l");
			herr[l]->SetLineStyle(2);
			herr[l]->Draw("hist,same");
		}
		// else if  (inputmode[l] == "Total"){ // haven't got this to work properly/ dont know what it means
		// 	herr[l]->SetLineColor(kGray);
		// 	herr[l]->SetLineWidth(2);
		// 	lfrac->AddEntry(herr[l], inputmode[l].c_str(), "l");
		// 	herr[l]->Draw("hist,same");
		// }
		
		// if (inputmode[l] == "ms_PPFX") herr[l]->Draw("hist,same");
		// if (inputmode[l] == "Master" || inputmode[l] == "Total"  || inputmode[l] == "ms_PPFX")herr[l]->Draw("hist,same");

		if (mode == "numu")  herr[l]->SetTitle("#nu_{#mu}");
		if (mode == "nue")   herr[l]->SetTitle("#nu_{e}");
		if (mode == "numubar")  herr[l]->SetTitle("#bar{#nu_{#mu}}");
		if (mode == "nuebar")   herr[l]->SetTitle("#bar{#nu_{e}}");
		
		herr[l]->GetYaxis()->SetRangeUser(0,0.35);
		// herr[l]->GetYaxis()->SetRangeUser(0,1.75);
		
	}

	f2->cd();
	herr2 = (TH1D*) f2->Get(Err_names);

	herr2->SetLineColor(kBlue);
	herr2->SetLineWidth(2);
	herr2->Draw("same");
	lfrac->AddEntry(herr2, "NOvA","l");

	// Draw the legend
	lfrac->Draw();

	// create plots folder if it does not exist
	gSystem->Exec("if [ ! -d \"plots\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir plots; fi"); 

	// Save the plots as pdfs in the plots folder
	
	// using mipp
	if (mipp == "mippon"){
		if (mode == "numu"){ 	
			c1->Print("plots/CV_Flux_Prediction_NuMu_MIPPOn.pdf");
			c2->Print("plots/Ratio_FLux_Prediction_NuMu_MIPPOn.pdf");
			c3->Print("plots/Correlation_Matrix_NuMu_MIPPOn.pdf");
			c4->Print("plots/Fractional_Uncertainties_NuMu_MIPPOn.pdf");
			std::cout << "\n"<< std::endl;
		}
		else if (mode == "nue") {
			c1->Print("plots/CV_Flux_Prediction_Nue_MIPPOn.pdf");
			c2->Print("plots/Ratio_FLux_Prediction_Nue_MIPPOn.pdf");
			c3->Print("plots/Correlation_Matrix_Nue_MIPPOn.pdf");
			c4->Print("plots/Fractional_Uncertainties_Nue_MIPPOn.pdf");
			std::cout << "\n"<< std::endl;

		}
		else if (mode == "numubar"){
			c1->Print("plots/CV_Flux_Prediction_Numubar_MIPPOn.pdf");
			c2->Print("plots/Ratio_FLux_Prediction_Numubar_MIPPOn.pdf");
			c3->Print("plots/Correlation_Matrix_Numubar_MIPPOn.pdf");
			c4->Print("plots/Fractional_Uncertainties_Numubar_MIPPOn.pdf");
			std::cout << "\n"<< std::endl;

		}
		else {
			c1->Print("plots/CV_Flux_Prediction_Nuebar_MIPPOn.pdf");
			c2->Print("plots/Ratio_FLux_Prediction_Nuebar_MIPPOn.pdf");
			c3->Print("plots/Correlation_Matrix_Nuebar_MIPPOn.pdf");
			c4->Print("plots/Fractional_Uncertainties_Nuebar_MIPPOn.pdf");
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
			std::cout << "\n"<< std::endl;
		}
		else if (mode == "nue"){
			c1->Print("plots/CV_Flux_Prediction_Nue_MIPPOff.pdf");
			c2->Print("plots/Ratio_FLux_Prediction_Nue_MIPPOff.pdf");
			c3->Print("plots/Correlation_Matrix_Nue_MIPPOff.pdf");
			c4->Print("plots/Fractional_Uncertainties_Nue_MIPPOff.pdf");
			std::cout << "\n"<< std::endl;

		}
		else if (mode == "numubar"){
			c1->Print("plots/CV_Flux_Prediction_Numubar_MIPPOff.pdf");
			c2->Print("plots/Ratio_FLux_Prediction_Numubar_MIPPOff.pdf");
			c3->Print("plots/Correlation_Matrix_Numubar_MIPPOff.pdf");
			c4->Print("plots/Fractional_Uncertainties_Numubar_MIPPOff.pdf");
			std::cout << "\n"<< std::endl;

		}
		else {
			c1->Print("plots/CV_Flux_Prediction_Nuebar_MIPPOff.pdf");
			c2->Print("plots/Ratio_FLux_Prediction_Nuebar_MIPPOff.pdf");
			c3->Print("plots/Correlation_Matrix_Nuebar_MIPPOff.pdf");
			c4->Print("plots/Fractional_Uncertainties_Nuebar_MIPPOff.pdf");
			std::cout << "\n"<< std::endl;

		}
	}

} // end of main




