/**
 * NuMI Flux at uboone Plotting
 *
 * Plots each individual weighting mode instead of a single one
 * does not plot the correlation matrix for each individual mode to speed up the time
 * 
 * A. Mastbaum <mastbaum@uchicago.edu> 2018/11
 * Modified by K. Mistry 12/18
 */

// ----------------------------------------------------------------------------
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
// ----------------------------------------------------------------------------
// function that makes a legend for multiple histograms and draws them to the canvas
void legDraw(TLegend *legend, TH1D *hist, TString prodmode, TString mipp, std::string inputmode, TString mode){
	
	if (inputmode == "PPFXMaster"){
	hist->SetLineColor(kBlack);
	hist->SetLineWidth(2);
	legend->AddEntry(hist, "PPFXMaster", "l");
	hist->Draw("hist,same");
	} 
	else if  (inputmode == "ms_PPFX"){
		hist->SetLineColor(kMagenta+2);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, prodmode, "l"); // prodmode is overridden in the terminal input
		hist->Draw("hist,same");
	}
	else if  (inputmode == "PPFXMIPPKaon" && mipp =="mippon"){ // only turn on if there is MIPP
		hist->SetLineColor(30);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "MIPPKaon", "l");
		hist->SetLineStyle(2);
		hist->Draw("hist,same");
	}
	else if  ( inputmode == "PPFXMIPPPion" && mipp =="mippon"){ // only turn on if there is MIPP
		hist->SetLineColor(38);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "MIPPPion", "l");
		hist->SetLineStyle(2);
		hist->Draw("hist,same");
	}
	else if  (inputmode == "PPFXOther"){
		hist->SetLineColor(28);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Other", "l");
		hist->SetLineStyle(2);
		hist->Draw("hist,same");
		
	}
	else if  (inputmode == "PPFXTargAtten" ){
		hist->SetLineColor(36);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "TargAtten", "l");
		hist->SetLineStyle(2);
		hist->Draw("hist,same");
	}
	else if  (inputmode == "PPFXThinKaon"){
		hist->SetLineColor(1001);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "pC #rightarrow KX", "l");
		hist->SetLineStyle(2);
		hist->Draw("hist,same");
	}
	else if  (inputmode == "PPFXThinMeson"){
		hist->SetLineColor(kBlue+1);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Meson Incident.", "l");
		hist->SetLineStyle(2);
		hist->Draw("hist,same");
	}
	else if  (inputmode == "PPFXThinNeutron"){
		hist->SetLineColor(42);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "nC #rightarrow #piX", "l");
		hist->SetLineStyle(2);
		hist->Draw("hist,same");
	}
	else if  (inputmode == "PPFXThinNucA"){
		hist->SetLineColor(50);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Nucleon-A", "l");
		hist->SetLineStyle(2);
		hist->Draw("hist,same");
	}
	else if  (inputmode == "PPFXThinNuc"){
		hist->SetLineColor(kOrange+10);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "pC #rightarrow NucleonX", "l");
		hist->SetLineStyle(2);
		hist->Draw("hist,same");
	}
	else if  (inputmode == "PPFXThinPion"){
		hist->SetLineColor(8);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "pC #rightarrow #piX", "l");
		hist->SetLineStyle(2);
		hist->Draw("hist,same");
	}
	else if  (inputmode == "PPFXTotAbsorp"){
		hist->SetLineColor(kMagenta-7);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "TotAbsorp", "l");
		hist->SetLineStyle(2);
		hist->Draw("hist,same");
	}
	else return;
	
	if (mode == "numu")  hist->SetTitle("#nu_{#mu}");
	if (mode == "nue")   hist->SetTitle("#nu_{e}");
	if (mode == "numubar")  hist->SetTitle("#bar{#nu_{#mu}}");
	if (mode == "nuebar")   hist->SetTitle("#bar{#nu_{e}}");


}
// ----------------------------------------------------------------------------
// Function to make the weight distribution plots
void weight_plots(TString mode, std::vector<std::string> &inputmode, TFile* f1, TString prodmode, TString mipp, TCanvas* c, TLegend * leg  ){
	TH1D* whist; // weight hist

	c->cd(); //canvas

	std::string mode_str;

	// Convert TString to a std::string
	if (mode == "numu") mode_str = "numu";
	else if (mode == "numubar") mode_str = "numubar";
	else if (mode == "nue") mode_str = "nue";
	else mode_str = "numubar";

	int loop = inputmode.size();
	
	// loop over labels
	for (int l=0; l<loop; l++){	
		if (inputmode[l] == "PPFXOther") continue;	// Skip this due to low stats
		char name[500];

		snprintf(name, 500, "%s/%s/%s_%s_wght_%s" ,mode_str.c_str(), inputmode[l].c_str(), mode_str.c_str(), "MS",  inputmode[l].c_str()); // MS
	
		whist = (TH1D*) f1->Get(name);
		
		// Check if sucessfully got histo
		if (whist == NULL) {
			std::cout << "\nfailed to get:\t" << name << "\tThis histogram might not exist in the file\n" << std::endl;
			return;
		}
		whist->Scale( 1./ (whist->GetEntries()) ); // Norm by num entries

		// Draw and customise the plot
		legDraw(leg, whist, prodmode, mipp, inputmode[l], mode);
		
		whist->GetXaxis()->SetTitle("Weight");
		
	}

	// CV -- needs to be separate from loop
	char name[500];
	snprintf(name, 500, "%s/%s_CV_wght",mode_str.c_str(), mode_str.c_str());  // CV
	whist = (TH1D*) f1->Get(name);

	// Check if sucessfully got histo
	if (whist == NULL) {
		std::cout << "\nfailed to get:\t" << name << "\tThis histogram might not exist in the file\n" << std::endl;
		return;
	}

	whist->Scale( 1./ (whist->GetEntries()) );
	whist->SetLineColor(kCyan);
	whist->SetLineWidth(2);
	leg->AddEntry(whist, "CV", "l");
	whist->Draw("hist,same");
	
	gPad->SetLogy();
	leg->Draw();

}
// ----------------------------------------------------------------------------
// Enumbers for the input mode 
enum e_mode{ enumu, enue, enumubar, enuebar};
// ----------------------------------------------------------------------------
// Function to retun enum from mode label
e_mode return_mode(TString mode){
		if (mode == "numu")    return enumu;
		if (mode == "nue")     return enue;
		if (mode == "numubar") return enumubar;
		if (mode == "nuebar")  return enuebar;
		else return enumu;

}
// ----------------------------------------------------------------------------
// Main
void plot_uboone_flux( TString mipp, TString inputfile, TString prodmode, TString wplot, TString mode) { // (mippon/mippoff, input, Product/noThinKaon etc. numu/nue)
	gStyle->SetOptStat(0); // say no to stats box

	std::vector<std::string> inputmode = loopdir(inputfile, mode); // Grab the names of the input reweighters

	TString Getmode;
	TString Gethist_TPC;
	TString Getflux;
	TString Cov_names;
	TString g_simp_names;

	// Select neutrino type to run with 
	switch (return_mode(mode)){
		case enumu:
			std::cout << "\nUsing NuMu Mode!\n" << std::endl;
			Getmode = "numu"; 												// Folder name
			Gethist_TPC = "numu_CV_AV_TPC";									// AV in TPC flux prediction
			Getflux = "flux_numu";											// CV flux from NOvA
			Cov_names = "numu/%s/Active_TPC_Volume/numu_%s_Uni_%i_AV_TPC";  // Covariance matrix names
			g_simp_names = "numuFluxHisto";									// G simple files
			break;

		case enue:
			std::cout << "\nUsing Nue Mode!\n" << std::endl;
			Getmode = "nue";
			Gethist_TPC = "nue_CV_AV_TPC";
			Getflux = "flux_nue";
			Cov_names = "nue/%s/Active_TPC_Volume/nue_%s_Uni_%i_AV_TPC";
			g_simp_names = "nueFluxHisto";
			break;

		case enumubar:
			std::cout << "\nUsing NuMubar Mode!\n" << std::endl;
			Getmode = "numubar";
			Gethist_TPC = "numubar_CV_AV_TPC";
			Getflux = "flux_numubar";
			Cov_names = "numubar/%s/Active_TPC_Volume/numubar_%s_Uni_%i_AV_TPC";
			g_simp_names = "anumuFluxHisto";
			break;

		case enuebar:
			std::cout << "\nUsing Nuebar Mode!\n" << std::endl;
			Getmode = "nuebar";
			Gethist_TPC = "nuebar_CV_AV_TPC";
			Getflux = "flux_nuebar";
			Cov_names = "nuebar/%s/Active_TPC_Volume/nuebar_%s_Uni_%i_AV_TPC";
			g_simp_names = "anueFluxHisto";
			break;

	}

	// Root is dumb and so need to pre-decalre  some stuff here
	TDirectory* d;
	TH1D* hCV_Flux;
	TH1D* hNOvA_CV_Flux;

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Overlay of output plot with official NOvA FHC numu flux
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
	TFile* f1 = TFile::Open(inputfile);
	TCanvas* c1 = new TCanvas();

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Get the POT in the file
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	TTree* TPOT = (TTree*) f1->Get("POT");
	if (TPOT == NULL) std::cout << "Error cant get POT info" << std::endl;

	double fPOT{0};
	TPOT->SetBranchAddress("POT", &fPOT); // Get the POT
	TPOT->GetEntry(0);
	std::cout << "TOTAL POT READ IN:\t" << fPOT << std::endl;
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	d = (TDirectory*)f1->Get(Getmode);

	// Check if sucessfully got DIR
	if (d == NULL) {
		std::cout << "\nfailed to get:\t" << Getmode << "\tThis directory might not exist in the file\n" << std::endl;
		return;
	}

	d->cd();
	hCV_Flux = (TH1D*) (gDirectory->Get(Gethist_TPC)->Clone("fx"));

	// Check if sucessfully got histo
	if (hCV_Flux == NULL) {
		std::cout << "\nfailed to get:\t" << Gethist_TPC << "\tThis histogram might not exist in the file\n" << std::endl;
		return;
	}

	hCV_Flux->SetDirectory(0);

	TH1D* horig = (TH1D*) hCV_Flux->Clone("horig"); // Clone for plotting so dont need to norm the ms histograms
	
	// Normalise flux by bin width (gives a flux/E [GeV])
	for (int i=1;i<hCV_Flux->GetNbinsX()+1;i++) {
		// std::cout << i<<"bc:\t"<< hCV_Flux->GetBinContent(i) <<"\tbw:\t" <<hCV_Flux->GetBinWidth(i)<<"\tbw norm:\t" <<hCV_Flux->GetBinWidth(i)/hCV_Flux->GetBinWidth(i)<< std::endl;
		hCV_Flux->SetBinContent(i, hCV_Flux->GetBinContent(i)/hCV_Flux->GetBinWidth(i));
		
	}
	
	// Get the gsimple flux given by colton
	TH1D* h_g_simp; // gsimple flux
	TFile* f_gsimple = TFile::Open("/uboone/data/users/kmistry/work/PPFX/uboone/NuMIFlux_update_morebins.root");

	h_g_simp = (TH1D*) (f_gsimple->Get(g_simp_names));
	// Check if sucessfully got histo
	if (h_g_simp == NULL) {
		std::cout << "\nfailed to get:\t" << g_simp_names << "\tThis histogram might not exist in the file\n" << std::endl;
		return;
	}
	

	hCV_Flux->Sumw2();
	hCV_Flux->SetLineColor(kRed);
	hCV_Flux->SetLineWidth(2);
	// Norm
	// 20 is to get the bins in 50 MeV from 1GeV pi is fudge factor gsimple comparison, POT counting done wrong becuase of >1 file per job
	hCV_Flux->Scale( (6.0e20)/ (100000*950*1.0e4) * (50./1000.) );  
	// hCV_Flux->Scale( 3.14159* (6.0e20)/ (100000*950*1.0e4*20) );  // 671.36 is the window area, 20 is to get the bins in 50 MeV from 1GeV pi is fudge factor
	hCV_Flux->SetTitle(";E_{#nu} (GeV);#nu / 6 #times 10^{20} POT / 5 MeV / cm^{2}");
	gPad->SetLogy();
	gPad->Update();
	h_g_simp->SetLineWidth(2);
	
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
	std::vector<TH2D*> cov;	// Covariance
	std::vector<TH2D*> cor; // Correlation
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
		herr[l] = new TH1D(Form("%s_herr",inputmode[l].c_str()),";E_{#nu} (GeV);Fractional Uncertainty", nbins, edges);
	}

	// Loop over all input modes, get cov matrix and then get fractional uncertainties
	for (unsigned int l = 0; l < inputmode.size(); l++){
		
		// Covariance matrix

		// Loop over universes
		for (int k=0; k<nuni; k++) {
			char name[500];
			snprintf(name, 500, Cov_names ,inputmode[l].c_str(),inputmode[l].c_str(), k); 
  
			hu = (TH1D*) f1->Get(name);

			// Check if sucessfully got histo
			if (hu == NULL) {
				std::cout << "\nfailed to get:\t" << name << "\tThis histogram might not exist in the file\n" << std::endl;
				return;
			}
			
			// Normalise new universe by bin width
			for (int m=1; m<nbins+1; m++) {
				hu->SetBinContent(m, hu->GetBinContent(m) / hu->GetBinWidth(m)); 
			}

			// Loop over rows
			for (int i=1; i<nbins+1; i++) {

				double cvi = horig->GetBinContent(i); // CV bin i
				double uvi = hu->GetBinContent(i);    // Univ bin i 

				// Loop over columns
				for (int j=1; j<nbins+1; j++) {
					
					double cvj = horig->GetBinContent(j); // CV bin j
					double uvj = hu->GetBinContent(j);    // Univ bin j 

					double c = (uvi - cvi) * (uvj - cvj);

					if (k != nuni - 1) cov[l]->SetBinContent(i, j, cov[l]->GetBinContent(i, j) + c ); // Fill with variance 
					else cov[l]->SetBinContent(i, j, (cov[l]->GetBinContent(i, j) + c) / nuni); // Fill with variance and divide by nuni
				}
			} // End cov calc for universe i

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

		// Make the plot
		legDraw(lfrac, herr[l], prodmode, mipp, inputmode[l], mode);
		
		herr[l]->GetYaxis()->SetRangeUser(0,0.35);
		// herr[l]->GetYaxis()->SetRangeUser(0,1.75);
		
	}


	// Draw the legend
	lfrac->Draw();

	// create plots folder if it does not exist
	gSystem->Exec("if [ ! -d \"plots\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir plots; fi"); 

	
	// ++++++++++++++++++++++++++++++++++
	// Make the weight histogram
	// ++++++++++++++++++++++++++++++++++

	TCanvas* c5;
	if (wplot == "wplot"){
		c5 = new TCanvas();
	
		TLegend* lwght = new TLegend(0.5, 0.65, 0.9, 0.9);
		lwght->SetNColumns(3);
		lwght->SetBorderSize(0);
		lwght->SetFillStyle(0);
		lwght->SetTextFont(62);

		weight_plots(mode, inputmode, f1, prodmode, mipp, c5, lwght);
	}

	double sigma{0}, stat{0}, sys{0};
	// ++++++++++++++++++++++++++++++++++
	// Update the CV flux prediction to include stat+sys errors
	// ++++++++++++++++++++++++++++++++++
	c1->cd();
	

	std::cout << "bins:\t" << hCV_Flux->GetNbinsX() << std::endl;
	std::cout << "bins:\t" << std::endl;

	// Loop over the bins
	for (int bin=1; bin<hCV_Flux->GetNbinsX()+1; bin++){
		
		// Get the bin error (stat)
		stat = hCV_Flux->GetBinError(bin);

		// Get the bin error (sys)
		sys = herr[11]->GetBinContent(bin);
		sys = hCV_Flux->GetBinContent(bin) * sys; 

		// add in quadrature
		sigma =std::sqrt( stat*stat + sys*sys );

		// std::cout << "stat:\t" << stat << "\t" << "sys:\t" <<sys << "\tsigma:\t"<< sigma<< "\tbin content:\t"<< hCV_Flux->GetBinContent(bin) <<  std::endl;

		// Update the error on the plot
		hCV_Flux->SetBinError(bin, sigma);
	}



	// Redraw the plot with the new errors
	// c1->Update();
	// c4->Update();


	// ++++++++++++++++++++++++++++++++++
	// Save the plots as pdfs in the plots folder
	// ++++++++++++++++++++++++++++++++++

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




