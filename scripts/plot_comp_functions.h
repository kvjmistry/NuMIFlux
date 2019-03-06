// ------------------------------------------------------------------------------------------------------------
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
// ------------------------------------------------------------------------------------------------------------
// function that makes a legend for multiple histograms and draws them to the canvas
void legDraw(TLegend * &legend, TH1D *&hist, TString prodmode, TString mipp, std::string inputmode, TString mode){
	
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
	// else if  (inputmode == "Total"){ // haven't got this to work properly/ dont know what it means
	// 	hist->SetLineColor(kGray);
	// 	hist->SetLineWidth(2);
	// 	legend->AddEntry(hist, "Total", "l");
	// 	hist->Draw("hist,same");
	// }
	
	// if (inputmode == "ms_PPFX") hist->Draw("hist,same");
	// if (inputmode == "Master" || inputmode == "Total"  || inputmode == "ms_PPFX")hist->Draw("hist,same");
	else return;
	
	if (mode == "numu")  hist->SetTitle("#nu_{#mu}");
	if (mode == "nue")   hist->SetTitle("#nu_{e}");
	if (mode == "numubar"){
		hist->GetYaxis()->SetRangeUser(0,0.45);
		hist->SetTitle("#bar{#nu_{#mu}}");
	}
	if (mode == "nuebar")   hist->SetTitle("#bar{#nu_{e}}");


}
// ------------------------------------------------------------------------------------------------------------
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
// ------------------------------------------------------------------------------------------------------------
// Enumbers for the input mode 
enum e_mode{ enumu, enue, enumubar, enuebar};
// ------------------------------------------------------------------------------------------------------------
// Function to retun enum from mode label
e_mode return_mode(TString mode){
		if (mode == "numu")    return enumu;
		if (mode == "nue")     return enue;
		if (mode == "numubar") return enumubar;
		if (mode == "nuebar")  return enuebar;
		else return enumu;

}
// ------------------------------------------------------------------------------------------------------------
bool GetDirectory(TFile* f, TDirectory* &d, TString string){
	d = (TDirectory*)f->Get(string);
	if (d == NULL) {
		std::cout << "\nfailed to get:\t" << string << "\tThis directory might not exist in the file\n" << std::endl;
		return false;
	}
	else {
		return true;
	}
}
// ------------------------------------------------------------------------------------------------------------
bool GetHist(TFile* f, TH1D* &h, TString string){
	h = (TH1D*) f->Get(string);
	if (h == NULL) {
		std::cout << "\nfailed to get:\t" << string << "\tThis histogram might not exist in the file\n" << std::endl;
		return false;
	}
	else {
		return true;
	}
}
// ------------------------------------------------------------------------------------------------------------
bool GetHist(TFile* f, TH2D* &h, TString string){
	h = (TH2D*) f->Get(string);
	if (h == NULL) {
		std::cout << "\nfailed to get:\t" << string << "\tThis histogram might not exist in the file\n" << std::endl;
		return false;
	}
	else {
		return true;
	}
}
// ------------------------------------------------------------------------------------------------------------
bool GetTree(TFile* f, TTree* &T, TString string){
	T = (TTree*) f->Get(string);
	if (T == NULL) {
		std::cout << "\nfailed to get:\t" << string << "\tThis tree might not exist in the file\n" << std::endl;
		return false;
	}
	else {
		return true;
	}
}
// ------------------------------------------------------------------------------------------------------------
bool GetFile(TFile* &f , TString string){
	f = TFile::Open(string);
	
	if (f == NULL) {
		std::cout << "failed to get:\t" << string << "\tThis file might not exist in the file" << std::endl;
		return false;
	}
	else {
		return true;
	}
}
// ------------------------------------------------------------------------------------------------------------
void CalcCovariance(std::string inputmode,TString Cov_names, TFile* f, TH2D* &cov, TH1D* horig, const int nbins  ){
	TH1D* hu; 					// Flux hist for each universe
	const int nuni = 100; // num universes

	// Loop over universes
	for (int k=0; k<nuni; k++) {
		char name[500];
		snprintf(name, 500, Cov_names ,inputmode.c_str(),inputmode.c_str(), k); 
		
		// Check if sucessfully got histo
		bool boolhist = GetHist(f, hu, name); if (boolhist == false) gSystem->Exit(0);
		
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

				if (k != nuni - 1) cov->SetBinContent(i, j, cov->GetBinContent(i, j) + c ); // Fill with variance 
				else cov->SetBinContent(i, j, (cov->GetBinContent(i, j) + c) / nuni);       // Fill with variance and divide by nuni
			}

		} // End cov calc for universe i

		hu->Reset();
	}
}
// ------------------------------------------------------------------------------------------------------------
// Function to calculate the correlation matrix
void CalcCorrelation(TH2D* &cor, TH2D* &cov, int nbins ){
	std::cout << "Calculating Correlation Matrix" << std::endl;
	double cor_bini;
	// loop over rows
	for (int i=1; i<nbins+1; i++) {
		double cii = cov->GetBinContent(i, i);

		// Loop over columns
		for (int j=1; j<nbins+1; j++) {
			double cjj = cov->GetBinContent(j, j);
			double n = sqrt(cii * cjj);

			// Catch Zeros, set to arbitary 1.0
			if (n == 0) cor_bini = 0;
			else cor_bini = cov->GetBinContent(i, j) / n;

			cor->SetBinContent(i, j, cor_bini );
		}
	}
}
// ------------------------------------------------------------------------------------------------------------
// Function to calculate the fractional covariance matrix
void CalcFracCovariance(TH1D* &hCV, TH2D* &frac_cov, int nbins ){
	std::cout << "Calculating Fractional Covariance Matrix" << std::endl;
	double setbin;

	// loop over rows
	for (int i=1; i<nbins+1; i++) {
		double cii = hCV->GetBinContent(i);

		// Loop over columns
		for (int j=1; j<nbins+1; j++) {
			double cjj = hCV->GetBinContent(j);
			double n = cii * cjj;

			// Catch Zeros, set to arbitary 1.0
			if (n == 0) setbin = 0;
			else setbin = frac_cov->GetBinContent(i, j) / n;

			frac_cov->SetBinContent(i, j, setbin );
		}
	}
}
// ------------------------------------------------------------------------------------------------------------
// Use this for covariance calculation in 4d
void CalcCovariance_4D(TH2D* &cov, TH1D* &hCV, TH1D* &hu, const int &nbinsX, const int &nbinsY, int k  ){

	int nuni{100};
	// Loop over rows
	for (int i=1; i<nbinsX+1; i++) {

		double cvi = hCV->GetBinContent(i);   // CV bin i
		double uvi = hu->GetBinContent(i);    // Univ bin i 

		// Loop over columns
		for (int j=1; j<nbinsY+1; j++) {
			
			double cvj = hCV->GetBinContent(j); // CV bin j
			double uvj = hu->GetBinContent(j);  // Univ bin j 
			double c = (uvi - cvi) * (uvj - cvj); 

			if (k != nuni - 1) cov->SetBinContent(i, j, cov->GetBinContent(i, j) + c ); // Fill with variance 
			else cov->SetBinContent(i, j, (cov->GetBinContent(i, j) + c) / nuni);       // Fill with variance and divide by nuni
		}

	} // End cov calc for universe i

}
// ------------------------------------------------------------------------------------------------------------
void BeamlineUncertainties(std::vector<TH1D*> &herr, TH1D* hCV_Flux, std::string beam_type, TString mode){
	
	// Convert TString to a std::string
	std::string mode_str;
	if (mode == "numu") mode_str = "numu";
	else if (mode == "numubar") mode_str = "numubar";
	else if (mode == "nue") mode_str = "nue";
	else mode_str = "numubar";

	TFile* f3;
	bool boolfile;
	// if (mode_str == "numu")   boolfile = GetFile(f3,Form("/uboone/data/users/kmistry/work/PPFX/nova/no_ThinKaon_Indiv/%s_Energy_ratios_quad.root", mode_str.c_str() )); 
	// else  boolfile = GetFile(f3,Form("/uboone/data/users/kmistry/work/PPFX/nova/no_ThinKaon_Indiv/%s_Energy_ratios.root", mode_str.c_str() ));
	boolfile = GetFile(f3,Form("/uboone/data/users/kmistry/work/PPFX/nova/no_ThinKaon_Indiv/beamline_uncertainties/%s_beam_errors.root", mode_str.c_str() )); 
	
	if ( boolfile == false) gSystem->Exit(0);
	
	
	// Declearations
	double sigma_beamline, sys, sigma;
	std::vector<double> vsqsum;
	std::vector<double> vquadsum;
	TH1D* h_beamline;

	if (beam_type == "file")  h_beamline = (TH1D*) f3->Get("herr_totalbeam");
	else {// Manually add errors using individual histos

		if (beam_type == "file" && h_beamline == NULL) std::cout << "Error, could not get hdev" << std::endl;
		else std::cout << "Adding in Beamline uncertainties" << std::endl;


		
		vsqsum.resize(hCV_Flux->GetNbinsX(),0);
		vquadsum.resize(hCV_Flux->GetNbinsX(),0);

		// Loop over each unverse
		for (int j=0; j < 23; j++){
			
			if (j == 7 ) j += 1; // skip run 15
			
			char name[500];
			if (j == 0 || j == 1) snprintf(name, 500, "hratio_run000%d" ,j+8);
			else snprintf(name, 500, "hratio_run00%d" , j+8);
			
			bool bool_hist = GetHist(f3, h_beamline,name); if (bool_hist == false) gSystem->Exit(0);

			// Loop over each bin in the universe
			for (int bin=1; bin < h_beamline->GetNbinsX()+1; bin++){
				vsqsum.at(bin - 1) +=  (h_beamline->GetBinContent(bin) - 1.0) * (h_beamline->GetBinContent(bin) - 1.0);
				vquadsum.at(bin - 1) += (h_beamline->GetBinContent(bin) - 1.0) * (h_beamline->GetBinContent(bin) - 1.0);
				
			}

		}

		// Calc std and quadrature
		for (int i = 0; i < vsqsum.size(); i++){
				vsqsum.at(i) = std::sqrt(vsqsum.at(i)/20);
				vquadsum.at(i) = std::sqrt(vquadsum.at(i)/2);
				// std::cout << vsqsum.at(i)*100 << "\t"<< vquadsum.at(i)*100 << std::endl; // Display errors
		}
	}

	// Loop over histograms
	for (int i=0; i < herr.size(); i++){
		
		// Loop over bins
		for (int bin=1; bin<herr[i]->GetNbinsX()+1; bin++){
			sys = herr[i]->GetBinContent(bin);
			
			// Choose beamline uncertainty to add in
			if (beam_type == "file") sigma_beamline = h_beamline->GetBinContent(bin); // errors from a file
			else if (beam_type == "stdev") sigma_beamline = vsqsum.at(bin - 1);		  // standard deviation
			else if (beam_type == "quad") sigma_beamline = vquadsum.at(bin - 1);	  // Quadrature
			else std::cout << "Unknown Option given" << std::endl;

			sigma = std::sqrt( sys*sys + sigma_beamline*sigma_beamline );
			herr[i]->SetBinContent(bin, sigma);
		}
	}
}
// ------------------------------------------------------------------------------------------------------------
// Returns a histogram with the mean in each universe, the error is the stddev
TH1D* getBand(std::vector<TH1D*> vhIn){

	int Nuniv = vhIn.size(); // Get num universes
	if(Nuniv==0) gSystem->Exit(1);  
  
	TH1D* hband = (TH1D*) vhIn[0]->Clone();
	int Nbins   = hband->GetXaxis()->GetNbins(); 

	// Loop over the bins
	for(int jj=1;jj<=Nbins;jj++){
	
		// Mean of each bin
		double mean = 0;

		for (int ii=0; ii < Nuniv; ii++){

			mean += vhIn[ii]->GetBinContent(jj);
		}

		mean /= double(Nuniv);
		
		// Error
		double err = 0.0;

		for(int ii=0; ii < Nuniv; ii++){
			
			err += pow(vhIn[ii]->GetBinContent(jj) - mean,2);
		}
		
		err /= double(Nuniv);
		err = sqrt(err);
		hband->SetBinContent(jj,mean);
		hband->SetBinError(jj,err);
	}
  
  return hband;
}
// ------------------------------------------------------------------------------------------------------------
// Calculates the fractional errors
TH1D* getFractionalError(TH1D* hIn){

	TH1D* hfe = (TH1D*)hIn->Clone();

	int Nbins   = hfe->GetXaxis()->GetNbins(); 
	for(int jj=1; jj <= Nbins; jj++){

		double cont = hIn->GetBinContent(jj); // means(CV)
		double err = hIn->GetBinError(jj);	  // Errors
	
		if ( cont > 0 ) hfe->SetBinContent(jj,err/cont);  
	}
	return hfe;
}
// ------------------------------------------------------------------------------------------------------------
// Leos method of calculating the HP uncertainties 
void HPUncertainties_Leo(TFile* fIn, TH1D* &herror, std::string inputmode, TString mode){

	// Convert TString to a std::string
	std::string mode_str;
	if (mode == "numu") mode_str = "numu";
	else if (mode == "numubar") mode_str = "numubar";
	else if (mode == "nue") mode_str = "nue";
	else mode_str = "numubar";

	std::vector<TH1D*> vhuniv;
	TH1D* hu;

	// Get the universe histograms
	TDirectory *udir = fIn->GetDirectory(Form("%s/%s/Active_TPC_Volume", mode_str.c_str(), inputmode.c_str()));
	TIter next(udir->GetListOfKeys());
	TKey* key;
	int i{0};

	// Loop over the directory and grab each histo from each uni
	while((key= (TKey*)next())){
		TClass* cl = gROOT->GetClass(key->GetClassName());
		if(!cl->InheritsFrom("TH1D"))continue; // if not a TH1D skip

		hu = (TH1D*)key->ReadObj(); // get the histogram
		i=vhuniv.size()-1;
		
		// Veto theta histograms
		std::string huname = hu->GetName();
		std::string thetaname = Form("Th_%s_PPFXMaster_Uni_%i_AV_TPC", mode_str.c_str(), i);
		if ( huname == thetaname ) continue;
		
		vhuniv.push_back(hu);
	}

	// Calculate the HP error
	TH1D* htemp = getBand(vhuniv); 
	TH1D* hcv  = (TH1D*)htemp->Clone(); // Contains the mean of all the universes
	herror  = getFractionalError(htemp);

}
// ------------------------------------------------------------------------------------------------------------
// Function to draw assymetric errorband with CV and mean
void DrawErrorBand(TFile* f, TString mode, TLegend* leg, std::string inputmode){

	// Convert TString to a std::string
	std::string mode_str;
	if (mode == "numu") mode_str = "numu";
	else if (mode == "numubar") mode_str = "numubar";
	else if (mode == "nue") mode_str = "nue";
	else mode_str = "numubar";
	TH1D* hcv;

	std::vector<TH1D*> vhuniv;
	TH1D* hu; 

	TDirectory *udir = f->GetDirectory(Form("%s/%s/Active_TPC_Volume", mode_str.c_str(), inputmode.c_str()));
	TIter next(udir->GetListOfKeys());
	TKey* key;
	int i{0};

	// Loop over the directory and grab each histo from each uni
	while((key= (TKey*)next())){
		TClass* cl = gROOT->GetClass(key->GetClassName());
		if(!cl->InheritsFrom("TH1D"))continue; // if not a TH1D skip

		hu = (TH1D*)key->ReadObj(); // get the histogram
		i=vhuniv.size()-1;

		// Veto theta histograms
		std::string huname = hu->GetName();
		std::string thetaname = Form("Th_%s_PPFXMaster_Uni_%i_AV_TPC",mode_str.c_str(), i);
		if ( huname == thetaname ) continue;
		
		vhuniv.push_back(hu);
	}

	// Calculate the HP error using leos method
	TH1D* htemp = getBand(vhuniv); 
	TH1D* htemp_clone = (TH1D*) htemp->Clone("htemp_clone");

	bool bool_hist = GetHist(f, hcv, Form("%s/%s_CV_AV_TPC", mode_str.c_str(), mode_str.c_str()) ); if (bool_hist == false) gSystem->Exit(0);

	htemp->SetFillColor(18);
	htemp->SetFillStyle(1001);
	htemp->SetLineColor(kRed);
	htemp->SetLineWidth(2);
	htemp_clone->SetLineColor(kRed);
	htemp_clone->SetLineWidth(2);

	hcv->SetLineWidth(2);
	hcv->SetMarkerSize(0);
	hcv->SetLineColor(kBlack);

	htemp->GetXaxis()->SetRangeUser(0,20);
	htemp->SetStats(0);
	htemp->SetMarkerSize(0);
	htemp->SetTitle(Form("Error band for %s",mode_str.c_str()));
	htemp->GetXaxis()->SetTitle("#nu energy (GeV)");

	htemp->GetYaxis()->SetTitle(Form("%s/m^{2}/1e7POT",mode_str.c_str()));

	htemp->Draw("e2");
	htemp_clone->Draw("same,hist");
	hcv->Draw("hist,same");

	leg->AddEntry(htemp,"Mean Flux","LF");
	leg->AddEntry(hcv,"CV Flux","l");
}
// ------------------------------------------------------------------------------------------------------------
// function to plot each neutrino flux on the same graph
void PlotFluxSame(TCanvas *c,TLegend *leg, TFile *f1, TString mode, double fPOT, TString Gethist_TPC ){
	TH1D *h_flux;
	c->cd();
	
	// Check if sucessfully got histo
	bool boolhist = GetHist(f1, h_flux, Gethist_TPC); if (boolhist == false) gSystem->Exit(0);
	h_flux->SetDirectory(0);
	
	// Normalise flux by bin width (gives a flux/E [GeV])
	for (int i=1;i<h_flux->GetNbinsX()+1;i++) {
		h_flux->SetBinContent(i, h_flux->GetBinContent(i)/h_flux->GetBinWidth(i));		
	}


	// 6e20 POT, 50/1000 for 50 MeV from 1GeV, 1e-4 for m2->cm2
	h_flux->Scale( (6.0e20)/ (fPOT*1.0e4) * (50./1000.) );  
	

	if (mode == "numu"){
		h_flux->SetLineColor(kRed+1);
		h_flux->SetLineWidth(2);
		h_flux->SetTitle(";Energy [GeV];#nu / 6 #times 10^{20} POT / 50 MeV / cm^{2}");
		leg->AddEntry(h_flux, "#nu_{#mu}","l");
		h_flux->Draw("hist,same");
	}

	else if (mode == "nue"){	
		h_flux->SetLineColor(kRed+1);
		h_flux->SetLineWidth(2);
		h_flux->SetLineStyle(2);
		h_flux->SetTitle(";Energy [GeV];#nu / 6 #times 10^{20} POT / 50 MeV / cm^{2}");
		leg->AddEntry(h_flux, "#nu_{e}","l");
		h_flux->Draw("hist,same");

	}
	
	else if (mode == "numubar"){
		h_flux->SetLineColor(kBlue+1);
		h_flux->SetLineWidth(2);
		h_flux->SetTitle(";Energy [GeV];#nu / 6 #times 10^{20} POT / 50 MeV / cm^{2}");
		leg->AddEntry(h_flux, "#bar{#nu_{#mu}}","l");
		h_flux->Draw("hist,same");
	}
	
	else if (mode == "nuebar"){	
		h_flux->SetLineColor(kBlue+1);
		h_flux->SetLineWidth(2);
		h_flux->SetLineStyle(2);
		h_flux->SetTitle(";Energy [GeV];#nu / 6 #times 10^{20} POT / 50 MeV / cm^{2}");
		leg->AddEntry(h_flux, "#bar{#nu_{e}} ","l");
		h_flux->Draw("hist,same");
	
	}
	else std::cout << "unknown neutrino type ...."<<std::endl;
	h_flux->GetXaxis()->SetRangeUser(0,6);
	
	gPad->SetLogy();
	// gPad->Update();
	

}
// ------------------------------------------------------------------------------------------------------------
// Function to calculate the mean 2d histogram from the multisim variations and retun the unwarapped version
void CalcMeanHist(TFile* fIn, TH1D* &hMean_unwrap, int nBinsEnu, int nBinsTh, TString mode){
	std::cout << "Calculating the Mean universe"<<std::endl;
	// Convert TString to a std::string
	std::string mode_str;
	if (mode == "numu") mode_str = "numu";
	else if (mode == "numubar") mode_str = "numubar";
	else if (mode == "nue") mode_str = "nue";
	else mode_str = "numubar";

	std::vector<TH1D*> vhuniv;
	TH2D* hu;
	TH1D* hu_unwrap;

	// Get the universe histograms
	TDirectory *udir = fIn->GetDirectory(Form("%s/%s/Active_TPC_Volume", mode_str.c_str(), "PPFXMaster"));
	TIter next(udir->GetListOfKeys());
	TKey* key;
	int i{0};
	double counter;

	// Loop over the directory and grab each histo from each uni
	while((key= (TKey*)next())){
		TClass* cl = gROOT->GetClass(key->GetClassName());
		if(!cl->InheritsFrom("TH2D"))continue; // if not a TH2D skip

		hu = (TH2D*)key->ReadObj(); // Get the histogram
		i=vhuniv.size()-1;
		
		// Veto theta histograms
		std::string huname = hu->GetName();
		std::string thetaname = Form("Th_%s_PPFXMaster_Uni_%i_AV_TPC", mode_str.c_str(), i);
		if ( huname == thetaname ) continue;

		// Normalise 2d hist by bin area / deg / GeV
		// Loop over rows
		for (int p=1; p<nBinsEnu+1; p++) {
			// Loop over columns
			for (int q=1; q<nBinsTh+1; q++) {
				hu->SetBinContent(p,q, hu->GetBinContent(p, q)/ ( hu->GetXaxis()->GetBinWidth(p) * hu->GetYaxis()->GetBinWidth(q) ));
			}
		} 
		hu->Scale((6.0e20)/ (2.5e8*1.0e4)); // scale to right POT and m2

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
		
		// Push back to vector
		vhuniv.push_back(hu_unwrap);
	}

	// Now get the mean universe
	hMean_unwrap = getBand(vhuniv); 

}
// ------------------------------------------------------------------------------------------------------------
// Function to caluclate the bias matrix from the mean and CV to correct for in the covariance matrix
void CalcBias(TH2D* &hBias, TH1D* hMean_unwrap, TH1D* hCV_unwrap, int nbins){
	std::cout << "Calculating the Bias Covariance Matrix"<<std::endl;
	// Loop over rows
	for (int i=1; i<nbins+1; i++) {
		double cvi = hCV_unwrap->GetBinContent(i);  // CV bin i 
		double mi  = hMean_unwrap->GetBinContent(i); // mean bin i
		
		// Loop over columns
		for (int j=1; j<nbins+1; j++) {
			double cvj = hCV_unwrap->GetBinContent(j);  // CV bin i 
			double mj  = hMean_unwrap->GetBinContent(j); // mean bin i
			double c   = (mi - cvi) * (mj - cvj); 
			
			hBias->SetBinContent(i, j, c);
		}

	}

}
// ------------------------------------------------------------------------------------------------------------
// Function that calculates the percentage differnece between the CV and mean
void CalcRatioMeanCV(TH1D* hCV, TH1D *hMean, TH1D* &hRatioCVMean){

	std::cout << "Calculating the Ratio of CV to Mean"<<std::endl;
	hRatioCVMean = (TH1D*) hCV->Clone("hratioCVMean");
	hRatioCVMean->Divide(hMean);
	hRatioCVMean->SetLineColor(kBlack);
	hRatioCVMean->SetMarkerStyle(7);
	// hRatioCVMean->SetMarkerSize(3);
	hRatioCVMean->SetTitle(";Bin i;Ratio of CV to Mean");
	// hRatioCVMean->GetYaxis()->SetRangeUser(0.85, 1.15);

}
// ------------------------------------------------------------------------------------------------------------