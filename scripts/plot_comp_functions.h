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
	if (mode == "numubar")  hist->SetTitle("#bar{#nu_{#mu}}");
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
void BeamlineUncertainties(std::vector<TH1D*> &herr, TH1D* hCV_Flux, std::string beam_type){
	TFile* f3 = TFile::Open("/uboone/data/users/kmistry/work/PPFX/nova/no_ThinKaon_Indiv/numu_Energy_ratios_quad.root");
		// Declearations
		double sigma_beamline, sys, sigma;
		std::vector<double> vsqsum;
		std::vector<double> vquadsum;
		TH1D* h_beamline;

		if (beam_type == "file")  h_beamline = (TH1D*) f3->Get("h_error2");
		
		// Manually add errors using individual histos
		else {

			if (h_beamline == NULL) std::cout << "Error, could not get hdev" << std::endl;
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

	int Nuniv = vhIn.size();
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
void HPUncertainties_Leo(TFile* fIn, TH1D* &herror, std::string inputmode){

	std::vector<TH1D*> vhuniv;

	// Get the universe histograms
	TDirectory *udir = fIn->GetDirectory(Form("numu/%s/Active_TPC_Volume", inputmode.c_str()));
	TIter next(udir->GetListOfKeys());
	TKey* key;

	// Loop over the directory and grab each histo from each uni
	while((key= (TKey*)next())){
		TClass* cl = gROOT->GetClass(key->GetClassName());
		if(!cl->InheritsFrom("TH1D"))continue;
		vhuniv.push_back((TH1D*)key->ReadObj());
	}

	// Calculate the HP error
	TH1D* htemp = getBand(vhuniv); 
	TH1D* hcv  = (TH1D*)htemp->Clone(); // Contains the mean of all the universes
	herror  = getFractionalError(htemp);

}
// ------------------------------------------------------------------------------------------------------------