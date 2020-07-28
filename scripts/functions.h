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
// Function to increase to axes labels
void IncreaseLabelSize(TH1D* h){

    // h->GetXaxis()->SetRangeUser(0,3.5);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.12);
}
void IncreaseLabelSize(TH2D* h){

    // h->GetXaxis()->SetRangeUser(0,3.5);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetZaxis()->SetLabelSize(0.05);
    h->GetZaxis()->SetTitleSize(0.05);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.2);
    gPad->SetBottomMargin(0.13);
    h->SetMarkerSize(1.8);
    // gPad->SetGridx(); 
}
// ------------------------------------------------------------------------------------------------------------
// function that makes a legend for multiple histograms and draws them to the canvas
void legDraw(TLegend * &legend, TH1D *&hist, std::string inputmode, TString mode){
    

    if (inputmode == "ppfx_ms_UBPPFX"){
    hist->SetLineColor(kBlack);
    hist->SetLineWidth(2);
    legend->AddEntry(hist, "All", "l");
    hist->Draw("hist,same");
    } 
    else if  (inputmode == "ms_PPFX"){
        hist->SetLineColor(kMagenta+2);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "ms_ppfx", "l");
        hist->Draw("hist,same");
    }
    else if  (inputmode == "ppfx_other_PPFXOther"){
        hist->SetLineColor(28);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "Other", "l");
        hist->SetLineStyle(2);
        hist->Draw("hist,same");
        
    }
    else if  (inputmode == "ppfx_targatt_PPFXTargAtten" ){
        hist->SetLineColor(36);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "TargAtten", "l");
        hist->SetLineStyle(2);
        hist->Draw("hist,same");
    }
    else if  (inputmode == "ppfx_think_PPFXThinKaon"){
        hist->SetLineColor(1001);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "pC #rightarrow KX", "l");
        hist->SetLineStyle(2);
        hist->Draw("hist,same");
    }
    else if  (inputmode == "ppfx_thinmes_PPFXThinMeson"){
        hist->SetLineColor(kBlue+1);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "Meson Incident.", "l");
        hist->SetLineStyle(2);
        hist->Draw("hist,same");
    }
    else if  (inputmode == "ppfx_thinnpi_PPFXThinNeutronPion"){
        hist->SetLineColor(42);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "nC #rightarrow #piX", "l");
        hist->SetLineStyle(2);
        hist->Draw("hist,same");
    }
    else if  (inputmode == "thinna_PPFXThinNucA"){
        hist->SetLineColor(50);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "Nucleon-A", "l");
        hist->SetLineStyle(2);
        hist->Draw("hist,same");
    }
    else if  (inputmode == "thinn_PPFXThinNuc"){
        hist->SetLineColor(kOrange+10);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "pC #rightarrow NucleonX", "l");
        hist->SetLineStyle(2);
        hist->Draw("hist,same");
    }
    else if  (inputmode == "ppfx_thinpi_PPFXThinPion"){
        hist->SetLineColor(8);
        hist->SetLineWidth(2);
        legend->AddEntry(hist, "pC #rightarrow #piX", "l");
        hist->SetLineStyle(2);
        hist->Draw("hist,same");
    }
    else if  (inputmode == "ppfx_totabs_PPFXTotAbsorp"){
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
// Enumbers for the input mode 
enum e_mode{ enumu, enue, enumubar, enuebar};
// ------------------------------------------------------------------------------------------------------------
// Function to retun enum from mode label
e_mode return_mode(const char* mode){
        if (strncmp("numu", mode, 4) == 0)    return enumu;
        if (strncmp("nue", mode, 3) == 0)     return enue;
        if (strncmp("numubar", mode, 7) == 0) return enumubar;
        if (strncmp("nuebar", mode, 6) == 0)  return enuebar;
        else return enumu;

}
// Function to retun enum from mode label
e_mode return_mode(TString mode){
        if (mode == "numu")    return enumu;
        if (mode == "nue")      return enue;
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
double GetPOT(TFile* f){
    TTree* TPOT = (TTree*) f->Get("POT");
    if (TPOT == NULL) std::cout << "Error cant get POT info" << std::endl;

    double fPOT{0};
    TPOT->SetBranchAddress("POT", &fPOT); // Get the POT
    
    double total_entries = TPOT->GetEntries(); // if using hadd, this will not be 1 equal to 1 anymore

    double POT = 0;

    // Loop over the entries and add to the total pot
    for (int k = 0; k < total_entries; k++){
        TPOT->GetEntry(k);
        POT += fPOT;
    }
    
    std::cout << "TOTAL POT READ IN:\t" << POT << std::endl;
    
    return POT;
}
// ------------------------------------------------------------------------------------------------------------
// Overload to override the POT branch name
double GetPOT(TFile* f, TString POT_name, TString var_name){
    TTree* TPOT = (TTree*) f->Get(POT_name);
    if (TPOT == NULL) std::cout << "Error cant get POT info" << std::endl;

    double fPOT{0};
    TPOT->SetBranchAddress(var_name, &fPOT); // Get the POT
    
    double total_entries = TPOT->GetEntries(); // if using hadd, this will not be 1 equal to 1 anymore

    double POT = 0;

    // Loop over the entries and add to the total pot
    for (int k = 0; k < total_entries; k++){
        TPOT->GetEntry(k);
        POT += fPOT;
    }
    
    std::cout << "TOTAL POT READ IN:\t" << POT << std::endl;
    
    return POT;
}
// ------------------------------------------------------------------------------------------------------------
// Overloaded to supress the POT message
double GetPOT(TFile* f, bool disp){
    TTree* TPOT = (TTree*) f->Get("POT");
    if (TPOT == NULL) std::cout << "Error cant get POT info" << std::endl;

    double fPOT{0};
    TPOT->SetBranchAddress("POT", &fPOT); // Get the POT
    TPOT->GetEntry(0);
    double total_entries = TPOT->GetEntries(); // if using hadd, this will not be 1 equal to 1 anymore
    fPOT*=total_entries;
    if (disp) std::cout << "TOTAL POT READ IN:\t" << fPOT << std::endl;

    return fPOT;
}
// ------------------------------------------------------------------------------------------------------------
// Normalise a 1D histogram
void Normalise(TH1D* &h){
    
    // Normalise flux by bin width (gives a flux/E [GeV])
    for (int i=1;i<h->GetNbinsX()+1;i++) {
        h->SetBinContent(i, h->GetBinContent(i)/h->GetBinWidth(i));
        h->SetBinError(i, h->GetBinError(i)/h->GetBinWidth(i));
        // std::cout << h->GetBinWidth(i) << std::endl;
    }
}
// ------------------------------------------------------------------------------------------------------------
// Normalise a 2D histogram
void Normalise(TH2D* &h){
    
    for (int p=1; p<h->GetXaxis()->GetNbins()+1; p++) {
            // Loop over columns
            for (int q=1; q<h->GetYaxis()->GetNbins()+1; q++) {
                h->SetBinContent(p,q, h->GetBinContent(p, q)/ ( h->GetXaxis()->GetBinWidth(p) * h->GetYaxis()->GetBinWidth(q) ));
                h->SetBinError(p,q, h->GetBinError(p, q)/ ( h->GetXaxis()->GetBinWidth(p) * h->GetYaxis()->GetBinWidth(q) ));
                
            }
        } 
}
// ------------------------------------------------------------------------------------------------------------
// Unwrap histogram
void UnwrapHist(TH2D* h2d, TH1D* &h_unwrap){
    
    int counter{0};
    for (int i=1; i<h2d->GetXaxis()->GetNbins()+1; i++) { // Loop over rows
        for (int j=1; j<h2d->GetYaxis()->GetNbins()+1; j++){// Loop over columns
            counter++;
            h_unwrap->SetBinContent(counter, h2d->GetBinContent(i , j)  );
        }
    }
}
// ------------------------------------------------------------------------------------------------------------
void CalcCovariance(std::string inputmode, const char* mode, TString Cov_names, TFile* f, TH2D* &cov, TH1D* horig, const int nbins  ){
    TH1D* hu; 					// Flux hist for each universe
    const int nuni = 200; // num universes

    // Loop over universes
    for (int k=0; k<nuni; k++) {
        char name[500];
        snprintf(name, 500, Cov_names ,mode, mode, inputmode.c_str(), k); 
        
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

    int nuni{200};
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
void BeamlineUncertainties(std::vector<TH1D*> &herr, TH1D* hCV_Flux, std::string beam_type, const char* mode){
    
    TFile* f3;
    bool boolfile;
    boolfile = GetFile(f3,Form("/uboone/data/users/kmistry/work/PPFX/nova/no_ThinKaon_Indiv/beamline_uncertainties/%s_beam_errors.root", mode )); 
    
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
    for(int jj=1; jj < Nbins+1; jj++){
    
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
void HPUncertainties_Leo(TFile* fIn, TH1D* &herror, std::string inputmode, const char* mode){

    std::vector<TH1D*> vhuniv;
    TH1D* hu;

    // Get the universe histograms
    TDirectory *udir = fIn->GetDirectory(Form("%s/Multisims", mode));
    TIter next(udir->GetListOfKeys());
    TKey* key;
    int i{0};

    // Loop over the directory and grab each histo from each uni
    while((key= (TKey*)next())){
        TClass* cl = gROOT->GetClass(key->GetClassName());
        if(!cl->InheritsFrom("TH1D"))continue; // if not a TH1D skip

        hu = (TH1D*)key->ReadObj(); // get the histogram
        i=vhuniv.size()-1;
        
        std::string huname = hu->GetName();
        // Find the label name 
        if (huname.find(inputmode) != std::string::npos) {
            // Veto theta and 2D histograms
            std::string thetaname = Form("Th_%s_%s_Uni_%i_AV_TPC",mode, inputmode.c_str() , i);
            if ( huname == thetaname ) continue;
            if (huname.find("2D") != std::string::npos) continue;
        }
        else continue;
        
        vhuniv.push_back(hu);
    }

    // Calculate the HP error
    TH1D* htemp = getBand(vhuniv); 
    TH1D* hcv  = (TH1D*)htemp->Clone(); // Contains the mean of all the universes
    herror  = getFractionalError(htemp);

}
// ------------------------------------------------------------------------------------------------------------
// Function to draw assymetric errorband with CV and mean
void DrawErrorBand(TFile* f, const char* mode, TLegend* leg, std::string inputmode){

    const char* mode_title;
    // Create title characters from input
    if (strncmp("numu", mode, 4) == 0)		mode_title = "#nu_{#mu}";
    if (strncmp("nue", mode, 3) == 0)		mode_title = "#nu_{e}";
    if (strncmp("numubar", mode, 7) == 0)	mode_title = "#bar{#nu_{#mu}}";
    if (strncmp("nuebar", mode, 6) == 0)	mode_title = "#bar{#nu_{e}}";

    TH1D* hcv;

    std::vector<TH1D*> vhuniv;
    TH1D* hu; 

    TDirectory *udir = f->GetDirectory(Form("%s/Multisims", mode));
    TIter next(udir->GetListOfKeys());
    TKey* key;
    int i{0};

    // Loop over the directory and grab each histo from each uni
    while((key= (TKey*)next())){
        TClass* cl = gROOT->GetClass(key->GetClassName());
        if(!cl->InheritsFrom("TH1D"))continue; // if not a TH1D skip

        hu = (TH1D*)key->ReadObj(); // get the histogram
        i=vhuniv.size()-1;

        std::string huname = hu->GetName();
        // Find the label name 
        if (huname.find(inputmode) != std::string::npos) {
            // Veto theta and 2D histograms
            std::string thetaname = Form("Th_%s_%s_Uni_%i_AV_TPC",mode, inputmode.c_str() , i);
            if ( huname == thetaname ) continue;
            if (huname.find("2D") != std::string::npos) continue;
        }
        else continue;
        
        vhuniv.push_back(hu);
    }

    // Calculate the HP error using leos method
    TH1D* htemp = getBand(vhuniv); 
    TH1D* htemp_clone = (TH1D*) htemp->Clone("htemp_clone");

    bool bool_hist = GetHist(f, hcv, Form("%s/Detsmear/%s_CV_AV_TPC", mode, mode) ); if (bool_hist == false) gSystem->Exit(0);

    Normalise(hcv);
    Normalise(htemp);
    Normalise(htemp_clone);

    htemp->Scale(6.0e20 / (GetPOT(f, false) * 1.0e4));
    htemp_clone->Scale(6.0e20 / (GetPOT(f, false) * 1.0e4));
    hcv->Scale(6.0e20 / (GetPOT(f, false) * 1.0e4));

    htemp->SetFillColor(18);
    htemp->SetFillStyle(1001);
    htemp->SetLineColor(kRed);
    htemp->SetLineWidth(2);
    htemp_clone->SetLineColor(kRed);
    htemp_clone->SetLineWidth(2);

    hcv->SetLineWidth(2);
    hcv->SetMarkerSize(0);
    hcv->SetLineColor(kBlack);

    // htemp_clone->GetXaxis()->SetRangeUser(0, 6);
    htemp->SetStats(0);
    htemp->SetMarkerSize(0);
    // htemp->SetTitle(Form("Error band for %s",mode));
    htemp->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
    htemp->GetYaxis()->SetTitle(Form("%s / 6 #times 10^{20} POT / cm^{2} / GeV", mode_title));
    htemp->SetTitle(Form("%s", mode_title));

    htemp->Draw("e2");
    htemp_clone->Draw("same,hist");
    hcv->Draw("hist,same");

    IncreaseLabelSize(htemp);

    leg->AddEntry(htemp,"Mean Flux","LF");
    leg->AddEntry(hcv,"CV Flux","l");

    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(62); 

    gPad->SetLogy();
}
// ------------------------------------------------------------------------------------------------------------
// function to plot each neutrino flux on the same graph
void PlotFluxSame(TCanvas *c,TLegend *leg, TFile *f1, TString mode, double fPOT, TString Gethist_TPC, double tot_flux ){
    TH1D *h_flux;
    c->cd();
    
    // Check if sucessfully got histo
    bool boolhist = GetHist(f1, h_flux, Gethist_TPC); if (boolhist == false) gSystem->Exit(0);
    h_flux->SetDirectory(0);

    // Get the flux for the neutrino flavour
    double flux_flav = h_flux->Integral(0,  h_flux->GetNbinsX()+1);

    // Convert the flux percentage to a char
    char flux_pcent[15];
    snprintf(flux_pcent, 15,"%2.1f" ,100 * flux_flav / tot_flux);
    
    // Rebon
    h_flux->Rebin(2);
    
    // Normalise flux by bin width (gives a flux/E [GeV])
    //for (int i=1;i<h_flux->GetNbinsX()+1;i++) {
    //	h_flux->SetBinContent(i, h_flux->GetBinContent(i)/h_flux->GetBinWidth(i));		
    //}
    Normalise(h_flux);


    // 50/1000 for 50 MeV from 1GeV, 1e-4 for m2->cm2

    double POT_Scale = 1.0; // The POT to scale to // was 6.0e20

    h_flux->Scale( (POT_Scale)/ (fPOT*1.0e4) );  

    IncreaseLabelSize(h_flux);	

    if (mode == "numu"){
        h_flux->SetLineColor(kRed+1);
        h_flux->SetLineWidth(2);
        h_flux->SetTitle(";Neutrino Energy [GeV];#nu / POT / GeV / cm^{2}");
        leg->AddEntry(h_flux, Form("#nu_{#mu} (%s%%)", flux_pcent),"l");
        h_flux->Draw("hist,same");
    }

    else if (mode == "nue"){	
        h_flux->SetLineColor(kRed+1);
        h_flux->SetLineWidth(2);
        h_flux->SetLineStyle(2);
        h_flux->SetTitle(";Neutrino Energy [GeV];#nu / POT / GeV / cm^{2}");
        leg->AddEntry(h_flux, Form("#nu_{e} (%s%%)", flux_pcent),"l");
        h_flux->Draw("hist,same");

    }
    
    else if (mode == "numubar"){
        h_flux->SetLineColor(kBlue+1);
        h_flux->SetLineWidth(2);
        h_flux->SetTitle(";Neutrino Energy [GeV];#nu / POT / GeV / cm^{2}");
        leg->AddEntry(h_flux, Form("#bar{#nu_{#mu}} (%s%%)", flux_pcent),"l");
        h_flux->Draw("hist,same");
    }
    
    else if (mode == "nuebar"){	
        h_flux->SetLineColor(kBlue+1);
        h_flux->SetLineWidth(2);
        h_flux->SetLineStyle(2);
        h_flux->SetTitle(";Neutrino Energy [GeV];#nu / POT / GeV / cm^{2}");
        leg->AddEntry(h_flux, Form("#bar{#nu_{e}} (%s%%)", flux_pcent),"l");
        h_flux->Draw("hist,same");
    
    }
    else std::cout << "unknown neutrino type ...."<<std::endl;
    h_flux->GetXaxis()->SetRangeUser(0,5);
    
    gPad->SetLogy();
    // gPad->Update();
    

}
// ------------------------------------------------------------------------------------------------------------
// Function to get the total flux for each neutrino flavour
double GetTotalFlux(TFile *f1){

    bool boolhist;
    TH1D *h;

    double flux_cv{0.0}; 

    boolhist = GetHist(f1, h, "numu/Detsmear/numu_UW_AV_TPC_5MeV_bin"); if (boolhist == false) gSystem->Exit(0);
    flux_cv += h->Integral(0,  h->GetNbinsX()+1);

    boolhist = GetHist(f1, h, "numubar/Detsmear/numubar_UW_AV_TPC_5MeV_bin"); if (boolhist == false) gSystem->Exit(0);
    flux_cv += h->Integral(0,  h->GetNbinsX()+1);

    boolhist = GetHist(f1, h, "nue/Detsmear/nue_UW_AV_TPC_5MeV_bin"); if (boolhist == false) gSystem->Exit(0);
    flux_cv += h->Integral(0,  h->GetNbinsX()+1);

    boolhist = GetHist(f1, h, "nuebar/Detsmear/nuebar_UW_AV_TPC_5MeV_bin"); if (boolhist == false) gSystem->Exit(0);
    flux_cv += h->Integral(0,  h->GetNbinsX()+1);

    return flux_cv;

}
// ------------------------------------------------------------------------------------------------------------
// Function to calculate the mean 2d histogram from the multisim variations and retun the unwarapped version
void CalcMeanHist(TFile* fIn, TH1D* &hMean_unwrap, int nBinsEnu, int nBinsTh, const char* mode){
    std::cout << "Calculating the Mean universe"<<std::endl;

    std::vector<TH1D*> vhuniv;
    TH2D* hu;
    TH1D* hu_unwrap;

    // Get the universe histograms
    TDirectory *udir = fIn->GetDirectory(Form("%s/Multisims", mode));
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
        
        // Veto theta histograms and 1D
        std::string huname = hu->GetName();
        if (huname.find("ppfx_ms_UBPPFX") != std::string::npos) {
            // Veto theta and 2D histograms
            std::string thetaname = Form("Th_%s_ppfx_ms_UBPPFX_Uni_%i_AV_TPC", mode, i);
            if ( huname == thetaname ) continue;
            if (huname.find("1D") != std::string::npos) continue;
        }
        else continue;

        // Normalise 2d hist by bin area / deg / GeV
        Normalise(hu);
        hu->Scale((6.0e20)/ (GetPOT(fIn, false) *1.0e4)); // scale to right POT and m2

        // Unwrap the histogram to binindex
        hu_unwrap = new TH1D("", "",nBinsEnu*nBinsTh, 0, nBinsEnu*nBinsTh );
        UnwrapHist( hu, hu_unwrap);
        
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
    // hRatioCVMean->SetMarkerSize(3);.q
    hRatioCVMean->SetTitle(";Bin i;Ratio of CV to Mean");
    // hRatioCVMean->GetYaxis()->SetRangeUser(0.85, 1.15);

}
// ------------------------------------------------------------------------------------------------------------
// Function that calculates the Pull
void CalcPull(TH1D* hCV, TH1D *hMean, TH1D* &hPull){

    std::cout << "Calculating the Pull"<<std::endl;
    
    for (unsigned int i = 1; i< hCV->GetNbinsX() + 1; i++){
        hPull->Fill( (hCV->GetBinContent(i) - hMean->GetBinContent(i) )/ (hMean->GetBinError(i)/sqrt(200))  );
        }
    
    hPull->SetLineWidth(2);
    hPull->SetLineColor(kBlack);
}
// ------------------------------------------------------------------------------------------------------------
// Function to caluclate the fractional error 
void CalcFractionalError(TH2D* cov4d, TH1D* hCV, TH1D* &hFracError4d ){

    std::cout << "Calculating the 4D Fractional Covariance Matrix"<<std::endl;

    int nbins = cov4d->GetXaxis()->GetNbins(); // Get the number of bins

    // Loop over rows
    for (int i=1; i<nbins+1; i++) {
        double cii = cov4d->GetBinContent(i, i); // Get diag

        if (hCV->GetBinContent(i) == 0) hFracError4d->SetBinContent(i, 0); // catch dividing by zero
        else hFracError4d->SetBinContent(i, sqrt(cii) / hCV->GetBinContent(i) );

    }
    hFracError4d->SetTitle("Covariance Matrix Fractional Error; Bin index; Fractional Error");
    hFracError4d->SetLineColor(kBlack);
    hFracError4d->SetLineWidth(2);
}
// ------------------------------------------------------------------------------------------------------------
double IntegrateHist1D(TH1D* h){

    // find the bin which has the energy of 0.75*0.2065
    double xbin_th = h->GetXaxis()->FindBin( 0.75*0.2065);
    // double xbin_th = h->GetXaxis()->FindBin( 0.75*0.6065);
    // double ybin_th = h->GetXaxis()->FindBin( 3);

    // std::cout << "xbin_th:\t" << xbin_th << std::endl;

    double flux  = h->Integral(xbin_th,  h->GetNbinsX()+1);
    return flux;
}
// ------------------------------------------------------------------------------------------------------------
// Overload to give it an energy threshold
double IntegrateHist1D(TH1D* h, double E_th){

    // find the bin which has the energy of threshold
    double xbin_th = h->GetXaxis()->FindBin( E_th);
    
    double flux  = h->Integral(xbin_th, h->GetNbinsX());
    return flux;
}
// ------------------------------------------------------------------------------------------------------------
double IntegrateHist2D(TH2D* h){

    // find the bin which has the energy of 0.75*0.2065
    // double xbin_th = h->GetXaxis()->FindBin( 0.75*0.2065);
    // double ybin_th = h->GetXaxis()->FindBin( 920);

    double xbin_th = h->GetXaxis()->FindBin( 0.75*0.2065);
    // double ybin_th = h->GetXaxis()->FindBin( 0.3);

    // std::cout << "xbin_th:\t" << xbin_th << std::endl;

    double flux  = h->Integral(xbin_th, h->GetNbinsX()+1, 0, h->GetNbinsY()+1);
    return flux;
}
// ------------------------------------------------------------------------------------------------------------
// Function to return a vector of tlines for a 2d hist. Contaons hardcoded bins for this analysis
std::vector<TLine*> MakeTLineVector(TString mode){

    std::vector<TLine*> vLine;
    std::vector<double> Ebins, Thbins;

    if (mode == "nue" || mode =="nuebar"){ // nue, nuebar
        Ebins = { 0.00 ,0.06, 0.125,  0.25, 0.5, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.50, 3.00, 3.50};
    } 
    else{ // numu or numubar
        Ebins = {0.00, 0.025, 0.03, 0.235 ,0.24, 0.50, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 3.00, 4.00};
        // Ebins = {0.00, 0.025, 0.03, 0.235 ,0.24, 0.50}; // for zoomed version

    }

    Thbins = {  0, 20, 110,  160 };

    double Th0 = Thbins[0]; 
    double Th1 = Thbins[Thbins.size()-1];

    double E0 = Ebins[0]; 
    double E1 = Ebins[Ebins.size()-1]; 
    
    for (unsigned int i = 0; i < Thbins.size() - 1; i++){

        double y0 = Thbins[i];
        TLine *l2  = new TLine(E1, y0, E0 , y0 );
        vLine.push_back(l2);

        for (unsigned int j = 0; j < Ebins.size() - 1 ; j++){
            double x0 = Ebins[j];
            TLine *l  = new TLine(x0, Th0, x0, Th1 );
            vLine.push_back(l);
        }
    }

    return vLine;
}
// ------------------------------------------------------------------------------------------------------------
// Function to draw whether neutrino mode or anti neutrino mode on a plot
void Draw_Nu_Mode(TCanvas* c, const char* horn){
    c->cd();

    TPaveText *pt;

    if (!strcmp(horn, "fhc")){
        pt = new TPaveText(0.115, 0.89, 0.315, 0.96,"NDC");
        pt->AddText("FHC Mode");
        pt->SetTextColor(kRed+2);
    }
    else {
        pt = new TPaveText(0.115, 0.89, 0.315, 0.96,"NDC");
        pt->AddText("RHC Mode");
        pt->SetTextColor(kBlue+2);
    }
    
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetTextSize(0.04);
    pt->Draw();
}
// ------------------------------------------------------------------------------------------------------------