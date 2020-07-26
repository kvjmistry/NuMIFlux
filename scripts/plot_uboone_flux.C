/**
 * NuMI Flux at uboone Plotting
 *
 * Script will generate all the correlation and covariance matrices of interest

 You can execute it with the commmand:
 root -l 'plot_uboone_flux.C("<horn_mode>", "<nu_flav>")'
 where <horn_mode> = fhc or rhc
 and <nu_flav> = nue or nuebar or numu or numubar 

 * 
 * A. Mastbaum <mastbaum@uchicago.edu> 2018/11
 * Modified by K. Mistry 22/06/19
 */

#include "functions.h"

// ----------------------------------------------------------------------------
// Main
void plot_uboone_flux(const char* horn,  const char* mode, const char* detector) { // (input, numu/nue)
    gStyle->SetOptStat(0); // say no to stats box

    // Declare variables
    TH1D *hCV_Flux, *hUW_Flux;
    TFile* f1;
    bool overwrite_errors{false};
    bool unweighted{false};
    bool boolfile;
    const char* mode_title;

    // Choose which file to get based on the detector
    if (!strcmp(detector,"uboone")) {
        // Load in the file, either fhc or rhc
        if (!strcmp(horn,"fhc")) {
            boolfile  = GetFile(f1 ,"/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold_v46/FHC/output_uboone_run0.root");
            if (boolfile == false) gSystem->Exit(0);
        }
        else {
            boolfile  = GetFile(f1 ,"/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold/RHC/output_uboone_run0.root");
            if (boolfile == false) gSystem->Exit(0);
        }
    }
    else if (!strcmp(detector,"nova")){
        boolfile  = GetFile(f1 ,"/uboone/data/users/kmistry/work/PPFX/nova/output.root");
        if (boolfile == false) gSystem->Exit(0);

    }
    else if (!strcmp(detector,"minerva")){
        boolfile  = GetFile(f1 ,"/uboone/data/users/kmistry/work/PPFX/minerva/output.root");
        if (boolfile == false) gSystem->Exit(0);

    }
    else {
        std::cout << "unknown detector type specified: " << detector << std::endl;
        gSystem->Exit(0);
    }

    // Create title characters from input
    if (strncmp("numu", mode, 4) == 0)		mode_title = "#nu_{#mu}";
    if (strncmp("nue", mode, 3) == 0)		mode_title = "#nu_{e}";
    if (strncmp("numubar", mode, 7) == 0)	mode_title = "#bar{#nu_{#mu}}";
    if (strncmp("nuebar", mode, 6) == 0)	mode_title = "#bar{#nu_{e}}";

    // PPFX Weight Modes
    std::vector<std::string> inputmode = {
        "ppfx_mippk_PPFXMIPPKaon",
        "ppfx_mipppi_PPFXMIPPPion",
        "ppfx_other_PPFXOther",
        "ppfx_targatt_PPFXTargAtten",
        "ppfx_think_PPFXThinKaon",
        "ppfx_thinmes_PPFXThinMeson",
        "ppfx_thinnpi_PPFXThinNeutronPion",
        "thinna_PPFXThinNucA",
        "thinn_PPFXThinNuc",
        "ppfx_thinpi_PPFXThinPion",
        "ppfx_totabs_PPFXTotAbsorp",
        "ppfx_ms_UBPPFX"};

    // Overwrite the updated names. This has only been done for nova and minerva so far
    if (!strcmp(detector,"nova") || !strcmp(detector,"minerva") ){
        inputmode[7] = "thinna_PPFXThinNucA";
        inputmode[8] = "thinn_PPFXThinNuc";
    }

    // ------------------------------------------------------------------------------------------------------------
    //                                                   CV Flux
    // ------------------------------------------------------------------------------------------------------------
    TCanvas* c1 = new TCanvas();
    TLegend* lFlux = new TLegend(0.5, 0.65, 0.9, 0.9);

    bool boolhist = GetHist(f1, hCV_Flux, Form("%s/Detsmear/%s_CV_AV_TPC", mode, mode)); if (boolhist == false) gSystem->Exit(0);
    hCV_Flux->SetDirectory(0);
    Normalise(hCV_Flux); 							// Normalise flux by bin width (gives a flux/E [GeV])
    TH1D* horig = (TH1D*) hCV_Flux->Clone("horig"); // Clone for plotting so dont need to norm the ms histograms

    // Scale and pretty
    double fPOT = GetPOT(f1); // POT
    hCV_Flux->Scale( (6.0e20)/ (fPOT*1.0e4) );
    hCV_Flux->SetLineColor(kRed+1);
    hCV_Flux->SetLineWidth(2);
    hCV_Flux->SetTitle(Form("%s; E_{#nu} [GeV];#nu / 6 #times 10^{20} POT / GeV / cm^{2}", mode_title));
    hCV_Flux->Draw("hist,same");
    gPad->SetLogy();
    
    // Legend
    lFlux->SetNColumns(1);
    lFlux->SetBorderSize(0);
    lFlux->SetFillStyle(0);
    lFlux->SetTextFont(62); 
    lFlux->AddEntry(hCV_Flux, "PPFX Flux","l");
    lFlux->Draw();
    IncreaseLabelSize(hCV_Flux);

    // IncreaseLabelSize(hCV_Flux);
    Draw_Nu_Mode(c1, horn); // Draw FHC Mode/RHC Mode Text
    
    // ------------------------------------------------------------------------------------------------------------
    //                                     Draw weighted flux vs unweighted flux
    // ------------------------------------------------------------------------------------------------------------
    TCanvas* c_uw_v_w = new TCanvas();
    TLegend* l_uw_v_w = new TLegend(0.6, 0.65, 0.75, 0.9);

    c_uw_v_w->cd();
    
    // Check if sucessfully got histo
    boolhist = GetHist(f1, hUW_Flux, Form("%s/Detsmear/%s_UW_AV_TPC", mode, mode )); if (boolhist == false) gSystem->Exit(0);
    Normalise(hUW_Flux); // Normalise flux by bin width (gives a flux/E [GeV])

    // 6e20 POT, 1e-4 for m2->cm2
    hUW_Flux->Scale( (6.0e20)/ (fPOT*1.0e4));  
    hUW_Flux->SetLineColor(kBlue+1);
    hUW_Flux->SetLineWidth(2);
    hUW_Flux->SetTitle(";Energy [GeV];#nu / 6 #times 10^{20} POT / GeV / cm^{2}");
    hCV_Flux->GetXaxis()->SetRangeUser(0,5);
    hCV_Flux->Draw("hist");
    hUW_Flux->Draw("hist,same");
    gPad->SetLogy();
    
    l_uw_v_w->AddEntry(hUW_Flux, "Unweighted","l");
    l_uw_v_w->AddEntry(hCV_Flux, "PPFX","l");
    l_uw_v_w->SetNColumns(1);
    l_uw_v_w->SetBorderSize(0);
    l_uw_v_w->SetFillStyle(0);
    l_uw_v_w->SetTextFont(62); 
    l_uw_v_w->Draw();

    // IncreaseLabelSize(hUW_Flux);
    IncreaseLabelSize(hCV_Flux);
    Draw_Nu_Mode(c_uw_v_w, horn); // Draw FHC Mode/RHC Mode Text

    // ------------------------------------------------------------------------------------------------------------
    //                                     Correlations, Covariance & uncertainties
    // ------------------------------------------------------------------------------------------------------------
    f1->cd();

    // Varables and histograms
    const int nbins = hCV_Flux->GetNbinsX();
    double* edges = new double[nbins+1];
    TCanvas* c3 = new TCanvas();
    TCanvas* c4 = new TCanvas();
    std::vector<TH2D*> cov;		// Covariance
    std::vector<TH2D*> cor; 	// Correlation
    std::vector<TH1D*> herr ;	// Fractional Uncertenties 

    // Set bin widths to be the same as NOvA
    for (int i=1; i<nbins+1; i++) edges[i-1] = hCV_Flux->GetBinLowEdge(i);
    
    // Get bin Edges
    edges[nbins] = hCV_Flux->GetBinLowEdge(nbins-1) + 2 * (hCV_Flux->GetBinWidth(nbins-1));
    
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

    // Legened
    TLegend* lfrac = new TLegend(0.5, 0.65, 0.9, 0.9);
    lfrac->SetNColumns(3);
    lfrac->SetBorderSize(0);
    lfrac->SetFillStyle(0);
    lfrac->SetTextFont(62); 

    // Loop over all input modes, get cov matrix and then get fractional uncertainties
    for (unsigned int l = 0; l < inputmode.size(); l++){
        
        // +++++++++++++++++
        // Covariance matrix
        // +++++++++++++++++
        CalcCovariance(inputmode[l], mode, "%s/Multisims/%s_%s_Uni_%i_AV_TPC" , f1, cov[l], horig, nbins);

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
        if (inputmode.at(l).find("ppfx_ms_UBPPFX") != std::string::npos) {
            std::cout << "Drawing correlation plot" << std::endl;
            cor[l]->SetTitle("Correlation PPFX All Weight");
            cor[l]->Draw("colz");
            IncreaseLabelSize(cor[l]);
        }

        // Draw_Nu_Mode(c3, horn); // Draw FHC Mode/RHC Mode Text


        gStyle->SetPalette(55); // kRainbow

        // Plot fractional errors overlaid with official NOvA plot
        c4->cd();

        // Make the plot
        if (overwrite_errors == false) legDraw(lfrac, herr[l], inputmode[l], mode);
        herr[l]->GetYaxis()->SetRangeUser(0,0.5);
        herr[l]->SetTitle(Form("%s; Energy [GeV];Fractional Uncertainty", mode_title));
        IncreaseLabelSize(herr[l]);
        // herr[l]->GetYaxis()->SetRangeUser(0,1.75);
        
    } // End loop over input labels

    // Draw the legend
    if (overwrite_errors == false) lfrac->Draw();
    Draw_Nu_Mode(c4, horn); // Draw FHC Mode/RHC Mode Text

    // ------------------------------------------------------------------------------------------------------------
    //                          Override the errors to use Leos method (mean instead of CV)
    // ------------------------------------------------------------------------------------------------------------
    // Decide if the errors need overwriting
    if (overwrite_errors == true ){
        c4->cd();
        std::cout << "Overwriting the errors" << std::endl;
        for (unsigned int l = 0; l < inputmode.size(); l++){
            HPUncertainties_Leo(f1, herr[l], inputmode[l], mode);
            legDraw(lfrac, herr[l], inputmode[l], mode);
            herr[l]->SetLineColor(kBlack);
            herr[l]->SetLineWidth(2);
            herr[l]->SetTitle(Form("%s; Energy [GeV];Fractional Uncertainty", mode_title));
            herr[l]->Draw("hist");
            lfrac->Draw();
        }
        c4->Update();
    }


    // ------------------------------------------------------------------------------------------------------------
    //                             Make a plot with the uncertainties with errorbands
    // ------------------------------------------------------------------------------------------------------------
    TCanvas* cband = new TCanvas();
    cband->cd();
    TLegend* leg = new TLegend(0.60,0.70,0.90,0.90);
    DrawErrorBand(f1, mode, leg, "ppfx_ms_UBPPFX"); // Plot for masterweight
    Draw_Nu_Mode(cband, horn); // Draw FHC Mode/RHC Mode Text
    leg->Draw();

    // ------------------------------------------------------------------------------------------------------------
    //                                Make a 4D Covariance matrix for re-weighing
    // ------------------------------------------------------------------------------------------------------------
    // Call caclulate covzriance matrix function. For now write function here.
    TH1D *hCV_unwrap, *hu_unwrap; 	// Flux hist for each universe
    TH2D *cov4d, *hCV2d, *hu;

    // Get the CV 2D matrix
    boolhist = GetHist(f1, hCV2d, Form("%s/Detsmear/%s_CV_AV_TPC_2D", mode, mode));   if (boolhist == false) gSystem->Exit(0);
    const int nBinsEnu = hCV2d->GetXaxis()->GetNbins(); // Enu
    const int nBinsTh  = hCV2d->GetYaxis()->GetNbins(); // Theta
    int nuni{200};
    
    //------------------------------
    // Normalise 2d hist by bin area / deg / GeV
    Normalise(hCV2d);
    double POT_2d = GetPOT(f1,false);
    hCV2d->Scale((6.0e20)/ (POT_2d * 1.0e4)); // scale to POT
    //------------------------------
    // Unwrap the histogram to binindex
    hCV_unwrap = new TH1D("", "",nBinsEnu*nBinsTh, 0, nBinsEnu*nBinsTh );
    UnwrapHist( hCV2d, hCV_unwrap);

    // Draw the CV Unwrapped
    TCanvas* c_CV_unwrap= new TCanvas();
    c_CV_unwrap->cd();
    hCV_unwrap->SetTitle("CV Unwrapped; Bin index; Flux ");
    hCV_unwrap->SetLineWidth(2);
    hCV_unwrap->SetLineColor(kBlack);
    hCV_unwrap->Draw("hist");
    gPad->SetLogy();
    IncreaseLabelSize(hCV_unwrap);
    Draw_Nu_Mode(c_CV_unwrap, horn); // Draw FHC Mode/RHC Mode Text

    // Draw the CV in 2D
    TCanvas* c_CV2d= new TCanvas();
    c_CV2d->cd();
    hCV2d->SetTitle(Form("%s CV 2D; Energy [GeV]; Theta [deg] ", mode_title));
    gPad->SetLogz();
    // gPad->Update();
    hCV2d->Draw("colz");
    hCV2d->GetZaxis()->SetTitle("#nu / 6 #times 10^{20} POT / GeV / cm^{2}");
    hCV2d->GetZaxis()->SetTitleOffset(+1.3);
    IncreaseLabelSize(hCV2d);
    Draw_Nu_Mode(c_CV2d, horn); // Draw FHC Mode/RHC Mode Text

     

    //------------------------------
    // 4D Covariance matrix	
    cov4d  = new TH2D(Form("%s_cov_4d", "PPFX_All"), ";Bin i; Bin j", nBinsEnu*nBinsTh, 0, nBinsEnu*nBinsTh , nBinsEnu*nBinsTh, 0, nBinsEnu*nBinsTh); // Bin i  Bin j

    // Loop over universes
    for (int k=0; k < nuni; k++) {
        
        char name[500];
        snprintf(name, 500,"%s/Multisims/%s_%s_Uni_%i_AV_TPC_2D" ,mode, mode, "ppfx_ms_UBPPFX", k); // Get uni i
        
        // Check if sucessfully got histo
        boolhist = GetHist(f1, hu, name); if (boolhist == false) gSystem->Exit(0);
        //------------------------------
        // Normalise 2d hist by bin area / deg / GeV
        Normalise(hu);
        hu->Scale((6.0e20)/ (POT_2d*1.0e4)); // scale to right POT and m2
        //------------------------------
        // Unwrap the histogram to binindex
        hu_unwrap = new TH1D("", "",nBinsEnu*nBinsTh, 0, nBinsEnu*nBinsTh );
        UnwrapHist(hu, hu_unwrap);
        //------------------------------
        // Now calculate the covariance matrix
        CalcCovariance_4D(cov4d, hCV_unwrap, hu_unwrap, nBinsEnu*nBinsTh, nBinsEnu*nBinsTh, k  );

    } // End cov calc for universe i
    
    TCanvas* c_cov= new TCanvas();
    c_cov->cd();
    cov4d->SetTitle(Form("%s Covariance Matrix ; Bin i; Bin j", mode_title));
    gStyle->SetPalette(kViridis);
    cov4d->Draw("colz");
    cov4d->GetZaxis()->SetTitle("#nu^{2} / 3.6 #times 10^{41} POT / GeV^{2} / cm^{4}");
    cov4d->GetZaxis()->SetTitleOffset(+1.3);
    IncreaseLabelSize(cov4d);
    Draw_Nu_Mode(c_cov, horn); // Draw FHC Mode/RHC Mode Text
    gPad->SetLogz();
    // gPad->Update();

    //------------------------------
    // Calculate the correlation matrix
    TCanvas *c_corr4d = new TCanvas();
    c_corr4d->cd();
    TH2D *hcorr4d = (TH2D*) cov4d->Clone("hCorr4d");
    CalcCorrelation(hcorr4d, cov4d, nBinsEnu*nBinsTh );
    hcorr4d->SetTitle(Form("%s Correlation Matrix ; Bin i; Bin j", mode_title));
    hcorr4d->Draw("colz");
    hcorr4d->GetZaxis()->SetTitle(" ");
    IncreaseLabelSize(hcorr4d);
    Draw_Nu_Mode(c_corr4d, horn); // Draw FHC Mode/RHC Mode Text
    // gPad->SetLogz();
    // gPad->Update();
    //------------------------------
    // Calculate the fractional covariance matrix
    TCanvas *c_fraccov4d = new TCanvas();
    c_fraccov4d->cd();
    TH2D *hfraccov4d = (TH2D*) cov4d->Clone("hfraccov4d");
    CalcFracCovariance(hCV_unwrap, hfraccov4d, nBinsEnu*nBinsTh );
    hfraccov4d->SetTitle(Form("%s Fractional Covariance Matrix ; Bin i; Bin j", mode_title));
    hfraccov4d->GetZaxis()->SetTitle(" ");
    hfraccov4d->Draw("colz");
    IncreaseLabelSize(hfraccov4d);
    Draw_Nu_Mode(c_fraccov4d, horn); // Draw FHC Mode/RHC Mode Text

    //------------------------------
    // Create a histogram with the bin indexes
    TCanvas *c_binidx = new TCanvas();
    c_binidx->cd();
    TH2D *hbinidx = (TH2D*) hCV2d->Clone("hCVClone");
    int counter = 0;
    // Loop over rows
    for (int i=1; i<nBinsEnu+1; i++) {  // Enu
        // Loop over columns
        for (int j=1; j<nBinsTh+1; j++){ // Theta
            counter ++;
            hbinidx->SetBinContent( i, j, counter );
        }
    }
    // hbinidx->SetTitle("Bin Indexes; Energy [GeV]; Theta [deg]");
    hbinidx->SetTitle(Form("%s Bin Indexes; Energy [GeV]; Theta [deg]", mode_title));
    hbinidx->Draw("text00");
    IncreaseLabelSize(hbinidx);

    std::vector<TLine*> vLine = MakeTLineVector(mode);
    for (unsigned int i =0; i < vLine.size(); i++){
       
        vLine[i]->SetLineColor(kBlue+1); // Line specifiers
        vLine[i]->SetLineStyle(3);
        vLine[i]->Draw("SAME");
    }
    
    hbinidx->GetXaxis()->SetRangeUser(0,0.5);
    hbinidx->GetXaxis()->SetRangeUser(0,3.5);
    gPad->SetRightMargin(0.1);

    // ------------------------------------------------------------------------------------------------------------
    // Compare the percentage difference between the mean and the CV in each unwrapped histogram bin right now we output the pull
    // ------------------------------------------------------------------------------------------------------------
    TCanvas *cMeanCV = new TCanvas();
    TH1D* hRatioCVMean;
    TH1D *hMean_unwrap = (TH1D*) hCV_unwrap->Clone("hMean_unwrap");
    CalcMeanHist(f1, hMean_unwrap, nBinsEnu, nBinsTh,mode);

    // Fill a histogram with the Ratio values to better visualise
    TH1D* hRatioMeanCVHist = new TH1D("hRatioMeanCVHist", "Ratio CV to Mean; Ratio; Entries", 30, 0.6, 1.4);
    TH1D* hPull = new TH1D("hPull", "CV - Mean / (stdev/#sqrt{n}); Pull; Entries", 25,-10, 10); // histogram with the pull 

    CalcPull(hCV_unwrap, hMean_unwrap, hPull);
    hPull->Draw("hist");
    // gPad->SetLogy();
    IncreaseLabelSize(hPull);
    Draw_Nu_Mode(cMeanCV, horn); // Draw FHC Mode/RHC Mode Text

    // // Now call a function that calcuates the ratio between them
    // CalcRatioMeanCV(hCV_unwrap, hMean_unwrap, hRatioCVMean);

    // // Get rid of zeros
    // for (unsigned int i =1; i <  hRatioCVMean->GetNbinsX()+1; i++ ){
    // 	if (hRatioCVMean->GetBinContent(i) == 0) continue;
    // 	// std::cout << hRatioCVMean->GetBinContent(i) << std::endl;
    // 	hRatioMeanCVHist->Fill(hRatioCVMean->GetBinContent(i));
    // }

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
    //                   Make a plot of the fractional uncertainties from the 4d covariance matrix
    // ------------------------------------------------------------------------------------------------------------
    TCanvas *c_FracError4d = new TCanvas();
    TH1D *hFracError4d = (TH1D*) hCV_unwrap->Clone("hFracError4d");

    // Function to caluclate the fractional uncertainties
    CalcFractionalError(cov4d, hCV_unwrap, hFracError4d );
    c_FracError4d->cd();
    hFracError4d->Draw("his");
    IncreaseLabelSize(hFracError4d);
    Draw_Nu_Mode(c_FracError4d, horn); // Draw FHC Mode/RHC Mode Text

    
    // ++++++++++++++++++++++++++++++++++
    // Save the plots as pdfs in the plots folder
    // ++++++++++++++++++++++++++++++++++
    std::cout << "-------------------------------------"<< std::endl;
    std::cout << "Now Saving Histograms.....\n"<< std::endl;

    // create plots folder if it does not exist
    gSystem->Exec("if [ ! -d \"plots/Hadron_Production\" ]; then echo \"\nplots/Hadron_Production folder does not exist... creating\"; mkdir -p plots/Hadron_Production; fi");
    gSystem->Exec("if [ ! -d \"plots/CV_Flux\" ]; then echo \"\nplots/CV_Flux folder does not exist... creating\"; mkdir -p plots/CV_Flux; fi"); 
    
    c1->Print(Form("plots/CV_Flux/CV_Flux_Prediction_Uneven_bins_%s_%s.pdf",horn, mode)); // CV 
    
    // c3->Print(Form("plots/Correlation_Matrix_%s.pdf", mode));
    
    c4->Print(Form("plots/Hadron_Production/Fractional_Uncertainties_%s_%s.pdf",horn, mode));
    
    c_uw_v_w->Print(Form("plots/CV_Flux/CV_Flux_unweighted_vs_ppfx_uneven_bins_%s_%s.pdf",horn, mode));
    
    c_cov->Print(Form("plots/Hadron_Production/Covariance_Marix_%s_%s.pdf", horn, mode));
    
    c_corr4d->Print(Form("plots/Hadron_Production/Correlation_Marix_%s_%s.pdf",horn,  mode));
    
    c_FracError4d->Print(Form("plots/Hadron_Production/Frac_Covariance_Marix_%s_%s.pdf",horn, mode));
    
    c_binidx->Print(Form("plots/Hadron_Production/BinIndex_%s.pdf", mode));
    
    c_fraccov4d->Print(Form("plots/Hadron_Production/Frac_Covariance_%s_%s.pdf",horn, mode));
    
    c_CV2d->Print(Form("plots/CV_Flux/CV_Flux_Energy_Theta_%s_%s.pdf",horn, mode));
    
    cband->Print(Form("plots/CV_Flux/CV_vs_Mean_Flux_%s_%s.pdf",horn, mode));

    cMeanCV->Print(Form("plots/Hadron_Production/CV_Pull_%s_%s.pdf",horn, mode));
    
    
} // end of main
