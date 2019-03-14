// Function to plot the beamline uncertainties compared to the CV
#include "plot_comp_functions.h"
// ------------------------------------------------------------------------------------------------------------
// Function to unwrap the histogram
void UnwrapTH2D(TH2D* &hFlux2d,TH1D* &h_unwrap, double POT){

    // Get bins
    const int nBinsEnu = hFlux2d->GetXaxis()->GetNbins(); // Enu
	const int nBinsTh = hFlux2d->GetYaxis()->GetNbins(); // Theta

    // Normalise 2d hist by bin area / deg / GeV
	// Loop over rows
	for (int i=1; i<nBinsEnu+1; i++) {

		// Loop over columns
		for (int j=1; j<nBinsTh+1; j++) {
			hFlux2d->SetBinContent(i,j, hFlux2d->GetBinContent(i, j)/ ( hFlux2d->GetXaxis()->GetBinWidth(i) * hFlux2d->GetYaxis()->GetBinWidth(j) ));
		}
	} 

	hFlux2d->Scale((6.0e20)/ (POT*1.0e4)); // scale to right POT and m2
    
    h_unwrap = new TH1D("", "",nBinsEnu*nBinsTh, 0, nBinsEnu*nBinsTh );
    h_unwrap->SetDirectory(0);
	int counter{0};
	for (int i=1; i<nBinsEnu+1; i++) { // Loop over rows
		for (int j=1; j<nBinsTh+1; j++){// Loop over columns
			counter++;
			h_unwrap->SetBinContent(counter, hFlux2d->GetBinContent(i , j)  );
		}
	}

}
// ------------------------------------------------------------------------------------------------------------
// Legend and draw specifiers for beamine variations
void DrawSpecifiers(TH1D* hist, TLegend* legend, std::string param, const char* specifier){
    // ----------------------
	//    Draw Specifiers
	// ----------------------
	if (param == "CV"){
	hist->SetLineColor(kBlack);
	hist->SetLineWidth(2);
	hist->SetLineStyle(1);
	legend->AddEntry(hist, "CV", "l");
	hist->Draw(specifier);
	} 
	else if  (param == "HP"){
		hist->SetLineColor(kMagenta+2);
		hist->SetLineWidth(2);
		hist->SetLineStyle(1);
		legend->AddEntry(hist, "Hadron Prod.", "l");
		hist->Draw(specifier);
	}
	else if  (param == "Horn_p2kA" ){ 
		hist->SetLineColor(30);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Horn +2kA", "l");
		hist->SetLineStyle(1);
		hist->Draw(specifier);
	}
	else if  (param == "Horn_m2kA"){ 
		hist->SetLineColor(30);
		hist->SetLineWidth(2);
		legend->AddEntry(hist,"Horn -2kA", "l");
		hist->SetLineStyle(2);
		hist->Draw(specifier);
	}
	else if  (param == "Horn1_x_p3mm"){
		hist->SetLineColor(28);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Horn1 x +3mm", "l");
		hist->SetLineStyle(1);
		hist->Draw(specifier);
		
	}
	else if  (param == "Horm1_x_m3mm" ){
		hist->SetLineColor(28);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Horm1 x -3mm", "l");
		hist->SetLineStyle(2);
		hist->Draw(specifier);
	}
	else if  (param == "Horn1_y_p3mm"){
		hist->SetLineColor(1001);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Horn1 y +3mm", "l");
		hist->SetLineStyle(1);
		hist->Draw(specifier);
	}
    else if  (param == "Horn1_y_m3mm"){
		hist->SetLineColor(1001);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Horn1 y -3mm", "l");
		hist->SetLineStyle(2);
		hist->Draw(specifier);
	}
	else if  (param == "Beam_spot_1_1mm"){
		hist->SetLineColor(kBlue+1);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Beam spot 1.1mm", "l");
		hist->SetLineStyle(1);
		hist->Draw(specifier);
	}
	else if  (param == "Nominal"){
		hist->SetLineColor(kBlue+1);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Nominal", "l");
		hist->SetLineStyle(3);
		hist->Draw(specifier);
	}
	else if  (param == "Beam_spot_1_5mm"){
		hist->SetLineColor(kBlue+1);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Beam spot 1.5mm", "l");
		hist->SetLineStyle(2);
		hist->Draw(specifier);
	}
	else if  (param == "Horn2_x_p3mm"){
		hist->SetLineColor(kOrange+1);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Horn2 x +3mm", "l");
		hist->SetLineStyle(1);
		hist->Draw(specifier);
	}
	else if  (param == "Horm2_x_m3mm"){
		hist->SetLineColor(kOrange+1);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Horm2 x -3mm", "l");
		hist->SetLineStyle(2);
		hist->Draw(specifier);
	}
	else if  (param == "Horn2_y_p3mm"){
		hist->SetLineColor(kMagenta);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Horn2 y +3mm", "l");
		hist->SetLineStyle(1);
		hist->Draw(specifier);
	}
	else if  (param == "Horn2_y_m3mm"){
		hist->SetLineColor(kMagenta);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Horn2 y -3mm", "l");
		hist->SetLineStyle(2);
		hist->Draw(specifier);
	}
	else if  (param == "Horns_0mm_water"){
		hist->SetLineColor(kSpring-7);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Horns 0mm water", "l");
		hist->SetLineStyle(1);
		hist->Draw(specifier);
	}
	else if  (param == "Horns_2mm_water"){
		hist->SetLineColor(kSpring-7);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Horns 2mm water", "l");
		hist->SetLineStyle(2);
		hist->Draw(specifier);
	}
    else if  (param == "Old_Horn" ){ 
		hist->SetLineColor(30);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Old Horn", "l");
		hist->SetLineStyle(1);
		hist->Draw(specifier);
	}
    else if  (param == "Beam_shift_x_p1mm" ){
		hist->SetLineColor(36);
		hist->SetLineWidth(2);
		legend->AddEntry(hist,"Beam shift x +1mm", "l");
		hist->SetLineStyle(1);
		hist->Draw(specifier);
	}
    else if  (param == "Beam_shift_x_m1mm" ){
		hist->SetLineColor(36);
		hist->SetLineWidth(2);
		legend->AddEntry(hist,"Beam shift x -1mm", "l");
		hist->SetLineStyle(2);
		hist->Draw(specifier);
	}
    else if  (param == "Beam_shift_y_p1mm"){
		hist->SetLineColor(42);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Beam shift y +1mm", "l");
		hist->SetLineStyle(1);
		hist->Draw(specifier);
	}
    else if  (param == "Beam_shift_y_m1mm"){
		hist->SetLineColor(42);
		hist->SetLineWidth(2);
		legend->AddEntry(hist,"Beam shift y -1mm", "l");
		hist->SetLineStyle(2);
		hist->Draw(specifier);
	}
    else if  (param == "Target_z_p7mm"){
		hist->SetLineColor(kOrange+10);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Target z +7mm", "l");
		hist->SetLineStyle(1);
		hist->Draw(specifier);
	}
    else if  (param == "Target_z_m7mm"){
		hist->SetLineColor(kOrange+10);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Target z -7mm", "l");
		hist->SetLineStyle(2);
		hist->Draw(specifier);
	}
    else if  (param == "Decay_pipe_Bfield"){
		hist->SetLineColor(50);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Decay pipe B field", "l");
		hist->SetLineStyle(1);
		hist->Draw(specifier);
	}
    else if  (param == "Horn1_refined_descr"){
		hist->SetLineColor(kMagenta-10);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Horn1 refined descr.", "l");
		hist->SetLineStyle(1);
		hist->Draw(specifier);
	}
    else if  (param == "Beam_divergence_54urad"){
		hist->SetLineColor(kTeal+6);
		hist->SetLineWidth(2);
		legend->AddEntry(hist, "Beam divergence 54#murad", "l");
		hist->SetLineStyle(1);
		hist->Draw(specifier);
	}
	else return;
}
// ------------------------------------------------------------------------------------------------------------
void DivideHists(TH1D* hCV, TH1D* hUniv, TH1D* &htemp){

    for (unsigned int i = 1; i < hCV->GetNbinsX()+1;i++){

        if (hUniv->GetBinContent(i) == 0) htemp->SetBinContent(i, 0);
        else htemp->SetBinContent(i, hCV->GetBinContent(i) / hUniv->GetBinContent(i));
        htemp->SetBinError(i, 0.000000001);

    }
}

// ------------------------------------------------------------------------------------------------------------
void plot_uboone_beamline(){
    gStyle->SetOptStat(0); // say no to stats box


    // Declearation of variables
	std::vector<std::string> params = { // A vector with the variations
		"CV", "Horn_p2kA","Horn_m2kA","Horn1_x_p3mm","Horm1_x_m3mm",
		"Horn1_y_p3mm","Horn1_y_m3mm","Beam_spot_1_1mm","Nominal",
		"Beam_spot_1_5mm","Horn2_x_p3mm","Horm2_x_m3mm","Horn2_y_p3mm",
		"Horn2_y_m3mm","Horns_0mm_water","Horns_2mm_water","Old_Horn",
		"Beam_shift_x_p1mm","Beam_shift_x_m1mm","Beam_shift_y_p1mm",
		"Beam_shift_y_m1mm","Target_z_p7mm","Target_z_m7mm",
		"Decay_pipe_Bfield","Horn1_refined_descr",
		"Beam_divergence_54urad" };


    bool boolfile, boolhist, useBeamlineCV{true};
    TFile *f;
    TH2D *htemp;
    TH1D *h_unwrap, *hCV_clone, *htemp_unwrap;
    double POT;
    TCanvas *c_beamline = new TCanvas();
    TCanvas *c_beamline_ratio = new TCanvas();
    // c_beamline->SetWindowSize(1000, 1000);
    TLegend *lFlux = new TLegend(0.8, 0.30, 1.0, 0.91);
    lFlux->SetNColumns(1);
	lFlux->SetBorderSize(0);
	lFlux->SetFillStyle(0);
	lFlux->SetTextFont(62);

    TLegend *lFlux_ratio = new TLegend(0.8, 0.30, 1.0, 0.91);
    lFlux_ratio->SetNColumns(1);
	lFlux_ratio->SetBorderSize(0);
	lFlux_ratio->SetFillStyle(0);
	lFlux_ratio->SetTextFont(62); 

    // ------------------------------------------------------------------------------------------------------------
    // CV
    // boolfile  = GetFile(f , "/uboone/data/users/kmistry/work/PPFX/uboone/with_tilt_2Dhists/output.root"); if (boolfile == false) gSystem->Exit(0); 
	if (useBeamlineCV) {boolfile  = GetFile(f , "/uboone/data/users/kmistry/work/PPFX/uboone/beamline/run15/output.root"); if (boolfile == false) gSystem->Exit(0);} // Beamline CV
	else boolfile  = GetFile(f , "/uboone/data/users/kmistry/work/PPFX/uboone/bugfix_release_novafiles/2dhists/output.root"); if (boolfile == false) gSystem->Exit(0); // Nova files CV
	
    boolhist = GetHist(f, htemp, "nue/nue_CV_AV_TPC"); if (boolhist == false) gSystem->Exit(0);
    
    // Get POT and Unwrap histogram
    POT = GetPOT(f);
    UnwrapTH2D( htemp, h_unwrap, POT );

    hCV_clone = (TH1D*) h_unwrap->Clone("hCV_clone");
    hCV_clone->SetDirectory(0);

    // Now Draw
	if (useBeamlineCV) h_unwrap->SetTitle("Beamline CV;Bin index; #nu / 6 #times 10^{20} POT / 50 MeV / deg / cm^{2}");
	else h_unwrap->SetTitle("Threshold File CV;Bin index; #nu / 6 #times 10^{20} POT / 50 MeV / deg / cm^{2}");;
    c_beamline->cd();
    DrawSpecifiers(h_unwrap, lFlux, "CV","hist,same");
    f->Close();

    // ------------------------------------------------------------------------------------------------------------
    // Loop over the beamline
    for (int i = 8; i < 27; i++){
		if (useBeamlineCV && i == 15) continue; // used to see ratio of nova file CV to Beamline
        boolfile  = GetFile(f , Form("/uboone/data/users/kmistry/work/PPFX/uboone/beamline/run%d/output.root",i)); if (boolfile == false) continue;
        boolhist = GetHist(f, htemp, "nue/nue_CV_AV_TPC"); if (boolhist == false) gSystem->Exit(0);
        
        // Get POT and Unwrap histogram
        POT = GetPOT(f);
        UnwrapTH2D( htemp, h_unwrap, POT );

        htemp_unwrap = (TH1D*) hCV_clone->Clone("htemp_unwrap");
		

        DivideHists(hCV_clone, h_unwrap, htemp_unwrap); // CV_clone, Beamline uni i, out

        // Now Draw
        c_beamline->cd();
		gPad->SetRightMargin(0.2);
        DrawSpecifiers(h_unwrap, lFlux, params.at(i - 7),"hist,same");
    
        c_beamline_ratio->cd();
        htemp_unwrap->SetMarkerStyle(7);
        htemp_unwrap->GetYaxis()->SetRangeUser(0.0, 1.5);
		gPad->SetRightMargin(0.2);
		if (useBeamlineCV) htemp_unwrap->SetTitle("Beamline CV;Bin index; Ratio Beamline/CV");
		else htemp_unwrap->SetTitle("Threshold File CV;Bin index; Ratio Beamline/CV");

		for (unsigned int i=1; i<htemp_unwrap->GetNbinsX()+1; i++ ) {if (htemp_unwrap->GetBinContent(i) == 0) htemp_unwrap->SetBinContent(i, 1);} // Get zero bins to 1 
        DrawSpecifiers(htemp_unwrap, lFlux_ratio, params.at(i - 7),"E1,same");
        
    
    }
    
    // ------------------------------------------------------------------------------------------------------------
    c_beamline->cd();
    lFlux->Draw();
    gPad->SetLogy();
	gPad->Update();
    
    c_beamline_ratio->cd();
    lFlux_ratio->Draw();
    TLine* flat = new TLine(0, 1, 42, 1);
	flat->SetLineStyle(7);
	flat->Draw();

	// ++++++++++++++++++++++++++++++++++
	// Save the plots as pdfs in the plots folder
	// ++++++++++++++++++++++++++++++++++
	// create plots folder if it does not exist
	gSystem->Exec("if [ ! -d \"plots\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir plots; fi"); 
	if (useBeamlineCV) c_beamline->Print("plots/nue_Beamline_2D_unwrapped_Flux_BeamlineCV.pdf");
	else c_beamline->Print("plots/nue_Beamline_2D_unwrapped_Flux_HPCV.pdf");
	if (useBeamlineCV) c_beamline_ratio->Print("plots/nue_Beamline_2D_unwrapped_Flux_ratio_BeamlineCV.pdf");
	else c_beamline_ratio->Print("plots/nue_Beamline_2D_unwrapped_Flux_ratio_HPCV.pdf");

} // End