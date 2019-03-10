// Script with finctions needed for fluxsystematics script to run
#include "plot_comp_functions.h"
// ------------------------------------------------------------------------------------------------------------
// Function to take the ratio of universe i with nominal and return a weight for an E, theta
double GetWeight(int universe, TH2D* hCV2d , int index, double Enu, double Theta)){
    
    double weight;
    TH2D *hBeamine2d, *hBeamine2dnue, *hBeamine2dnuebar, *hHP2d, *hHP2dnue ,*hHP2dnuebar;
    TFile *fBeamline, *fCV;

    // if index is 0 then we are running over HP
    // get histogram from fCV for Masterweight PPFX universe i and divide out  by the cv
    if (index == 1){
        bool boolfile  = GetFile(fCV , "/uboone/data/users/kmistry/work/PPFX/uboone/with_tilt_2Dhists/output.root"); if (boolfile == false) gSystem->Exit(0);
        bool boolhist = GetHist(fCV, hHP2dnue, Form("nue/PPFXMaster/Active_TPC_Volume/nue_PPFXMaster_Uni_%i_AV_TPC"), universe); if (boolhist == false) gSystem->Exit(0);
        boolhist = GetHist(fCV, hHP2dnuebar, Form("nuebar/PPFXMaster/Active_TPC_Volume/nuebar_PPFXMaster_Uni_%i_AV_TPC"), universe); if (boolhist == false) gSystem->Exit(0);
        
        hHP2d = (TH2D*) hHP2dnue->Clone("hHP2d");
	    hHP2d->Add(hHP2dnuebar); // Combine the fluxes

        TH2D* hRatio = (TH2D*) hHP2d->Clone("hRatio");
        hRatio->Divide(hCV2d); // Divide the histogram
        fCV->Close();
    }
    // else we have a beamline variations indexes from 1 to N
    else if (index > 1) {

        bool boolfile  = GetFile(fBeamline , Form("/uboone/data/users/kmistry/work/PPFX/uboone/with_tilt_2Dhists/beamline/output%i.root", index+7)); if (boolfile == false) gSystem->Exit(0);
        bool boolhist = GetHist(fBeamline, hBeamline2dnue, "nue/nue_CV_AV_TPC"); if (boolhist == false) gSystem->Exit(0);
        boolhist = GetHist(fCV, hBeamline2dnuebar, "nuebar/nuebar_CV_AV_TPC"); if (boolhist == false) gSystem->Exit(0);
        
        hBeamline2d = (TH2D*) hBeamline2dnue->Clone("hBeamline2d");
	    hBeamline2d->Add(hBeamline2dnuebar); // Combine the fluxes

        TH2D* hRatio = (TH2D*) hBeamline2d->Clone("hRatio");
        hRatio->Divide(hCV2d); // Divide the histogram
        fBeamline->Close();
    }
    else return weight = 1; // CV 

    // Now pick the value at Enu and Theta
    double xbin = hRatio->GetXaxis()->FindBin(Enu);
    double ybin = hRatio->GetXaxis()->FindBin(Theta);
    weight = hRatio->GetBinContent(xbin, ybin);

    return weight;
}
// ------------------------------------------------------------------------------------------------------------
// Function to calculate the interated flux
double IntegrateFlux(int universe, TFile* fCV, int index){

    double flux{0};
    TH2D *hBeamine2d, *hBeamine2dnue, *hBeamine2dnuebar, *hHP2d, *hHP2dnue ,*hHP2dnuebar, *hCV2dnue, *hCV2dnuebar, *hCV2d;

    // if index is 0 then we are running over HP
    // get histogram from fCV for Masterweight PPFX universe i and divide out  by the cv
    if (index == 1){

        bool boolhist = GetHist(fCV, hHP2dnue, Form("nue/PPFXMaster/Active_TPC_Volume/nue_PPFXMaster_Uni_%i_AV_TPC"), universe); if (boolhist == false) gSystem->Exit(0);
        boolhist = GetHist(fCV, hHP2dnuebar, Form("nuebar/PPFXMaster/Active_TPC_Volume/nuebar_PPFXMaster_Uni_%i_AV_TPC"), universe); if (boolhist == false) gSystem->Exit(0);
        boolhist = GetHist(fCV, hHP2dnuebar, Form("nuebar/PPFXMaster/Active_TPC_Volume/nuebar_PPFXMaster_Uni_%i_AV_TPC"), universe); if (boolhist == false) gSystem->Exit(0);
        
        hHP2d = (TH2D*) hHP2dnue->Clone("hHP2d");
	    hHP2d->Add(hHP2dnuebar); // Combine the fluxes

        flux  = hHP2d->Integral(0, 20, 0, 180); // Integrate over whole phase space
    }
    // else we have a beamline variations indexes from 1 to N
    else if (index > 1){

        bool boolfile  = GetFile(fBeamline , Form("/uboone/data/users/kmistry/work/PPFX/uboone/bugfix_release/beamline/output%i.root", index+7)); if (boolfile == false) gSystem->Exit(0);

        // Get the CV histogram in 2D
        bool boolhist = GetHist(fBeamline, hBeamline2d, "nue/nue_CV_AV_TPC"); if (boolhist == false) gSystem->Exit(0);
        boolhist = GetHist(fCV, hBeamline2dnuebar, "nuebar/nuebar_CV_AV_TPC"); if (boolhist == false) gSystem->Exit(0);
        
        hBeamline2d = (TH2D*) hBeamline2dnue->Clone("hBeamline2d");
	    hBeamline2d->Add(hBeamline2dnuebar); // Combine the fluxes

        flux  = hBeamline2d->Integral(0, 20, 0, 180); // Integrate over whole phase space

        fBeamline->Close();
    }
    else { // CV
        bool boolhist = GetHist(fCV, hCV2dnue, "nue/nue_CV_AV_TPC"); if (boolhist == false) gSystem->Exit(0);
        boolhist = GetHist(fCV, hCV2dnuebar, "nuebar/nuebar_CV_AV_TPC"); if (boolhist == false) gSystem->Exit(0);

        hCV2d = (TH2D*) hCV2dnue->Clone("hCV2d");
        hCV2d->Add(hCV2dnuebar); // Combine the fluxes

        flux  = hCV2d->Integral(0, 20, 0, 180); // Integrate over whole phase space

    }

    return flux;

}
// ------------------------------------------------------------------------------------------------------------
// Function to calculate the MC cross section
double CalcMCXSec(sel, bkg, sig, gen, flux, targets){
	return (sel  - bkg) / ( (sig / gen) * flux * targets ); 
}
// ------------------------------------------------------------------------------------------------------------
// Function to add the weights for each universe -- will want to overload this with another function to handle unisims
AddWeights(std::vector<double> &N, int Universes, TH2D* hCV2d, TFile* fCV, int index){
    
    // Call GetWeight()

	// Initialise the size of the counter if it is the first event loop. 
	if (N.size() == 0 )  N.resize( Universes );  // Resize to number of Universes. 

	// Loop over each universe
	for (unsigned int i = 0; i < Universes; i++){ 

		// Get the weight
		// -- need the energy and theta of the event

		N.at(i) += weight; // Add weight to vector of counters.

	} // loop over each universe

}
// ------------------------------------------------------------------------------------------------------------
// Function to read in the event lists
ReadEvents(const char *filename, std::vector<int> &N_evtnum ){
	
	std::ifstream fileIN; 

	fileIN.open(filename); // Open the file
	
	if (!fileIN.good()) {  // Check if the file opened correctly
		std::cerr << "Error: file:\t" << filename <<"\tcould not be opened" << std::endl;
		exit(1);
	}

	double temp{0}; // Use a temp var to get the values and push back

	if (fileIN.is_open()) { 
		
		while ( !fileIN.eof()) {   // loop over lines in file
			
			fileIN >> temp;        // Add number to temp var
			N_evtnum.push_back(temp);
		}
		
		fileIN.close();
	}
}
// ------------------------------------------------------------------------------------------------------------
// Function to convert detector co-ords to beam coords
// ------------------------------------------------------------------------------------------------------------