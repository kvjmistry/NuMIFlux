// Script with finctions needed for fluxsystematics script to run
#include "plot_comp_functions.h"
// ------------------------------------------------------------------------------------------------------------
class event {
	public:
		std::string type; // Gen, Sig, Bkg, Sel, dirt
		double E; // Energy
		double Theta; // Theta
		int nu_flav; // nuetrino flavour type for picking the right histogram to reweight
};
// ------------------------------------------------------------------------------------------------------------
// Class to store the ratios of the CV to new universe for each nu flavour
class HistWeights {
	public:
		std::vector<TH2D*> HP; // Vector of histograms for each universe/cv
		std::vector<TH2D*> Beamline; // As above but for the beamline variations
}; 
		
// ------------------------------------------------------------------------------------------------------------
// Function to retun the right neutrino flavour ro weight by
std::string GetMode(int nu_flav){
	if  (nu_flav == 1 || nu_flav == 3) return "nue";
	else if  (nu_flav == 5 || nu_flav == 7) return "nuebar";
	else if  (nu_flav == 2 || nu_flav == 4) return "numu";
	else if  (nu_flav == 6 || nu_flav == 8) return "numubar";
	else return " ";
}
// ------------------------------------------------------------------------------------------------------------
void DivideHists(TH2D* hCV, TH2D* hUniv, TH2D* &htemp){

	// Loop over rows
	for (unsigned int i=1; i<hCV->GetNbinsX()+1; i++) { 
		// Loop over columns
		for (unsigned int j=1; j<hCV->GetNbinsY()+1; j++){ 
			if (hCV->GetBinContent(i, j) == 0) htemp->SetBinContent(i, j, 0);
			else htemp->SetBinContent(i, j, hUniv->GetBinContent(i, j) / hCV->GetBinContent(i, j) );
			//std::cout << htemp->GetBinContent(i,j) << std::endl;
		}
	}
}
// ------------------------------------------------------------------------------------------------------------
// Function to integrate a histogram
double IntegrateHist(TH2D* h){
	double integral{0};
	// Loop over rows
	for (unsigned int i=1; i<h->GetNbinsX()+1; i++) { 
		// Loop over columns
		for (unsigned int j=1; j<h->GetNbinsY()+1; j++){ 
			integral+= h->GetBinContent(i, j);
		}
	}
	return integral;
}
// ------------------------------------------------------------------------------------------------------------
// Function to take the ratio of universe i with nominal and return a weight for an E, theta
double GetWeight(int universe, TH2D* &hCV2d , int index, event event, TFile* fCV){
	std::string mode = GetMode(event.nu_flav); // get the neutrino flavour to reweigh the event type
	double weight;
	TH2D *hBeamline2d , *hHP2d, *hRatio;
	TFile *fBeamline;

	// if index is 0 then we are running over HP
	// get histogram from fCV for Masterweight PPFX universe i and divide out  by the cv
	if (index == 1){
		// bool boolfile  = GetFile(fCV , "/uboone/data/users/kmistry/work/PPFX/uboone/with_tilt_2Dhists/output.root"); if (boolfile == false) gSystem->Exit(0);
		bool boolhist = GetHist(fCV, hHP2d, Form("%s/PPFXMaster/Active_TPC_Volume/%s_PPFXMaster_Uni_%i_AV_TPC",mode.c_str(), mode.c_str(), universe)); if (boolhist == false) gSystem->Exit(0);
		
		hRatio = (TH2D*) hHP2d->Clone("hRatio");
		DivideHists(hCV2d, hHP2d, hRatio);
		// hRatio->Divide(hCV2d); // Divide the histogram
		// fCV->Close();
		// delete fCV;
	}
	// else we have a beamline variations indexes from 1 to N
	else if (index > 1) {
		bool boolfile  = GetFile(fBeamline , Form("/uboone/data/users/kmistry/work/PPFX/uboone/beamline/run%i/output.root", index+7)); if (boolfile == false) gSystem->Exit(0);
		bool boolhist = GetHist(fBeamline, hBeamline2d, Form("%s/%s_CV_AV_TPC",mode.c_str(), mode.c_str())); if (boolhist == false) gSystem->Exit(0);
	
		hRatio = (TH2D*) hBeamline2d->Clone("hRatio");
		// hRatio->Divide(hCV2d); // Divide the histogram
		DivideHists(hCV2d, hBeamline2d, hRatio);
		// fBeamline->Close();
	}
	else return weight = 1; // CV 

	// Now pick the value at Enu and Theta
	double xbin = hRatio->GetXaxis()->FindBin(event.E);
	double ybin = hRatio->GetYaxis()->FindBin(event.Theta);
	weight = hRatio->GetBinContent(xbin, ybin);

	// std::cout << weight << std::endl;

	return weight;
}
// ------------------------------------------------------------------------------------------------------------
// Function to calculate the interated flux
double IntegrateFlux(int universe, TFile* fCV, int index, double POTScale){

	double flux{0};
	TH2D *hBeamline2d, *hBeamline2dnue, *hBeamline2dnuebar, *hHP2d, *hHP2dnue ,*hHP2dnuebar, *hCV2dnue, *hCV2dnuebar, *hCV2d;
	TFile *fBeamline;

	// if index is 0 then we are running over HP
	// get histogram from fCV for Masterweight PPFX universe i and divide out  by the cv
	if (index == 1){

		bool boolhist = GetHist(fCV, hHP2dnue, Form("nue/PPFXMaster/Active_TPC_Volume/nue_PPFXMaster_Uni_%i_AV_TPC", universe)); if (boolhist == false) gSystem->Exit(0);
		boolhist = GetHist(fCV, hHP2dnuebar, Form("nuebar/PPFXMaster/Active_TPC_Volume/nuebar_PPFXMaster_Uni_%i_AV_TPC", universe)); if (boolhist == false) gSystem->Exit(0);
		
		hHP2d = (TH2D*) hHP2dnue->Clone("hHP2d");
		hHP2d->Add(hHP2dnuebar); // Combine the fluxes

		flux  = hHP2d->Integral(0, hHP2d->GetNbinsX(), 0, hHP2d->GetNbinsY()); // Integrate over whole phase space
		flux*= (POTScale / (GetPOT(fCV)*1.0e4));
	}
	// else we have a beamline variations indexes from 1 to N
	else if (index > 1){

		bool boolfile  = GetFile(fBeamline , Form("/uboone/data/users/kmistry/work/PPFX/uboone/bugfix_release/beamline/output%i.root", index+7)); if (boolfile == false) gSystem->Exit(0);

		// Get the CV histogram in 2D
		bool boolhist = GetHist(fBeamline, hBeamline2d, "nue/nue_CV_AV_TPC"); if (boolhist == false) gSystem->Exit(0);
		boolhist = GetHist(fBeamline, hBeamline2dnuebar, "nuebar/nuebar_CV_AV_TPC"); if (boolhist == false) gSystem->Exit(0);
		
		hBeamline2d = (TH2D*) hBeamline2dnue->Clone("hBeamline2d");
		hBeamline2d->Add(hBeamline2dnuebar); // Combine the fluxes

		flux  = hBeamline2d->Integral(0, 20, 0, 180); // Integrate over whole phase space
		
		flux*= (POTScale / (GetPOT(fBeamline)*1.0e4));

		fBeamline->Close();
	}
	else { // CV where index = 0
		bool boolhist = GetHist(fCV, hCV2dnue, "nue/nue_CV_AV_TPC"); if (boolhist == false) gSystem->Exit(0);
		boolhist = GetHist(fCV, hCV2dnuebar, "nuebar/nuebar_CV_AV_TPC"); if (boolhist == false) gSystem->Exit(0);

		hCV2d = (TH2D*) hCV2dnue->Clone("hCV2d");
		hCV2d->Add(hCV2dnuebar); // Combine the fluxes

		// Normalise 2d hist by bin area / deg / GeV
		// Loop over rows
		for (int p=1; p<hCV2d->GetNbinsX()+1; p++) {
			// Loop over columns
			for (int q=1; q<hCV2d->GetNbinsY()+1; q++) {
				hCV2d->SetBinContent(p,q, hCV2d->GetBinContent(p, q)/ ( hCV2d->GetXaxis()->GetBinWidth(p) * hCV2d->GetYaxis()->GetBinWidth(q) ));
			}
		} 
		
		flux  = hCV2d->Integral(0, hCV2d->GetNbinsX(), 0, hCV2d->GetNbinsY()), "width"; // Integrate over whole phase space
		// flux = IntegrateHist(hCV2d);

		flux*= (POTScale / (GetPOT(fCV)*1.0e4));
	}

	return flux;

}
// ------------------------------------------------------------------------------------------------------------
// Function to calculate the MC cross section -- deprecated for now
double CalcMCXSec(double sel, double bkg, double sig, double gen, double flux, double targets){
	return (sel  - bkg) / ( (sig / gen) * flux * targets ); 
}
// ------------------------------------------------------------------------------------------------------------
// Function to calculate the data cross section
double CalcDataXSec(double sel, double bkg , double flux,
					double targets, double intime_cosmics_bkg, double intime_cosmic_scale_factor,
					double dirt, double dirt_scale_factor, double mc_scale_factor, double efficiency ){

	return (sel - (intime_cosmics_bkg * intime_cosmic_scale_factor) - (dirt * dirt_scale_factor) - (bkg * mc_scale_factor)) / (efficiency * targets * flux); 
}
// ------------------------------------------------------------------------------------------------------------
// Function to add the weights for each universe -- will want to overload this with another function to handle unisims
void AddWeights(std::vector<double> &N, int Universes, TH2D* hCV2d, TFile* fCV, int index, event event){
	double weight{0};

	// Initialise the size of the counter if it is the first event loop. 
	if (N.size() == 0 )  N.resize( Universes );  // Resize to number of Universes.
	// std::cout << "Size of N:\t" << N.size() << " " << event.type<< std::endl; 

	// Loop over each universe
	for (unsigned int i = 0; i < Universes; i++){ 

		// Get the weight
		weight = GetWeight(i, hCV2d ,index, event, fCV);

		N.at(i) += weight; // Add weight to vector of counters.

	} 

}
// ------------------------------------------------------------------------------------------------------------
// Function to read in the event lists
void ReadEvents(const char *filename, std::vector<int> &N_evtnum ){
	
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
// Function to convert detector co-ords to a beam coordinate theta
double GetTheta(double detx, double dety, double detz){

	// Variables
	TRotation RotDet2Beam;     // Rotations
	TVector3 TransDetCoord, TransBeamCoord, detxyz; // Translations
	std::vector<double> rotmatrix;  // Inputs

	detxyz = {detx, dety, detz}; // input detector coordinates to translate
	TVector3 TransDettoTarg = {-31387.58422, -3316.402543 , -60100.2414}; // Translation from detector to target in detector coords

	TVector3 BeamCoords = detxyz; //+TransDettoTarg;

	// std::cout << BeamCoords.X() << "  "<< BeamCoords.Y() << "  " <<BeamCoords.Z() << std::endl;

	// From detector to beam coords -- need to check
	rotmatrix = {
		0.921038538,	4.625400126e-05,	-0.3894714486,
   		0.0227135048,	0.9982916247,		0.05383241394,
   		0.3888085752,	-0.05842798945, 	0.9194640079};

	// Return the TRotation
	TVector3 newX, newY, newZ;
	newX = TVector3(rotmatrix[0], rotmatrix[1], rotmatrix[2]);
	newY = TVector3(rotmatrix[3], rotmatrix[4], rotmatrix[5]);
	newZ = TVector3(rotmatrix[6], rotmatrix[7], rotmatrix[8]);

	RotDet2Beam.RotateAxes(newX, newY, newZ); // Return the TRotation
	
	BeamCoords = RotDet2Beam * BeamCoords;

	detx = BeamCoords.X();
	dety = BeamCoords.Y();
	detz = BeamCoords.Z();

	double costheta = detz / std::sqrt(detx*detx + dety*dety + detz*detz);
	double theta = std::acos(costheta) * 180 / 3.14159265;

	// Override because the angle is coming out weird from the above calculation...
	TVector3 temp= {0,0,1};
	theta = detxyz.Angle(temp) * 180 / 3.1415926;

	return theta;
}
// ------------------------------------------------------------------------------------------------------------
// Updated function for flux lists
void ReadEventList(const char *filename, std::vector<int> &N_evtnum,
					std::vector<std::string> &class_type, std::vector<int> &mc_nu_id,
					std::vector<double> &mc_nu_dir_x, std::vector<double> &mc_nu_dir_y, std::vector<double> &mc_nu_dir_z, std::vector<double> &mc_nu_energy    ){
	
	// event number, classifier type, mc_nu_id, mc_nu_dir_x, mc_nu_dir_y, mc_nu_dir_z, mc_nu_energy


	std::ifstream fileIN; 

	fileIN.open(filename); // Open the file
	
	if (!fileIN.good()) {  // Check if the file opened correctly
		std::cerr << "Error: file:\t" << filename <<"\tcould not be opened" << std::endl;
		exit(1);
	}

	int temp_evtnum,  temp_mc_nu_id;
	double temp_mc_nu_dir_x, temp_mc_nu_dir_y, temp_mc_nu_dir_z, temp_mc_nu_energy;
	std::string temp_class_type;

	if (fileIN.is_open()) { 
		
		// loop over lines in file
		while ( fileIN >> temp_evtnum >> temp_class_type >> temp_mc_nu_id >> temp_mc_nu_dir_x >> temp_mc_nu_dir_y >> temp_mc_nu_dir_z >> temp_mc_nu_energy) {  
			
			N_evtnum.push_back(temp_evtnum);
			class_type.push_back(temp_class_type);
			mc_nu_id.push_back(temp_mc_nu_id);
			mc_nu_dir_x.push_back(temp_mc_nu_dir_x);
			mc_nu_dir_y.push_back(temp_mc_nu_dir_y);
			mc_nu_dir_z.push_back(temp_mc_nu_dir_z);
			mc_nu_energy.push_back(temp_mc_nu_energy);
		}
		
		fileIN.close();
	}
}
// ------------------------------------------------------------------------------------------------------------
//  Function to calculate histogram ratios and return them as a vector of Th2D
void PrecalcHistRatio(HistWeights &flav, const char* mode){
	TH2D *hBeamline2d , *hHP2d, *hCV2d, *hCV2d_Beam;
	TFile *fBeamline, *fCV, *fCV_Beam;
	bool boolhist, boolfile;

	std::cout << "Calculating Hist ratios for flavour:\t" << mode << std::endl;

	// --- ---- ----- HP ---- --- --- ----- //
	// File with HP
	boolfile  = GetFile(fCV , "/uboone/data/users/kmistry/work/PPFX/uboone/with_tilt_2Dhists/output.root"); if (boolfile == false) gSystem->Exit(0);
	boolhist = GetHist(fCV, hCV2d, Form("%s/%s_CV_AV_TPC", mode, mode)); if (boolhist == false) gSystem->Exit(0); // Get the CV

	for (unsigned int i=0; i<100; i++){
		boolhist = GetHist(fCV, hHP2d, Form("%s/PPFXMaster/Active_TPC_Volume/%s_PPFXMaster_Uni_%i_AV_TPC",mode, mode, i)); if (boolhist == false) gSystem->Exit(0);
		hHP2d->Divide(hCV2d); // Divide hists
		// DivideHists(hCV2d, hHP2d, hRatio);
		flav.HP.push_back(hHP2d); // push back to vector
	}
	
	// --- ---- ----- Beamline ---- --- --- ----- //
	// Get the CV from the nova files
	boolfile  = GetFile(fCV_Beam , "/uboone/data/users/kmistry/work/PPFX/uboone/bugfix_release_novafiles/2dhists/output.root"); if (boolfile == false) gSystem->Exit(0); 
	boolhist = GetHist(fCV_Beam, hCV2d_Beam, Form("%s/%s_CV_AV_TPC",mode, mode)); if (boolhist == false) gSystem->Exit(0);
	
	for (unsigned int i=8; i<27; i++){
		boolfile  = GetFile(fBeamline , Form("/uboone/data/users/kmistry/work/PPFX/uboone/beamline/run%i/output.root", i)); if (boolfile == false) gSystem->Exit(0);
		boolhist = GetHist(fBeamline, hBeamline2d, Form("%s/%s_CV_AV_TPC",mode, mode)); if (boolhist == false) gSystem->Exit(0);
		hBeamline2d->Divide(hCV2d_Beam); // Divide hists
		// DivideHists(hCV2d, hHP2d, hRatio);
		flav.Beamline.push_back(hBeamline2d); // push back to vector
	}
	
}
// ------------------------------------------------------------------------------------------------------------