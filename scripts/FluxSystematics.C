/* 
This Macro will handle the calculation of the systematic error
om the cross section in a similar way to genie systematics.
Don't need to use the original Event_Weight_reader.cc module
since we dont need to get any dataproducts e.g eventweight. 

Simply read in a file list and loop over it on an event by event basis.
Match the tyoe of event, generate a weight from the histogram and recalculate the cross section
The intenstion is to treat the hadron production and unisim beamline variations separately
and then add them all together at the end. 
*/
#include "FluxSyst_functions.h"
// ------------------------------------------------------------------------------------------------------------
void FluxSystematics(){
	
	// Declearation of variables
	std::vector<std::string> params = { // A vector with the variations
		"CV","HP", "Horn_p2kA","Horn_m2kA","Horn1_x_p3mm","Horm1_x_m3mm",
		"Horn1_y_p3mm","Horn1_y_m3mm","Beam_spot_1_1mm","Nominal",
		"Beam_spot_1_5mm","Horn2_x_p3mm","Horm2_x_m3mm","Horn2_y_p3mm",
		"Horn2_y_m3mm","Horns_0mm_water","Horns_2mm_water","Old_Horn",
		"Beam_shift_x_p1mm","Beam_shift_x_m1mm","Beam_shift_y_p1mm",
		"Beam_shift_y_m1mm","Target_z_p7mm","Target_z_m7mm",
		"Decay_pipe_Bfield","Horn1_refined_descr",
		"Beam_divergence_54urad" };


	// Declare member data here.
	int run, subrun, evt;
	int Universes;        		// The number of universes simulated

	double total_in{0}; 		// Total number of matched events 
	double tot_gen{0};  		// counter for the total number of gen events read in
	double tot_sel{0};  		// counter for the total number of sel events read in
	double tot_sig{0};  		// counter for the total number of sig events read in
	double tot_bkg{0};  		// counter for the total number of bkg events read in
	double tot_filt{0};  		// counter for the total number of filtered events read in
	double tot_unmatched{0};  	// counter for the total number of filtered events read in

	// Vectors for cross section calculation. 
	std::vector<double> N_gen, N_sig, N_bkg, N_sel, N_dirt ,Data_x_sec, Efficiency, Flux;	// Vectors of num of events with new weights for generated, signal, background, selected
	double N_gen_CV{0}, N_sig_CV{0}, N_bkg_CV{0}, N_sel_CV{0}, N_dirt_CV{0};					// Input counters for num sel, gen, sig, bkg from the input lists + dirt as POT scaling diff
	std::vector<std::vector<double>> vec_Data_x_sec , vec_Efficiency, vec_gen, vec_sig, vec_dirt, vec_bkg , vec_flux;      				// vectors to store the differnt versions of the cross section

	// Flux
	const double flux_mc{4.19844e+10};  // POT Scaled divide these to get scale factor
	const double flux_gsimple_data{5383349994.0}; //-- gsimple flux --
	const double POTScale_MC{1.82027e+21}; // MC POT Scale
	const double POTScale{2.334e+20};

	// Num Targets
	const double targets_mc{3.50191e+31};
	const double targets_data{3.4723e+31};

	// DATA
	const double intime_cosmics_bkg{81};                // Number of intime cosmics for background
	const double num_selected_data{203};                // The number of selected events in data
	const double intime_cosmic_scale_factor{1.0154};    // Scale factor to apply to the intime cosimic background
	const double mc_scale_factor{0.1301};               // Scale factor to apply to the mc background
	const double dirt_scale_factor = 0.16411;

	// DEBUG
	bool DEBUG{true};
	bool UseGSimpleFlux{false};
	if (UseGSimpleFlux) std::cout << "\033[1;37mUsing GSimple Flux in the CV calcualtion not dk2nu\033[0m"<< std::endl; 

	// Histograms
	TGraph* gModelXsec;             		          // Graph of the Model vs cross section
	TH2D *hCV2dnue, *hCV2dnuebar, *hCV2d;

	// TTree
	TTree *DataTree;

	// TFile
	TFile *fCV;

	// Precalculate the weight histograms
	HistWeights nue("nue");
	HistWeights nuebar("nuebar");
	HistWeights numu("numu");
	HistWeights numubar("numubar");

	PrecalcHistRatio( nue,         "nue");
	PrecalcHistRatio( nuebar,   "nuebar");
	PrecalcHistRatio( numu,       "numu");
	PrecalcHistRatio( numubar, "numubar");

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Read in the events
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Load in the file containing the event number for Signal_Generated, N_gen or Signal_Selected N_sig
	if (DEBUG) std::cout << "\nNow reading in event files!" << std::endl;

	std::vector<std::string> class_type;
	std::vector<int> evtnum;
	std::vector<int> mc_nu_id;
	std::vector<double> mc_nu_dir_x;
	std::vector<double> mc_nu_dir_y;
	std::vector<double> mc_nu_dir_z;
	std::vector<double> mc_nu_energy;    
	std::vector<event> event;
	
	// Read in the selected list
	ReadEventList("filelists/selected_events_spaced.txt", evtnum, class_type, mc_nu_id, mc_nu_dir_x, mc_nu_dir_y, mc_nu_dir_z, mc_nu_energy    );
	
	// Read in the generated list
	ReadEventList("filelists/neutrino_in_tpc_list_spaced.txt", evtnum, class_type, mc_nu_id, mc_nu_dir_x, mc_nu_dir_y, mc_nu_dir_z, mc_nu_energy    );
	for (unsigned int i = 0; i < class_type.size(); i++){
	
		class event temp; // create a temp event class

		// Calculate theta
		double theta = GetTheta(mc_nu_dir_x[i], mc_nu_dir_y[i], mc_nu_dir_z[i]);

		// std::cout << theta << std::endl;

		// got a generated event
		if (evtnum[i] == 0 && (mc_nu_id[i] == 1 || mc_nu_id[i] == 5)){ // generated events have flag of zero as the event number, and we want only nue and nuebar cc
			// Push back generated events
			temp.E     =  mc_nu_energy[i];
			temp.Theta = theta;
			temp.type  = "gen";
			temp.nu_flav = mc_nu_id[i];
			event.push_back(temp);
			N_gen_CV++;
		}
		// got a selected event
		else {
		
			// Push back selected events
			temp.E     =  mc_nu_energy[i];
			temp.Theta = theta;
			temp.type  = "sel";
			temp.nu_flav = mc_nu_id[i];
			event.push_back(temp);
			N_sel_CV++;

			// Signal
			if ((class_type[i].compare(0,6,"nue_cc") == 0 && class_type[i] != "nue_cc_out_fv" && class_type[i] != "nue_cc_mixed") || (class_type[i].compare(0,10,"nue_bar_cc") == 0  && class_type[i] != "nue_bar_mixed") ) {
				temp.E     =  mc_nu_energy[i];
				temp.Theta = theta;
				temp.type  = "sig";
				temp.nu_flav = mc_nu_id[i];
				event.push_back(temp);
				N_sig_CV++;
			}
			// Dirt
			else if (class_type[i] == "Dirt"){
				
				temp.E     =  mc_nu_energy[i];
				temp.Theta = theta;
				temp.type  = "dirt";
				temp.nu_flav = mc_nu_id[i];
				event.push_back(temp);
				N_dirt_CV++;

			}
			// Background
			else {
				temp.E     =  mc_nu_energy[i];
				temp.Theta = theta;
				temp.type  = "bkg";
				temp.nu_flav = mc_nu_id[i];
				event.push_back(temp);
				N_bkg_CV++;
			
			}
		}
	}

	if (DEBUG) std::cout << "Total Generated in:  \t"<< N_gen_CV << std::endl;
	if (DEBUG) std::cout << "Total signal in:     \t"<< N_sig_CV << std::endl;
	if (DEBUG) std::cout << "Total selected in:   \t"<< N_sel_CV << std::endl;
	if (DEBUG) std::cout << "Total background in: \t"<< N_bkg_CV << std::endl;
	if (DEBUG) std::cout << "Total dirt in: \t\t"<< N_dirt_CV << std::endl;
		
	
	// Get the file with the CV and Hadron Production uncertainties
	bool boolfile  = GetFile(fCV , "/uboone/data/users/kmistry/work/PPFX/uboone/with_tilt_2Dhists/output.root"); if (boolfile == false) gSystem->Exit(0);

	// Loop over the parameters
	for (unsigned int i = 0; i < params.size() - 7; i++){
		if (i > 1) continue; // skip the beamline uncertainties for now...

		std::string param = params.at(i); // Get the parameter
		if (DEBUG) std::cout << "\n++++++++++++++++++++++++++" << std::endl;
		std::cout << "\033[1;31mRunning over Parameter:\033[0m\t\033[1;29m" << param << "\033[0m"<< std::endl;
		
		
		// Set the number of universes
		if (i == 1 ) Universes = 100; // change this to 100 for all the universes!
		else Universes = 1;

		// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// Event loop
		// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		for (unsigned int j = 0; j < event.size(); j++){ 
			// Alert the user
    		if (j % 1000 == 0 && i == 1) std::cout << "On entry " << j/1000 << "k"<< std::endl;
			// std::cout << event[j].type<< std::endl;
			
			// Call the Add weights function
			if      (event[j].type == "gen")  AddWeights(N_gen,  Universes, i, event[j],  nue, nuebar, numu, numubar);  // Gen
			else if (event[j].type == "sig")  AddWeights(N_sig,  Universes, i, event[j],  nue, nuebar, numu, numubar);  // Sig
			else if (event[j].type == "sel")  AddWeights(N_sel,  Universes, i, event[j],  nue, nuebar, numu, numubar);  // Sel
			else if (event[j].type == "bkg")  AddWeights(N_bkg,  Universes, i, event[j],  nue, nuebar, numu, numubar);  // Bkg
			else if (event[j].type == "dirt") AddWeights(N_dirt, Universes, i, event[j],  nue, nuebar, numu, numubar);  // Dirt
			else std::cout << "unknown type :("<< std::endl;
		
		} // END of event loop
		

		// ++++++++++++++++++++ Now calculate the new cross section in parameter iniverse i +++++++++++++++++++++++++++++

		Data_x_sec.resize(N_gen.size()); // Resize
		Efficiency.resize(N_gen.size()); 
		Flux.resize(N_gen.size()); 

		for (unsigned int k{0}; k < N_gen.size(); k++){
			if (i == 1) std::cout << "\033[1;32mUniverse:\t" << k << "\n\033[0m"<< std::endl;

			double flux = IntegrateFlux(k, fCV, i, POTScale);
			Flux[k] = flux; // save the flux
			std::cout << "\nflux: " << flux << std::endl;

			// Recalculations due to not weighting non MC genie stuff and other bugs
			N_gen[k] = N_gen[k];  
			// N_sel[k] = N_sel_CV; // ANDY F: do not reweight the selected events in MC so it is like data
			N_sel[k] = num_selected_data; // Use this for data x-section calcuation
			N_sig[k] = N_sig[k];
			N_bkg[k] = N_bkg[k];

			Efficiency[k] = N_sig[k] / N_gen[k];  // 0.0884133 CV efficiency
			if (DEBUG) std::cout << "\nEfficiency\t" << Efficiency[k] << std::endl;

			std::cout << "\nN_gen:\t" << N_gen[k] << std::endl;
			std::cout << "N_sel:\t" << N_sel[k] << std::endl;
			std::cout << "N_sig:\t" << N_sig[k] << std::endl;
			std::cout << "N_bkg:\t" << N_bkg[k] << std::endl;
			std::cout << "N_dirt:\t" << N_dirt[k] << "\n"<< std::endl;

			// Make flux calculation -- choice of gsimple or dk2nu flux controlled by bool at top of code
			if (UseGSimpleFlux ) Data_x_sec[k] = CalcDataXSec(num_selected_data, N_bkg[k], flux_gsimple_data, targets_data, intime_cosmics_bkg, intime_cosmic_scale_factor, N_dirt[k], dirt_scale_factor, mc_scale_factor, Efficiency[k] );
			else Data_x_sec[k] = CalcDataXSec(num_selected_data, N_bkg[k], flux, targets_data, intime_cosmics_bkg, intime_cosmic_scale_factor, N_dirt[k], dirt_scale_factor, mc_scale_factor, Efficiency[k] );
			if (DEBUG) std::cout << "\033[1;33mNew DATA X-section [10^-39 cm^2]\t\t" << Data_x_sec[k]/1.0e-39 << "\033[0m"<<std::endl;
			if (DEBUG) std::cout << "\n++++++++++++++++++++++++++++++++++++\n" << std::endl;
			
		}

		// Add the x section and efficiencies to a vector
		vec_Data_x_sec.push_back(Data_x_sec);
		vec_Efficiency.push_back(Efficiency);
		vec_gen.push_back( N_gen );
		vec_sig.push_back( N_sig );
		vec_dirt.push_back( N_dirt );
		vec_bkg.push_back( N_bkg );
		vec_flux.push_back(Flux);

		// Clear and do it again for next paramter
		N_gen.clear(); N_sig.clear(); N_bkg.clear(); N_sel.clear(); Data_x_sec.clear(), N_dirt.clear(); Efficiency.clear(); Flux.clear();
	
	} // END loop over all paramters

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Now we want to calucalte the uncertainties
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// First get the standard deviation of the efficiency and the Data cross section from the HP
	double HP_XSEC_STD = STD_Calculator(vec_Data_x_sec[1], vec_Data_x_sec[0][0]); // vec<double>, double
	double Eff_STD     = STD_Calculator(vec_Efficiency[1], vec_Efficiency[0][0]); 
	double HP_Gen_STD     = STD_Calculator(vec_gen[1], vec_gen[0][0]);
	double HP_Sig_STD     = STD_Calculator(vec_sig[1], vec_sig[0][0]);
	double HP_Dirt_STD    = STD_Calculator(vec_dirt[1], vec_dirt[0][0]);
	double HP_Bkg_STD     = STD_Calculator(vec_bkg[1], vec_bkg[0][0]);
	double HP_Flux_STD     = STD_Calculator(vec_flux[1], vec_flux[0][0]);

	std::cout << "HP_XSEC_STD:\t" <<  HP_XSEC_STD << "\tHP XSEC Err:\t" << (100 * HP_XSEC_STD) / vec_Data_x_sec[0][0] << " \%"<<std::endl;
	std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	std::cout << "Eff_STD:\t" <<  Eff_STD << "\tEfficiency Err:\t" << (100 * Eff_STD) / vec_Efficiency[0][0] << " \%" << std::endl;
	std::cout << "HP_Gen_STD:\t" <<  HP_Gen_STD << "\t\tHP Gen Err:\t" << (100 * HP_Gen_STD) / vec_gen[0][0] << " \%"<<std::endl;
	std::cout << "HP_Sig_STD:\t" <<  HP_Sig_STD << "\t\tHP Sig Err:\t" << (100 * HP_Sig_STD) / vec_sig[0][0] << " \%"<<std::endl;
	std::cout << "HP_Dirt_STD:\t" <<  HP_Dirt_STD << "\t\tHP Dirt Err:\t" << (100 * HP_Dirt_STD) / vec_dirt[0][0] << " \%"<<std::endl;
	std::cout << "HP_Bkg_STD:\t" <<  HP_Bkg_STD << "\t\tHP Bkg Err:\t" << (100 * HP_Bkg_STD) / vec_bkg[0][0] << " \%"<<std::endl;
	std::cout << "HP_Flux_STD:\t" <<  HP_Flux_STD << "\tHP Flux Err:\t" << (100 * HP_Flux_STD) / vec_flux[0][0] << " \%"<<std::endl;


	gSystem->Exit(0);
} // END
// ------------------------------------------------------------------------------------------------------------