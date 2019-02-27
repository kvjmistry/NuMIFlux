/*
A script to read in a flugg file and get the number of Kaon events in the file
it will then compare thid to the number of events in a dk file.

*/


#include "plot_comp_functions.h"

void flugg_comparison(){

	double KaonEvents_Flugg;
	double POT_Flugg;
	double KaonEvents_dk;
	double POT_dk;

	// ------------------------------------------------------------------------------------------------------------
	// FLUGG
	// ------------------------------------------------------------------------------------------------------------
	// File loop
	for (int i = 0 ; i < 0; i++){
	
		// ------------------------------------------------------------------------------------------------------------
		// Grab the flugg files
		// ------------------------------------------------------------------------------------------------------------
		char name[1000];
		if (i < 10) snprintf(name, 1000, "/uboone/data/flux/numi/v2/flugg_mn000z200i_rp11_lowth_pnut_f112c0f093bbird/flugg_mn000z200i_rp11_bs1.1_pnut_lowth_f112c0f093bbird_000%d.root" ,i);
		else if  (i >= 10 && i < 100) snprintf(name, 1000, "/uboone/data/flux/numi/v2/flugg_mn000z200i_rp11_lowth_pnut_f112c0f093bbird/flugg_mn000z200i_rp11_bs1.1_pnut_lowth_f112c0f093bbird_00%d.root" , i);
		else snprintf(name, 1000, "/uboone/data/flux/numi/v2/flugg_mn000z200i_rp11_lowth_pnut_f112c0f093bbird/flugg_mn000z200i_rp11_bs1.1_pnut_lowth_f112c0f093bbird_0%d.root" , i);

		TFile* f_flugg;
		bool bool_file = GetFile(f_flugg , name); 
		if (bool_file == false) continue; // Get the file, if it doesnt exist, then try the next one
		else std::cout << "Got the file:\t"<< name <<std::endl;
		// ------------------------------------------------------------------------------------------------------------

		TTree* flugg_tree = (TTree*)f_flugg->Get("h10");
		int ptype_fg, evtno;
		double Nimpwt_fg;
		TBranch *b_ptype_fg, *b_Nimpwt_fg, *b_evtno;
		flugg_tree->SetBranchAddress("ptype", &ptype_fg, &b_ptype_fg);
		flugg_tree->SetBranchAddress("Nimpwt", &Nimpwt_fg, &b_Nimpwt_fg);
		flugg_tree->SetBranchAddress("evtno", &evtno, &b_evtno);

		// ------------------------------------------------------------------------------------------------------------
		// Event loop in file
		// ------------------------------------------------------------------------------------------------------------
		int i_entry{0};
		while(flugg_tree->GetEntry(i_entry)) {
			++i_entry;

			if (ptype_fg == 11 || ptype_fg == 12){ // Select Kaons
			KaonEvents_Flugg+=Nimpwt_fg; // Count Kaons

			}

		}
		POT_Flugg+=evtno; // Increment the POT_Flugg

		f_flugg->Close();

	}
	std::cout << "\n-------------------------------------------" << std::endl;
	std::cout << "Total Kaon Events in Flugg File:\t" << KaonEvents_Flugg << std::endl;
	std::cout << "Total POT in Flugg File:\t\t" << POT_Flugg << std::endl;
	std::cout << "-------------------------------------------\n" << std::endl;


}