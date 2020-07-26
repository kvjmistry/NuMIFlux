#ifndef FUNCTIONS_MAKEHIST_H
#define FUNCTIONS_MAKEHIST_H

// Functions file for make hist script
// Contains the predetermined detector geometries too

// Std Includes
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <stdlib.h>
#include <string>
#include <vector>
#include <chrono> 

// ROOT Includes
#include "TDirectory.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TFile.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "canvas/Utilities/InputTag.h"

// LArsoft Includes
#include "geo/GeoVector.h"
#include "geo/GeoAABox.h"
#include "geo/GeoHalfLine.h"
#include "geo/GeoAlgo.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"
#include "dk2nu/tree/calcLocationWeights.h"

//___________________________________________________________________________
TVector3 FromDetToBeam( const TVector3 det, bool rotate_only ) {

    TVector3 beam;
    TRotation R;
    bool debug{false};

    // Rotation matrix using the 0,0,0 position for MicroBooNE (beam to det input)
    TVector3 Rot_row_x = { 0.92103853804025681562, 0.022713504803924120662, 0.38880857519374290021  };
    TVector3 Rot_row_y = { 4.6254001262154668408e-05, 0.99829162468141474651, -0.058427989452906302359 };
    TVector3 Rot_row_z = { -0.38947144863934973769, 0.053832413938664107345, 0.91946400794392302291   };

    R.RotateAxes(Rot_row_x, Rot_row_y, Rot_row_z); // Also inverts so now to det to beam
    R.Invert(); // go back to beam to det
    if (debug) {
        std::cout << "R_{beam to det} = " << std::endl;
        std::cout << " [ " << R.XX() << " " << R.XY() << " " << R.XZ() << " ] " << std::endl;
        std::cout << " [ " << R.YX() << " " << R.YY() << " " << R.YZ() << " ] " << std::endl;
        std::cout << " [ " << R.ZX() << " " << R.ZY() << " " << R.ZZ() << " ] " << std::endl;
        std::cout << std::endl;
    }
    R.Invert(); // R is now the inverse
    if (debug) {
        std::cout << "R_{det to beam} = " << std::endl;
        std::cout << " [ " << R.XX() << " " << R.XY() << " " << R.XZ() << " ] " << std::endl;
        std::cout << " [ " << R.YX() << " " << R.YY() << " " << R.YZ() << " ] " << std::endl;
        std::cout << " [ " << R.ZX() << " " << R.ZY() << " " << R.ZZ() << " ] " << std::endl;
        std::cout << std::endl;
    }

    beam = R * det;               // Only rotate the vector
    
    return beam;
}
//___________________________________________________________________________
// Function to check the status of weights, if they are bad, set to zero
void check_weight(float &weight){
    
    if (weight < 0) weight = 0; // get rid of them pesky negative weights
    
    if (std::isnan(weight) == 1) { // catch NaN values
        std::cout << "got a nan:\t" << weight << std::endl;
        weight = 0;
    }

    // Catch infinate values
    if (std::isinf(weight)) { // catch NaN values
        std::cout << "got a inf:\t" << weight << std::endl;
        weight = 0;
    }
}
//___________________________________________________________________________
// Get the PPFX weights
void GetPPFXWeights(std::vector< std::vector< double > > &Weights, int shift, std::map<std::string, std::vector<double>> evtwgt_map, std::vector<std::string> labels){


    // Loop over all options
    for (unsigned l=0; l<labels.size(); l++) { 

        // Look for string name wishing to push back
        if (evtwgt_map.find(labels[l]) != evtwgt_map.end()) {

            // Loop over ms universes
            for (unsigned i=0; i < evtwgt_map.find(labels[l])->second.size(); i++) { 

                // Update the weight if it is negative or greater than 50
                if (evtwgt_map.find(labels[l])->second.at(i) < 0 || evtwgt_map.find(labels[l])->second.at(i) > 50){ 
                    Weights[l][i+shift] *= 1; 
                }     
                else {
                    Weights[l][i+shift] *= evtwgt_map.find(labels[l])->second.at(i);
                    // std::cout <<evtwgt_map.find(labels[l])->second.at(i) << std::endl;
                    if (std::isnan(Weights[l][i+shift]) == 1) { // catch NaN values
                        std::cout << "got a nan:\t" << Weights[l][i+shift] << std::endl;
                        Weights[l][i+shift] = 0;
                    }

                    // Catch infinate values
                    if (std::isinf(Weights[l][i+shift])) { // catch NaN values
                        std::cout << "got a inf:\t" << Weights[l][i+shift] << std::endl;
                        Weights[l][i+shift] = 0;
                    }
                }

            } // End loop over universes

        }

    }

    return;
}
//___________________________________________________________________________






#endif