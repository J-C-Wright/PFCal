#include <string>
#include <iostream>
#include "HGCSSGenParticle.hh"

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TVector3.h"

#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"

#include "PositionTruth.h"


int main(int argc, char** argv){

    //////
    //Test
    //////
    int nLayers = 51;
    unsigned int versionNumber = 33;
    PositionTruth positionTruth(nLayers);
    std::cout << "Hello world! Here are the HGC layer z positions" << std::endl;
    positionTruth.getLayerZPositions(versionNumber);

    ////////////
    //Input Test
    ////////////






    return 0;
}














