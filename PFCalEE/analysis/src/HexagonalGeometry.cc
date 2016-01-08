
#define _USE_MATH_DEFINES

#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include <math.h>

#include "HGCSSGenParticle.hh"
#include "HGCSSRecoHit.hh"
#include "TrackTruthProducer.hh"

#include "HexagonalGeometry.hh"

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

#define PI 3.14159265

    std::pair<unsigned,unsigned> Hexagon::getNearestEdges_XY(ROOT::Math::XYPoint hit) {

        std::pair<unsigned,unsigned> closestEdges(999,999);
        //First is closest in x, second is closest in y
        //See personal notes, 7/1/16 pg 2

        //Is it closer to the top or bottom edge?
        bool closerToTop = hit.Y() > centre_.Y();
        //Closer to the left or right?
        bool closerToRight = hit.X() > centre_.X();
        //Is it in the central region?
        bool inCentralRegion = fabs( hit.X() - centre_.X() ) < edgeSize_/2.0;

        //Associate to nearest edges
        if (closerToTop && !closerToRight && !inCentralRegion) {
            closestEdges.first  = 0;
            closestEdges.second = 0;
        }else if (closerToTop && !closerToRight && inCentralRegion) {
            closestEdges.first  = 0;
            closestEdges.second = 1;
        }else if (closerToTop && closerToRight && inCentralRegion) {
            closestEdges.first  = 2;
            closestEdges.second = 1;
        }else if (closerToTop && closerToRight && !inCentralRegion) {
            closestEdges.first  = 2;
            closestEdges.second = 2;
        }else if (!closerToTop && closerToRight && !inCentralRegion) {
            closestEdges.first  = 3;
            closestEdges.second = 3;
        }else if (!closerToTop && closerToRight && inCentralRegion) {
            closestEdges.first  = 3;
            closestEdges.second = 4;
        }else if (!closerToTop && !closerToRight && inCentralRegion) {
            closestEdges.first  = 5;
            closestEdges.second = 4;
        }else if (!closerToTop && closerToRight && inCentralRegion) {
            closestEdges.first  = 5;
            closestEdges.second = 5;
        }

        return closestEdges;

    }

    std::pair<float,float> Hexagon::getDistanceToEdges_XY(ROOT::Math::XYPoint hit) {

        std::pair<unsigned,unsigned> closestEdges = getNearestEdges_XY(hit); 
        std::pair<float,float> centreDisplacement = getDisplacementFromCentre_XY(hit);
        std::pair<float,float> distanceToEdges(-999.,-999.);

        //Is it in the central region of the hex?
        if (closestEdges.first != closestEdges.second) {

            //Y
            distanceToEdges.second = centreDisplacement.second - edgeSize_*sqrt(3);
            if ( centreDisplacement.second < 0 ) distanceToEdges.second *= -1.0;
            //X
            distanceToEdges.first = centreDisplacement.first - (edgeSize_/2.0) - fabs(distanceToEdges.second)/sqrt(3);
            if ( centreDisplacement.first < 0 ) distanceToEdges.first *= -1.0;

        }else if (closestEdges.first == closestEdges.second && closestEdges.first < 6) {
            
            //Calculate the displacements to edge
            float delta = atan( fabs(edgeSize_ - centreDisplacement.first / centreDisplacement.second) );
            float d = (edgeSize_*sqrt(3)/4) * 1.0/sin(5*PI/6 - delta);
            distanceToEdges.second = 2.0 * (d - sqrt( pow(edgeSize_/2.0 - fabs(centreDisplacement.first), 2) + pow(centreDisplacement.second,2)))
                                               * sin(5*PI/6 - delta);
            distanceToEdges.first  = distanceToEdges.second/sqrt(3);

            //Set the minus signs according to the associated edge
            if (closestEdges.first == 0) {
                distanceToEdges.second *= -1.0;
            }else if (closestEdges.first == 2) {
                distanceToEdges.first *= -1.0;
                distanceToEdges.second *= -1.0;
            }else if (closestEdges.first == 3) {
                distanceToEdges.first *= -1.0;
            }

        }
        return distanceToEdges;    

    }


    std::pair<float,float> Hexagon::getDisplacementFromCentre_XY(ROOT::Math::XYPoint hit) {

        std::pair<float,float> displacementFromCentre(-999.,-999.);

        displacementFromCentre.first  = hit.X() - centre_.X();
        displacementFromCentre.second = hit.Y() - centre_.Y();

        return displacementFromCentre;

    }








