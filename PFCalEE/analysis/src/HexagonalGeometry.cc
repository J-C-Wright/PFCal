
#define _USE_MATH_DEFINES

#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include <math.h>

#include "HGCSSGenParticle.hh"
#include "HGCSSRecoHit.hh"

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

        std::pair<unsigned,unsigned> closestEdges(9999,9999);
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
        std::pair<float,float> distanceToEdges(9999.,9999.);

        float dY = fabs(fabs(centreDisplacement.second) - edgeSize_*sqrt(3)/2);
        float b = fabs(edgeSize_*0.5 + (1/sqrt(3))*dY - fabs(centreDisplacement.first));
        float c = dY;
        if (fabs(centreDisplacement.first) > edgeSize_/2) {
            c = b*sqrt(3);
        }

        if (centreDisplacement.first < 0 && centreDisplacement.second > 0){
            c *= -1.0;
        }else if (centreDisplacement.first > 0 && centreDisplacement.second > 0){
            b *= -1.0;
            c *= -1.0;
        }else if (centreDisplacement.first > 0 && centreDisplacement.second < 0){
            b *= -1.0; 
        }

        distanceToEdges.first  = b;
        distanceToEdges.second = c;

        return distanceToEdges;    
    }
    std::pair<float,float> Hexagon::getDisplacementFromCentre_XY(ROOT::Math::XYPoint hit) {

        std::pair<float,float> displacementFromCentre(9999.,9999.);

        displacementFromCentre.first  = hit.X() - centre_.X();
        displacementFromCentre.second = hit.Y() - centre_.Y();

        return displacementFromCentre;
    }

    std::vector<float> Hexagon::getDistanceToEdges_UVW(ROOT::Math::XYPoint hit) {

        UVWPoint centreUVW(centre_);
        UVWPoint hitUVW(hit);
        UVWPoint distFromCentre = hitUVW - centreUVW;

        float dU = fabs(distFromCentre.getU() - edgeSize_*sqrt(3.0)/2.0);
        if (distFromCentre.getU() > 0) dU *= -1.0;
        float dV = fabs(distFromCentre.getV() - edgeSize_*sqrt(3.0)/2.0);
        if (distFromCentre.getV() > 0) dV *= -1.0;
        float dW = fabs(distFromCentre.getW() - edgeSize_*sqrt(3.0)/2.0);
        if (distFromCentre.getW() > 0) dW *= -1.0;
            
        std::vector<float> edgeDists(3);
        edgeDists[0] = dU;
        edgeDists[1] = dV;
        edgeDists[2] = dW;

        return edgeDists;
    }






