
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

    ROOT::Math::XYPoint Hexagon::getDistanceToEdges_XY(ROOT::Math::XYPoint hit) {

        std::pair<float,float> centreDisplacement = getDisplacementFromCentre_XY(hit);

        float dA,dB,A,B;
        float x,y;
        if (isEdgeUp_) {
            A = centreDisplacement.second;
            B = centreDisplacement.first;
        }else{
            A = centreDisplacement.first;
            B = centreDisplacement.second;
        }
        dA = fabs(edgeSize_ - A);
        dB = fabs(B - fabs( edgeSize_*0.5 + dA/sqrt(3) ));
        if (fabs(B) > edgeSize_*0.5) {
            dA = dB*sqrt(3);
        }

        if (isEdgeUp_) {
            if (A > 0 && B > 0) {
                dA *= -1.0;
                dB *= -1.0;
            }else if (A > 0 && B < 0) {
                dA *= -1.0;
            }else if (A < 0 && B > 0) {
                dB *= -1.0;
            }
            x = B;
            y = A;
        }else{
            if (A > 0 && B > 0) {
                dA *= -1.0;
                dB *= -1.0;
            }else if (A > 0 && B < 0) {
                dB *= -1.0;
            }else if (A < 0 && B > 0) {
                dA *= -1.0;
            }
            x = A;
            y = B;
        }

        ROOT::Math::XYPoint edgeDists(x,y);
        return edgeDists;
    }


    std::pair<float,float> Hexagon::getDisplacementFromCentre_XY(ROOT::Math::XYPoint hit) {

        std::pair<float,float> displacementFromCentre(9999.,9999.);

        displacementFromCentre.first  = hit.X() - centre_.X();
        displacementFromCentre.second = hit.Y() - centre_.Y();

        return displacementFromCentre;
    }

    std::vector<float> Hexagon::getDistanceToEdges_UVW(ROOT::Math::XYPoint hit) {

        UVWPoint hitUVW(hit);
        UVWPoint difference = uvwCentre_ - hitUVW;

        float dU = fabs(difference.getU() - edgeSize_*sqrt(3.0)/2.0);
        if (difference.getU() > 0) dU *= -1.0;

        float dV = fabs(difference.getV() - edgeSize_*sqrt(3.0)/2.0);
        if (difference.getV() > 0) dV *= -1.0;

        float dW = fabs(difference.getW() - edgeSize_*sqrt(3.0)/2.0);
        if (difference.getW() > 0) dW *= -1.0;
            
        std::vector<float> edgeDists(3);
        edgeDists[0] = dU;
        edgeDists[1] = dV;
        edgeDists[2] = dW;

        return edgeDists;
    }






