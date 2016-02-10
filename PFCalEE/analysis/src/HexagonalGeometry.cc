
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

        ROOT::Math::XYPoint disp = getDisplacementFromCentre_XY(hit);

        float dA,dB,A,B;
        float a = edgeSize_;
        float b = edgeSize_*sqrt(3)*0.5;

        if (isEdgeUp_) {
            A = disp.Y();
            B = disp.X();
        }else{
            A = disp.X();
            B = disp.Y();
        }
        dA = fabs(b - fabs(A));
        dB = fabs(a*0.5 + fabs(dA)/sqrt(3) - fabs(B));
        if (fabs(B) > a*0.5) {
            dA = dB*sqrt(3);
        }
        /*
        if (A > 0 && B > 0) {
            dA *= -1.0;
            dB *= -1.0;
        }else if (A > 0 && B < 0) {
            dB *= -1.0;
        }else if (A < 0 && B > 0) {
            dA *= -1.0;
        }
        */

        float x,y;
        if (isEdgeUp_) {
            x = dB;
            y = dA;
        }else{
            x = dA;
            y = dB;
        }

        if (disp.X() > 0 && disp.Y() > 0) {
            x *= -1.0;
            y *= -1.0;
        }else if (disp.X() > 0 && disp.Y() < 0) {
            x *= -1.0;
        }else if (disp.X() < 0 && disp.Y() > 0) {
            y *= -1.0;
        }

        ROOT::Math::XYPoint edgeDists(x,y);
        return edgeDists;
    }


    ROOT::Math::XYPoint Hexagon::getDisplacementFromCentre_XY(ROOT::Math::XYPoint hit) {

        ROOT::Math::XYPoint displacement(hit.X() - centre_.X(), hit.Y() - centre_.Y());
        return displacement;

    }

    UVWEdgeDisplacements Hexagon::getDistanceToEdges_UVW(UVWPoint hitUVW){

        UVWPoint difference = hitUVW - uvwCentre_;

        float dU = fabs(fabs(difference.getU()) - edgeSize_*sqrt(3.0)/2.0);
        if (difference.getU() > 0) dU *= -1.0;

        float dV = fabs(fabs(difference.getV()) - edgeSize_*sqrt(3.0)/2.0);
        if (difference.getV() > 0) dV *= -1.0;

        float dW = fabs(fabs(difference.getW()) - edgeSize_*sqrt(3.0)/2.0);
        if (difference.getW() > 0) dW *= -1.0;
            
        UVWEdgeDisplacements edgeDists;
        edgeDists.dU = dU;
        edgeDists.dV = dV;
        edgeDists.dW = dW;

        return edgeDists;
    }






