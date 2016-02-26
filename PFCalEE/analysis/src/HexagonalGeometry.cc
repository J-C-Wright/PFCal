
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

        float a = edgeSize_;
        float b = edgeSize_*sqrt(3)*0.5;

        float A = isEdgeUp_ ? disp.Y() : disp.X();
        float B = isEdgeUp_ ? disp.X() : disp.Y();
        float dA,dB,eA,eB;

        eB = a/2.0 + fabs(b-fabs(A))/sqrt(3);
        eA = fabs(B) > a/2.0 ? fabs(eB - fabs(B))*sqrt(3) : b;
        
        dA = A > 0 ? A - eA : A + eA;
        dB = B > 0 ? B - eB : B + eB;

        float x,y;
        x = isEdgeUp_ ? dB : dA;
        y = isEdgeUp_ ? dA : dB;
        ROOT::Math::XYPoint edgeDists(x,y);

        return edgeDists;
    }


    ROOT::Math::XYPoint Hexagon::getDisplacementFromCentre_XY(ROOT::Math::XYPoint hit) {

        ROOT::Math::XYPoint displacement(hit.X() - centre_.X(), hit.Y() - centre_.Y());
        return displacement;

    }

    UVWEdgeDisplacements Hexagon::getDistanceToEdges_UVW(UVWPoint hitUVW){

        UVWPoint difference = hitUVW - uvwCentre_;
        float b = 0.5*edgeSize_*sqrt(3);

        float dU = difference.getU();
        float dV = difference.getV();
        float dW = difference.getW();

        UVWEdgeDisplacements edgeDists;
        edgeDists.dU = dU > 0 ? dU - b : dU + b;
        edgeDists.dV = dV > 0 ? dV - b : dV + b;
        edgeDists.dW = dW > 0 ? dW - b : dW + b;

        return edgeDists;
    }






