#include <string>
#include <iostream>

#include "HGCSSGenParticle.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGeometryConversion.hh"

#include <utility>
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

using namespace std;

class Hexagon {

    private:
        float edgeSize_;
        bool  isEdgeUp_;
        ROOT::Math::XYPoint centre_;

    public:
        //Creator and destructor
        Hexagon(float edgeSize, bool isEdgeUp, ROOT::Math::XYPoint centre) {
            edgeSize_ = edgeSize;
            isEdgeUp_ = isEdgeUp;
            centre_ = centre;
        };
        ~Hexagon(){}

        //Setters 
        void setEdgeSize(float edgeSize) { edgeSize_ = edgeSize; }
        void setIsEdgeUp(bool isEdgeUp) { isEdgeUp_ = isEdgeUp; }
        void setCentre(ROOT::Math::XYPoint centre) { centre_ = centre; }

        //Getters
        float getEdgeSize() { return edgeSize_; } 
        bool isEdgeUp() { return isEdgeUp_; }
        ROOT::Math::XYPoint getCentre() { return centre_; }

        //Derived
        std::pair<unsigned,unsigned> getNearestEdges_XY(ROOT::Math::XYPoint hit);
        std::pair<float,float> getDistanceToEdges_XY(ROOT::Math::XYPoint hit);
        std::pair<float,float> getDisplacementFromCentre_XY(ROOT::Math::XYPoint hit);


};

