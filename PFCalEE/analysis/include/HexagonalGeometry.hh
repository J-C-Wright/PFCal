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

#define PI 3.14159265

using namespace std;

class UVWPoint {

    private: 
        float u_, v_, w_;
        float uPerp_, vPerp_, wPerp_;
    public:
        UVWPoint(){};
        //Insert creator from xy
        UVWPoint(ROOT::Math::XYPoint XY) {
            float phi = XY.Phi();
            float r   = XY.R();
            u_ = r*cos(phi-5.0*PI/6.0);
            v_ = r*cos(phi-PI/6.0);
            w_ = -1.0*XY.Y();
            uPerp_ = r*sin(phi-5.0*PI/6.0);
            vPerp_ = r*sin(phi-PI/6.0);
            wPerp_ = -1.0*XY.X();
        };

        float getU(){return u_;}
        float getV(){return v_;}
        float getW(){return w_;}

        float getUPerp(){return uPerp_;}
        float getVPerp(){return vPerp_;}
        float getWPerp(){return wPerp_;}

        void setUVW(ROOT::Math::XYPoint XY){
            float phi = XY.Phi();
            float r   = XY.R();
            u_ = r*cos(phi-5.0*PI/6.0);
            v_ = r*cos(phi-PI/6.0);
            w_ = -1.0*XY.Y();
            uPerp_ = r*sin(phi-5.0*PI/6.0);
            vPerp_ = r*sin(phi-PI/6.0);
            wPerp_ = -1.0*XY.X();
        }


        void print() {
            std::cout << setw(12) << u_ << setw(12) << v_ << setw(12) << w_ << std::endl;
            std::cout << setw(12) << uPerp_ << setw(12) << vPerp_ << setw(12) << wPerp_ << std::endl;
        }
        ROOT::Math::XYPoint getXYPoint() {
            ROOT::Math::XYPoint point( -1.0*wPerp_, w_*-1.0);
            return point;
        }

        //Operator overloading
        UVWPoint operator+(UVWPoint point) {

            float x = getXYPoint().X() + point.getXYPoint().X();
            float y = getXYPoint().Y() + point.getXYPoint().Y();
            ROOT::Math::XYPoint XY(x,y);
            UVWPoint newPoint(XY);

            return newPoint;
        }
        UVWPoint operator-(UVWPoint point) {

            float x = getXYPoint().X() - point.getXYPoint().X();
            float y = getXYPoint().Y() - point.getXYPoint().Y();
            ROOT::Math::XYPoint XY(x,y);
            UVWPoint newPoint(XY);

            return newPoint;
        }
};

        

            

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
        std::vector<float>     getDistanceToEdges_UVW(ROOT::Math::XYPoint hit);


};

