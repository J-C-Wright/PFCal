#include <string>
#include <iostream>

#include "HGCSSGenParticle.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGeometryConversion.hh"

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

struct PositionTruth {
    unsigned int eventNumber;
    unsigned int particleNumber;
    double x0;
    double y0;
    double z0;
    double eta;
    double phi;
    std::vector<ROOT::Math::XYPoint> truthPositions;
    std::vector<std::pair<float,float>> distsFromHitCentre;
}; 

class PositionTruthProducer{

    public: 
        PositionTruthProducer(unsigned int nLayers, bool debug){
            layerZsLoaded_ = false;
            nLayers_ = nLayers;
            debug_ = debug;
        };
    
        ~PositionTruthProducer(){};

        void calcTruthPositions(std::vector<HGCSSGenParticle> *genvec, unsigned int eventNumber);

        void getLayerZPositions(const unsigned versionNumber);

        PositionTruth getPosition(unsigned int entry){
            return positionTruthEntries_[entry];
        };

        std::vector<PositionTruth> getAllPositions(){
            return positionTruthEntries_;
        };
   
        std::vector<std::vector<ROOT::Math::XYPoint>> getAllEnergyWeightedXY() {
            return energyWeightedXY_; 
        } 

        std::vector<ROOT::Math::XYPoint> getEnergyWeightedXY(unsigned int index) {
            return energyWeightedXY_[index]; 
        } 
        
        bool hasLayerZs(){
            return layerZsLoaded_;
        };

        unsigned int numberOfPhotons(std::vector<HGCSSGenParticle> *genvec);

        void calcEnergyWeightedXYMax(std::vector<HGCSSRecoHit> *recoHitVec, const HGCSSGeometryConversion & geomConv, const unsigned nSR);

        void calcEnergyWeightedXYTruth(std::vector<HGCSSRecoHit> *recoHitVec, const HGCSSGeometryConversion & geomConv, const unsigned nSR);

        void clearEntries() {
            positionTruthEntries_.clear();
            energyWeightedXY_.clear();
        };

    private:
        bool layerZsLoaded_;
        bool debug_;
        unsigned int nLayers_;
        std::vector<double> layerZPositions_;
        std::vector<PositionTruth> positionTruthEntries_;
        std::vector<std::vector<ROOT::Math::XYPoint>> energyWeightedXY_;

};
