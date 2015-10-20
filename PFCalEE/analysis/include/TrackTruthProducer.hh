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

struct TrackInfo {
    unsigned showerStart;
    Float_t energyWeightedX[28];
    Float_t energyWeightedY[28];
    Float_t truthX[28];
    Float_t truthY[28];
    Float_t centralE[28];
    Float_t totalE[28];
};

class TrackTruth {

    private:
        HGCSSGenParticle particleInfo_;
        std::vector<ROOT::Math::XYPoint> truthPositions_;
        std::vector<ROOT::Math::XYPoint> energyWeightedXY_;
        std::vector<std::vector<HGCSSRecoHit>> hitsByLayer3x3_;
        std::vector<ROOT::Math::XYPoint> distsFromHitCentre_;

    public:
    //Setters
        void setParticleInfo(HGCSSGenParticle particleInfo) {particleInfo_ = particleInfo;}
        void setTruthPositions(std::vector<ROOT::Math::XYPoint> truthPositions) {truthPositions_ = truthPositions;}
        void setEnergyWeightedXY(std::vector<ROOT::Math::XYPoint> energyWeightedXY) {energyWeightedXY_ = energyWeightedXY;}
        void setHitsByLayer(std::vector<std::vector<HGCSSRecoHit>> hitsByLayer3x3) {hitsByLayer3x3_ = hitsByLayer3x3;}
        void setDistsFromHitCentre(std::vector<ROOT::Math::XYPoint> distsFromHitCentre) {distsFromHitCentre_ = distsFromHitCentre;}
    //Getters
        //Member vars
        HGCSSGenParticle getParticleInfo() {return particleInfo_;}
        std::vector<ROOT::Math::XYPoint> getTruthPositions() {return truthPositions_;}
        std::vector<ROOT::Math::XYPoint> getEnergyWeightedXY() {return energyWeightedXY_;}
        std::vector<std::vector<HGCSSRecoHit>> getHitsByLayer3x3() {return hitsByLayer3x3_;}
        std::vector<ROOT::Math::XYPoint> getDistsFromHitCentre() {return DistsFromHitCentre_;}
        //Derived
        unsigned getShowerStart() {
            for (unsigned layer(0);layer<energyWeightedXY_.size();layer++) {
                if (energyWeightedXY_[layer].X() != 0 && energyWeightedXY_[layer].Y()) {return layer;}
            }
        }
                
};

class TrackTruthProducer{

    private:
        bool debug_;
        bool layerZsLoaded_;
        unsigned nLayers_;
        unsigned versionNumber_;
        std::vector<std::vector<HGCSSRecoHit>> hitsByLayer_;
        std::vector<double> layerZPositions_;
        std::vector<TrackTruth> tracks_;

    public:
        TrackTruthProducer( bool debug, 
                            unsigned nLayers, 
                            unsigned versionNumber );

        void produce(std::vector<HGCSSGenParticle> *genvec,
                     std::vector<HGCSSRecoHit> *recoHitVec,
                     const HGCSSGeometryConversion & geomConv,
                     int mipCut);

        ~TrackTruthProducer(){};

        TrackTruth getTrack(unsigned index) {
            return tracks_[index];
        };

        std::vector<TrackTruth> getAllTracks() {
            return tracks_;
        };

        TrackInfo trackStruct(unsigned index);

        void clear(){
            hitsByLayer_.clear();
            tracks_.clear();
        };

};
