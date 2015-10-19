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

struct TrackTruth {

    HGCSSGenParticle particleInfo;
    std::vector<ROOT::Math::XYPoint> truthPositions;
    std::vector<ROOT::Math::XYPoint> energyWeightedXY;
    std::vector<std::vector<HGCSSRecoHit>> hitsByLayer3x3;
    std::vector<ROOT::Math::XYPoint> distsFromHitCentre;
    std::vector<std::pair<float,float>> energyWeightedrPhi;
    std::vector<std::pair<float,float>> truthPositionsrPhi;
    std::vector<std::pair<float,float>> distsFromHitCentrerPhi;
    
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

        void clear(){
            hitsByLayer_.clear();
            tracks_.clear();
        };

};
