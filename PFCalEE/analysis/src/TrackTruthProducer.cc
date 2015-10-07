
#define _USE_MATH_DEFINES

#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include "HGCSSGenParticle.hh"
#include "HGCSSRecoHit.hh"
#include "TrackTruthProducer.hh"

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
#include "MyFunctions.h"

using namespace std;

    TrackTruthProducer::TrackTruthProducer(bool debug, 
                                           unsigned nLayers, 
                                           unsigned versionNumber){
        debug_ = debug;
        nLayers_ = nLayers;
        versionNumber_ = versionNumber;

        //Find the layer Zs
        std::ifstream inputFile;
        std::ostringstream inputFileName;
        inputFileName << "data/zPositions_v"<< versionNumber_ << ".dat";
        inputFile.open(inputFileName.str());
        if (!inputFile.is_open()){
            std::cout << "Can't open input file " << inputFileName.str() << std::endl;
            layerZsLoaded_ = false;
        }
        while (!inputFile.eof()) {
            unsigned l=nLayers_;
            double z=0;
            inputFile >> l >> z;
            if (l < nLayers_) {
                layerZPositions_.push_back(z);
                if (debug_) {std::cout << "Layer " << l << ", z = " << z << std::endl;}
            }
        }
        if (layerZPositions_.size() != nLayers_) {
            std::cout << "Warning! Problem in extracting z positions, did not find one value per layer. Please check input file: ";
            std::cout << inputFileName.str() << std::endl;
            layerZsLoaded_ = false;
            return;
        }
        inputFile.close();
        layerZsLoaded_ = true;
   
    }

    void TrackTruthProducer::produce(std::vector<HGCSSGenParticle> *genvec,
                                        std::vector<HGCSSRecoHit> *recoHitVec,
                                        const HGCSSGeometryConversion & geomConv){
        
        //Load hit info
        std::vector<std::vector<HGCSSRecoHit>> hitsByLayer_(nLayers_);
        for (unsigned int hitLoop(0);hitLoop<(*recoHitVec).size();hitLoop++) {
            HGCSSRecoHit hit = (*recoHitVec)[hitLoop];
            hitsByLayer_[hit.layer()].push_back(hit);
        }

        //Loop over genParticles
        std::vector<TrackTruth> trackVec((*genvec).size());
        for (unsigned int trackLoop(0);trackLoop<(*genvec).size();trackLoop++) {
            
            //Calculate the trackposition truth
            trackVec[trackLoop].particleInfo = (*genvec)[trackLoop]; 
            for (unsigned int layer(0);layer<nLayers_;layer++) {
                double x = trackVec[trackLoop].particleInfo.x() 
                            + (layerZPositions_[layer]-trackVec[trackLoop].particleInfo.z())*((*genvec)[trackLoop].px())/((*genvec)[trackLoop].pz());
                double y = trackVec[trackLoop].particleInfo.y() 
                            + (layerZPositions_[layer]-trackVec[trackLoop].particleInfo.z())*((*genvec)[trackLoop].py())/((*genvec)[trackLoop].pz());
                trackVec[trackLoop].truthPositions.push_back(ROOT::Math::XYPoint(x,y));
                std::pair<float,float> rPhi = rPhiCoordinates(trackVec[trackLoop].truthPositions[layer].X(),trackVec[trackLoop].truthPositions[layer].Y());
                trackVec[trackLoop].truthPositionsrPhi.push_back(rPhi);
            }

            //Energy-weighted positions, hits,  and position within cell
            for (unsigned layerLoop(0);layerLoop<nLayers_;layerLoop++) {

                if (hitsByLayer_[layerLoop].size() != 0) {

                    //Find closest hit to truth track
                    float distToTruth(9999.);
                    float closestCellIndex(0);
                    for (unsigned int hitLoop(0);hitLoop<hitsByLayer_[layerLoop].size();hitLoop++) {
                        float dx = hitsByLayer_[layerLoop][hitLoop].get_x() - trackVec[trackLoop].truthPositions[layerLoop].X();
                        float dy = hitsByLayer_[layerLoop][hitLoop].get_y() - trackVec[trackLoop].truthPositions[layerLoop].Y();
                        float dr = sqrt(pow(dx,2)+pow(dy,2));
                        if (dr < distToTruth) {
                            distToTruth = dr;
                            closestCellIndex = hitLoop;
                        }
                    }     

                    if (hitsByLayer_[layerLoop][closestCellIndex].energy() < 2) {

                        //Treat hit as empty, fill with dummy info
                        std::cout << "---- " << "In layer " << layerLoop << " Energy of central hit is < 2 MIPs ----" << std::endl;
                        std::vector<HGCSSRecoHit> empty(0);
                        trackVec[trackLoop].hitsByLayer3x3.push_back(empty);
                        //XY
                        ROOT::Math::XYPoint point(0,0);
                        trackVec[trackLoop].energyWeightedXY.push_back(point);
                        ROOT::Math::XYPoint dispFromCentre(9999.,9999.);
                        trackVec[trackLoop].distsFromHitCentre.push_back(dispFromCentre);
                        //rPhi
                        std::pair<float,float> displacementrPhiEmpty(9999.,9999.);
                        std::pair<float,float> energyWeightedrPhiEmpty(0,0);
                        trackVec[trackLoop].energyWeightedrPhi.push_back(energyWeightedrPhiEmpty);
                        trackVec[trackLoop].distsFromHitCentrerPhi.push_back(displacementrPhiEmpty);

                    }else{

                        //Hit is good, do EW calcs
                        double radialDisplacement = sqrt(pow(hitsByLayer_[layerLoop][closestCellIndex].get_x(),2)+pow(hitsByLayer_[layerLoop][closestCellIndex].get_y(),2));
                        double step = geomConv.cellSize(layerLoop,radialDisplacement)+0.1;

                        //Calculate the distance from the truth to the centre of the closest cell
                        //XY
                        float dx = (trackVec[trackLoop].truthPositions[layerLoop].X()
                                                    - hitsByLayer_[layerLoop][closestCellIndex].get_x())/geomConv.cellSize(layerLoop,radialDisplacement);
                        float dy = (trackVec[trackLoop].truthPositions[layerLoop].Y()
                                                    - hitsByLayer_[layerLoop][closestCellIndex].get_y())/geomConv.cellSize(layerLoop,radialDisplacement);
                        ROOT::Math::XYPoint displacement(dx,dy);
                        trackVec[trackLoop].distsFromHitCentre.push_back(displacement);
                        //r rPhi
                        std::pair<float,float> rPhi = rPhiCoordinates(hitsByLayer_[layerLoop][closestCellIndex].get_x(),hitsByLayer_[layerLoop][closestCellIndex].get_y());
                        std::pair<float,float> displacementrPhi(trackVec[trackLoop].truthPositionsrPhi[layerLoop].first-rPhi.first,
                                                                trackVec[trackLoop].truthPositionsrPhi[layerLoop].second-rPhi.second);
                        trackVec[trackLoop].distsFromHitCentrerPhi.push_back(displacementrPhi);

                        //find the hits that belong to the 3x3 grid centred on the closest
                        std::vector<HGCSSRecoHit> retainedHits;
                        for (unsigned int hitLoop(0);hitLoop<hitsByLayer_[layerLoop].size();hitLoop++) {
                            if (fabs(hitsByLayer_[layerLoop][hitLoop].get_x() - hitsByLayer_[layerLoop][closestCellIndex].get_x() ) < step &&
                                fabs(hitsByLayer_[layerLoop][hitLoop].get_y() - hitsByLayer_[layerLoop][closestCellIndex].get_y() ) < step) {
                                retainedHits.push_back(hitsByLayer_[layerLoop][hitLoop]);
                            }
                        } 
                        trackVec[trackLoop].hitsByLayer3x3.push_back(retainedHits);

                        //Calculate energy-weighted position
                        //X Y position
                        double energyWeightedX(0.0), energyWeightedY(0.0);
                        double totalEnergy(0.0);
                        for (unsigned int hitLoop(0);hitLoop<retainedHits.size();hitLoop++) {
                            energyWeightedX += retainedHits[hitLoop].get_x()*retainedHits[hitLoop].energy();
                            energyWeightedY += retainedHits[hitLoop].get_y()*retainedHits[hitLoop].energy();
                            totalEnergy += retainedHits[hitLoop].energy();
                        }
                        energyWeightedX /= totalEnergy;
                        energyWeightedY /= totalEnergy;
                        ROOT::Math::XYPoint point(energyWeightedX,energyWeightedY);
                        trackVec[trackLoop].energyWeightedXY.push_back(point);
                        //r rPhi
                        double energyWeightedr(0.0), energyWeightedPhi(0.0);
                        for (unsigned int hitLoop(0);hitLoop<retainedHits.size();hitLoop++) {
                            std::pair<float,float> tilerPhi = rPhiCoordinates(retainedHits[hitLoop].get_x(),retainedHits[hitLoop].get_y());
                            energyWeightedr   += tilerPhi.first*retainedHits[hitLoop].energy(); 
                            energyWeightedPhi += tilerPhi.second*retainedHits[hitLoop].energy(); 
                        }
                        energyWeightedr /= totalEnergy;
                        energyWeightedPhi /= totalEnergy;
                        std::pair<float,float> pointrPhi(energyWeightedr,energyWeightedPhi);
                        trackVec[trackLoop].energyWeightedrPhi.push_back(pointrPhi);

                    }
                }else {

                    std::cout << "---- " << "In layer " << layerLoop << " there are no hits ----" << std::endl;
                    std::vector<HGCSSRecoHit> empty(0);
                    trackVec[trackLoop].hitsByLayer3x3.push_back(empty);
                    //XY
                    ROOT::Math::XYPoint point(0,0);
                    trackVec[trackLoop].energyWeightedXY.push_back(point);
                    ROOT::Math::XYPoint dispFromCentre(9999.,9999.);
                    trackVec[trackLoop].distsFromHitCentre.push_back(dispFromCentre);
                    //rPhi
                    std::pair<float,float> displacementrPhiEmpty(9999.,9999.);
                    std::pair<float,float> energyWeightedrPhiEmpty(0,0);
                    trackVec[trackLoop].energyWeightedrPhi.push_back(energyWeightedrPhiEmpty);
                    trackVec[trackLoop].distsFromHitCentrerPhi.push_back(displacementrPhiEmpty);
                    
                }
            }

            //Print info to screen
            if (debug_) {
                std::cout << "Track " << trackLoop << std::endl;
                std::cout << setw(12) << "x0"  << setw(12) << "y0" << setw(12) << "z0" << setw(12) << "PDG Id" ;
                std::cout << setw(12) << "Eta" << setw(12) << "Phi" << std::endl;
                std::cout << setw(12) << trackVec[trackLoop].particleInfo.x();
                std::cout << setw(12) << trackVec[trackLoop].particleInfo.y();
                std::cout << setw(12) << trackVec[trackLoop].particleInfo.z();
                std::cout << setw(12) << trackVec[trackLoop].particleInfo.pdgid();
                std::cout << setw(12) << trackVec[trackLoop].particleInfo.eta();
                std::cout << setw(12) << trackVec[trackLoop].particleInfo.phi() << std::endl;
                std::cout << setw(12) << "Hit locations over " << nLayers_ << " layers" << std::endl;
                std::cout << setw(24) << "Truth" << setw(24) << "E-Weighted" << setw(24) << "D from centre" << setw(24) << "Truth (rPhi)";
                std::cout << setw(24) << "E-Weighted  (rPhi)" << setw(24) << "Dist from centre (rPhi)" << std::endl;
                std::cout << setw(12) << "X" << setw(12) << "Y" << setw(12) << "X" << setw(12) << "Y";
                std::cout << setw(12) << "X" << setw(12) << "Y";
                std::cout << setw(12) << "r" << setw(12) << "rPhi"; 
                std::cout << setw(12) << "r" << setw(12) << "rPhi" << setw(12) << "r" << setw(12) << "rPhi";
                std::cout << setw(12) << "Layer";
                std::cout << std::endl;
                for (unsigned layerLoop(1);layerLoop<nLayers_;layerLoop++) {
                    std::cout << setw(12) << trackVec[trackLoop].truthPositions[layerLoop].X();
                    std::cout << setw(12) << trackVec[trackLoop].truthPositions[layerLoop].Y();
                    std::cout << setw(12) << trackVec[trackLoop].energyWeightedXY[layerLoop].X();
                    std::cout << setw(12) << trackVec[trackLoop].energyWeightedXY[layerLoop].Y();
                    std::cout << setw(12) << trackVec[trackLoop].distsFromHitCentre[layerLoop].X();
                    std::cout << setw(12) << trackVec[trackLoop].distsFromHitCentre[layerLoop].Y();

                    std::cout << setw(12) << trackVec[trackLoop].truthPositionsrPhi[layerLoop].first;
                    std::cout << setw(12) << trackVec[trackLoop].truthPositionsrPhi[layerLoop].second;
                    std::cout << setw(12) << trackVec[trackLoop].energyWeightedrPhi[layerLoop].first;
                    std::cout << setw(12) << trackVec[trackLoop].energyWeightedrPhi[layerLoop].second;
                    std::cout << setw(12) << trackVec[trackLoop].distsFromHitCentrerPhi[layerLoop].first;
                    std::cout << setw(12) << trackVec[trackLoop].distsFromHitCentrerPhi[layerLoop].second;
                    std::cout << setw(12) << layerLoop;
                    
                    std::cout << std::endl;
                }
            }
        }
        tracks_ = trackVec;
    }



