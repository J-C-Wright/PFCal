
#include <string>
#include <iostream>
#include <fstream>
#include "HGCSSGenParticle.hh"
#include "HGCSSRecoHit.hh"
#include "PositionTruth.h"

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

    unsigned int PositionTruthProducer::numberOfPhotons(std::vector<HGCSSGenParticle> *genvec) {
        unsigned int numPhotons(0);
        for (unsigned int genLoop(0);genLoop<(*genvec).size();genLoop++) {
            if ((*genvec)[genLoop].pdgid() == 22) {
                numPhotons++;
            } 
        }
        return numPhotons;
    }

    void PositionTruthProducer::calcTruthPositions(std::vector<HGCSSGenParticle> *genvec, unsigned int eventNumber) {

        bool foundMatch;

        if (!layerZsLoaded_) {
            std::cout << "Error! Z positions of the HGCal layers not loaded!" << std::endl;
            return; 
        }
        std::vector<PositionTruth> positionEntries;
        for (unsigned int genLoop(0);genLoop<(*genvec).size();genLoop++) {
        
            if ((*genvec)[genLoop].pdgid() != 22) {
                std::cout << "Error! Particle trackID " << genLoop << " is not a photon ! pdgid = " << (*genvec)[genLoop].pdgid() << std::endl;
            }else{
                foundMatch = true;

                double x0 = (*genvec)[genLoop].x();
                double y0 = (*genvec)[genLoop].y();
                double z0 = (*genvec)[genLoop].z();

                if (debug_) {
                    std::cout << "Particle " << genLoop << std::endl;
                    std::cout << setw(12) << x0 << setw(12) << y0 << setw(12) << z0 << std::endl; 
                    std::cout << setw(12) << "Truth hit locations" << std::endl;
                }

                PositionTruth entry;
                entry.eventNumber = eventNumber;
                entry.particleNumber = genLoop;
                entry.x0 = x0;
                entry.y0 = y0;
                entry.z0 = z0;
                entry.eta = (*genvec)[genLoop].eta();
                entry.phi = (*genvec)[genLoop].phi();

                for (unsigned int layer(0);layer<layerZPositions_.size();layer++) {
                    double x = x0 + (layerZPositions_[layer]-z0)*((*genvec)[genLoop].px())/((*genvec)[genLoop].pz()); 
                    double y = y0 + (layerZPositions_[layer]-z0)*((*genvec)[genLoop].py())/((*genvec)[genLoop].pz()); 
                    entry.truthPositions.push_back(ROOT::Math::XYPoint(x,y));                    
                    if (debug_) {std::cout << setw(12) << x << setw(12) << y << setw(12) << layerZPositions_[layer] << std::endl;}
                }

                positionTruthEntries_.push_back(entry);
            }
        }
    }

    void PositionTruthProducer::getLayerZPositions(const unsigned versionNumber){

        std::ifstream inputFile;
        std::ostringstream inputFileName;
        inputFileName << "data/zPositions_v"<< versionNumber << ".dat";
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
        }

        inputFile.close();
        layerZsLoaded_ = true;

    }

     void PositionTruthProducer::calcEnergyWeightedXYMax(std::vector<HGCSSRecoHit> *recoHitVec, const HGCSSGeometryConversion & geomConv, const unsigned nSR){

        std::vector<std::vector<HGCSSRecoHit>> hitsByLayer(layerZPositions_.size());
        std::vector<ROOT::Math::XYPoint> energyWeightedXY;

        unsigned int layerNum(0);
        for (unsigned int hitLoop(0);hitLoop<(*recoHitVec).size();hitLoop++) {
            HGCSSRecoHit hit = (*recoHitVec)[hitLoop];
            hitsByLayer[hit.layer()].push_back(hit);
        }

        for (unsigned int layerLoop(0);layerLoop<layerZPositions_.size();layerLoop++) {

            float maxEnergy(0.0);
            unsigned int maxCellIndex(0);
            bool emptyLayer = false;

            if (hitsByLayer[layerLoop].size() == 0) {
                std::cout << "\tLayer " << layerLoop << " doesn't have any hits" << std::endl;
                emptyLayer = true;
            }

            if (!emptyLayer) {

                for (unsigned int cellLoop(0);cellLoop<hitsByLayer[layerLoop].size();cellLoop++) {
                    if (hitsByLayer[layerLoop][cellLoop].energy() > maxEnergy){
                        maxEnergy = hitsByLayer[layerLoop][cellLoop].energy();
                        maxCellIndex = cellLoop;
                    } 
                }

                double radialDisplacement = sqrt(pow(hitsByLayer[layerLoop][maxCellIndex].get_x(),2)+pow(hitsByLayer[layerLoop][maxCellIndex].get_y(),2));
                double step = geomConv.cellSize(layerLoop,radialDisplacement)*sqrt(2)+0.1;

                //dist from cell centre goes here
                std::pair<float,float> dispFromCentre;
                dispFromCentre.first  = (positionTruthEntries_[0].truthPositions[layerLoop].X() 
                                            - hitsByLayer[layerLoop][maxCellIndex].get_x())/geomConv.cellSize(layerLoop,radialDisplacement);
                dispFromCentre.second = (positionTruthEntries_[0].truthPositions[layerLoop].Y() 
                                            - hitsByLayer[layerLoop][maxCellIndex].get_y())/geomConv.cellSize(layerLoop,radialDisplacement);
 
                positionTruthEntries_[0].distsFromHitCentre.push_back(dispFromCentre);

                std::vector<HGCSSRecoHit> retainedHits;
                for (unsigned int cellLoop(0);cellLoop<hitsByLayer[layerLoop].size();cellLoop++) {
                    if (fabs(hitsByLayer[layerLoop][cellLoop].get_x() - hitsByLayer[layerLoop][maxCellIndex].get_x() ) < step && 
                        fabs(hitsByLayer[layerLoop][cellLoop].get_y() - hitsByLayer[layerLoop][maxCellIndex].get_y() ) < step) {
                        retainedHits.push_back(hitsByLayer[layerLoop][cellLoop]);
                    }            
                }

                double energyWeightedX(0.0), energyWeightedY(0.0);
                double totalEnergy(0.0);
                for (unsigned int cellLoop(0);cellLoop<retainedHits.size();cellLoop++) {
                    energyWeightedX += retainedHits[cellLoop].get_x()*retainedHits[cellLoop].energy();
                    energyWeightedY += retainedHits[cellLoop].get_y()*retainedHits[cellLoop].energy();
                    totalEnergy += retainedHits[cellLoop].energy();
                } 

                energyWeightedX /= totalEnergy;
                energyWeightedY /= totalEnergy;
                ROOT::Math::XYPoint point(energyWeightedX,energyWeightedY);
                energyWeightedXY.push_back(point);

            }else{

                ROOT::Math::XYPoint point(0,0);
                energyWeightedXY.push_back(point);

                std::pair<float,float> dispFromCentre;
                dispFromCentre.first = 9999.;
                dispFromCentre.second = 9999.;
                positionTruthEntries_[0].distsFromHitCentre.push_back(dispFromCentre);

            }
        }
        energyWeightedXY_.push_back(energyWeightedXY);
    }


    void PositionTruthProducer::calcEnergyWeightedXYTruth(std::vector<HGCSSRecoHit> *recoHitVec, const HGCSSGeometryConversion & geomConv, const unsigned nSR) {

        if (positionTruthEntries_.size() == 0) {
            std::cout << "Position truth information is unavailable. Probably needs calculating" << std::endl;
            return;
        }

        std::vector<std::vector<HGCSSRecoHit>> hitsByLayer(layerZPositions_.size());
        for (unsigned int hitLoop(0);hitLoop<(*recoHitVec).size();hitLoop++) {
            HGCSSRecoHit hit = (*recoHitVec)[hitLoop];
            hitsByLayer[hit.layer()].push_back(hit);
        }

        for (unsigned int trackLoop(0);trackLoop<positionTruthEntries_.size();trackLoop++) {

            std::vector<ROOT::Math::XYPoint> energyWeightedXY;
            for (unsigned int layerLoop(0);layerLoop<nLayers_;layerLoop++){

                if (hitsByLayer[layerLoop].size() != 0) {
                    
                    float distToTruth(9999.);
                    float maxCellIndex(0);
                    for (unsigned int hitLoop(0);hitLoop<hitsByLayer[layerLoop].size();hitLoop++) {

                        float dx = hitsByLayer[layerLoop][hitLoop].get_x() - positionTruthEntries_[trackLoop].truthPositions[layerLoop].X();        
                        float dy = hitsByLayer[layerLoop][hitLoop].get_y() - positionTruthEntries_[trackLoop].truthPositions[layerLoop].Y();        
                        float dr = sqrt(pow(dx,2)+pow(dy,2));

                        if (dr < distToTruth) {
                            distToTruth = dr;
                            maxCellIndex = hitLoop;
                        }
                    } 

                    double radialDisplacement = sqrt(pow(hitsByLayer[layerLoop][maxCellIndex].get_x(),2)+pow(hitsByLayer[layerLoop][maxCellIndex].get_y(),2));
                    double step = geomConv.cellSize(layerLoop,radialDisplacement)*sqrt(2)+0.1;

                    std::pair<float,float> dispFromCentre;
                    dispFromCentre.first  = (positionTruthEntries_[0].truthPositions[layerLoop].X() 
                                                - hitsByLayer[layerLoop][maxCellIndex].get_x())/geomConv.cellSize(layerLoop,radialDisplacement);
                    dispFromCentre.second = (positionTruthEntries_[0].truthPositions[layerLoop].Y() 
                                                - hitsByLayer[layerLoop][maxCellIndex].get_y())/geomConv.cellSize(layerLoop,radialDisplacement);
                    positionTruthEntries_[trackLoop].distsFromHitCentre.push_back(dispFromCentre);

                    std::vector<HGCSSRecoHit> retainedHits;
                    for (unsigned int hitLoop(0);hitLoop<hitsByLayer[layerLoop].size();hitLoop++) {
                        if (fabs(hitsByLayer[layerLoop][hitLoop].get_x() - hitsByLayer[layerLoop][maxCellIndex].get_x() ) < step && 
                            fabs(hitsByLayer[layerLoop][hitLoop].get_y() - hitsByLayer[layerLoop][maxCellIndex].get_y() ) < step) {
                            retainedHits.push_back(hitsByLayer[layerLoop][hitLoop]);
                        }            
                    }

                    double energyWeightedX(0.0), energyWeightedY(0.0);
                    double totalEnergy(0.0);
                    for (unsigned int cellLoop(0);cellLoop<retainedHits.size();cellLoop++) {
                        energyWeightedX += retainedHits[cellLoop].get_x()*retainedHits[cellLoop].energy();
                        energyWeightedY += retainedHits[cellLoop].get_y()*retainedHits[cellLoop].energy();
                        totalEnergy += retainedHits[cellLoop].energy();
                    } 

                    energyWeightedX /= totalEnergy;
                    energyWeightedY /= totalEnergy;
                    ROOT::Math::XYPoint point(energyWeightedX,energyWeightedY);
                    energyWeightedXY.push_back(point);

                }else{

                    ROOT::Math::XYPoint point(0,0);
                    energyWeightedXY.push_back(point);

                    std::pair<float,float> dispFromCentre;
                    dispFromCentre.first = 9999.;
                    dispFromCentre.second = 9999.;
                    positionTruthEntries_[trackLoop].distsFromHitCentre.push_back(dispFromCentre);

                    std::cout << "\tLayer " << layerLoop << " doesn't have any hits" << std::endl;
                }
            }

            energyWeightedXY_.push_back(energyWeightedXY);

        }    
    }












