
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
        }else{
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
    }

    void TrackTruthProducer::getZpositions(const unsigned versionNumber,
                                           TTree *aSimTree,
                                           const unsigned nEvts,
                                           const unsigned numSiLayers){

        std::vector<HGCSSSimHit> * simhitvec = 0;
        aSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);

        std::ofstream fout;
        std::ostringstream foutname;
        foutname << "data/zPositions_v" << versionNumber << ".dat";
        fout.open(foutname.str());
        if (!fout.is_open()){
            std::cout << " Cannot open outfile " << foutname.str() << " for writing ! Exiting..." << std::endl;
            exit(1);
        }

        std::cout << "--- Filling z positions:" << std::endl
        << "- Processing = " << nEvts  << " events out of " << aSimTree->GetEntries() << std::endl;

        layerZPositions_.resize(nLayers_,0);

        for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
            if (debug_) std::cout << "... Processing entry: " << ievt << " (Layer location)" << std::endl;
            else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << " (Layer location)" << std::endl;
            aSimTree->GetEntry(ievt);
            for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on rechits
                const HGCSSSimHit & lHit = (*simhitvec)[iH];
                unsigned layer = lHit.layer();
                if (layer >= nLayers_) {
                continue;
                }

                //discard some si layers...
                if (lHit.silayer() >= numSiLayers) continue; 

                double posz = lHit.get_z();
                //get z position of hits
                if (layerZPositions_[layer]<posz || layerZPositions_[layer] == 0) layerZPositions_[layer]=posz;
            }
        }

        std::cout << " --- Z positions of layers: " << std::endl;
        for (unsigned iL(0); iL<nLayers_;++iL){
            std::cout << " Layer " << iL << ", z = " << layerZPositions_[iL] << std::endl;
            fout << iL << " " << layerZPositions_[iL] << std::endl;

        }

        fout.close();
        layerZsLoaded_ = true;
    }

    void TrackTruthProducer::produce(std::vector<HGCSSGenParticle> *genvec,
                                        std::vector<HGCSSRecoHit> *recoHitVec,
                                        const HGCSSGeometryConversion & geomConv,
                                        float mipCut, float centralMipCut, float adjacentMipCut){
        
        //Load hit info
        std::vector<std::vector<HGCSSRecoHit>> hitsByLayer_(nLayers_);
        for (unsigned int hitLoop(0);hitLoop<(*recoHitVec).size();hitLoop++) {
            HGCSSRecoHit hit = (*recoHitVec)[hitLoop];
            if (hit.energy() > mipCut) {hitsByLayer_[hit.layer()].push_back(hit);}
        }
        
        //Loop over genParticles
        std::vector<TrackTruth> trackVec((*genvec).size());
        for (unsigned int trackLoop(0);trackLoop<(*genvec).size();trackLoop++) {
            
            //Calculate the trackposition truth
            trackVec[trackLoop].setParticleInfo( (*genvec)[trackLoop] ); 
            std::vector<ROOT::Math::XYPoint> truthPositions;
            std::vector<ROOT::Math::XYPoint> truthDistsFromEdge;
            for (unsigned int layer(0);layer<nLayers_;layer++) {
                double x = trackVec[trackLoop].getParticleInfo().x() 
                            + (layerZPositions_[layer]-trackVec[trackLoop].getParticleInfo().z())*((*genvec)[trackLoop].px())/((*genvec)[trackLoop].pz());
                double y = trackVec[trackLoop].getParticleInfo().y() 
                            + (layerZPositions_[layer]-trackVec[trackLoop].getParticleInfo().z())*((*genvec)[trackLoop].py())/((*genvec)[trackLoop].pz());
                truthPositions.push_back(ROOT::Math::XYPoint(x,y));

                //Dists from edge
                float nearestCentreX = 10.0*round(x/10);
                float nearestCentreY = 10.0*round(y/10);
                float distToEdgeX;
                float distToEdgeY;
                if (x < nearestCentreX) {
                    distToEdgeX = x - (nearestCentreX - 5);
                }else{
                    distToEdgeX = x - (nearestCentreX + 5);
                }
                if (y < nearestCentreY) {
                    distToEdgeY = y - (nearestCentreY - 5);
                }else{
                    distToEdgeY = y - (nearestCentreY + 5);
                }
                truthDistsFromEdge.push_back(ROOT::Math::XYPoint(distToEdgeX,distToEdgeY));

            }
            trackVec[trackLoop].setTruthPositions(truthPositions);
            trackVec[trackLoop].setTruthDistsFromEdge(truthDistsFromEdge);

            //Energy-weighted positions, hits,  and position within cell
            std::vector<std::vector<HGCSSRecoHit>> hitsByLayer3x3(nLayers_);
            std::vector<std::vector<HGCSSRecoHit>> hitsByLayer5x5(nLayers_);
            std::vector<HGCSSRecoHit> centralHitsByLayer(nLayers_);
            std::vector<ROOT::Math::XYPoint> energyWeightedXY(nLayers_);
            std::vector<ROOT::Math::XYPoint> energyWeighted5x5XY(nLayers_);
            std::vector<ROOT::Math::XYPoint> distsFromHitCentre(nLayers_);
            std::vector<ROOT::Math::XYPoint> distsFromTileEdges(nLayers_);
            std::vector<ROOT::Math::XYPoint> distsFromTileEdges5x5(nLayers_);
            std::vector<std::vector<bool>>   adjacentCutsApplied(nLayers_);

            for (unsigned layerLoop(0);layerLoop<nLayers_;layerLoop++) {

                std::vector<HGCSSRecoHit> hit3x3(0);
                std::vector<HGCSSRecoHit> hit5x5(0);
                ROOT::Math::XYPoint energyWeightedPoint(9999,9999);
                ROOT::Math::XYPoint energyWeighted5x5Point(9999,9999);
                ROOT::Math::XYPoint distFromHitCentre(9999,9999);
                ROOT::Math::XYPoint distsFromTileEdge(9999,9999);
                ROOT::Math::XYPoint distsFromTileEdge5x5(9999,9999);
                std::vector<bool> adjacentCuts(4);

                bool mipCutsPass = false;
                float closestCellIndex(0);
                if (hitsByLayer_[layerLoop].size() != 0) {

                    //Find closest hit to truth track
                    float distToTruth(9999.);
                    for (unsigned int hitLoop(0);hitLoop<hitsByLayer_[layerLoop].size();hitLoop++) {
                        float dx = hitsByLayer_[layerLoop][hitLoop].get_x() - trackVec[trackLoop].getTruthPosition(layerLoop).X();
                        float dy = hitsByLayer_[layerLoop][hitLoop].get_y() - trackVec[trackLoop].getTruthPosition(layerLoop).Y();
                        float dr = sqrt(pow(dx,2)+pow(dy,2));
                        if (dr < distToTruth) {
                            distToTruth = dr;
                            closestCellIndex = hitLoop;
                        }
                    }     
                
                    //Does it pass the central Mip cut?
                    if (hitsByLayer_[layerLoop][closestCellIndex].energy() < centralMipCut && debug_) {
                        std::cout << "Layer " << layerLoop << " Central hit fails the central hit mip cut" << std::endl;
                    }
                    mipCutsPass = hitsByLayer_[layerLoop][closestCellIndex].energy() > centralMipCut;


                }
                
                if (mipCutsPass) {
                    //Hit is good, do EW calcs
                    double radialDisplacement = sqrt(pow(hitsByLayer_[layerLoop][closestCellIndex].get_x(),2)+pow(hitsByLayer_[layerLoop][closestCellIndex].get_y(),2));
                    double step = geomConv.cellSize(layerLoop,radialDisplacement)+0.1;

                    //Calculate the distance from the truth to the centre of the closest cell
                    //XY
                    float dx = (trackVec[trackLoop].getTruthPosition(layerLoop).X()
                                                - hitsByLayer_[layerLoop][closestCellIndex].get_x())/geomConv.cellSize(layerLoop,radialDisplacement);
                    float dy = (trackVec[trackLoop].getTruthPosition(layerLoop).Y()
                                                - hitsByLayer_[layerLoop][closestCellIndex].get_y())/geomConv.cellSize(layerLoop,radialDisplacement);
                    distFromHitCentre.SetX(dx);
                    distFromHitCentre.SetY(dy);

                    //Store central hit
                    centralHitsByLayer[layerLoop] = hitsByLayer_[layerLoop][closestCellIndex];

                    //find the hits that belong to the 3x3 grid centred on the closest
                    for (unsigned int hitLoop(0);hitLoop<hitsByLayer_[layerLoop].size();hitLoop++) {
                        if (fabs(hitsByLayer_[layerLoop][hitLoop].get_x() - hitsByLayer_[layerLoop][closestCellIndex].get_x() ) < step &&
                            fabs(hitsByLayer_[layerLoop][hitLoop].get_y() - hitsByLayer_[layerLoop][closestCellIndex].get_y() ) < step) {
                            hit3x3.push_back(hitsByLayer_[layerLoop][hitLoop]);
                        }
                    } 

                    //Find the hits that belong to the 5x5 grid centred on the closest
                    for (unsigned int hitLoop(0);hitLoop<hitsByLayer_[layerLoop].size();hitLoop++) {
                        if (fabs(hitsByLayer_[layerLoop][hitLoop].get_x() - hitsByLayer_[layerLoop][closestCellIndex].get_x() ) < step*2.0 &&
                            fabs(hitsByLayer_[layerLoop][hitLoop].get_y() - hitsByLayer_[layerLoop][closestCellIndex].get_y() ) < step*2.0) {
                            hit5x5.push_back(hitsByLayer_[layerLoop][hitLoop]);
                        }
                    } 

                    //Flag when the adjacent cells are less than the adjacent cell MIP cut
                    adjacentCuts[0] = false;
                    adjacentCuts[1] = false;
                    adjacentCuts[2] = false;
                    adjacentCuts[3] = false;
                    for (unsigned hit(0);hit<hit3x3.size();hit++) {
                        // - B -
                        // A x C
                        // - D -
                        //A
                        bool aFail = hit3x3[hit].get_x()-centralHitsByLayer[layerLoop].get_x()-0.1 < -1.0*geomConv.cellSize(layerLoop,radialDisplacement) &&
                                     hit3x3[hit].get_y() == centralHitsByLayer[layerLoop].get_y() && hit3x3[hit].energy() < adjacentMipCut;
                        //B
                        bool bFail = hit3x3[hit].get_y()-centralHitsByLayer[layerLoop].get_y()+0.1 > geomConv.cellSize(layerLoop,radialDisplacement) &&
                                     hit3x3[hit].get_x() == centralHitsByLayer[layerLoop].get_x() && hit3x3[hit].energy() < adjacentMipCut;
                        //C
                        bool cFail = hit3x3[hit].get_x()-centralHitsByLayer[layerLoop].get_x()+0.1 > geomConv.cellSize(layerLoop,radialDisplacement) &&
                                     hit3x3[hit].get_y() == centralHitsByLayer[layerLoop].get_y() && hit3x3[hit].energy() < adjacentMipCut;
                        //D
                        bool dFail = hit3x3[hit].get_y()-centralHitsByLayer[layerLoop].get_y()-0.1 < -1.0*geomConv.cellSize(layerLoop,radialDisplacement) &&
                                     hit3x3[hit].get_x() == centralHitsByLayer[layerLoop].get_x() && hit3x3[hit].energy() < adjacentMipCut;
                        if (aFail) adjacentCuts[0] = true;
                        if (bFail) adjacentCuts[1] = true;
                        if (cFail) adjacentCuts[2] = true;
                        if (dFail) adjacentCuts[3] = true;
                    }

                    //Calculate energy-weighted position
                    //XY position
                    //3x3
                    double energyWeightedX(0.0), energyWeightedY(0.0);
                    double totalEnergy(0.0);
                    for (unsigned int hitLoop(0);hitLoop<hit3x3.size();hitLoop++) {
                        energyWeightedX += hit3x3[hitLoop].get_x()*hit3x3[hitLoop].energy();
                        energyWeightedY += hit3x3[hitLoop].get_y()*hit3x3[hitLoop].energy();
                        totalEnergy += hit3x3[hitLoop].energy();
                    }
                    energyWeightedPoint.SetX(energyWeightedX/totalEnergy);
                    energyWeightedPoint.SetY(energyWeightedY/totalEnergy);
                    //5x5
                    energyWeightedX = 0.0; energyWeightedY = 0.0;
                    totalEnergy = 0.0;
                    for (unsigned int hitLoop(0);hitLoop<hit5x5.size();hitLoop++) {
                        energyWeightedX += hit5x5[hitLoop].get_x()*hit5x5[hitLoop].energy();
                        energyWeightedY += hit5x5[hitLoop].get_y()*hit5x5[hitLoop].energy();
                        totalEnergy += hit5x5[hitLoop].energy();
                    }
                    energyWeighted5x5Point.SetX(energyWeightedX/totalEnergy);
                    energyWeighted5x5Point.SetY(energyWeightedY/totalEnergy);


                    //Calculate the displacement of the energy-weighted position from the edges
                    //3x3
                    //X
                    float positionEdgeX1 = hitsByLayer_[layerLoop][closestCellIndex].get_x() - geomConv.cellSize(layerLoop,radialDisplacement)/2.0;
                    float positionEdgeX2 = hitsByLayer_[layerLoop][closestCellIndex].get_x() + geomConv.cellSize(layerLoop,radialDisplacement)/2.0;
                    if (fabs(positionEdgeX1 - energyWeightedPoint.X()) < fabs(positionEdgeX2 - energyWeightedPoint.X())) {
                        distsFromTileEdge.SetX(energyWeightedPoint.X()-positionEdgeX1);
                    }else{
                        distsFromTileEdge.SetX(energyWeightedPoint.X()-positionEdgeX2);
                    } 
                    //Y
                    float positionEdgeY1 = hitsByLayer_[layerLoop][closestCellIndex].get_y() - geomConv.cellSize(layerLoop,radialDisplacement)/2.0;
                    float positionEdgeY2 = hitsByLayer_[layerLoop][closestCellIndex].get_y() + geomConv.cellSize(layerLoop,radialDisplacement)/2.0;
                    if (fabs(positionEdgeY1 - energyWeightedPoint.Y()) < fabs(positionEdgeY2 - energyWeightedPoint.Y())) {
                        distsFromTileEdge.SetY(energyWeightedPoint.Y()-positionEdgeY1);
                    }else{
                        distsFromTileEdge.SetY(energyWeightedPoint.Y()-positionEdgeY2);
                    }    

                    //5x5
                    //X
                    if (fabs(positionEdgeX1 - energyWeighted5x5Point.X()) < fabs(positionEdgeX2 - energyWeighted5x5Point.X())) {
                        distsFromTileEdge5x5.SetX(energyWeighted5x5Point.X()-positionEdgeX1);
                    }else{
                        distsFromTileEdge5x5.SetX(energyWeighted5x5Point.X()-positionEdgeX2);
                    } 
                    //Y
                    if (fabs(positionEdgeY1 - energyWeighted5x5Point.Y()) < fabs(positionEdgeY2 - energyWeighted5x5Point.Y())) {
                        distsFromTileEdge5x5.SetY(energyWeighted5x5Point.Y()-positionEdgeY1);
                    }else{
                        distsFromTileEdge5x5.SetY(energyWeighted5x5Point.Y()-positionEdgeY2);
                    }    


                }else {
                    if (debug_) {std::cout << "---- " << "In layer " << layerLoop << " there are no hits above " << mipCut << " MIPs ----" << std::endl;}
                }

                hitsByLayer3x3[layerLoop] = hit3x3;
                hitsByLayer5x5[layerLoop] = hit5x5;
                energyWeightedXY[layerLoop] = energyWeightedPoint;
                energyWeighted5x5XY[layerLoop] = energyWeighted5x5Point;
                distsFromHitCentre[layerLoop] = distFromHitCentre;
                distsFromTileEdges[layerLoop] = distsFromTileEdge;
                distsFromTileEdges5x5[layerLoop] = distsFromTileEdge5x5;
                adjacentCutsApplied[layerLoop] = adjacentCuts;

            }
            //setters
            trackVec[trackLoop].setHitsByLayer3x3(hitsByLayer3x3);
            trackVec[trackLoop].setHitsByLayer5x5(hitsByLayer5x5);
            trackVec[trackLoop].setEnergyWeightedXY(energyWeightedXY);
            trackVec[trackLoop].setEnergyWeighted5x5XY(energyWeighted5x5XY);
            trackVec[trackLoop].setDistsFromHitCentre(distsFromHitCentre);
            trackVec[trackLoop].setDistsFromTileEdges(distsFromTileEdges);
            trackVec[trackLoop].setDistsFromTileEdges5x5(distsFromTileEdges5x5);
            trackVec[trackLoop].setCentralHitsByLayer(centralHitsByLayer);
            trackVec[trackLoop].setAdjacentCutsStatus(adjacentCutsApplied);

            //Print info to screen
            if (debug_) {
                std::cout << "Track " << trackLoop << std::endl;
                std::cout << "Shower start: " << trackVec[trackLoop].getShowerStart() << std::endl;
                std::cout << setw(12) << "x0"  << setw(12) << "y0" << setw(12) << "z0" << setw(12) << "PDG Id" ;
                std::cout << setw(12) << "Eta" << setw(12) << "Phi" << std::endl;
                std::cout << setw(12) << trackVec[trackLoop].getParticleInfo().x();
                std::cout << setw(12) << trackVec[trackLoop].getParticleInfo().y();
                std::cout << setw(12) << trackVec[trackLoop].getParticleInfo().z();
                std::cout << setw(12) << trackVec[trackLoop].getParticleInfo().pdgid();
                std::cout << setw(12) << trackVec[trackLoop].getParticleInfo().eta();
                std::cout << setw(12) << trackVec[trackLoop].getParticleInfo().phi() << std::endl;
                std::cout << setw(12) << "Hit locations over " << nLayers_ << " layers" << std::endl;
                std::cout << setw(24) << "Truth" << setw(24) << "E-Weighted" << setw(24) << "D from centre" << std::endl;
                std::cout << setw(12) << "X" << setw(12) << "Y" << setw(12) << "X" << setw(12) << "Y";
                std::cout << setw(12) << "X" << setw(12) << "Y" << setw(12) << "X" << setw(12) << "Y";
                std::cout << setw(12) << "Num Hits" << setw(12) << "Layer";
                std::cout << std::endl;
                for (unsigned layerLoop(1);layerLoop<nLayers_;layerLoop++) {
                    std::cout << setw(12) << trackVec[trackLoop].getTruthPosition(layerLoop).X();
                    std::cout << setw(12) << trackVec[trackLoop].getTruthPosition(layerLoop).Y();
                    std::cout << setw(12) << trackVec[trackLoop].getEnergyWeightedXYAtLayer(layerLoop).X();
                    std::cout << setw(12) << trackVec[trackLoop].getEnergyWeightedXYAtLayer(layerLoop).Y();
                    std::cout << setw(12) << trackVec[trackLoop].getDistsFromHitCentreAtLayer(layerLoop).X();
                    std::cout << setw(12) << trackVec[trackLoop].getDistsFromHitCentreAtLayer(layerLoop).Y();
                    std::cout << setw(12) << trackVec[trackLoop].getDistsFromTileEdgesAtLayer(layerLoop).X();
                    std::cout << setw(12) << trackVec[trackLoop].getDistsFromTileEdgesAtLayer(layerLoop).Y();
                    std::cout << setw(12) << trackVec[trackLoop].numberOfHitsInLayer3x3(layerLoop);
                    std::cout << setw(12) << layerLoop;
                    std::cout << std::endl;
                }
                std::cout << "Shower starts at layer " << trackVec[trackLoop].getShowerStart() << std::endl;
            }
        }
        tracks_ = trackVec;
    }

    TrackInfo TrackTruthProducer::trackStruct(unsigned index){

        TrackInfo outStruct;        
    
        std::vector<ROOT::Math::XYPoint> truthPositions = tracks_[index].getTruthPositions();
        std::vector<ROOT::Math::XYPoint> energyWeightedXY = tracks_[index].getEnergyWeightedXY();
        std::vector<ROOT::Math::XYPoint> distsFromHitCentre = tracks_[index].getDistsFromHitCentre();
        std::vector<ROOT::Math::XYPoint> distsFromTileEdges = tracks_[index].getDistsFromTileEdges();
        std::vector<HGCSSRecoHit> centralHitsByLayer = tracks_[index].getCentralHitsByLayer();
        std::vector<std::vector<bool>> adjacentCutsStatus = tracks_[index].getAdjacentCutsStatus();

        outStruct.showerStart = tracks_[index].getShowerStart();
        for (unsigned layer(0);layer<28;layer++) {
            outStruct.truthX[layer] = truthPositions[layer].X();
            outStruct.truthY[layer] = truthPositions[layer].Y();
            outStruct.energyWeightedX[layer] = energyWeightedXY[layer].X();
            outStruct.energyWeightedY[layer] = energyWeightedXY[layer].Y();
            outStruct.centralE[layer] = centralHitsByLayer[layer].energy();
            outStruct.totalE[layer] = tracks_[index].totalEnergyOf3x3Hit(layer);
            outStruct.numHitsInLayer[layer] = tracks_[index].numberOfHitsInLayer3x3(layer);
            outStruct.distsFromHitCentreX[layer]  = tracks_[index].getDistsFromHitCentreAtLayer(layer).X();
            outStruct.distsFromHitCentreY[layer]  = tracks_[index].getDistsFromHitCentreAtLayer(layer).Y();
            outStruct.distsFromTileEdgesX[layer]  = tracks_[index].getDistsFromTileEdgesAtLayer(layer).X();
            outStruct.distsFromTileEdgesY[layer]  = tracks_[index].getDistsFromTileEdgesAtLayer(layer).Y();
            outStruct.aAdjacentCut[layer] = adjacentCutsStatus[layer][0];
            outStruct.bAdjacentCut[layer] = adjacentCutsStatus[layer][1];
            outStruct.cAdjacentCut[layer] = adjacentCutsStatus[layer][2];
            outStruct.dAdjacentCut[layer] = adjacentCutsStatus[layer][3];
            outStruct.truthDistsFromEdgeX[layer] = tracks_[index].getTruthDistFromEdge(layer).X();
            outStruct.truthDistsFromEdgeY[layer] = tracks_[index].getTruthDistFromEdge(layer).Y();
            outStruct.energyWeighted5x5X[layer]  = tracks_[index].getEnergyWeighted5x5XYAtLayer(layer).X();
            outStruct.energyWeighted5x5Y[layer]  = tracks_[index].getEnergyWeighted5x5XYAtLayer(layer).Y();
            outStruct.distsFromTileEdges5x5X[layer] = tracks_[index].getDistsFromTileEdges5x5AtLayer(layer).X();
            outStruct.distsFromTileEdges5x5Y[layer] = tracks_[index].getDistsFromTileEdges5x5AtLayer(layer).Y();
        }
        
        return outStruct;
    }

