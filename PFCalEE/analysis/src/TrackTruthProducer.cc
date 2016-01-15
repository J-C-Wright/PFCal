
#define _USE_MATH_DEFINES

#include <string>
#include <iostream>
#include <fstream>
#include <utility>
#include "HGCSSGenParticle.hh"
#include "HGCSSRecoHit.hh"
#include "TrackTruthProducer.hh"
//#include "HexagonalGeometry.hh"

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
                                        int mipCut, int centralMipCut){
        
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
            std::vector<UVWPoint> truthPositionsUVW;

            for (unsigned int layer(0);layer<nLayers_;layer++) {
                double x = trackVec[trackLoop].getParticleInfo().x() 
                            + (layerZPositions_[layer]-trackVec[trackLoop].getParticleInfo().z())*((*genvec)[trackLoop].px())/((*genvec)[trackLoop].pz());
                double y = trackVec[trackLoop].getParticleInfo().y() 
                            + (layerZPositions_[layer]-trackVec[trackLoop].getParticleInfo().z())*((*genvec)[trackLoop].py())/((*genvec)[trackLoop].pz());

                truthPositions.push_back(ROOT::Math::XYPoint(x,y));
                UVWPoint uvwPoint(ROOT::Math::XYPoint(x,y));
                truthPositionsUVW.push_back(uvwPoint);

            }

            trackVec[trackLoop].setTruthPositions(truthPositions);
            trackVec[trackLoop].setTruthPositionsUVW(truthPositionsUVW);

            //Energy-weighted positions, hits,  and position within cell
            std::vector<std::vector<HGCSSRecoHit>> hitsByLayer3x3(nLayers_);
            std::vector<HGCSSRecoHit> centralHitsByLayer(nLayers_);

            //XY
            std::vector<ROOT::Math::XYPoint> energyWeightedXY(nLayers_);
            std::vector<ROOT::Math::XYPoint> distsFromHitCentre(nLayers_);
            std::vector<ROOT::Math::XYPoint> distsFromTileEdges(nLayers_);
            //UVW
            std::vector<UVWPoint> energyWeightedUVW(nLayers_);
            std::vector<UVWPoint> distsFromHitCentreUVW(nLayers_);
            std::vector<std::vector<float>> distsFromTileEdgesUVW(nLayers_);

            for (unsigned layerLoop(0);layerLoop<nLayers_;layerLoop++) {

                std::vector<HGCSSRecoHit> hit3x3(0);
                ROOT::Math::XYPoint dummy(9999,9999);
                //XY
                ROOT::Math::XYPoint energyWeightedPoint = dummy;
                ROOT::Math::XYPoint distFromHitCentre   = dummy;
                ROOT::Math::XYPoint distFromTileEdge   = dummy;
                //UVW
                UVWPoint energyWeightedPointUVW(dummy);
                UVWPoint distFromHitCentreUVW(dummy);
                std::vector<float> distFromTileEdgeUVW(3,9999);


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
                    mipCutsPass = hitsByLayer_[layerLoop][closestCellIndex].energy() > centralMipCut;
                    if (!mipCutsPass && debug_) {
                        std::cout << "Layer " << layerLoop << " Central hit fails the central hit mip cut" << std::endl;
                    }

                }
                
                if (mipCutsPass) {
                    //Hit is good, do EW calcs
                    double radialDisplacement = sqrt(pow(hitsByLayer_[layerLoop][closestCellIndex].get_x(),2)+pow(hitsByLayer_[layerLoop][closestCellIndex].get_y(),2));
                    double step = sqrt(3.0)*geomConv.cellSize(layerLoop,radialDisplacement)+0.1;

                    //Calculate the distance from the truth to the centre of the closest cell
                    //find the hits that belong to the heptad centred on the closest
                    for (unsigned int hitLoop(0);hitLoop<hitsByLayer_[layerLoop].size();hitLoop++) {
                        if (fabs(hitsByLayer_[layerLoop][hitLoop].get_x() - hitsByLayer_[layerLoop][closestCellIndex].get_x() ) < step &&
                            fabs(hitsByLayer_[layerLoop][hitLoop].get_y() - hitsByLayer_[layerLoop][closestCellIndex].get_y() ) < step) {
                            hit3x3.push_back(hitsByLayer_[layerLoop][hitLoop]);
                        }
                    } 
                    centralHitsByLayer[layerLoop] = hitsByLayer_[layerLoop][closestCellIndex];

                    //Calculate energy-weighted position
                    //XY position
                    double energyWeightedX(0.0), energyWeightedY(0.0);
                    double totalEnergy(0.0);
                    for (unsigned int hitLoop(0);hitLoop<hit3x3.size();hitLoop++) {
                        energyWeightedX += hit3x3[hitLoop].get_x()*hit3x3[hitLoop].energy();
                        energyWeightedY += hit3x3[hitLoop].get_y()*hit3x3[hitLoop].energy();
                        totalEnergy += hit3x3[hitLoop].energy();
                    }
                    energyWeightedPoint.SetX(energyWeightedX/totalEnergy);
                    energyWeightedPoint.SetY(energyWeightedY/totalEnergy);

                    //Make hexagon object for hex geometry stuff
                    ROOT::Math::XYPoint centre(hitsByLayer_[layerLoop][closestCellIndex].get_x(),hitsByLayer_[layerLoop][closestCellIndex].get_y());
                    Hexagon tile(geomConv.cellSize(layerLoop,radialDisplacement), true, centre);

                    //Set energy-weighted position from centre
                    std::pair<float,float> displacementFromCentre = tile.getDisplacementFromCentre_XY( energyWeightedPoint );
                    distFromHitCentre.SetX(displacementFromCentre.first);
                    distFromHitCentre.SetY(displacementFromCentre.second);

                    //Set energy-weighted position from edge
                    std::pair<float,float> distanceToEdges_EWP = tile.getDistanceToEdges_XY(energyWeightedPoint); 
                    distFromTileEdge.SetX(distanceToEdges_EWP.first);
                    distFromTileEdge.SetY(distanceToEdges_EWP.second);
                    
                    //UVW Conversions
                    energyWeightedPointUVW.setUVW(energyWeightedPoint);
                    distFromHitCentreUVW.setUVW(distFromHitCentre);
                    distFromTileEdgeUVW = tile.getDistanceToEdges_UVW(distFromHitCentre);

                }else {
                    if (debug_) {std::cout << "---- " << "In layer " << layerLoop << " there are no hits above " << mipCut << " MIPs ----" << std::endl;}
                }

                hitsByLayer3x3[layerLoop] = hit3x3;
                //XY
                energyWeightedXY[layerLoop]   = energyWeightedPoint;
                distsFromHitCentre[layerLoop] = distFromHitCentre;
                distsFromTileEdges[layerLoop] = distFromTileEdge;
                //UVW
                energyWeightedUVW[layerLoop]     = energyWeightedPointUVW;
                distsFromHitCentreUVW[layerLoop] = distFromHitCentreUVW;
                distsFromTileEdgesUVW[layerLoop] = distFromTileEdgeUVW; 

            }
            //setters
            trackVec[trackLoop].setHitsByLayer(hitsByLayer3x3);
            trackVec[trackLoop].setCentralHitsByLayer(centralHitsByLayer);
            //XY
            trackVec[trackLoop].setEnergyWeightedXY(energyWeightedXY);
            trackVec[trackLoop].setDistsFromHitCentre(distsFromHitCentre);
            trackVec[trackLoop].setDistsFromTileEdges(distsFromTileEdges);
            //UVW
            trackVec[trackLoop].setEnergyWeightedUVW(energyWeightedUVW);
            trackVec[trackLoop].setDistsFromHitCentreUVW(distsFromHitCentreUVW);
            trackVec[trackLoop].setDistsFromTileEdgesUVW(distsFromTileEdgesUVW);

            //Print info to screen
            if (debug_) {
                std::cout << "Track " << trackLoop << std::endl;
                std::cout << "Shower start: " << trackVec[trackLoop].getShowerStart(0.1) << std::endl;
                std::cout << setw(12) << "x0"  << setw(12) << "y0" << setw(12) << "z0" << setw(12) << "PDG Id" ;
                std::cout << setw(12) << "Eta" << setw(12) << "Phi" << std::endl;
                std::cout << setw(12) << trackVec[trackLoop].getParticleInfo().x();
                std::cout << setw(12) << trackVec[trackLoop].getParticleInfo().y();
                std::cout << setw(12) << trackVec[trackLoop].getParticleInfo().z();
                std::cout << setw(12) << trackVec[trackLoop].getParticleInfo().pdgid();
                std::cout << setw(12) << trackVec[trackLoop].getParticleInfo().eta();
                std::cout << setw(12) << trackVec[trackLoop].getParticleInfo().phi() << std::endl;
                std::cout << setw(12) << "Hit locations over " << nLayers_ << " layers" << std::endl;
                std::cout << setw(24) << "Truth" << setw(24) << "E-Weighted" << setw(24) << "D from centre" << setw(24) << "D from edge" << std::endl;
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
                    std::cout << setw(12) << trackVec[trackLoop].numberOfHitsInLayer(layerLoop);
                    std::cout << setw(12) << layerLoop;
                    std::cout << std::endl;
                }
                std::cout << "Shower starts at layer " << trackVec[trackLoop].getShowerStart(0.1) << std::endl;
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

        outStruct.showerStart = tracks_[index].getShowerStart(0.1);
        for (unsigned layer(0);layer<28;layer++) {
            outStruct.truthX[layer] = truthPositions[layer].X();
            outStruct.truthY[layer] = truthPositions[layer].Y();
            outStruct.energyWeightedX[layer] = energyWeightedXY[layer].X();
            outStruct.energyWeightedY[layer] = energyWeightedXY[layer].Y();
            outStruct.centralE[layer] = centralHitsByLayer[layer].energy();
            outStruct.totalE[layer] = tracks_[index].totalEnergyOf3x3Hit(layer);
            outStruct.numHitsInLayer[layer] = tracks_[index].numberOfHitsInLayer(layer);
            outStruct.distsFromHitCentreX[layer]  = tracks_[index].getDistsFromHitCentreAtLayer(layer).X();
            outStruct.distsFromHitCentreY[layer]  = tracks_[index].getDistsFromHitCentreAtLayer(layer).Y();
            outStruct.distsFromTileEdgesX[layer]  = tracks_[index].getDistsFromTileEdgesAtLayer(layer).X();
            outStruct.distsFromTileEdgesY[layer]  = tracks_[index].getDistsFromTileEdgesAtLayer(layer).Y();
            //UVW
            outStruct.energyWeightedU[layer] = tracks_[index].getEnergyWeightedUVWAtLayer(layer).getU();
            outStruct.energyWeightedV[layer] = tracks_[index].getEnergyWeightedUVWAtLayer(layer).getV();
            outStruct.energyWeightedW[layer] = tracks_[index].getEnergyWeightedUVWAtLayer(layer).getW();
            outStruct.truthU[layer] = tracks_[index].getTruthPositionUVWAtLayer(layer).getU();
            outStruct.truthV[layer] = tracks_[index].getTruthPositionUVWAtLayer(layer).getV();
            outStruct.truthW[layer] = tracks_[index].getTruthPositionUVWAtLayer(layer).getW();
            outStruct.distsFromHitCentreU[layer] = tracks_[index].getDistsFromHitCentreUVWAtLayer(layer).getU();
            outStruct.distsFromHitCentreV[layer] = tracks_[index].getDistsFromHitCentreUVWAtLayer(layer).getV();
            outStruct.distsFromHitCentreW[layer] = tracks_[index].getDistsFromHitCentreUVWAtLayer(layer).getW();
            outStruct.distsFromHitCentreUPerp[layer] = tracks_[index].getDistsFromHitCentreUVWAtLayer(layer).getUPerp();
            outStruct.distsFromHitCentreVPerp[layer] = tracks_[index].getDistsFromHitCentreUVWAtLayer(layer).getVPerp();
            outStruct.distsFromHitCentreWPerp[layer] = tracks_[index].getDistsFromHitCentreUVWAtLayer(layer).getWPerp();
        }
        
        return outStruct;
    }
/*
    UVWPoint getTruthPositionUVWAtLayer(unsigned layer) {return truthPositionsUVW_[layer];}
    UVWPoint getEnergyWeightedUVWAtLayer(unsigned layer) {return energyWeightedUVW_[layer];}
    UVWPoint getDistsFromHitCentreUVWAtLayer(unsigned layer) {return distsFromHitCentreUVW_[layer];}
    std::vector<float> getDistsFromTileEdgesUVWAtLayer(unsigned layer) {return distsFromTileEdgesUVW_[layer];}

    Float_t energyWeightedU[28];
    Float_t energyWeightedV[28];
    Float_t energyWeightedW[28];
    Float_t truthU[28];
    Float_t truthV[28];
    Float_t truthW[28];
    Float_t distsFromHitCentreU[28];
    Float_t distsFromHitCentreV[28];
    Float_t distsFromHitCentreW[28];
    Float_t distsFromHitCentreUPerp[28];
    Float_t distsFromHitCentreVPerp[28];
    Float_t distsFromHitCentreWPerp[28];
*/


