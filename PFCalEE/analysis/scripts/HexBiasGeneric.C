#include <vector>

void HexBiasGeneric() {

    gROOT->ProcessLine(".L scripts/effSigmaMacro.C");

    float eRatioCut(1.0);
   
    TFile *file = TFile::Open("out_V100_100k.root"); 
    TFile *outputFile = new TFile("Plots/HexBiases.root","RECREATE");
    outputFile->cd();

    TTree *tree = (TTree*)file->Get("tracks");

    //Bias curves
    TH2F* layerBiasesXEdge[28];
    TH2F* layerBiasesXCentre[28];
    TH2F* layerBiasesUEdge[28];
    TH2F* layerBiasesUCentre[28];
    TH2F* layerBiasesVEdge[28];
    TH2F* layerBiasesVCentre[28];
    TH2F* layerBiasesWEdge[28];
    TH2F* layerBiasesWCentre[28];

    unsigned layerNum(28);
    unsigned bins(200);
 
    for (unsigned layer(0);layer<layerNum;layer++) {
        
        std::cout << "Processing layer " << layer << std::endl;

        float minX(-6.4),maxX(6.4);
        float minY(-6.4),maxY(6.4);
        //Bias curve
        //X
        layerBiasesXEdge[layer] = getBiasCurve(layer,bins,tree,1,3,7);
        layerBiasesXEdge[layer]->Write();
        layerBiasesXCentre[layer] = getBiasCurve(layer,bins,tree,1,3,5);
        layerBiasesXCentre[layer]->Write();
        //UVW
        layerBiasesUCentre[layer] = getBiasCurve(layer,bins,tree,12,15,18);
        layerBiasesUCentre[layer]->Write();
        layerBiasesVCentre[layer] = getBiasCurve(layer,bins,tree,13,16,19);
        layerBiasesVCentre[layer]->Write();
        layerBiasesWCentre[layer] = getBiasCurve(layer,bins,tree,14,17,20);
        layerBiasesWCentre[layer]->Write();

    }

    outputFile->Close();
}
TH2F* getBiasCurve(unsigned layer, unsigned numBins, TTree* tree, 
                   unsigned ewCoord, unsigned truthCoord, unsigned edgeCoord
                   ) {

    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 
    std::ostringstream name;
    name << "Layer"<< layer << "_BiasVsEdgeX_" << ewCoord << "_" << truthCoord << "_" << edgeCoord;
    TString histName = name.str(); 
    TH2F *biasCurve = new TH2F(histName,histName,numBins,-6.4,6.4,numBins,-6.4,-6.4); 

    unsigned skipCount(0);

    for (unsigned event(0);event<tree->GetEntries();event++) {
        tree->GetEntry(event);
        int start  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue();
        if (start + layer > 27) continue;

        //Within hex?
        /*
        float u = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(18)->GetName())->GetValue(start + layer);
        float v = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(19)->GetName())->GetValue(start + layer);
        float w = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(20)->GetName())->GetValue(start + layer);
        if (fabs(u) > 5.54 || fabs(v) > 5.54 || fabs(w) > 5.54) {
            skipCount++; 
            continue;
        }
        */

        float x_ew = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(ewCoord)->GetName())->GetValue(start + layer);
        float x_t  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(truthCoord)->GetName())->GetValue(start + layer);
        float x_e  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(edgeCoord)->GetName())->GetValue(start + layer);
        if (x_ew < 9999 && fabs(x_e) < 6.4 && fabs(x_ew - x_t) < 6.4 ) {
            biasCurve->Fill(x_e,x_ew - x_t);
        }
    }
    biasCurve->SetOption("COLZ");
    biasCurve->SetStats(kFALSE);

    //std::cout << "There were " << skipCount << " events outside the central hexagon" << std::endl;

    return biasCurve;
} 

struct TrackInfo {

    unsigned showerStart;               //0
    Float_t energyWeightedX[28];        //1
    Float_t energyWeightedY[28];        //2
    Float_t truthX[28];                 //3
    Float_t truthY[28];                 //4
    Float_t distsFromHitCentreX[28];    //5
    Float_t distsFromHitCentreY[28];    //6
    Float_t distsFromTileEdgesX[28];    //7
    Float_t distsFromTileEdgesY[28];    //8
    Float_t centralE[28];               //9
    Float_t totalE[28];                 //10
    UInt_t  numHitsInLayer[28];         //11
    //UVWs
    Float_t energyWeightedU[28];        //12
    Float_t energyWeightedV[28];        //13
    Float_t energyWeightedW[28];        //14
    Float_t truthU[28];                 //15
    Float_t truthV[28];                 //16
    Float_t truthW[28];                 //17
    Float_t distsFromHitCentreU[28];    //18
    Float_t distsFromHitCentreV[28];    //19
    Float_t distsFromHitCentreW[28];    //20
    Float_t distsFromHitCentreUPerp[28];//21
    Float_t distsFromHitCentreVPerp[28];//22
    Float_t distsFromHitCentreWPerp[28];//23

};

