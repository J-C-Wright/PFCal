#include <vector>

void HexBiases() {

    gROOT->ProcessLine(".L scripts/effSigmaMacro.C");

    float eRatioCut(1.0);
   
    TFile *file = TFile::Open("out_V100_100k.root"); 
    TFile *outputFile = new TFile("Plots/HexBiases.root","RECREATE");
    outputFile->cd();

    TTree *tree = (TTree*)file->Get("tracks");

    //Bias curves
    TH2F* layerBiasesXEdge[28];
    TH2F* layerBiasesXCentre[28];
    TH2F* layerBiasesYEdge[28];
    TH2F* layerBiasesYCentre[28];

    unsigned layerNum(28);
    unsigned bins(200);
 
    for (unsigned layer(0);layer<layerNum;layer++) {
        
        std::cout << "Processing layer " << layer << std::endl;

        float minX(-6.4),maxX(6.4);
        float minY(-6.4),maxY(6.4);
        //Bias curve
        //X
        layerBiasesXEdge[layer] = getEdgeBiasCurveX(layer,bins,tree,minX,maxX,minY,maxY);
        layerBiasesXEdge[layer]->Write();
        layerBiasesXCentre[layer] = getCentreBiasCurveX(layer,bins,tree,minX,maxX,minY,maxY);
        layerBiasesXCentre[layer]->Write();
        //Y
        layerBiasesYEdge[layer] = getEdgeBiasCurveY(layer,bins,tree,minX,maxX,minY,maxY);
        layerBiasesYEdge[layer]->Write();
        layerBiasesYCentre[layer] = getCentreBiasCurveY(layer,bins,tree,minX,maxX,minY,maxY);
        layerBiasesYCentre[layer]->Write();
    }

    outputFile->Close();
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

TH2F* getEdgeBiasCurveX(unsigned layer, unsigned numBins, TTree* tree, float minX, float maxX, float minY, float maxY) {

    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 
    std::ostringstream name;
    name << "Layer"<< layer << "_BiasVsEdgeX";
    TString histName = name.str(); 
    TH2F *biasCurve = new TH2F(histName,histName,numBins,-6.4,6.4,numBins,-6.4,-6.4); 

    for (unsigned event(0);event<tree->GetEntries();event++) {
        tree->GetEntry(event);
        int start  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue();
        if (start + layer > 27) continue;

        float x_c  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(5)->GetName())->GetValue(start + layer);
        float y_c  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(6)->GetName())->GetValue(start + layer);
        if (x_c > maxX || x_c < minX) continue;
        if (y_c > maxY || y_c < minY) continue;

        float x_ew = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(1)->GetName())->GetValue(start + layer);
        float x_t  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(3)->GetName())->GetValue(start + layer);
        float x_e  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(7)->GetName())->GetValue(start + layer);
        if (x_ew < 9999 && fabs(x_e) < 6.4 && fabs(x_ew - x_t) < 6.4 ) {
            biasCurve->Fill(x_e,x_ew - x_t);
        }
    }
    biasCurve->SetOption("COLZ");
    biasCurve->SetStats(kFALSE);

    return biasCurve;
} 

TH2F* getEdgeBiasCurveY(unsigned layer, unsigned numBins, TTree* tree, float minX, float maxX, float minY, float maxY) {

    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 
    std::ostringstream name;
    name << "Layer"<< layer << "_BiasVsEdgeY";
    TString histName = name.str(); 
    TH2F *biasCurve = new TH2F(histName,histName,numBins,-6.4,6.4,numBins,-6.4,-6.4); 

    for (unsigned event(0);event<tree->GetEntries();event++) {
        tree->GetEntry(event);
        int start  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue();
        if (start + layer > 27) continue;

        float x_c  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(5)->GetName())->GetValue(start + layer);
        float y_c  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(6)->GetName())->GetValue(start + layer);
        if (x_c > maxX || x_c < minX) continue;
        if (y_c > maxY || y_c < minY) continue;

        float y_ew = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(2)->GetName())->GetValue(start + layer);
        float y_t  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(4)->GetName())->GetValue(start + layer);
        float y_e  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(8)->GetName())->GetValue(start + layer);
        if (y_ew < 9999 && fabs(y_e) < 6.4 && fabs(y_ew - y_t) < 6.4 ) {
            biasCurve->Fill(y_e,y_ew - y_t);
        }
    }
    biasCurve->SetOption("COLZ");
    biasCurve->SetStats(kFALSE);

    return biasCurve;
} 

TH2F* getCentreBiasCurveX(unsigned layer, unsigned numBins, TTree* tree, float minX, float maxX, float minY, float maxY) {

    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 
    std::ostringstream name;
    name << "Layer"<< layer << "_BiasVsCentreX";
    TString histName = name.str(); 
    TH2F *biasCurve = new TH2F(histName,histName,numBins,-6.4,6.4,numBins,-6.4,-6.4); 

    for (unsigned event(0);event<tree->GetEntries();event++) {
        tree->GetEntry(event);
        int start  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue();
        if (start + layer > 27) continue;

        float x_c  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(5)->GetName())->GetValue(start + layer);
        float y_c  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(6)->GetName())->GetValue(start + layer);
        if (x_c > maxX || x_c < minX) continue;
        if (y_c > maxY || y_c < minY) continue;

        float x_ew = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(1)->GetName())->GetValue(start + layer);
        float x_t  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(3)->GetName())->GetValue(start + layer);
        if (x_ew < 9999 && fabs(x_c) < 6.4 && fabs(x_ew - x_t) < 6.4 ) {
            biasCurve->Fill(x_c,x_ew - x_t);
        }
    }
    biasCurve->SetOption("COLZ");
    biasCurve->SetStats(kFALSE);

    return biasCurve;
}

TH2F* getCentreBiasCurveY(unsigned layer, unsigned numBins, TTree* tree, float minX, float maxX, float minY, float maxY) {

    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 
    std::ostringstream name;
    name << "Layer"<< layer << "_BiasVsCentreY";
    TString histName = name.str(); 
    TH2F *biasCurve = new TH2F(histName,histName,numBins,-6.4,6.4,numBins,-6.4,-6.4); 

    for (unsigned event(0);event<tree->GetEntries();event++) {
        tree->GetEntry(event);
        int start  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue();
        if (start + layer > 27) continue;

        float x_c  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(5)->GetName())->GetValue(start + layer);
        float y_c  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(6)->GetName())->GetValue(start + layer);
        if (x_c > maxX || x_c < minX) continue;
        if (y_c > maxY || y_c < minY) continue;

        float y_ew = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(2)->GetName())->GetValue(start + layer);
        float y_t  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(4)->GetName())->GetValue(start + layer);
        if (y_ew < 9999 && fabs(y_c) < 6.4 && fabs(y_ew - y_t) < 6.4 ) {
            biasCurve->Fill(y_c,y_ew - y_t);
        }
    }
    biasCurve->SetOption("COLZ");
    biasCurve->SetStats(kFALSE);

    return biasCurve;
}


