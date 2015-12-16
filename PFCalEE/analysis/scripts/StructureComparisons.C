
struct TrackInfo {
    int showerStart;                    //0
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
};

void StructureComparisons() {

    TCanvas C1("c1");

    TFile *file100 = TFile::Open("out_V100.root");
    TFile *file101 = TFile::Open("out_V101.root");
    TFile *file102 = TFile::Open("out_V102.root");
    TFile *file103 = TFile::Open("out_V103.root");

    TFile *outputFile = new TFile("Plots/StructurePlots.root","RECREATE");
    outputFile->cd();

    TTree *tree100 = (TTree*)file100->Get("tracks");
    TTree *tree101 = (TTree*)file101->Get("tracks");
    TTree *tree102 = (TTree*)file102->Get("tracks");
    TTree *tree103 = (TTree*)file103->Get("tracks");

    TGraph* totalE_100 = getMeanTotalEVsLayer(tree100);
    totalE_100->SetLineColor(1);
    TGraph* totalE_101 = getMeanTotalEVsLayer(tree101);
    totalE_101->SetLineColor(2);
    TGraph* totalE_102 = getMeanTotalEVsLayer(tree102);
    totalE_102->SetLineColor(3);
    TGraph* totalE_103 = getMeanTotalEVsLayer(tree103);
    totalE_103->SetLineColor(4);
    totalE_101->Draw();
    totalE_100->Draw("same");
    totalE_102->Draw("same");
    totalE_103->Draw("same");
    C1.Print("Plots/TotalEnergyByLayer.pdf");
    C1.Clear();

    TGraph* eRatio_100 = getMeanERatioVsLayer(tree100);
    eRatio_100->SetLineColor(1);
    TGraph* eRatio_101 = getMeanERatioVsLayer(tree101);
    eRatio_101->SetLineColor(2);
    TGraph* eRatio_102 = getMeanERatioVsLayer(tree102);
    eRatio_102->SetLineColor(3);
    TGraph* eRatio_103 = getMeanERatioVsLayer(tree103);
    eRatio_103->SetLineColor(4);
    eRatio_101->Draw();
    eRatio_100->Draw("same");
    eRatio_102->Draw("same");
    eRatio_103->Draw("same");
    C1.Print("Plots/ERatioByLayer.pdf");
    C1.Clear();

    TGraph* drift_100 = getDriftVsLayer(tree100);
    drift_100->SetLineColor(1);
    TGraph* drift_101 = getDriftVsLayer(tree101);
    drift_101->SetLineColor(2);
    TGraph* drift_102 = getDriftVsLayer(tree102);
    drift_102->SetLineColor(3);
    TGraph* drift_103 = getDriftVsLayer(tree103);
    drift_103->SetLineColor(4);
    drift_103->Draw();
    drift_101->Draw("same");
    drift_100->Draw("same");
    drift_102->Draw("same");
    C1.Print("Plots/DriftByLayer.pdf");
    C1.Clear();
}


TGraph* getMeanTotalEVsLayer(TTree* tree) {
    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 
    Double_t totalE[28];
    Double_t layerN[28];
    for (unsigned layer(0);layer<28;layer++) {
        unsigned count(0);
        layerN[layer] = layer;
        float tempTotal(0.0);
        for (unsigned event(0);event<tree->GetEntries();event++) {
            tree->GetEntry(event);
            tempTotal += tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(10)->GetName())->GetValue(layer);
            count++;
        }
        totalE[layer] = tempTotal/(float)count;
    }
    TGraph * graph = new TGraph(28,layerN,totalE);
    graph->SetLineWidth(2);
    graph->SetTitle("Total Energy by Layer");
    graph->GetXaxis()->SetTitle("Layer");
    graph->GetYaxis()->SetTitle("Total Energy in 3x3");
    return graph;    
}

TGraph* getMeanERatioVsLayer(TTree* tree) {
    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 
    Double_t eRatio[28];
    Double_t layerN[28];
    for (unsigned layer(0);layer<28;layer++) {
        unsigned count(0);
        layerN[layer] = layer;
        float sumRatio(0.0);
        for (unsigned event(0);event<tree->GetEntries();event++) {
            tree->GetEntry(event);
            float centralE = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(9)->GetName())->GetValue(layer);
            float totalE   = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(10)->GetName())->GetValue(layer);
            if (totalE == 0) continue;
            sumRatio += centralE/totalE;
            count++;
        }
        eRatio[layer] = sumRatio/(float)count;
    }
    TGraph * graph = new TGraph(28,layerN,eRatio);
    graph->SetLineWidth(2);
    graph->SetTitle("Ratio of central energy to total by Layer");
    graph->GetXaxis()->SetTitle("Layer");
    graph->GetYaxis()->SetTitle("E_{centre}/E_{total}");
    return graph;    
}

TGraph* getDriftVsLayer(TTree* tree) {
    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 
    Double_t bias[28];
    Double_t layerN[28];
    for (unsigned layer(0);layer<28;layer++) {
        unsigned count(0);
        layerN[layer] = layer;
        float sumBias(0.0);
        for (unsigned event(0);event<tree->GetEntries();event++) {
            tree->GetEntry(event);
            float y_ew = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(2)->GetName())->GetValue(layer);
            float y_t  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(4)->GetName())->GetValue(layer);
            if (y_ew < 9999 && fabs(y_ew - y_t) < 5){sumBias += y_ew - y_t;}
            count++;
        }
        bias[layer] = sumBias/(float)count;
    }
    TGraph * graph = new TGraph(28,layerN,bias);
    graph->SetLineWidth(2);
    graph->SetTitle("Mean bias by Layer");
    graph->GetXaxis()->SetTitle("Layer");
    graph->GetYaxis()->SetTitle("Mean bias");
    return graph;    
}









/*

TH2F* getBiasCurve(unsigned layer, unsigned numBins, TTree* tree, float eRatioCut) {

    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 
    std::ostringstream name;
    name << "Layer"<< layer << "_Bias";
    TString histName = name.str(); 
    TH2F *biasCurve = new TH2F(histName,histName,numBins,-5,5,numBins,-5,-5); 

    for (unsigned event(0);event<tree->GetEntries();event++) {
        tree->GetEntry(event);
        int start  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue();
        if (start + layer > 27) continue;
        float x_ew = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(1)->GetName())->GetValue(start + layer);
        float x_t  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(3)->GetName())->GetValue(start + layer);
        float x_e  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(7)->GetName())->GetValue(start + layer);
        float centralE = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(9)->GetName())->GetValue(start + layer);
        float totalE   = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(10)->GetName())->GetValue(start + layer);
        if (x_ew < 9999 && fabs(x_e) < 5 && fabs(x_ew - x_t) < 5 && centralE/totalE < eRatioCut) { 
            biasCurve->Fill(x_e,x_ew - x_t);
        }
    }
    biasCurve->SetOption("COLZ");
    biasCurve->SetStats(kFALSE);

    return biasCurve;

}

*/


