#include <vector>

void Slices() {

    float eRatioCut(1.0);
   
    TFile *file = TFile::Open("out_V100.root"); 
    TFile *outputFile = new TFile("Plots/AllSlices.root","RECREATE");
    outputFile->cd();

    TrackInfo test;

    TTree *tree = (TTree*)file->Get("tracks");
    tree->SetBranchAddress("truthInfo",&test.showerStart);
    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 

    unsigned layerNum(28);

    std::vector<std::vector<TH1F*>> slicesByLayer(layerNum); 
    for (unsigned layer(0);layer<layerNum;layer++) {
        
        std::cout << "Processing layer " << layer << std::endl;
        unsigned const numSlices = 30, numBins = 50;
        std::vector<TH1F*> slices =  getBiasSlices(layer,numSlices,numBins,tree,eRatioCut);
        slicesByLayer[layer] = slices;

    }
   
    for (unsigned layer(0);layer<slicesByLayer.size();layer++) {
        for (unsigned slice(0);slice<slicesByLayer[0].size();slice++) {
            slicesByLayer[layer][slice]->Write(); 
        }
    }    

    TCanvas c1("c1");

    for (unsigned int layer(0);layer<layerNum;layer++) {
        c1.Clear();
        std::ostringstream name;
        name << "Plots/gifs/LayerSlices/" << "Layer"<< layer << "_Slices.gif+20";
        TString fileName = name.str(); 
        for (unsigned int slice(0);slice<slicesByLayer[0].size();slice++) {
            slicesByLayer[layer][slice]->Draw(); 
            c1.Print(fileName);
        }
        name.str("");
        name.clear();
        name << "Plots/gifs/LayerSlices/" << "Layer"<< layer << "_Slices.gif++";
        fileName = name.str(); 
        c1.Print(fileName);
    }

    for (unsigned int slice(0);slice<slicesByLayer[0].size();slice++) {
        c1.Clear();
        std::ostringstream name;
        name << "Plots/gifs/SliceLayers/" << "Slice"<< slice << "_Layers.gif+20";
        TString fileName = name.str(); 
        for (unsigned int layer(0);layer<layerNum;layer++) {
            slicesByLayer[layer][slice]->Draw(); 
            c1.Print(fileName);
        } 
        name.str("");
        name.clear();
        name << "Plots/gifs/SliceLayers/" << "Slice"<< slice << "_Layers.gif++";
        fileName = name.str(); 
        c1.Print(fileName);
    }

    outputFile->Close();
}









TH1F* getCorrectedSliceProjection(unsigned layer, unsigned numBins, TTree* tree, float min, float max, TF1 *fit) {
    
    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 
    std::ostringstream name;
    name << "Layer"<< layer << "_CorrEdgeProj";
    TString histName = name.str(); 
    TH1F *hist = new TH1F(histName,histName,numBins,-5,5); 

    for (unsigned event(0);event<tree->GetEntries();event++) {
        tree->GetEntry(event);
        int start  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue();
        if (start + layer > 27) continue;
        float x_ew = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(1)->GetName())->GetValue(start + layer);
        float x_t  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(3)->GetName())->GetValue(start + layer);
        float x_e  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(7)->GetName())->GetValue(start + layer);
        if (x_e > max || x_e < min) continue;
        if (x_ew < 9999 && fabs(x_e) < 5 && fabs(x_ew - x_t) < 5 ) {
            hist->Fill(x_ew - x_t - fit(x_e));
        }
    }       
    return hist;
}


TH1F* getBiasedSliceProjection(unsigned layer, unsigned numBins, TTree* tree, float min, float max) {
    
    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 
    std::ostringstream name;
    name << "Layer"<< layer << "_BiasEdgeProj";
    TString histName = name.str(); 
    TH1F *hist = new TH1F(histName,histName,numBins,-5,5); 

    for (unsigned event(0);event<tree->GetEntries();event++) {
        tree->GetEntry(event);
        int start  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue();
        if (start + layer > 27) continue;
        float x_ew = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(1)->GetName())->GetValue(start + layer);
        float x_t  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(3)->GetName())->GetValue(start + layer);
        float x_e  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(7)->GetName())->GetValue(start + layer);
        if (x_e > max || x_e < min) continue;
        if (x_ew < 9999 && fabs(x_e) < 5 && fabs(x_ew - x_t) < 5 ) {
            hist->Fill(x_ew - x_t);
        }
    }       
    return hist;
}



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


TH2F* getBiasVsEratio(unsigned layer, unsigned numBins, TTree* tree, float sliceMin,float  sliceMax) {
    
    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 
    std::ostringstream name;
    name << "Layer"<< layer << "_BiasVsEratio_" << sliceMin << "_" << sliceMax;
    TString histName = name.str();
    TH2F* hist = new TH2F(histName,histName,numBins,-5,5,numBins,0,1); 
    
    for (unsigned event(0);event<tree->GetEntries();event++) {
        tree->GetEntry(event);
        int start  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue();
        if (start + layer > 27) continue;
        float x_ew = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(1)->GetName())->GetValue(start + layer);
        float x_t  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(3)->GetName())->GetValue(start + layer);
        float x_e  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(7)->GetName())->GetValue(start + layer);
        float centralE = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(9)->GetName())->GetValue(start + layer);
        float totalE   = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(10)->GetName())->GetValue(start + layer);
        if (x_e > sliceMax || x_e < sliceMin) continue;
        if (x_ew < 9999 && fabs(x_e) < 5 && fabs(x_ew - x_t) < 5 ) {
            hist->Fill( x_ew-x_t, centralE/totalE );
        }
    }
    hist->SetOption("COLZ");
    hist->SetStats(kFALSE);

    return hist;
}


TH2F* getCorrectedBiasCurve(unsigned layer, unsigned numBins, TTree* tree, float eRatioCut, TF1* fit) {

    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 
    std::ostringstream name;
    name << "Layer"<< layer << "_Corrected";
    TString histName = name.str(); 
    TH2F* hist = new TH2F(histName,histName,200,-5,5,200,-5,-5); 

    for (unsigned event(0);event<tree->GetEntries();event++) {
        tree->GetEntry(event);
        int start  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue();
        if (start + layer > 27) continue;
        float x_ew = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(1)->GetName())->GetValue(start+layer);
        float x_t  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(3)->GetName())->GetValue(start+layer);
        float x_e  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(7)->GetName())->GetValue(start+layer);
        float centralE = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(9)->GetName())->GetValue(start+layer);
        float totalE   = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(10)->GetName())->GetValue(start+layer);
        if (x_ew < 9999 && fabs(x_e) < 5 && fabs(x_ew - x_t) < 5 && centralE/totalE < eRatioCut) { 
            hist->Fill(x_e,x_ew - x_t - fit(x_e));
        }
    }
    hist->SetOption("COLZ");
    hist->SetStats(kFALSE);

    return hist;
}

std::vector<TH1F*> getBiasSlices(unsigned const layer, unsigned const numSlices, unsigned const numBins, TTree* tree, float eRatioCut) {

    std::vector<TH1F*> sliceHists;
    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 

    float width = 10/(float)numSlices;
    float startVal(-5.0);

    for (unsigned slice(0);slice<numSlices;slice++) {

        float min = startVal + width*(float)slice;
        float max = startVal + width*( (float)slice + 1);
        
        std::ostringstream name;
        name << "Layer"<< layer << "_Bias" << "_slice" << slice;
        TString histName = name.str();
        TH1F* sliceHist = new TH1F(histName,histName,50,-5,5);
        for (unsigned event(0);event<tree->GetEntries();event++) {
            tree->GetEntry(event);
            int start  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue();
            if (start+layer > 27) continue;
            float x_e  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(7)->GetName())->GetValue(start+layer);
            if (x_e > max || x_e < min) continue;
            float x_ew = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(1)->GetName())->GetValue(start+layer);
            float x_t  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(3)->GetName())->GetValue(start+layer);
            float centralE = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(9)->GetName())->GetValue(start+layer);
            float totalE   = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(10)->GetName())->GetValue(start+layer);
            if (x_ew < 9999 && fabs(x_e) < 5 && fabs(x_ew-x_t) < 5 && centralE/totalE < eRatioCut) {
                sliceHist->Fill(x_ew - x_t);       
            }
        }
        sliceHists.push_back(sliceHist);
    }   
    return sliceHists;
}


float slicePeak3(TH1F * hist) {

    int histMin(-5),histMax(5);
    float peakX,peakY;
    peakX = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
    peakY = hist->GetYaxis()->GetBinCenter(hist->GetMaximumBin());

    TF1 *poly   = new TF1(hist->GetName() + TString("_poly"),
                          "[0] + [1]*x + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4) + [5]*pow(x,5)",
                          histMin, histMax);
    hist->Fit(poly,"QNR");

    TF1 *gaus   = new TF1(hist->GetName() + TString("_gaus"),
                          "gaus(0)", 
                          histMin, histMax);
    gaus->SetParameters( peakY, peakX, 1.0);
    hist->Fit(gaus,"QNR");
    
    if (gaus->GetChisquare() < poly->GetChisquare()) {
        return gaus->GetMaximumX(-5,5);
    }else{
        return poly->GetMaximumX(-5,5);
    }
}

float slicePeak2(TH1F * hist) {

    int histMin(-5),histMax(5);
    float peakX,peakY;
    peakX = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
    peakY = hist->GetYaxis()->GetBinCenter(hist->GetMaximumBin());

    TF1 *gaus   = new TF1(hist->GetName() + TString("_gaus"),
                         "gaus(0)", 
                         histMin, histMax);

    //Inital fit
    gaus->SetParameters( peakY, peakX, 1.5);
    hist->Fit(gaus,"QNR");

    float mean = gaus->GetParameter(1); 
    float rms  = gaus->GetParameter(2);
    gaus   = new TF1( hist->GetName() + TString("_gaus"),
                      "gaus(0)", mean - rms*2, mean + rms*2);

    hist->Fit(gaus, "QNR");
    return gaus->GetMaximumX(-5,5);
}

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

std::pair<float,float> errorBars(TH1F * hist) {
    
    unsigned peakBin = hist->GetMaximumBin();
    float leftIntegral  = hist->Integral(-1,peakBin);
    float rightIntegral = hist->Integral(peakBin,-1);

    std::pair<float,float> errorBars;

    for (unsigned i(0);i<peakBin;i++) {
        errorBars.first = fabs( hist->GetXaxis()->GetBinCenter( peakBin ) 
                                - hist->GetXaxis()->GetBinCenter((peakBin - i)) );
        float areaFraction = hist->Integral(peakBin - i, peakBin)/leftIntegral;
        if (areaFraction > 0.66) break;
    }
    for (unsigned i(0); i< hist->GetNbinsX() - peakBin; i++) {
        errorBars.second = fabs( hist->GetXaxis()->GetBinCenter( peakBin ) 
                                - hist->GetXaxis()->GetBinCenter((peakBin + i)) );
        float areaFraction = hist->Integral(peakBin, peakBin + i)/rightIntegral;
        if (areaFraction > 0.66) break;
    }

    return errorBars;

}

float slicePeak1(TH1F* hist) {
    int binmax = hist->GetMaximumBin();
    float x    = hist->GetXaxis()->GetBinCenter(binmax);
    return x;
}

Double_t skewNormal(Double_t *x, Double_t *par) {
    Float_t xx = x[0];
    Double_t f = (2.0/par[0]) * pow(2*3.14159,-0.5) * TMath::Exp((xx-par[1])/par[0]) * 0.5*(1 + TMath::Erf(pow(2,-0.5)*par[2]*(xx-par[1])/par[0]));
    return f;
}

