#include <vector>

Double_t skewNormal(Double_t *x, Double_t *par) {
    Float_t xx = x[0];
    Double_t f = (2.0/par[0]) * pow(2*3.14159,-0.5) * TMath::Exp((xx-par[1])/par[0]) * 0.5*(1 + TMath::Erf(pow(2,-0.5)*par[2]*(xx-par[1])/par[0]));
    return f;
}

void Slices() {

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

    TCanvas c1("c1");
    c1.cd();

    TFile *file = TFile::Open("out.root"); 
    TFile *outputFile = new TFile("Plots.root","RECREATE");

    TrackInfo test;

    TTree *testTree = (TTree*)file->Get("tracks");
    testTree->SetBranchAddress("truthInfo",&test.showerStart);
    TObjArray *leafNames = testTree->GetBranch("truthInfo")->GetListOfLeaves(); 

    //Bias curves
    std::cout << "Making the bias curves" << std::endl;
    TH2F* layerBiases[28];
    for (unsigned layer(0);layer<28;layer++) {
        std::ostringstream name;
        name << "Layer"<< layer << "_bias";
        TString histName = name.str(); 
        TH2F *biasCurve = new TH2F(histName,histName,200,-5,5,200,-5,-5); 

        for (unsigned event(0);event<testTree->GetEntries();event++) {
            testTree->GetEntry(event);
            float x_ew = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(1)->GetName())->GetValue(layer);
            float x_t  = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(3)->GetName())->GetValue(layer);
            float x_e  = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(7)->GetName())->GetValue(layer);
            if (x_ew < 9999 && fabs(x_e) < 5 && fabs(x_ew - x_t) < 5) { 
                biasCurve->Fill(x_e,x_ew - x_t);
            }
        }
        layerBiases[layer] = biasCurve;
        layerBiases[layer]->SetOption("COLZ");
        layerBiases[layer]->SetStats(kFALSE);
    }

    //Fitting to slices
    std::cout << "Calculating slices and fitting them by layer" << std::endl;
    std::vector<TF1*> layerFits;
    for (unsigned layer(0);layer<28;layer++) {

        std::cout << "Processing layer " << layer << std::endl;        

        //Collection of slices for one layer
        unsigned const sliceNum(30);
        float width = 10/(float)sliceNum;
        float start(-5.0);

        std::vector<TH1F*> slices;
        Double_t sliceModes[sliceNum + 2];
        Double_t sliceCentres[sliceNum + 2];

        for (unsigned slice(0);slice<sliceNum;slice++) {

            float min = start + width*(float)slice;
            float max = start + width*( (float)slice + 1);
            
            std::ostringstream name;
            name << "layer"<< layer << "_Slice" << slice;
            TString histName = name.str(); 
            TH1F* sliceHist = new TH1F(histName,histName,50,-5,5);
            for (unsigned event(0);event<testTree->GetEntries();event++) {
                testTree->GetEntry(event);
                float x_e  = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(7)->GetName())->GetValue(layer);
                if (x_e > max || x_e < min) continue;
                float x_ew = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(1)->GetName())->GetValue(layer);
                float x_t  = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(3)->GetName())->GetValue(layer);
                if (x_ew < 9999 && fabs(x_e) < 5 && fabs(x_ew-x_t) < 5) {
                    sliceHist->Fill(x_ew - x_t);       
                }
            }
            
            int binmax = sliceHist->GetMaximumBin();
            double x   = sliceHist->GetXaxis()->GetBinCenter(binmax);
            std::cout << setw(12) << slice << setw(12) << min << setw(12) << max << setw(12) << binmax << setw(12) << x << std::endl;
            sliceModes[slice+1] = x;
            sliceCentres[slice+1] = min+0.5*width;
            slices.push_back(sliceHist);

        }
        sliceModes[0] = 0.0;
        sliceCentres[0] = -5.0;
        sliceModes[sliceNum + 1] = 0.0;
        sliceCentres[sliceNum + 1] = 5.0;
        
        TGraph *layerGraph = new  TGraph(sliceNum+2,sliceCentres,sliceModes);

        std::ostringstream fitname;
        fitname << "layer" << layer << "_fit";
        TF1 *fit   = new TF1(TString(fitname.str()),"[0]*x + [1]*pow(x,3) + [2]*pow(x,5) + [3]*pow(x,7)",-5,5);
        layerGraph->Fit(fit);
        layerFits.push_back(fit);
    }

    //Corrections
    std::cout << "Making the bias-corrected curve" << std::endl;
    TH2F* layerBiasCorrection[28];
    TH2F* layerTrueVsMeasured[28];
    for (unsigned layer(0);layer<28;layer++) {
        
        std::ostringstream name;
        name << "Layer"<< layer << "_BiasCorrection";
        TString histName = name.str(); 
        TH2F *biasCurve  = new TH2F(histName,histName,200,-5,5,200,-5,-5); 
        name.str(std::string());
        name << "Layer"<< layer << "_TrueVsMeasured";
        histName = name.str();
        TH2F *trueVsMeasured = new TH2F(histName,histName,200,-5,5,200,-5,-5); 

        float p0 = layerFits[layer]->GetParameter(0);
        float p1 = layerFits[layer]->GetParameter(1);
        float p2 = layerFits[layer]->GetParameter(2);
        float p3 = layerFits[layer]->GetParameter(3);

        for (unsigned event(0);event<testTree->GetEntries();event++) {

            testTree->GetEntry(event);

            float x_ew = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(1)->GetName())->GetValue(layer);
            float x_t  = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(3)->GetName())->GetValue(layer);
            float x_e  = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(7)->GetName())->GetValue(layer);
            float correction = p0*x_e + p1*pow(x_e,3) + p2*pow(x_e,5) + p3*pow(x_e,7);

            if (x_ew < 9999 && fabs(x_e) < 5 && fabs(x_ew - x_t) < 5) { 
                biasCurve->Fill(x_e,x_ew - x_t - correction);
                if (x_ew - correction < -5) {
                    trueVsMeasured->Fill(x_t + 10, x_ew - correction + 10);
                }else if (x_ew - correction > 5) {
                    trueVsMeasured->Fill(x_t - 10, x_ew - correction - 10);
                }else{
                    trueVsMeasured->Fill(x_t, x_ew - correction);
                }

            }
        }

        layerBiasCorrection[layer] = biasCurve;
        layerBiasCorrection[layer]->SetOption("COLZ");
        layerBiasCorrection[layer]->SetStats(kFALSE);
        layerTrueVsMeasured[layer] = trueVsMeasured;
        layerTrueVsMeasured[layer]->SetOption("COLZ");
        layerTrueVsMeasured[layer]->SetStats(kFALSE);
    }

    c1.Clear();
    for (unsigned layer(0);layer<28;layer++) {
        layerBiases[layer]->Draw();
        layerFits[layer]->Draw("same");
        c1.Print("Plots/gifs/BiasAtEdge_Fitted.gif+20");
    }
    c1.Print("Plots/gifs/BiasAtEdge_Fitted.gif++");

    outputFile->cd();
    for (unsigned layer(0);layer<28;layer++) {
        layerFits[layer]->Write();
        layerBiases[layer]->Write();
        layerBiasCorrection[layer]->Write();
        layerTrueVsMeasured[layer]->Write();
    }
    outputFile->Close(); 

}


