#include <vector>

void SquarePlots3x3() {

    gROOT->ProcessLine(".L scripts/effSigmaMacro.C");

    float eRatioCut(1.0);
   
    TFile *file = TFile::Open("RootFiles/out_Sq_V100.root"); 
    TFile *outputFile = new TFile("Plots_Squares/Plots2.root","RECREATE");
    outputFile->cd();

    TTree *tree = (TTree*)file->Get("tracks");
    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 

    int startLayer(-1);

    std::ostringstream fileName;
    fileName << "scripts/EffSigmaSquares/EffSigma_Truth3x3_LS" << startLayer << ".csv";
    string name = fileName.str();
    ofstream effsigma;
    effsigma.open(name.c_str());

    //Bias curves
    TH2F* layerBiases[28];
    TH2F* layerCorrection[28];
    TGraph * layerCurve[28];
    TF1* layerFits[28];
    TH2F* layerBiasVsEratio[28];
    TH1D* projectionCorr[28];
    TH1D* projectionBias[28];
    TH1F* projectionCorrEdge[28];
    TH1F* projectionBiasEdge[28];

    unsigned layerNum(28);
    TCanvas c1("c1");
 
    for (unsigned layer(5);layer<layerNum;layer++) {
        
        std::cout << "Processing layer " << layer << std::endl;

        //Bias curve
        layerBiases[layer] = getBiasCurve(layer,200,tree,eRatioCut,startLayer);

        //Divide into slices and fit a polynomial for correction
        unsigned const numSlices = 30, numBins = 50;
        std::vector<TH1F*> slices =  getBiasSlices(layer,numSlices,numBins,tree,eRatioCut,startLayer);
        float width = 10.0/(float)numSlices;

        Double_t sliceModes[numSlices + 2];
        Double_t sliceCentres[numSlices + 2];
        Double_t sliceErr1[numSlices + 2];
        Double_t sliceErr2[numSlices + 2];
        Double_t zeroes[numSlices + 2];

        for (unsigned slice(0);slice<slices.size();slice++) {
            sliceModes[slice+1] = slicePeak(slices[slice]);
            sliceCentres[slice+1] = -5.0 + width*(float)slice + 0.5*width;
        }

        sliceModes[0] = 0.0;
        sliceCentres[0] = -5.0;
        sliceErr1[0] = 0.01;
        sliceErr2[0] = 0.01;
        zeroes[0] = 0.0;
        sliceModes[numSlices + 1] = 0.0;
        sliceCentres[numSlices + 1] = 5.0;
        sliceErr1[numSlices + 1] = 0.01;
        sliceErr2[numSlices + 1] = 0.01;
        zeroes[numSlices + 1] = 0.0;

        for (unsigned slice(0);slice<slices.size();slice++) {
            std::pair<float,float> errBars = errorBars(slices[slice]);    
            sliceErr1[slice+1] = errBars.first;
            sliceErr2[slice+1] = errBars.second;
            zeroes[slice+1] = 0.0;
        }

        TGraphAsymmErrors *graph = new TGraphAsymmErrors(numSlices+2,sliceCentres,sliceModes,zeroes,zeroes,sliceErr1,sliceErr2);
        TF1 *fit   = new TF1("fit" ,"[0]*x + [1]*pow(x,3) + [2]*pow(x,5) + [3]*pow(x,7)",-5,5);
        graph->Fit(fit,"QNR");
        graph->SetLineColor(kBlue);
        graph->SetLineWidth(2);
        layerCurve[layer] = graph;
        layerFits[layer] = fit;

        //Correction
        layerCorrection[layer] = getCorrectedBiasCurve(layer, numBins, tree, eRatioCut, fit, startLayer);

        //Projections onto Y
        projectionCorr[layer] = layerCorrection[layer]->ProjectionY();
        projectionBias[layer] = layerBiases[layer]->ProjectionY();

        //Correction in range
        projectionCorrEdge[layer] = getCorrectedSliceProjection(layer, 100, tree, -1.0, 1.0, layerFits[layer], startLayer);
        projectionBiasEdge[layer] = getBiasedSliceProjection(layer, 100, tree, -1.0, 1.0, startLayer);


        if (layer == 6) {
            //Print example plots
            layerBiases[layer]->GetXaxis()->SetTitle("x_{measured} distance to edge (mm)");
            layerBiases[layer]->GetYaxis()->SetTitle("x_{measured} - x_{true} (mm)");
            layerBiases[layer]->Draw();
            c1.Print("Plots/LayerBiasExample.pdf");
            c1.Clear();
            graph->GetXaxis()->SetTitle("x_{measured} distance to edge (mm)");
            graph->GetYaxis()->SetTitle("x_{measured} - x_{true} (mm)");
            graph->Draw();
            fit->Draw("same");
            c1.Print("Plots/BiasSlicesExample.pdf");
            c1.Clear();
            layerCorrection[layer]->GetXaxis()->SetTitle("x_{measured} distance to edge (mm)");
            layerCorrection[layer]->GetYaxis()->SetTitle("x_{measured} - x_{true} (mm)");
            layerCorrection[layer]->Draw();
            c1.Print("Plots/BiasCorrectedExample.pdf");
            c1.Clear();
            projectionCorrEdge[layer]->GetXaxis()->SetTitle("x_{measured} - x_{true} (mm)");
            projectionCorrEdge[layer]->SetLineWidth(2);
            projectionCorrEdge[layer]->Draw();
            c1.Print("Plots/ProjectionCorrEdgeExample.pdf");
            c1.Clear();
            projectionBiasEdge[layer]->GetXaxis()->SetTitle("x_{measured} - x_{true} (mm)");
            projectionBiasEdge[layer]->SetLineWidth(2);
            projectionBiasEdge[layer]->Draw();
            c1.Print("Plots/ProjectionBiasEdgeExample.pdf");
            c1.Clear();
            projectionCorr[layer]->GetXaxis()->SetTitle("x_{measured} - x_{true} (mm)");
            projectionCorr[layer]->SetLineWidth(2);
            projectionCorr[layer]->Draw();
            c1.Print("Plots/ProjectionCorrExample.pdf");
            c1.Clear();
            projectionBias[layer]->GetXaxis()->SetTitle("x_{measured} - x_{true} (mm)");
            projectionBias[layer]->SetLineWidth(2);
            projectionBias[layer]->Draw();
            c1.Print("Plots/ProjectionBiasExample.pdf");
            c1.Clear();

        }
            



    }

    outputFile->cd();
    for (unsigned layer(0);layer<28;layer++) {
        layerBiases[layer]->Write();
        layerCorrection[layer]->Write();
        layerFits[layer]->Write();
    }

    c1.Print(Form("LayerBiases2_S%d.pdf(",startLayer));
    for (unsigned layer(0);layer<28;layer++) {
        layerBiases[layer]->Draw();
        c1.Print(Form("LayerBiases2_S%d.pdf",startLayer));
    }
    c1.Print(Form("LayerBiases2_S%d.pdf)",startLayer));

    c1.Print(Form("LayerCorrections2_S%d.pdf(",startLayer));
    for (unsigned layer(0);layer<28;layer++) {
        layerCorrection[layer]->Draw();
        c1.Print(Form("LayerCorrections2_S%d.pdf",startLayer));
    }
    c1.Print(Form("LayerCorrections2_S%d.pdf)",startLayer));



    Double_t effSigmaBiasEdge[28];
    Double_t effSigmaCorrEdge[28];
    Double_t layerNumbers[28];
    std::cout << setw(24) << "Eff Sigmas: " << std::endl;
    std::cout << setw(12) << "Bias" << setw(12) << "Corr" << std::endl;
    for (unsigned layer(0);layer<28;layer++) {
        layerNumbers[layer] = layer;
        effSigmaBiasEdge[layer] = effSigmaMacro(projectionBiasEdge[layer]);
        effSigmaCorrEdge[layer] = effSigmaMacro(projectionCorrEdge[layer]);
    }
   
    for (unsigned layer(0);layer<28;layer++) {
        effsigma << effSigmaBiasEdge[layer] << ",";
        effsigma << effSigmaCorrEdge[layer] << ",";
        effsigma << layerNumbers[layer] << ",";
        effsigma << std::endl;

        std::cout << setw(12) << effSigmaBiasEdge[layer];
        std::cout << setw(12) << effSigmaCorrEdge[layer];
        std::cout << setw(12) << layerNumbers[layer];
        std::cout << std::endl;
    }
    effsigma.close();

    outputFile->Close();
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
    Bool_t  aAdjacentCut[28];           //12
    Bool_t  bAdjacentCut[28];           //13
    Bool_t  cAdjacentCut[28];           //14
    Bool_t  dAdjacentCut[28];           //15
    Float_t truthDistsFromEdgeX[28];    //16
    Float_t truthDistsFromEdgeY[28];    //17

};

TH1F* getCorrectedSliceProjection(unsigned layer, unsigned numBins, TTree* tree, float min, float max, TF1 *fit, int startPick) {
    
    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 
    std::ostringstream name;
    name << "Layer"<< layer << "_CorrEdgeProj";
    TString histName = name.str(); 
    TH1F *hist = new TH1F(histName,histName,numBins,-5,5); 

    for (unsigned event(0);event<tree->GetEntries();event++) {
        tree->GetEntry(event);
        int start  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue();
        if (start != startPick && startPick >= 0) continue;
        if (start + layer > 27) continue;
        float x_ew = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(1)->GetName())->GetValue(start + layer);
        float x_t  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(3)->GetName())->GetValue(start + layer);
        float x_e  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(7)->GetName())->GetValue(start + layer);
        float x_et = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(16)->GetName())->GetValue(start + layer);
        bool  aCut = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(12)->GetName())->GetValue(start + layer);
        bool  cCut = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(14)->GetName())->GetValue(start + layer);
        if (x_e > 0 && cCut) continue;
        if (x_e < 0 && aCut) continue;
        if (x_e > max || x_e < min) continue;
        if (x_et > max || x_et < min) continue;
        if (x_ew < 9999 && fabs(x_e) < 5 && fabs(x_ew - x_t) < 5 ) {
            hist->Fill(x_ew - x_t - fit(x_e));
        }
    }       
    return hist;
}

TH1F* getBiasedSliceProjection(unsigned layer, unsigned numBins, TTree* tree, float min, float max, int startPick) {
    
    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 
    std::ostringstream name;
    name << "Layer"<< layer << "_BiasEdgeProj";
    TString histName = name.str(); 
    TH1F *hist = new TH1F(histName,histName,numBins,-5,5); 

    for (unsigned event(0);event<tree->GetEntries();event++) {
        tree->GetEntry(event);
        int start  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue();
        if (start != startPick && startPick >= 0) continue;
        if (start + layer > 27) continue;
        float x_ew = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(1)->GetName())->GetValue(start + layer);
        float x_t  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(3)->GetName())->GetValue(start + layer);
        float x_e  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(7)->GetName())->GetValue(start + layer);
        bool  aCut = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(12)->GetName())->GetValue(start + layer);
        bool  cCut = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(14)->GetName())->GetValue(start + layer);
        float x_et = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(16)->GetName())->GetValue(start + layer);
        if (x_et > max || x_et < min) continue;
        if (x_e > 0 && cCut) continue;
        if (x_e < 0 && aCut) continue;
        if (x_e > max || x_e < min) continue;
        if (x_ew < 9999 && fabs(x_e) < 5 && fabs(x_ew - x_t) < 5 ) {
            hist->Fill(x_ew - x_t);
        }
    }       
    return hist;
}

TH2F* getBiasCurve(unsigned layer, unsigned numBins, TTree* tree, float eRatioCut, int startPick) {

    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 
    std::ostringstream name;
    name << "Layer"<< layer << "_Bias";
    TString histName = name.str(); 
    TH2F *biasCurve = new TH2F(histName,histName,numBins,-5,5,numBins,-5,-5); 

    for (unsigned event(0);event<tree->GetEntries();event++) {
        tree->GetEntry(event);
        int start  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue();
        if (start != startPick && startPick >= 0) continue;
        if (start + layer > 27) continue;
        float x_ew = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(1)->GetName())->GetValue(start + layer);
        float x_t  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(3)->GetName())->GetValue(start + layer);
        float x_e  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(7)->GetName())->GetValue(start + layer);
        float x_et = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(16)->GetName())->GetValue(start + layer);
        float centralE = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(9)->GetName())->GetValue(start + layer);
        float totalE   = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(10)->GetName())->GetValue(start + layer);
        bool  aCut = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(12)->GetName())->GetValue(start + layer);
        bool  cCut = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(14)->GetName())->GetValue(start + layer);
        if (x_e > 0 && cCut) continue;
        if (x_e < 0 && aCut) continue;
        if (fabs(x_et) >= 5) continue;
        if (x_ew < 9999 && fabs(x_e) < 5 && fabs(x_ew - x_t) < 5 && centralE/totalE < eRatioCut) { 
            biasCurve->Fill(x_e,x_ew - x_t);
        }
    }
    biasCurve->SetOption("COLZ");
    biasCurve->SetStats(kFALSE);

    return biasCurve;
} 

TH2F* getCorrectedBiasCurve(unsigned layer, unsigned numBins, TTree* tree, float eRatioCut, TF1* fit, int startPick) {

    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 
    std::ostringstream name;
    name << "Layer"<< layer << "_Corrected";
    TString histName = name.str(); 
    TH2F* hist = new TH2F(histName,histName,200,-5,5,200,-5,-5); 

    for (unsigned event(0);event<tree->GetEntries();event++) {
        tree->GetEntry(event);
        int start  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue();
        if (start != startPick && startPick >= 0) continue;
        if (start + layer > 27) continue;
        float x_ew = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(1)->GetName())->GetValue(start+layer);
        float x_t  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(3)->GetName())->GetValue(start+layer);
        float x_e  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(7)->GetName())->GetValue(start+layer);
        float x_et = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(16)->GetName())->GetValue(start+layer);
        float centralE = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(9)->GetName())->GetValue(start+layer);
        float totalE   = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(10)->GetName())->GetValue(start+layer);
        bool  aCut = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(12)->GetName())->GetValue(start + layer);
        bool  cCut = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(14)->GetName())->GetValue(start + layer);
        if (x_e > 0 && cCut) continue;
        if (x_e < 0 && aCut) continue;
        if (fabs(x_et) >= 5) continue;
        if (x_ew < 9999 && fabs(x_e) < 5 && fabs(x_ew - x_t) < 5 && centralE/totalE < eRatioCut) { 
            hist->Fill(x_e,x_ew - x_t - fit(x_e));
        }
    }
    hist->SetOption("COLZ");
    hist->SetStats(kFALSE);

    return hist;
}

std::vector<TH1F*> getBiasSlices(unsigned const layer, unsigned const numSlices, unsigned const numBins, TTree* tree, float eRatioCut, int startPick) {

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
            if (start != startPick && startPick >= 0) continue;
            if (start+layer > 27) continue;
            float x_e  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(7)->GetName())->GetValue(start+layer);
            float x_et = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(16)->GetName())->GetValue(start+layer);
            if (x_e > max || x_e < min) continue;
            if (fabs(x_et) >= 5) continue;
            float x_ew = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(1)->GetName())->GetValue(start+layer);
            float x_t  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(3)->GetName())->GetValue(start+layer);
            float centralE = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(9)->GetName())->GetValue(start+layer);
            float totalE   = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(10)->GetName())->GetValue(start+layer);
            bool  aCut = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(12)->GetName())->GetValue(start + layer);
            bool  cCut = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(14)->GetName())->GetValue(start + layer);
            if (x_e > 0 && cCut) continue;
            if (x_e < 0 && aCut) continue;
            if (x_ew < 9999 && fabs(x_e) < 5 && fabs(x_ew-x_t) < 5 && centralE/totalE < eRatioCut) {
                sliceHist->Fill(x_ew - x_t);       
            }
        }
        sliceHists.push_back(sliceHist);
    }   
    return sliceHists;
}

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

float slicePeak(TH1F* hist) {
    int binmax = hist->GetMaximumBin();
    float x    = hist->GetXaxis()->GetBinCenter(binmax);
    return x;
}

