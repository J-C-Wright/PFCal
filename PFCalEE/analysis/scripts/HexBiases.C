#include <vector>

void HexBiases() {

    gROOT->ProcessLine(".L scripts/effSigmaMacro.C");

    TCanvas c1("c1");

    TFile *file = TFile::Open("RootFiles/out_Hex_V100_100k_ZS0p5.root");
    TTree *tree = (TTree*)file->Get("tracks");

    float a = 6.4;
    float b = 0.5*a*sqrt(3);
    unsigned bins(200);
    unsigned numSlices(20),sliceBins(40);
    unsigned numFilled(0);

    for (int showerStart(1);showerStart<2;showerStart++) {

        std::ofstream csv;
        std::ostringstream fileName;
        if (showerStart < 0) {
            fileName << "Plots_Hexes/EffSigmaAtEdge1mm_All.csv";
        }else{
            fileName << "Plots_Hexes/EffSigmaAtEdge1mm_S" << showerStart << ".csv";
        }
        csv.open(fileName.str().c_str());

        for (unsigned layer(0);layer<28;layer++) {

            std::cout << "Processing layer " << layer << std::endl;

            TF1* dummyFit = new TF1("noCorrection","0",b,b);

            TH2F * X  = getBiasCurve(layer,bins,tree,1,3,7,-1,b,0.0,a,dummyFit,numFilled,showerStart);
            TH2F * XR = getBiasCurve(layer,bins,tree,1,3,7,6,b,0.0,a/2.0,dummyFit,numFilled,showerStart);
            TH2F * UR = getBiasCurve(layer,bins,tree,12,15,24,21,b,0.0,a/2.0,dummyFit,numFilled,showerStart);
            TH2F * VR = getBiasCurve(layer,bins,tree,13,16,25,22,b,0.0,a/2.0,dummyFit,numFilled,showerStart);
            TH2F * WR = getBiasCurve(layer,bins,tree,14,17,26,23,b,0.0,a/2.0,dummyFit,numFilled,showerStart);

            X->Draw();
            c1.Print(Form("Plots_Hexes/HexBiasCurveXL%d.pdf",layer));
            XR->Draw();
            c1.Print(Form("Plots_Hexes/HexBiasCurveXRL%d.pdf",layer));
            UR->Draw();
            c1.Print(Form("Plots_Hexes/HexBiasCurveURL%d.pdf",layer));
            VR->Draw();
            c1.Print(Form("Plots_Hexes/HexBiasCurveVRL%d.pdf",layer));
            WR->Draw();
            c1.Print(Form("Plots_Hexes/HexBiasCurveWRL%d.pdf",layer));

            //Get slices of the bias curves
            std::vector<TH1F*> xrSlices = getBiasSlices(layer,numSlices,sliceBins,tree,1,3,7,6,b,0.0,a/2.0,numFilled,showerStart);
            std::vector<TH1F*> urSlices = getBiasSlices(layer,numSlices,sliceBins,tree,12,15,24,21,b,0.0,a/2.0,numFilled,showerStart);
            std::vector<TH1F*> vrSlices = getBiasSlices(layer,numSlices,sliceBins,tree,13,16,25,22,b,0.0,a/2.0,numFilled,showerStart);
            std::vector<TH1F*> wrSlices = getBiasSlices(layer,numSlices,sliceBins,tree,14,17,26,23,b,0.0,a/2.0,numFilled,showerStart);

            TF1 *corrFitX   = new TF1("corrFitX" ,"[0]*x + [1]*pow(x,3) + [2]*pow(x,5) + [3]*pow(x,7)",-b,b);
            TF1 *corrFitU   = new TF1("corrFitU" ,"[0]*x + [1]*pow(x,3) + [2]*pow(x,5) + [3]*pow(x,7)",-b,b);
            TF1 *corrFitV   = new TF1("corrFitV" ,"[0]*x + [1]*pow(x,3) + [2]*pow(x,5) + [3]*pow(x,7)",-b,b);
            TF1 *corrFitW   = new TF1("corrFitW" ,"[0]*x + [1]*pow(x,3) + [2]*pow(x,5) + [3]*pow(x,7)",-b,b);

            //Fit to the slices
            TGraphAsymmErrors* xrSlicesGraph = getSliceGraph(xrSlices,b);
            xrSlicesGraph->Fit("corrFitX","QNR");
            TGraphAsymmErrors* urSlicesGraph = getSliceGraph(urSlices,b);
            urSlicesGraph->Fit("corrFitU","QNR");
            TGraphAsymmErrors* vrSlicesGraph = getSliceGraph(vrSlices,b);
            vrSlicesGraph->Fit("corrFitV","QNR");
            TGraphAsymmErrors* wrSlicesGraph = getSliceGraph(wrSlices,b);
            wrSlicesGraph->Fit("corrFitW","QNR");

            //Make corrected bias plots
            TH2F* correctedXR = getBiasCurve(layer,bins,tree,1,3,7,6,b,0.0,a/2.0,corrFitX,numFilled,showerStart);
            correctedXR->Draw();
            c1.Print(Form("Plots_Hexes/XRCorrectedL%d.pdf",layer));
            TH2F* correctedUR = getBiasCurve(layer,bins,tree,12,15,24,21,b,0.0,a/2.0,corrFitU,numFilled,showerStart);
            correctedUR->Draw();
            c1.Print(Form("Plots_Hexes/URCorrectedL%d.pdf",layer));
            TH2F* correctedVR = getBiasCurve(layer,bins,tree,13,16,25,22,b,0.0,a/2.0,corrFitV,numFilled,showerStart);
            correctedVR->Draw();
            c1.Print(Form("Plots_Hexes/VRCorrectedL%d.pdf",layer));
            TH2F* correctedWR = getBiasCurve(layer,bins,tree,14,17,26,23,b,0.0,a/2.0,corrFitW,numFilled,showerStart);
            correctedWR->Draw();
            c1.Print(Form("Plots_Hexes/WRCorrectedL%d.pdf",layer));

            //Effective sigma at +/- 1mm
            //Corrected 
            TH1F* biasProjCorrEdgeXR = getBiasProjection(layer,bins,tree,1,3,7,6,1,0.0,a/2.0,corrFitX,showerStart);
            biasProjCorrEdgeXR->Draw();
            c1.Print(Form("Plots_Hexes/BiasProjectionAtEdgeXRL%d.pdf",layer));
            TH1F* biasProjCorrEdgeUR = getBiasProjection(layer,bins,tree,12,15,24,21,1,0.0,a/2.0,corrFitU,showerStart);
            biasProjCorrEdgeUR->Draw();
            c1.Print(Form("Plots_Hexes/BiasProjectionAtEdgeURL%d.pdf",layer));
            TH1F* biasProjCorrEdgeVR = getBiasProjection(layer,bins,tree,13,16,25,22,1,0.0,a/2.0,corrFitV,showerStart);
            biasProjCorrEdgeVR->Draw();
            c1.Print(Form("Plots_Hexes/BiasProjectionAtEdgeVRL%d.pdf",layer));
            TH1F* biasProjCorrEdgeWR = getBiasProjection(layer,bins,tree,14,17,26,23,1,0.0,a/2.0,corrFitW,showerStart);
            biasProjCorrEdgeWR->Draw();
            c1.Print(Form("Plots_Hexes/BiasProjectionAtEdgeWRL%d.pdf",layer));
            
            std::cout << "Effective sigma test XR-Edge" << effSigmaMacro(biasProjCorrEdgeXR) << std::endl;
            std::cout << "Effective sigma test UR-Edge" << effSigmaMacro(biasProjCorrEdgeUR) << std::endl;
            std::cout << "Effective sigma test VR-Edge" << effSigmaMacro(biasProjCorrEdgeVR) << std::endl;
            std::cout << "Effective sigma test WR-Edge" << effSigmaMacro(biasProjCorrEdgeWR) << std::endl;

            csv << layer << ",";
            csv << effSigmaMacro(biasProjCorrEdgeXR) << ",";
            csv << effSigmaMacro(biasProjCorrEdgeUR) << ",";
            csv << effSigmaMacro(biasProjCorrEdgeVR) << ",";
            csv << effSigmaMacro(biasProjCorrEdgeWR) << ",";

            //Uncorrected
            TH1F* biasProjEdgeXR = getBiasProjection(layer,bins,tree,1,3,7,6,1,0.0,a/2.0,dummyFit,showerStart);
            TH1F* biasProjEdgeUR = getBiasProjection(layer,bins,tree,12,15,24,21,1,0.0,a/2.0,dummyFit,showerStart);
            TH1F* biasProjEdgeVR = getBiasProjection(layer,bins,tree,13,16,25,22,1,0.0,a/2.0,dummyFit,showerStart);
            TH1F* biasProjEdgeWR = getBiasProjection(layer,bins,tree,14,17,26,23,1,0.0,a/2.0,dummyFit,showerStart);

            biasProjEdgeXR->Draw();
            c1.Print(Form("Plots_Hexes/BiasProjectionAtEdgeXRL%d.pdf",layer));
            biasProjEdgeUR->Draw();
            c1.Print(Form("Plots_Hexes/BiasProjectionAtEdgeURL%d.pdf",layer));
            biasProjEdgeVR->Draw();
            c1.Print(Form("Plots_Hexes/BiasProjectionAtEdgeVRL%d.pdf",layer));
            biasProjEdgeWR->Draw();
            c1.Print(Form("Plots_Hexes/BiasProjectionAtEdgeWRL%d.pdf",layer));

            std::cout << "Effective sigma test XR-Edge" << effSigmaMacro(biasProjEdgeXR) << std::endl;
            std::cout << "Effective sigma test UR-Edge" << effSigmaMacro(biasProjEdgeUR) << std::endl;
            std::cout << "Effective sigma test VR-Edge" << effSigmaMacro(biasProjEdgeVR) << std::endl;
            std::cout << "Effective sigma test WR-Edge" << effSigmaMacro(biasProjEdgeWR) << std::endl;

            csv << effSigmaMacro(biasProjEdgeXR) << ",";
            csv << effSigmaMacro(biasProjEdgeUR) << ",";
            csv << effSigmaMacro(biasProjEdgeVR) << ",";
            csv << effSigmaMacro(biasProjEdgeWR) << ",";
            csv << std::endl;

        }
        csv.close();
    }
}

//------------------------------------------------------------
//------------------------------------------------------------
//------------------------------------------------------------

TH2F* getBiasCurve(unsigned layer, unsigned numBins, TTree* tree, 
                   unsigned measured, unsigned truth, unsigned edge,
                   int perpCentreDisp, float limit, float perpMin, float perpMax, TF1* fit, unsigned numNonEmpty, int startPick) {

    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 
    std::ostringstream name;
    name << "Layer"<< layer << "_BiasVsEdge_" << measured << "_" << truth << "_";
    name << edge << "_" << perpCentreDisp << "_" << fit->GetName() << "_" << perpMin << "-" << perpMax << "_" << startPick;
    TString histName = name.str(); 
    std::cout << "Making bias curve " << histName << std::endl;
    TH2F *biasCurve = new TH2F(histName,histName,numBins,-limit,limit,numBins,-limit,-limit); 

    for (unsigned event(0);event<tree->GetEntries();event++) {

        tree->GetEntry(event);

        int start  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue();
        if (start + layer > 27) continue;
        if (start != startPick && startPick != -1) continue;
        float x_ew = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(measured)->GetName())->GetValue(start + layer);
        float x_t  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(truth)->GetName())->GetValue(start + layer);
        float x_e  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(edge)->GetName())->GetValue(start + layer);
        float x_cPerp;
        unsigned numInLayer = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(11)->GetName())->GetValue(start + layer);
        if (numInLayer < numNonEmpty) continue;
        if (perpCentreDisp > 0) {
            x_cPerp = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(perpCentreDisp)->GetName())->GetValue(start + layer);
        }else{
            x_cPerp = 0;
        }
        if (x_ew < 990 && fabs(x_e) < limit && fabs(x_ew - x_t) < limit && fabs(x_cPerp) < perpMax && fabs(x_cPerp) >= perpMin) {
            biasCurve->Fill(x_e,x_ew - x_t - fit(x_e));
        }

    }

    biasCurve->SetOption("COLZ");
    //biasCurve->SetStats(kFALSE);
    biasCurve->GetXaxis()->SetTitle("Distance from edge (mm)");
    biasCurve->GetYaxis()->SetTitle("Bias (x_{meas.} - x_{true}) (mm)");

    return biasCurve;
} 

std::vector<TH1F*> getBiasSlices(unsigned layer, unsigned numSlices, unsigned sliceBins, TTree* tree,
                                 unsigned measured, unsigned truth, unsigned edge,
                                 int perpCentreDisp, float limit, float perpMin, float perpMax, unsigned numNonEmpty, int startPick) {
    std::vector<TH1F*> sliceHists;
    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 

    float width = limit*2.0/(float)numSlices;
    float startVal(-limit);

    for (unsigned slice(0);slice<numSlices;slice++) {

        float xmin = startVal + width*(float)slice;
        float xmax = startVal + width*( (float)slice + 1);

        //Get slice limits
        float ymin = limit;
        float ymax = -limit;
        for (unsigned event(1);event<tree->GetEntries();event += 10) {
            tree->GetEntry(event);
            int start  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue();
            if (start+layer > 27) continue;
            if (start != startPick && startPick != -1) {continue;}
            float x_e  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(edge)->GetName())->GetValue(start+layer);
            if (x_e > xmax || x_e < xmin) continue;
            float x_ew = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(measured)->GetName())->GetValue(start+layer);
            float x_t  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(truth)->GetName())->GetValue(start+layer);
            unsigned numInLayer = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(11)->GetName())->GetValue(start + layer);
            if (numInLayer < numNonEmpty) continue;
            float x_cPerp;
            if (perpCentreDisp > 0) {
                x_cPerp = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(perpCentreDisp)->GetName())->GetValue(start + layer);
            }else{
                x_cPerp = 0;
            }
            if (x_ew < 9999 && fabs(x_e) < limit && fabs(x_ew-x_t) < limit && fabs(x_cPerp) < perpMax && fabs(x_cPerp) >= perpMin) {
                if (x_ew - x_t < ymin) ymin = x_ew - x_t;
                if (x_ew - x_t > ymax) ymax = x_ew - x_t;
            }
        }

        std::ostringstream name;
        name << "Layer"<< layer << "_Bias" << "_slice" << slice << "_" << measured << "_" << truth << "_" << edge << "_" << perpCentreDisp << "_" << startPick;
        TString histName = name.str();
        TH1F* sliceHist = new TH1F(histName,histName,sliceBins,ymin,ymax);
        std::cout << "Making slice histo " << histName << " with limits " << setw(12) << ymin << setw(12) << ymax << std::endl;
        for (unsigned event(0);event<tree->GetEntries();event++) {
            tree->GetEntry(event);
            int start  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue();
            if (start+layer > 27) continue;
            if (start != startPick && startPick != -1) {continue;}
            float x_e  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(edge)->GetName())->GetValue(start+layer);
            if (x_e > xmax || x_e < xmin) continue;
            float x_ew = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(measured)->GetName())->GetValue(start+layer);
            float x_t  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(truth)->GetName())->GetValue(start+layer);
            unsigned numInLayer = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(11)->GetName())->GetValue(start + layer);
            if (numInLayer < numNonEmpty) continue;
            float x_cPerp;
            if (perpCentreDisp > 0) {
                x_cPerp = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(perpCentreDisp)->GetName())->GetValue(start + layer);
            }else{
                x_cPerp = 0;
            }
            if (x_ew < 9999 && fabs(x_e) < limit && fabs(x_ew-x_t) < limit && fabs(x_cPerp) < perpMax && fabs(x_cPerp) >= perpMin) {
                sliceHist->Fill(x_ew - x_t);       
            }
        }
        sliceHists.push_back(sliceHist);
    }   
    return sliceHists;
}

std::pair<float,float> getErrorBars(TH1F * hist) {
    
    unsigned peakBin = hist->GetMaximumBin();
    float leftIntegral  = hist->Integral(-1,peakBin);
    float rightIntegral = hist->Integral(peakBin,-1);

    std::pair<float,float> errorBars;
    for (unsigned i(0);i<peakBin;i++) {
        errorBars.first = fabs( hist->GetXaxis()->GetBinCenter( peakBin ) 
                                - hist->GetXaxis()->GetBinCenter((peakBin - i)) )/leftIntegral;
        float areaFraction = hist->Integral(peakBin - i, peakBin)/leftIntegral;
        if (areaFraction > 0.66) break;
    }
    for (unsigned i(0); i< hist->GetNbinsX() - peakBin; i++) {
        errorBars.second = fabs( hist->GetXaxis()->GetBinCenter( peakBin ) 
                                - hist->GetXaxis()->GetBinCenter((peakBin + i)) )/rightIntegral;
        float areaFraction = hist->Integral(peakBin, peakBin + i)/rightIntegral;
        if (areaFraction > 0.66) break;
    }
    return errorBars;
}

TGraphAsymmErrors* getSliceGraph(std::vector<TH1F*> slices, float limit) {

    unsigned const numSlices = slices.size();

    Double_t sliceModes[numSlices + 2];
    Double_t sliceCentres[numSlices + 2];
    Double_t sliceErr1[numSlices + 2];
    Double_t sliceErr2[numSlices + 2];
    Double_t zeroes[numSlices + 2];

    //Modes and slice centres
    float width = limit*2.0/(float)numSlices;
    for (unsigned slice(0);slice<slices.size();slice++) {
        //sliceModes[slice+1] = slices[slice]->GetXaxis()->GetBinCenter(slices[slice]->GetMaximumBin());
        sliceModes[slice+1] = getPeak(slices[slice]);
        sliceCentres[slice+1] = -limit + width*(float)slice + 0.5*width;
    }

    //Zero points at either end of the range (cell centres)
    sliceModes[0] = 0.0;
    sliceCentres[0] = -limit;
    sliceErr1[0] = 0.01;
    sliceErr2[0] = 0.01;
    zeroes[0] = 0.0;
    sliceModes[numSlices + 1] = 0.0;
    sliceCentres[numSlices + 1] = limit;
    sliceErr1[numSlices + 1] = 0.01;
    sliceErr2[numSlices + 1] = 0.01;
    zeroes[numSlices + 1] = 0.0;

    //error bars
    for (unsigned slice(0);slice<slices.size();slice++) {
        std::pair<float,float> errorBars = getErrorBars(slices[slice]);    
        sliceErr1[slice+1] = errorBars.first;
        sliceErr2[slice+1] = errorBars.second;
        zeroes[slice+1] = 0.0;
    }

    TGraphAsymmErrors *graph = new TGraphAsymmErrors(numSlices+2,sliceCentres,sliceModes,zeroes,zeroes,sliceErr1,sliceErr2);
    graph->SetLineWidth(2);

    return graph;
}

float getPeak(TH1F* hist) {

    float mode = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
    float mean = hist->GetMean();
    float peak = (mode+mean)/2.0;

    std::pair<float,float> err = getErrorBars(hist);
    TF1 * f1 = new TF1("f1","gaus",peak - err.first, peak + err.second);
    hist->Fit("f1","QNR");
    return f1->GetMaximumX();

}
      
TH1F* getBiasProjection(unsigned layer, unsigned numBins, TTree* tree, 
                        unsigned measured, unsigned truth, unsigned edge,
                        int perpCentreDisp, float limit, float perpMin, float perpMax, TF1* fit, int startPick) {

    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 
    std::ostringstream name;
    name << "Layer"<< layer << "_BiasProjection_" << measured << "_" << truth << "_";
    name << edge << "_" << perpCentreDisp << "_" << fit->GetName() << "_" << perpMin << "-" << perpMax << "_" << limit << "_" << startPick;
    TString histName = name.str(); 
    std::cout << "Making bias curve projection " << histName << std::endl;
    TH1F *projection = new TH1F(histName,histName,numBins,-limit,limit);

    for (unsigned event(0);event<tree->GetEntries();event++) {

        tree->GetEntry(event);

        int start  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue();
        if (start + layer > 27) continue;
        if (start != startPick && startPick != -1) continue;
        float x_ew = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(measured)->GetName())->GetValue(start + layer);
        float x_t  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(truth)->GetName())->GetValue(start + layer);
        float x_e  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(edge)->GetName())->GetValue(start + layer);
        float x_cPerp;
        unsigned numInLayer = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(11)->GetName())->GetValue(start + layer);
        if (perpCentreDisp > 0) {
            x_cPerp = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(perpCentreDisp)->GetName())->GetValue(start + layer);
        }else{
            x_cPerp = 0;
        }
        if (x_ew < 990 && fabs(x_e) < limit && fabs(x_ew - x_t) < limit && fabs(x_cPerp) < perpMax && fabs(x_cPerp) >= perpMin) {
            projection->Fill(x_ew - x_t - fit(x_e));
        }

    }

    projection->GetXaxis()->SetTitle("Bias (x_{meas.} - x_{true}) (mm)");
    return projection;
}
