{

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
    
//Shower start hisogram
    std::cout << "Making the shower start event histogram...";
    TH1I *showerStart = new TH1I("ShowerStart","Shower Starting Layer",10,0,10); 
    for (unsigned event(0);event<testTree->GetEntries();event++) {
        testTree->GetEntry(event);
        showerStart->Fill(testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue()); 
    }
    std::cout << " done" << std::endl;

//Number of empty layers after shower start
    TH1I *emptyLayers = new TH1I("empties","empties",10,0,10);
    for (unsigned event(0);event<testTree->GetEntries();event++) {
        unsigned empties(0);
        unsigned showerStartLayer = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue();
        for (unsigned layer(showerStartLayer);layer<28;layer++) {
            if (testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(11)->GetName())->GetValue(layer) == 0) {
                 empties++;
            }
        }
        emptyLayers->Fill(empties);
    }

//E_grid
    std::cout << "Making the Etotal histograms...";
    TH1F* energyTotals[28];
    for (unsigned layer(0);layer<28;layer++) {
        std::ostringstream name;
        name << "E_{total}_layer"<< layer;
        TString histName = name.str(); 
        TH1F* energyTotal = new TH1F(histName,histName,50,0,300);
        for (unsigned event(0);event<testTree->GetEntries();event++) {
            testTree->GetEntry(event);
            float totalE   = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(10)->GetName())->GetValue(layer); 
            if (totalE > 0.0) {energyTotal->Fill(totalE);}
        }
        energyTotals[layer] = energyTotal;
    }
    std::cout << " done" << std::endl;
    Float_t energyTotalVal[28];
    Float_t layerNumber[28];
    for (unsigned layer(0);layer<28;layer++) {
        layerNumber[layer] = layer;
        energyTotalVal[layer] = energyTotals[layer]->GetMean();
    }
    TGraph *energyTotalVsLayer = new TGraph(28,layerNumber,energyTotalVal);
    energyTotalVsLayer->SetTitle("Energy total vs layer");
        

//E_hit/E_grid
    std::cout << "Making the Ehit/Etotal histograms...";
    TH1F* eRatios[28];
    for (unsigned layer(0);layer<28;layer++) {
        std::ostringstream name;
        name << "E_{hit}/E_{total}_layer"<< layer;
        TString histName = name.str(); 
        TH1F* eRatio = new TH1F(histName,histName,50,0,1);
        for (unsigned event(0);event<testTree->GetEntries();event++) {
            testTree->GetEntry(event);
            float centralE = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(9)->GetName())->GetValue(layer);
            float totalE   = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(10)->GetName())->GetValue(layer); 
            eRatio->Fill(centralE/totalE);
        }
        eRatios[layer] = eRatio;
    }
    std::cout << " done" << std::endl;
    Float_t energyRatio[28];
    for (unsigned layer(0);layer<28;layer++) {
        energyRatio[layer] = eRatios[layer]->GetMean();
    }
    TGraph *eRatioVsLayer = new TGraph(28,layerNumber,energyRatio);
    eRatioVsLayer->SetTitle("Energy ratio vs layer");

//Bias curves (2D Hists)
    std::cout << "Making the bias curve (X)...";
    TH2F* biasCurveX[28];
    for (unsigned layer(0);layer<28;layer++) {
        std::ostringstream name;
        name << "x_{EW}-x_{T}_layer"<< layer;
        TString histName = name.str(); 
        TH2F *biasCurve = new TH2F(histName,histName,200,-.5,.5,200,-5,-5); 

        for (unsigned event(0);event<testTree->GetEntries();event++) {
            testTree->GetEntry(event);
            float x_ew = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(1)->GetName())->GetValue(layer);
            float x_t  = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(3)->GetName())->GetValue(layer);
            float x_c  = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(5)->GetName())->GetValue(layer);
            if (x_ew < 9999 && fabs(x_c) < 0.5 && fabs(x_ew - x_t) < 5) { 
                biasCurve->Fill(x_c,x_ew - x_t);
            }
        }
        biasCurveX[layer] = biasCurve;
    }
    std::cout << " done" << std::endl;
    std::cout << "Making the bias curve (X measured)...";
    TH2F* biasCurveXM[28];
    for (unsigned layer(0);layer<28;layer++) {
        std::ostringstream name;
        name << "x_{EW}-x_{T}_XM_layer"<< layer;
        TString histName = name.str(); 
        TH2F *biasCurve = new TH2F(histName,histName,200,-.5,.5,200,-5,-5); 

        for (unsigned event(0);event<testTree->GetEntries();event++) {
            testTree->GetEntry(event);
            float x_ew = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(1)->GetName())->GetValue(layer);
            float x_t  = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(3)->GetName())->GetValue(layer);
            if (x_ew < 9999 && fabs(x_ew) < 5 && fabs(x_ew - x_t) < 5) { 
                biasCurve->Fill(x_ew,x_ew - x_t);
            }
        }
        biasCurveXM[layer] = biasCurve;
    }
    std::cout << " done" << std::endl;



    std::cout << "Making the bias curve (X measured vs dist from edge)...";
    TH2F* biasCurveXM_Edge[28];
    for (unsigned layer(0);layer<28;layer++) {
        std::ostringstream name;
        name << "x_{EW}-x_{T}_XM_Edge_layer"<< layer;
        TString histName = name.str(); 
        TH2F *biasCurve = new TH2F(histName,histName,200,-.5,.5,200,-5,-5); 

        for (unsigned event(0);event<testTree->GetEntries();event++) {
            testTree->GetEntry(event);
            float x_ew = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(1)->GetName())->GetValue(layer);
            float x_t  = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(3)->GetName())->GetValue(layer);
            float x_e  = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(7)->GetName())->GetValue(layer);
            if (x_ew < 9999 && fabs(x_e) < 5 && fabs(x_ew - x_t) < 5) { 
                biasCurve->Fill(x_e,x_ew - x_t);
            }
        }
        biasCurveXM_Edge[layer] = biasCurve;
    }
    std::cout << " done" << std::endl;




//Bias by shower start
    TH2F *biasCurveX_SS[28];
    for (unsigned layer(0);layer<28;layer++) {
        std::ostringstream name;
        name << "x_{EW}-x_{T}_layerAfterSS"<< layer;
        TString histName = name.str(); 
        biasCurveX_SS[layer] = new TH2F(histName,histName,200,-.5,.5,200,-5,-5); 
    }
    for (unsigned event(0);event<testTree->GetEntries();event++) {
        testTree->GetEntry(event);
        unsigned startLayer = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue();
        for (unsigned layer(startLayer);layer<28;layer++) {
            float x_ew = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(1)->GetName())->GetValue(layer);
            float x_t  = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(3)->GetName())->GetValue(layer);
            float x_c  = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(5)->GetName())->GetValue(layer);
            if (x_ew < 9999 && fabs(x_c) < 0.5 && fabs(x_ew - x_t) < 5) { 
                biasCurveX_SS[layer-startLayer]->Fill(x_c,x_ew - x_t);
            }
        }
    } 

//Bias curve profiles
    std::cout << "Making the bias curve profile (X)...";
    TProfile* biasProfileX[28];
    for (unsigned layer(0);layer<28;layer++) {
        std::ostringstream name;
        name << "x_profile_layer"<< layer;
        TString histName = name.str(); 
        biasProfileX[layer] = biasCurveX[layer]->ProfileX(histName,1,-1,"s");
    }
    std::cout << " done" << std::endl;
    std::cout << "Making the bias curve profile (X measured)...";
    TProfile* biasProfileXM[28];
    for (unsigned layer(0);layer<28;layer++) {
        std::ostringstream name;
        name << "x_profile_XM_layer_"<< layer;
        TString histName = name.str(); 
        biasProfileXM[layer] = biasCurveXM[layer]->ProfileX(histName,1,-1,"s");
    }
    std::cout << " done" << std::endl;
    std::cout << "Making the bias curve profile (X measured)...";
    TProfile* biasProfileXM_Edge[28];
    for (unsigned layer(0);layer<28;layer++) {
        std::ostringstream name;
        name << "x_profile_XM_Edge_layer_"<< layer;
        TString histName = name.str(); 
        biasProfileXM_Edge[layer] = biasCurveXM_Edge[layer]->ProfileX(histName,1,-1,"s");
    }
    std::cout << " done" << std::endl;

//Superimposed dots and profile
    c1.Clear();
    gStyle->SetPalette(1);
    gStyle->SetNumberContours(250);
    for (unsigned layer(0);layer<28;layer++) {
        biasCurveX[layer]->SetOption("COLZ");
        biasCurveX_SS[layer]->SetOption("COLZ");
        biasCurveXM[layer]->SetOption("COLZ");
        biasCurveXM_Edge[layer]->SetOption("COLZ");
        biasCurveX[layer]->SetStats(kFALSE);
        biasCurveX_SS[layer]->SetStats(kFALSE);
        biasCurveXM[layer]->SetStats(kFALSE);
        biasCurveXM_Edge[layer]->SetStats(kFALSE);
    }

    biasProfileX[11]->SetLineColor(0);
    biasCurveX[11]->Draw();
    biasProfileX[11]->Draw("same");
    c1.Print("Plots/superimposeTest.pdf"); 

//Correction

    std::cout << "Making the bias curve correction plots (X measured vs dist from edge)...";
    TH2F* corrected_XM_Edge[28];
    TH2F* truePosition_XM_Edge[28];
    for (unsigned layer(0);layer<28;layer++) {

        std::ostringstream name;
        name << "Correction_XM_Edge_layer"<< layer;
        TString histName = name.str(); 
        TH2F *corrected  = new TH2F(histName,histName,200,-.5,.5,200,-5,-5);
        name << "TrueVsCorrected_layer" << layer;
        histName = name.str();
        TH2F* truePosition = new TH2F(histName,histName,200,-15,15,200,-15,15);

        corrected->SetOption("COLZ");
        corrected->SetStats(kFALSE);

        truePosition->SetOption("COLZ");
        truePosition->SetStats(kFALSE);

        std::cout << "Fitting layer " << layer << std::endl;
        TF1 *fit   = new TF1("fit","[0] + [1]*x + [2]*pow(x,2) + [3]*pow(x,3) + [4]*pow(x,4)",-5,5);
        biasCurveXM_Edge[layer];
        fit->SetParameters(0.0,1.0,0.0,0.0,0.0);
        fit->SetLineColor(kRed);
        biasCurveXM_Edge[layer]->Fit(fit);

        float p0 = fit->GetParameter(0);
        float p1 = fit->GetParameter(1);
        float p2 = fit->GetParameter(2);
        float p3 = fit->GetParameter(3);
        float p4 = fit->GetParameter(4);

        for (unsigned event(0);event<testTree->GetEntries();event++) {

            testTree->GetEntry(event);
            float x_ew = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(1)->GetName())->GetValue(layer);
            float x_t  = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(3)->GetName())->GetValue(layer);
            float x_e  = testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(7)->GetName())->GetValue(layer);
            float correction = p0 + p1*x_e + p2*pow(x_e,2) + p3*pow(x_e,3) + p4*pow(x_e,4);
            if (x_ew < 9999 && fabs(x_e) < 5) { 
                corrected->Fill(x_e,x_ew-x_t-correction);
            }
            truePosition->Fill(x_t,x_ew - correction);
           
        }
        corrected_XM_Edge[layer] = corrected;
        truePosition_XM_Edge[layer] = truePosition;
    }

//Alternative fit test



//gifs
    c1.Clear();
    for (unsigned int layer(0);layer<28;layer++) {
        biasProfileX[layer]->SetLineColor(1);
        biasProfileX[layer]->Draw(); 
        c1.Print("Plots/gifs/BiasProfileX.gif+20");
    }
    c1.Print("Plots/gifs/BiasProfileX.gif++");
    c1.Clear();
    for (unsigned int layer(0);layer<28;layer++) {
        biasProfileXM[layer]->SetLineColor(1);
        biasProfileXM[layer]->Draw(); 
        c1.Print("Plots/gifs/BiasProfileXM.gif+20");
    }
    c1.Print("Plots/gifs/BiasProfileXM.gif++");
    c1.Clear();
    for (unsigned int layer(0);layer<28;layer++) {
        biasProfileXM_Edge[layer]->SetLineColor(1);
        biasProfileXM_Edge[layer]->Draw(); 
        c1.Print("Plots/gifs/BiasProfileXM_Edge.gif+20");
    }
    c1.Print("Plots/gifs/BiasProfileXM_Edge.gif++");
    c1.Clear();
    for (unsigned int layer(0);layer<28;layer++) {
        biasCurveX[layer]->Draw(); 
        c1.Print("Plots/gifs/BiasCurveX.gif+20");
    }
    c1.Print("Plots/gifs/BiasCurveX.gif++");
    c1.Clear();
    for (unsigned int layer(0);layer<28;layer++) {
        biasCurveXM[layer]->Draw(); 
        c1.Print("Plots/gifs/BiasCurveXM.gif+20");
    }
    c1.Print("Plots/gifs/BiasCurveXM.gif++");
    c1.Clear();
    for (unsigned int layer(0);layer<28;layer++) {
        biasCurveXM_Edge[layer]->Draw(); 
        c1.Print("Plots/gifs/BiasCurveXM_Edge.gif+20");
    }
    c1.Print("Plots/gifs/BiasCurveXM_Edge.gif++");
    c1.Clear();
    for (unsigned int layer(0);layer<28;layer++) {
        biasCurveX_SS[layer]->Draw(); 
        c1.Print("Plots/gifs/BiasCurveX_SS.gif+20");
    }
    c1.Print("Plots/gifs/BiasCurveX_SS.gif++");

    c1.Clear();
    for (unsigned int layer(0);layer<28;layer++) {
        eRatios[layer]->Draw(); 
        c1.Print("Plots/gifs/ERatios.gif+20");
    }
    c1.Print("Plots/gifs/ERatios.gif++");
    c1.Clear();
    for (unsigned int layer(0);layer<28;layer++) {
        corrected_XM_Edge[layer]->Draw();
        c1.Print("Plots/gifs/corrected_XM_Edge.gif+20");
    }
    c1.Print("Plots/gifs/corrected_XM_Edge.gif++");
    c1.Clear();
    for (unsigned int layer(0);layer<28;layer++) {
        truePosition_XM_Edge[layer]->Draw();
        c1.Print("Plots/gifs/truePosition_XM_Edge.gif+20");
    }
    c1.Print("Plots/gifs/truePosition_XM_Edge.gif++");
        

//.root output
    outputFile->cd();
    for (unsigned layer(0);layer<28;layer++) {eRatios[layer]->Write();}
    for (unsigned layer(0);layer<28;layer++) {energyTotals[layer]->Write();}
    for (unsigned layer(0);layer<28;layer++) {biasProfileX[layer]->Write();}
    for (unsigned layer(0);layer<28;layer++) {biasProfileXM[layer]->Write();}
    for (unsigned layer(0);layer<28;layer++) {biasProfileXM_Edge[layer]->Write();}
    for (unsigned layer(0);layer<28;layer++) {biasCurveX[layer]->Write();}
    for (unsigned layer(0);layer<28;layer++) {biasCurveX_SS[layer]->Write();}
    for (unsigned layer(0);layer<28;layer++) {biasCurveXM[layer]->Write();}
    for (unsigned layer(0);layer<28;layer++) {biasCurveXM_Edge[layer]->Write();}
    showerStart->Write();
    eRatioVsLayer->Write();
    energyTotalVsLayer->Write();
    emptyLayers->Write();
    truePosition->Write();
    
    outputFile->Close();

}
