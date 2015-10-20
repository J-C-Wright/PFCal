{


    struct TrackInfo {
        int showerStart;
        Float_t energyWeightedX[28];
        Float_t energyWeightedY[28];
        Float_t truthX[28];
        Float_t truthY[28];
        Float_t centralE[28];
        Float_t totalE[28];
        UInt_t  numHitsInLayer[28];
    };

    TFile *file = TFile::Open("out.root"); 
    TFile *outputFile = new TFile("Plots.root","RECREATE");

    TrackInfo test;

    TTree *testTree = (TTree*)file->Get("tracks");
    testTree->SetBranchAddress("truthInfo",&test.showerStart);
    TObjArray *branchNames = testTree->GetListOfBranches();
    
    TH1I *showerStart = new TH1I("ShowerStart","Shower Starting Layer",10,0,10); 
    for (unsigned event(0);event<testTree->GetEntries();event++) {
        testTree->GetEntry(event);
        TObjArray *leafNames = testTree->GetBranch("truthInfo")->GetListOfLeaves(); 
        showerStart->Fill(testTree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue()); 
    }

    outputFile->cd();
    showerStart->Write();
    outputFile->Close();

}
