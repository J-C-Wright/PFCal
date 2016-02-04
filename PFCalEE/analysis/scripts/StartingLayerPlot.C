
void StartingLayerPlot() {

    TCanvas c1("c1");

    TFile *file = TFile::Open("RootFiles/out_Sq_V100.root"); 
    TTree *tree = (TTree*)file->Get("tracks");

    TH1I *startingLayer = new TH1I("StartLayer","Shower Start Layers",20,0,20);

    for (unsigned event(0);event<tree->GetEntries();event++) {
        tree->GetEntry(event);
        unsigned start = tree->GetBranch("truthInfo")->GetLeaf("showerStart")->GetValue();
        startingLayer->Fill(start);
    }
    startingLayer->GetXaxis()->SetTitle("Starting Layer Number");
    startingLayer->Draw();
    c1.Print("StartingLayerPlot.pdf");

}
