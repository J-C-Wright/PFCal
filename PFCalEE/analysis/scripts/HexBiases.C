#include <vector>

void HexBiases() {

    gROOT->ProcessLine(".L scripts/effSigmaMacro.C");

    float eRatioCut(1.0);
   
    TFile *file = TFile::Open("RootFiles/out_Hex_V100_100k.root");

    TTree *tree = (TTree*)file->Get("tracks");

    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 
    for (unsigned branch(0);branch<leafNames->GetEntries();branch++) {
        std::cout << setw(6) << branch << "    " << leafNames->At(branch)->GetName() << std::endl;
    }


    TH2F * testX = getBiasCurve(6,200,tree,1,3,7,6.4*sqrt(3)*0.5);
    TH2F * testY = getBiasCurve(6,200,tree,2,4,8,6.4);
    
    testY->Draw();

}

TH2F* getBiasCurve(unsigned layer, unsigned numBins, TTree* tree, 
                   unsigned measured, unsigned truth, unsigned edge, float limit
                   ) {

    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 
    std::ostringstream name;
    name << "Layer"<< layer << "_BiasVsEdge_" << measured << "_" << truth << "_" << edge;
    TString histName = name.str(); 
    std::cout << "Making histo " << histName << std::endl;
    TH2F *biasCurve = new TH2F(histName,histName,numBins,-limit,limit,numBins,-limit,-limit); 

    unsigned skipCount(0);

    for (unsigned event(0);event<tree->GetEntries();event++) {
        tree->GetEntry(event);
        int start  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue();
        if (start + layer > 27) continue;
        float x_ew = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(measured)->GetName())->GetValue(start + layer);
        float x_t  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(truth)->GetName())->GetValue(start + layer);
        float x_e  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(edge)->GetName())->GetValue(start + layer);
        if (x_ew < 9999 && fabs(x_e) < limit && fabs(x_ew - x_t) < limit ) {
            biasCurve->Fill(x_e,x_ew - x_t);
        }
    }
    biasCurve->SetOption("COLZ");
    biasCurve->SetStats(kFALSE);

    return biasCurve;
} 
