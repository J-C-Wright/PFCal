#include <vector>

void HexCSVExport() {

    TCanvas c1("c1");

    TFile *file = TFile::Open("RootFiles/out_Hex_V100_100k_ZS0p5.root");
    TTree *tree = (TTree*)file->Get("tracks");

    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 
    for (unsigned branch(0);branch<leafNames->GetEntries();branch++) {
        std::cout << setw(6) << branch << "    " << leafNames->At(branch)->GetName() << std::endl;
    }

    float a = 6.4;
    float b = 0.5*a*sqrt(3);
    unsigned layer(6);
    unsigned bins(200);
    unsigned numSlices(20),sliceBins(40);
    unsigned numFilled(0);

    //.csv output
    outputToCSV(tree, layer, -1, b, a);
}
    




void outputToCSV(TTree* tree, unsigned layer, int size, float edgeDist, float cornerDist) {

    unsigned limit = size > 0 ? size : tree->GetEntries();

    TObjArray *leafNames = tree->GetBranch("truthInfo")->GetListOfLeaves(); 
    //Find normalisation factors
    std::pair<float,float> u_meas_maxmin(0,0);
    std::pair<float,float> v_meas_maxmin(0,0);
    std::pair<float,float> w_meas_maxmin(0,0);
    float maxE(0.0);
    std::pair<float,float> u_true_maxmin(0,0);
    std::pair<float,float> v_true_maxmin(0,0);
    std::pair<float,float> w_true_maxmin(0,0);
    float startMax(0.0);

    for (unsigned event(0);event<limit;event++) {

        tree->GetEntry(event);

        float u_meas = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(12)->GetName())->GetValue(layer);
        float v_meas = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(13)->GetName())->GetValue(layer);
        float w_meas = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(14)->GetName())->GetValue(layer);

        float totalE = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(10)->GetName())->GetValue(layer);

        float u_true = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(15)->GetName())->GetValue(layer);
        float v_true = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(16)->GetName())->GetValue(layer);
        float w_true = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(17)->GetName())->GetValue(layer);

        float uEdge_meas = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(24)->GetName())->GetValue(layer);
        float vEdge_meas = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(25)->GetName())->GetValue(layer);
        float wEdge_meas = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(26)->GetName())->GetValue(layer);

        if ( u_meas < 999 && fabs(uEdge_meas) < edgeDist &&
             v_meas < 999 && fabs(vEdge_meas) < edgeDist &&
             w_meas < 999 && fabs(wEdge_meas) < edgeDist &&
             fabs(u_meas-u_true) < edgeDist && 
             fabs(v_meas-v_true) < edgeDist && 
             fabs(w_meas-w_true) < edgeDist ) {

            if (u_meas < u_meas_maxmin.first)  u_meas_maxmin.first  = u_meas;
            if (u_meas > u_meas_maxmin.second) u_meas_maxmin.second = u_meas;
            if (v_meas < v_meas_maxmin.first)  v_meas_maxmin.first  = v_meas;
            if (v_meas > v_meas_maxmin.second) v_meas_maxmin.second = v_meas;
            if (w_meas < w_meas_maxmin.first)  w_meas_maxmin.first  = w_meas;
            if (w_meas > w_meas_maxmin.second) w_meas_maxmin.second = w_meas;

            if (totalE > maxE) maxE = totalE;

            if (u_true < u_true_maxmin.first)  u_true_maxmin.first  = u_true;
            if (u_true > u_true_maxmin.second) u_true_maxmin.second = u_true;
            if (v_true < v_true_maxmin.first)  v_true_maxmin.first  = v_true;
            if (v_true > v_true_maxmin.second) v_true_maxmin.second = v_true;
            if (w_true < w_true_maxmin.first)  w_true_maxmin.first  = w_true;
            if (w_true > w_true_maxmin.second) w_true_maxmin.second = w_true;

        }

    }
        
    float u_meas_norm = (u_meas_maxmin.second - u_meas_maxmin.first)/2.0;
    float v_meas_norm = (v_meas_maxmin.second - v_meas_maxmin.first)/2.0;
    float w_meas_norm = (w_meas_maxmin.second - w_meas_maxmin.first)/2.0;
    maxE = maxE/2.0;
    float u_true_norm = (u_true_maxmin.second - u_true_maxmin.first)/2.0;
    float v_true_norm = (v_true_maxmin.second - v_true_maxmin.first)/2.0;
    float w_true_norm = (w_true_maxmin.second - w_true_maxmin.first)/2.0;

    std::cout << "Normalization info: " << std::endl;
    std::cout << "u_meas " << setw(12) << u_meas_maxmin.first << setw(12) << u_meas_maxmin.second << setw(12) << u_meas_norm << std::endl;
    std::cout << setw(12)  << setw(12) << (u_meas_maxmin.first + u_meas_norm)/u_meas_norm << setw(12) << (u_meas_maxmin.second + u_meas_norm)/u_meas_norm << std::endl; 
    std::cout << "v_meas " << setw(12) << v_meas_maxmin.first << setw(12) << v_meas_maxmin.second << setw(12) << v_meas_norm << std::endl;
    std::cout << setw(12)  << setw(12) << (v_meas_maxmin.first - v_meas_norm)/v_meas_norm << setw(12) << (v_meas_maxmin.second - v_meas_norm)/v_meas_norm << std::endl; 
    std::cout << "w_meas " << setw(12) << w_meas_maxmin.first << setw(12) << w_meas_maxmin.second << setw(12) << w_meas_norm << std::endl;
    std::cout << setw(12)  << setw(12) << (w_meas_maxmin.first)/w_meas_norm << setw(12) << (w_meas_maxmin.second)/w_meas_norm << std::endl; 
    std::cout << "E max  " << setw(12) << maxE << std::endl;
    std::cout << "u_true " << setw(12) << u_true_maxmin.first << setw(12) << u_true_maxmin.second << setw(12) << u_true_norm << std::endl;
    std::cout << setw(12)  << setw(12) << (u_true_maxmin.first + u_true_norm)/u_true_norm << setw(12) << (u_true_maxmin.second + u_true_norm)/u_true_norm << std::endl; 
    std::cout << "v_true " << setw(12) << v_true_maxmin.first << setw(12) << v_true_maxmin.second << setw(12) << v_true_norm << std::endl;
    std::cout << setw(12)  << setw(12) << (v_true_maxmin.first - v_true_norm)/v_true_norm << setw(12) << (v_true_maxmin.second - v_true_norm)/v_true_norm << std::endl; 
    std::cout << "w_true " << setw(12) << w_true_maxmin.first << setw(12) << w_true_maxmin.second << setw(12) << w_true_norm << std::endl;
    std::cout << setw(12)  << setw(12) << (w_true_maxmin.first)/w_true_norm << setw(12) << (w_true_maxmin.second)/w_true_norm << std::endl; 

    std::ofstream csvNormalised;
    std::ofstream csvOriginal;
    std::ostringstream fileNameNormalised;
    std::ostringstream fileNameOriginal;
    fileNameNormalised << "Plots_Hexes/Layer" << layer << "_" << limit/1000 << "kSampleNormalised.csv";
    fileNameOriginal << "Plots_Hexes/Layer" << layer << "_" << limit/1000 << "kSampleOriginal.csv";
    csvNormalised.open(fileNameNormalised.str().c_str());
    csvOriginal.open(fileNameOriginal.str().c_str());
    std::cout << "Outputing hex info to " << fileNameNormalised.str();
    std::cout << "Outputing hex info to " << fileNameOriginal.str();
    
    for (unsigned event(0);event<limit;event++) {

        tree->GetEntry(event);

        float start  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(0)->GetName())->GetValue();

        float u_meas = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(12)->GetName())->GetValue(layer);
        float v_meas = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(13)->GetName())->GetValue(layer);
        float w_meas = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(14)->GetName())->GetValue(layer);

        float uPerp_meas = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(21)->GetName())->GetValue(layer);
        float vPerp_meas = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(22)->GetName())->GetValue(layer);
        float wPerp_meas = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(23)->GetName())->GetValue(layer);

        float uEdge_meas = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(24)->GetName())->GetValue(layer);
        float vEdge_meas = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(25)->GetName())->GetValue(layer);
        float wEdge_meas = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(26)->GetName())->GetValue(layer);

        float numIn7  = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(11)->GetName())->GetValue(layer);

        float totalE = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(10)->GetName())->GetValue(layer);

        float eRatio = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(9)->GetName())->GetValue(layer)/totalE;

        float u_true = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(15)->GetName())->GetValue(layer);
        float v_true = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(16)->GetName())->GetValue(layer);
        float w_true = tree->GetBranch("truthInfo")->GetLeaf(leafNames->At(17)->GetName())->GetValue(layer);

        if ( u_meas < 999 && fabs(uEdge_meas) < edgeDist &&
             v_meas < 999 && fabs(vEdge_meas) < edgeDist &&
             w_meas < 999 && fabs(wEdge_meas) < edgeDist &&
             fabs(u_meas-u_true) < edgeDist && 
             fabs(v_meas-v_true) < edgeDist && 
             fabs(w_meas-w_true) < edgeDist ) {

            csvNormalised << (start-14)/14  << ",";
            csvNormalised << (u_meas+u_meas_norm)/u_meas_norm << "," << (v_meas-v_meas_norm)/v_meas_norm << "," << w_meas/w_meas_norm << ",";
            csvNormalised << uPerp_meas/cornerDist << "," << vPerp_meas/cornerDist << "," << wPerp_meas/cornerDist << ",";
            csvNormalised << uEdge_meas/edgeDist << "," << vEdge_meas/edgeDist << "," << wEdge_meas/edgeDist << ",";
            csvNormalised << (numIn7-3.5)/3.5 << "," << (totalE-maxE)/maxE << "," << 2.0*(eRatio-0.5) << ",";
            csvNormalised << (u_true+u_true_norm)/u_true_norm << "," << (v_true-v_true_norm)/v_true_norm << "," << w_true/w_true_norm << std::endl;

            csvOriginal << start << ",";
            csvOriginal << u_meas << "," << v_meas << "," << w_meas << ",";
            csvOriginal << uPerp_meas << "," << vPerp_meas << "," << wPerp_meas << ",";
            csvOriginal << uEdge_meas << "," << vEdge_meas << "," << wEdge_meas << ",";
            csvOriginal << numIn7 << "," << totalE << "," << eRatio << ",";
            csvOriginal << u_true << "," << v_true << "," << w_true << std::endl;
            /*
            std::cout << (start-14.0)/14.0 << ",";
            std::cout << (u_meas+u_meas_norm)/u_meas_norm << "," << (v_meas-v_meas_norm)/v_meas_norm << "," << w_meas/w_meas_norm << ",";
            std::cout << uPerp_meas/cornerDist << "," << vPerp_meas/cornerDist << "," << wPerp_meas/cornerDist << ",";
            std::cout << uEdge_meas/edgeDist << "," << vEdge_meas/edgeDist << "," << wEdge_meas/edgeDist << ",";
            std::cout << (numIn7-3.5)/3.5 << "," << (totalE-maxE)/maxE << "," << 2.0*(eRatio-0.5) << ",";
            std::cout << (u_true+u_true_norm)/u_true_norm << "," << (v_true-v_true_norm)/v_true_norm << "," << w_true/w_true_norm << std::endl;
            */
        }
    }

    csvNormalised.close();
    csvOriginal.close();
    
    std::cout << " - Done." << std::endl;
    return;
}
       

