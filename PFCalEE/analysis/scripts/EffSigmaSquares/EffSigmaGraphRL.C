#include <fstream>
#include <sstream>
using namespace std;

void EffSigmaGraphRL() {

        std::ostringstream fileName;
        fileName << "EffSigma_Truth3x3_RL.csv";
        string name = fileName.str();
        ifstream file (name.c_str());

        float data[26][3];

        for (unsigned row(0);row<26;row++) {
            std::string line;
            std::getline(file,line);
            if ( !file.good() ) break;

            std::stringstream iss(line);

            for (unsigned col(0);col<3;col++) {
                std::string value;
                std::getline(iss,value,',');
                if (!iss.good()) break;

                std::stringstream convertor(value);
                convertor >> data[row][col];
           }
        }

        Float_t biased[26];
        Float_t corrected[26];
        Float_t layerNum[26];

        for (unsigned row(0);row<26;row++) {
            biased[row] = data[row][0];
            corrected[row] = data[row][1];
            layerNum[row] = data[row][2];
            std::cout << setw(12) << layerNum[row] << setw(12) << biased[row] << setw(12) << corrected[row] << std::endl;
        }

        TCanvas c1("c1");

        Int_t n = 26;
        TGraph *graph1 = new TGraph(n,layerNum,biased);
        graph1->SetLineColor(kRed);
        graph1->SetLineWidth(2);
        TGraph *graph2 = new TGraph(n,layerNum,corrected);
        graph2->SetLineColor(kBlue);
        graph2->SetLineWidth(2);
        
        graph1->Draw();
        graph1->GetYaxis()->SetRangeUser(0,1.0);
        graph1->GetXaxis()->SetTitle("Layer Number");
        graph1->GetYaxis()->SetTitle("x_{measured} - x_{true} (mm)");
        graph1->SetTitle("#sigma_{eff} within #pm1mm of edge before/after correction by layer");

        graph2->Draw("SAME");

        leg = new TLegend(0.1,0.7,0.4,0.9);
        leg->AddEntry(graph1,"Biased","l");
        leg->AddEntry(graph2,"Corrected","l");
        leg->Draw();

        c1.Print("EffectiveSigmaByLayerRL.pdf");

}
