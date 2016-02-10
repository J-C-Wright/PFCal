#include <fstream>
#include <sstream>
using namespace std;

void EffSigmaGraph() {

        TCanvas c1("c1");

        TGraph *biasGraph[5];
        TGraph *corrGraph[5];
        leg = new TLegend(0.1,0.5,0.4,0.9);

        for (unsigned start(0);start<5;start++) {
            std::cout << "Layer start " << start << std::endl;
            std::ostringstream fileName;
            fileName << "EffSigma_Truth_LS" << start << ".csv";
            string name = fileName.str();
            ifstream file (name.c_str());
            std::cout << name.c_str() << std::endl;

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

            Float_t biased[24];
            Float_t corrected[24];
            Float_t layerNum[24];

            for (unsigned row(0);row<24;row++) {
                biased[row] = data[row][0];
                corrected[row] = data[row][1];
                layerNum[row] = data[row][2];
                /*
                std::cout << setw(12) << layerNum[row];
                std::cout << setw(12) << biased[row];
                std::cout << setw(12) << corrected[row];
                std::cout << std::endl;
                */
            }

            Int_t n = 24;
            TGraph *graph1 = new TGraph(n,layerNum,biased);
            graph1->SetLineColor(start+1);
            graph1->SetLineWidth(2);
            graph1->SetLineStyle(2);
            TGraph *graph2 = new TGraph(n,layerNum,corrected);
            graph2->SetLineColor(start+1);
            graph2->SetLineWidth(2);
            
            graph1->Draw();
            graph1->GetYaxis()->SetRangeUser(0,1.0);
            graph1->GetXaxis()->SetTitle("Number of layers after shower start");
            graph1->GetYaxis()->SetTitle("#sigma_{eff} of x_{measured} - x_{true} (mm)");

            graph2->Draw("SAME");

            biasGraph[start] = graph1;
            corrGraph[start] = graph2;
        }
        
        leg->AddEntry(biasGraph[0],"Biased SL0","l");
        leg->AddEntry(corrGraph[0],"Corrected SL0","l");
        biasGraph[4]->Draw();
        corrGraph[4]->Draw();
        for (unsigned start(1);start<5;start++) {
            leg->AddEntry(biasGraph[start],Form("Biased SL%d",start),"l");
            leg->AddEntry(corrGraph[start],Form("Corrected SL%d",start),"l");
            biasGraph[4-start]->Draw("same");
            corrGraph[4-start]->Draw("same");
        }
        leg->Draw();

        c1.Print("EffectiveSigmaByLayer0To4.pdf");

}
