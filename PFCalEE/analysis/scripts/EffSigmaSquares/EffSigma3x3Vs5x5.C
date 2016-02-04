#include <fstream>
#include <sstream>
using namespace std;

void EffSigma3x3Vs5x5() {

        int start(-1);
    //Data in
        //3x3
        std::ostringstream fileName3x3;
        fileName3x3 << "EffSigma_Truth2_LS" << start << ".csv";
        string name3x3 = fileName3x3.str();
        ifstream file3x3(name3x3.c_str());

        float data3x3[26][3];

        for (unsigned row(0);row<26;row++) {
            std::string line;
            std::getline(file3x3,line);
            if ( !file3x3.good() ) break;

            std::stringstream iss(line);

            for (unsigned col(0);col<3;col++) {
                std::string value;
                std::getline(iss,value,',');
                if (!iss.good()) break;

                std::stringstream convertor(value);
                convertor >> data3x3[row][col];
           }
        }
        Float_t biased3x3[26];
        Float_t corrected3x3[26];
        Float_t layerNum[26];
        for (unsigned row(0);row<26;row++) {
            biased3x3[row] = data3x3[row][0];
            corrected3x3[row] = data3x3[row][1];
            layerNum[row] = data3x3[row][2];
        }

        //5x5
        std::ostringstream fileName5x5;
        fileName5x5 << "EffSigma_Truth5x5_LS" << start << ".csv";
        string name5x5 = fileName5x5.str();
        ifstream file5x5(name5x5.c_str());

        float data5x5[26][3];

        for (unsigned row(0);row<26;row++) {
            std::string line;
            std::getline(file5x5,line);
            if ( !file5x5.good() ) break;

            std::stringstream iss(line);

            for (unsigned col(0);col<3;col++) {
                std::string value;
                std::getline(iss,value,',');
                if (!iss.good()) break;

                std::stringstream convertor(value);
                convertor >> data5x5[row][col];
           }
        }
        Float_t biased5x5[26];
        Float_t corrected5x5[26];
        Float_t layerNum[26];
        for (unsigned row(0);row<26;row++) {
            biased5x5[row] = data5x5[row][0];
            corrected5x5[row] = data5x5[row][1];
            layerNum[row] = data5x5[row][2];
        }

    //Make graphs

        TCanvas c1("c1");

        Int_t n = 26;
        //3x3
        TGraph *graph1_3x3 = new TGraph(n,layerNum,biased3x3);
        graph1_3x3->SetLineColor(kRed);
        graph1_3x3->SetLineWidth(2);
        TGraph *graph2_3x3 = new TGraph(n,layerNum,corrected3x3);
        graph2_3x3->SetLineColor(kRed);
        graph2_3x3->SetLineWidth(2);
        graph2_3x3->SetLineStyle(2);
        //5x5
        TGraph *graph1_5x5 = new TGraph(n,layerNum,biased5x5);
        graph1_5x5->SetLineColor(kBlue);
        graph1_5x5->SetLineWidth(2);
        TGraph *graph2_5x5 = new TGraph(n,layerNum,corrected5x5);
        graph2_5x5->SetLineColor(kBlue);
        graph2_5x5->SetLineWidth(2);
        graph2_5x5->SetLineStyle(2);
        
        
        graph1_3x3->Draw();
        graph1_3x3->GetYaxis()->SetRangeUser(0,1.0);
        graph1_3x3->GetXaxis()->SetTitle("Number of layers after shower start");
        graph1_3x3->GetYaxis()->SetTitle("x_{measured} - x_{true} (mm)");
        if (start >= 0) {
            graph1_3x3->SetTitle(Form("#sigma_{eff} within #pm1mm of edge before/after correction (layer%d)",start));
        }else{
            graph1_3x3->SetTitle("#sigma_{eff} within #pm1mm of edge before/after correction");
        }
        graph2_3x3->Draw("SAME");
        graph1_5x5->Draw("SAME");
        graph2_5x5->Draw("SAME");

        leg = new TLegend(0.1,0.7,0.4,0.9);
        leg->AddEntry(graph1_3x3,"3x3 Biased","l");
        leg->AddEntry(graph2_3x3,"3x3 Corrected","l");
        leg->AddEntry(graph1_5x5,"5x5 Biased","l");
        leg->AddEntry(graph2_5x5,"5x5 Corrected","l");
        leg->Draw();

        if (start >= 0) {
            c1.Print(Form("EffSigmaByLayerComparison_%d.pdf",start));
        }else{
            c1.Print("EffSigmaByLayerComparisonAll.pdf");
        }

}
