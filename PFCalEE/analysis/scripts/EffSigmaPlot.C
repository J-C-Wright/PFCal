
{
    
    Float_t biased[26];
    Float_t corrected[26];
    Float_t layerNum[26];
    
    for (unsigned i(0);i<26;i++){
        layerNum[i] = i;
    }

    biased[0]  = 0.667665;
    biased[1]  = 0.607636;
    biased[2]  = 0.54963;
    biased[3]  = 0.579656;
    biased[4]  = 0.548637;
    biased[5]  = 0.547993; 
    biased[6]  = 0.543895;
    biased[7]  = 0.543769;
    biased[8]  = 0.511194;
    biased[9]  = 0.511893;
    biased[10] = 0.514434;
    biased[11] = 0.534443;
    biased[12] = 0.637665;
    biased[13] = 0.748733;
    biased[14] = 0.964178;
    biased[15] = 1.38924;
    biased[16] = 1.73794;
    biased[17] = 2.08476;
    biased[18] = 2.37447;
    biased[19] = 2.63418;
    biased[20] = 2.77303;
    biased[21] = 2.94005;
    biased[22] = 2.99536;
    biased[23] = 2.44778;
    biased[24] = 2.4611;
    biased[25] = 3.2647;
    
    corrected[0]  = 0.351831;
    corrected[1]  = 0.183066;
    corrected[2]  = 0.173379;
    corrected[3]  = 0.159082;
    corrected[4]  = 0.168923;
    corrected[5]  = 0.165208;
    corrected[6]  = 0.181237;
    corrected[7]  = 0.194054;
    corrected[8]  = 0.234673;
    corrected[9]  = 0.266813;
    corrected[10] = 0.327392;
    corrected[11] = 0.40381;
    corrected[12] = 0.623627;
    corrected[13] = 0.759603;
    corrected[14] = 1.05974;
    corrected[15] = 1.45927;
    corrected[16] = 1.76632;
    corrected[17] = 2.0998;
    corrected[18] = 2.38488;
    corrected[19] = 2.65896;
    corrected[20] = 2.76003;
    corrected[21] = 2.90864;
    corrected[22] = 2.92629;
    corrected[23] = 2.49889;
    corrected[24] = 2.4111;
    corrected[25] = 2.8647;
    
    TGraph graph1 = new TGraph(26,layerNum,biased);
    graph1->SetLineColor(kRed);
    TGraph graph2 = new TGraph(26,layerNum,corrected);
    graph2->SetLineColor(kBlue);
    
    graph1->Draw();
    graph2->Draw("same");
    
}



