#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/function.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFitResult.h"

#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSParameters.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"
#include "HGCSSDetector.hh"
#include "HGCSSGeometryConversion.hh"
#include "HGCSSPUenergy.hh"

#include "TrackTruthProducer.hh"


#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"

using namespace std;
using boost::lexical_cast;
namespace po=boost::program_options;

bool testInputFile(std::string input, TFile* & file){
  file = TFile::Open(input.c_str());
  
  if (!file) {
    std::cout << " -- Error, input file " << input.c_str() << " cannot be opened. Skipping..." << std::endl;
    return false;
  }
  else std::cout << " -- input file " << file->GetName() << " successfully opened." << std::endl;
  return true;
};


int main(int argc, char** argv){//main  

  //Input output and config options
  std::string cfg;
  bool concept;
  //size of signal region to perform Chi2 position fit.
  //in units of 2.5mm cells to accomodate different granularities
  unsigned nSR;
  //maximum value of residuals to use in error matrix: discard positions that are too far away 
  double residualMax;//mm
  unsigned pNevts;
  std::string filePath;
  std::string digifilePath;
  unsigned nRuns;
  std::string simFileName;
  std::string recoFileName;
  std::string outPath;
  unsigned nSiLayers;
  //0:do just the energies, 1:do fit+energies, 2: do zpos+fit+energies
  unsigned redoStep;
  unsigned debug;
  bool applyPuMixFix;

  po::options_description preconfig("Configuration"); 
  preconfig.add_options()("cfg,c",po::value<std::string>(&cfg)->required());
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
  po::notify(vm);
  po::options_description config("Configuration");
  config.add_options()
    //Input output and config options //->required()
    ("concept",        po::value<bool>(&concept)->default_value(true))
    ("nSR",            po::value<unsigned>(&nSR)->default_value(3))
    ("residualMax",    po::value<double>(&residualMax)->default_value(25))
    ("pNevts,n",       po::value<unsigned>(&pNevts)->default_value(0))
    ("filePath,i",     po::value<std::string>(&filePath)->required())
    ("digifilePath",   po::value<std::string>(&digifilePath)->default_value(""))
    ("nRuns",          po::value<unsigned>(&nRuns)->default_value(0))
    ("simFileName,s",  po::value<std::string>(&simFileName)->required())
    ("recoFileName,r", po::value<std::string>(&recoFileName)->required())
    ("outPath,o",      po::value<std::string>(&outPath)->required())
    ("nSiLayers",      po::value<unsigned>(&nSiLayers)->default_value(2))
    ("redoStep",       po::value<unsigned>(&redoStep)->default_value(0))
    ("debug,d",        po::value<unsigned>(&debug)->default_value(0))
    ("applyPuMixFix",  po::value<bool>(&applyPuMixFix)->default_value(false))
    ;

  // ("output_name,o",            po::value<std::string>(&outputname)->default_value("tmp.root"))

  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);

  std::string inFilePath = filePath+simFileName;

  size_t end=outPath.find_last_of(".");
  std::string outFolder = outPath.substr(0,end);

  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << filePath << std::endl
	    << " -- Digi Input file path: " << digifilePath << std::endl
	    << " -- Output file path: " << outPath << std::endl
	    << " -- Output folder: " << outFolder << std::endl
	    << " -- Requiring " << nSiLayers << " si layers." << std::endl
	    << " -- Number cells in signal region for fit: " << nSR << " cells" << std::endl
	    << " -- Residual max considered for filling matrix and fitting: " << residualMax << " mm" << std::endl
	    << " -- Apply PUMix fix? " << applyPuMixFix << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events." << std::endl;

  TRandom3 lRndm(1);
  std::cout << " -- Random number seed: " << lRndm.GetSeed() << std::endl;

  /////////////////////////////////////////////////////////////
  //input
  /////////////////////////////////////////////////////////////

  std::ostringstream inputsim;
  inputsim << filePath << "/" << simFileName;
  std::ostringstream inputrec;
  if (digifilePath.size()==0)
    inputrec << filePath << "/" << recoFileName;
  else 
    inputrec << digifilePath << "/" << recoFileName;

  //std::cout << inputsim.str() << " " << inputrec.str() << std::endl;

  HGCSSInfo * info;

  TChain *lSimTree = new TChain("HGCSSTree");
  TChain *lRecTree = 0;
  
  TFile * simFile = 0;
  TFile * recFile = 0;

  if (recoFileName.find("Digi") != recoFileName.npos) 
    lRecTree = new TChain("RecoTree");
  else lRecTree = new TChain("PUTree");

  if (nRuns == 0){
    if (!testInputFile(inputsim.str(),simFile)) return 1;
    lSimTree->AddFile(inputsim.str().c_str());
    if (!testInputFile(inputrec.str(),recFile)) return 1;
    lRecTree->AddFile(inputrec.str().c_str());
  }
  else {
    for (unsigned i(0);i<nRuns;++i){
      std::ostringstream lstr;
      lstr << inputsim.str() << "_run" << i << ".root";
      if (testInputFile(lstr.str(),simFile)){  
        if (simFile) info =(HGCSSInfo*)simFile->Get("Info");
        else {
        std::cout << " -- Error in getting information from simfile!" << std::endl;
        return 1;
        }
      }
      else continue;
      lSimTree->AddFile(lstr.str().c_str());
      lstr.str("");
      lstr << inputrec.str() << "_run" << i << ".root";
      if (!testInputFile(lstr.str(),recFile)) continue;
      lRecTree->AddFile(lstr.str().c_str());
    }
  }

  if (!lSimTree){
    std::cout << " -- Error, tree HGCSSTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  if (!lRecTree){
    std::cout << " -- Error, tree RecoTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }else{
    std::cout << "RecoTree opened successfully" << std::endl;
  }

  /////////////////////////////////////////////////////////////
  //Info
  /////////////////////////////////////////////////////////////

  const double cellSize = info->cellSize();
  const unsigned versionNumber = info->version();
  const unsigned model = info->model();
  
  //models 0,1 or 3.
  //bool isTBsetup = (model != 2);
  bool isCaliceHcal = versionNumber==23;//inFilePath.find("version23")!=inFilePath.npos || inFilePath.find("version_23")!=inFilePath.npos;

  //extract input energy

  std::cout << " -- Version number is : " << versionNumber 
	    << ", model = " << model
	    << ", cellSize = " << cellSize
	    << std::endl;


  //initialise detector
  HGCSSDetector & myDetector = theDetector();
 
  myDetector.buildDetector(versionNumber,concept,isCaliceHcal);

  const unsigned nLayers = myDetector.nLayers();
  const unsigned nSections = myDetector.nSections();

  std::cout << " -- N layers = " << nLayers << std::endl
	    << " -- N sections = " << nSections << std::endl;


  HGCSSGeometryConversion geomConv(model,cellSize);
  //set granularity to get cellsize for PU subtraction
  std::vector<unsigned> granularity;
  granularity.resize(nLayers,4);
  geomConv.setGranularity(granularity);

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  ///////// Output File // /////////////////////////
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////

  TFile *outputFile = TFile::Open(outPath.c_str(),"RECREATE");
  
  if (!outputFile) {
    std::cout << " -- Error, output file " << outPath << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }
  outputFile->cd();


  ///initialise PU density object

  HGCSSPUenergy puDensity("data/EnergyDensity.dat");

    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
    ///////// PositionTruth //////////////////////////
    //////////////////////////////////////////////////
    //////////////////////////////////////////////////

    const unsigned nEvts = ((pNevts > lSimTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lSimTree->GetEntries()) : pNevts) ;

    //loop on events
    HGCSSEvent * event = 0;
    std::vector<HGCSSSamplingSection> * ssvec = 0;
    std::vector<HGCSSSimHit> * simhitvec = 0;
    std::vector<HGCSSRecoHit> * rechitvec = 0;
    std::vector<HGCSSGenParticle> * genvec = 0;
    unsigned nPuVtx = 0;
  
    lSimTree->SetBranchAddress("HGCSSEvent",&event);
    lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
    lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
    lSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);
  
    lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
    if (lRecTree->GetBranch("nPuVtx")) lRecTree->SetBranchAddress("nPuVtx",&nPuVtx);

//Calculating track truth
    std::vector<TrackTruth> tracks;
    bool ttpDebug = true;
    TrackTruthProducer trackTruthProducer(ttpDebug,nLayers,versionNumber);
    unsigned photonCount(0);
    for (unsigned ievt(0); ievt<nEvts; ++ievt){

        std::cout << "\n... Processing entry: " << ievt << std::endl;

        lSimTree->GetEntry(ievt);
        lRecTree->GetEntry(ievt);

        bool include = false;
        for (unsigned i(0);i<(*genvec).size();i++) {
            for (unsigned j(0);j<6;j++) {
                float edge = ((float)j-2.5)*10.0;
                if ((*genvec)[i].x() > edge - 0.5 && (*genvec)[i].x() < edge + 0.5) {
                    include  = true;
                }
            }
        }


        if ((*genvec).size() != 1) include = false;

        if (include){
            photonCount++;
            trackTruthProducer.produce(genvec,rechitvec,geomConv,2,5);
            tracks.push_back(trackTruthProducer.getTrack(0));
        }else{
            std::cout << "Not a single photon, or not in edge range -- Skipping event." << std::endl;
        }

    }
    std::cout << "There were " << photonCount << " photons" << std::endl;
    std::cout << "Track truth info calculated. Now generating plots" << std::endl;

/*
    std::vector<TH2F*> layerEWXHistos(nLayers);
    std::vector<TH2F*> layerEWYHistos(nLayers);

    outputFile->cd();
    for (unsigned layerLoop(0);layerLoop<nLayers;layerLoop++) {
       
        //Set up histos 
        std::ostringstream name;
        //XY
        name << "x_{T}-x_{EW}_layer" << layerLoop;
        TString histName = name.str(); 
        layerEWXHistos[layerLoop] = new TH2F(histName,histName,50,0.45,0.5,50,-1,5);   
    
        for (unsigned eventLoop(0);eventLoop<tracks.size();eventLoop++) {

            //Truth position vs truth minus energy weighted histos
            TrackTruth track = tracks[eventLoop];       
            //Check it's not a dummy entry and fill the histos
            if (track.distsFromHitCentre[layerLoop].X() < 9999 && track.distsFromHitCentre[layerLoop].Y() < 9999) {
                if (track.distsFromHitCentre[layerLoop].X() > 0.45) {
                    layerEWXHistos[layerLoop]->Fill(track.distsFromHitCentre[layerLoop].X(),track.truthPositions[layerLoop].X() - track.energyWeightedXY[layerLoop].X());
                } else {
                    layerEWXHistos[layerLoop]->Fill(-1.0*(track.distsFromHitCentre[layerLoop].X()), -1.0*(track.truthPositions[layerLoop].X()-track.energyWeightedXY[layerLoop].X()));
                }
            }
        }

        layerEWXHistos[layerLoop]->Write();
    }
*/    
    outputFile->Close();
  
    return 0;

}//main



    /*
//Apply restriction in phi here, remove entries in tracks with phi outside range


//Building plots
    std::vector<TH2F*> layerEWXHistos(nLayers);
    std::vector<TH2F*> layerEWYHistos(nLayers);
    std::vector<TH2F*> layerEWrHistos(nLayers);
    std::vector<TH2F*> layerEWrPhiHistos(nLayers);
    Double_t errorEWX[nLayers];
    Double_t errorEWY[nLayers];
    Double_t meanEnergyDeposited[nLayers];
    Double_t meanNumEmptyIn3x3[nLayers];
    Double_t fracNoHits[nLayers];
    Double_t fracAll3x3Energy[nLayers];    
    Double_t layerNumbers[nLayers];
    for (unsigned i(0);i<nLayers;i++) {layerNumbers[i] = i;}

    //Mean energy deposited in each layer
    for (unsigned int hitLoop(0);hitLoop<(*rechitvec).size();hitLoop++) {
        meanEnergyDeposited[(*rechitvec)[hitLoop].layer()] += (*rechitvec)[hitLoop].energy();
    }

    for (unsigned layerLoop(0);layerLoop<nLayers;layerLoop++) {
       
        //Set up histos 
        std::ostringstream name;
        //XY
        name << "x_{T}-x_{EW}_layer"<< layerLoop;
        TString histName = name.str(); 
        layerEWXHistos[layerLoop] = new TH2F(histName,histName,50,-0.5,0.5,50,-5,5);   
        name.str(std::string());
        name << "y_{T}-y_{EW}_layer"<< layerLoop;
        histName = name.str(); 
        layerEWYHistos[layerLoop] = new TH2F(histName,histName,50,-0.5,0.5,50,-5,5);   
        //rPhi
        name.str(std::string());
        name << "r_{T}-r_{EW}_layer"<< layerLoop;
        histName = name.str(); 
        layerEWrHistos[layerLoop] = new TH2F(histName,histName,50,-10,10,50,-10,10);   
        name.str(std::string());
        name << "rPhi_{T}-rPhi_{EW}_layer"<< layerLoop;
        histName = name.str(); 
        layerEWrPhiHistos[layerLoop] = new TH2F(histName,histName,50,-10,10,50,-15,15);   
        


        errorEWX[layerLoop] = 0;
        errorEWY[layerLoop] = 0;
        meanNumEmptyIn3x3[layerLoop] = 0;
        fracNoHits[layerLoop] = 0;
        fracAll3x3Energy[layerLoop] = 0;

        for (unsigned eventLoop(0);eventLoop<tracks.size();eventLoop++) {

            //Truth position vs truth minus energy weighted histos
            TrackTruth track = tracks[eventLoop];       

            //Check it's not a dummy entry and fill the histos
            if (track.distsFromHitCentre[layerLoop].X() < 9999 && track.distsFromHitCentre[layerLoop].Y() < 9999) {
                layerEWXHistos[layerLoop]->Fill(track.distsFromHitCentre[layerLoop].X(),track.truthPositions[layerLoop].X() - track.energyWeightedXY[layerLoop].X());
                layerEWYHistos[layerLoop]->Fill(track.distsFromHitCentre[layerLoop].Y(),track.truthPositions[layerLoop].Y() - track.energyWeightedXY[layerLoop].Y());
                layerEWrHistos[layerLoop]->Fill(track.distsFromHitCentrerPhi[layerLoop].first,
                                                track.truthPositionsrPhi[layerLoop].first - track.energyWeightedrPhi[layerLoop].first);
                layerEWrPhiHistos[layerLoop]->Fill(track.distsFromHitCentrerPhi[layerLoop].second,
                                                   track.truthPositionsrPhi[layerLoop].second - track.energyWeightedrPhi[layerLoop].second);
            }

            //Error in position by layer
            errorEWX[layerLoop] += track.truthPositions[layerLoop].X() - track.energyWeightedXY[layerLoop].X();       
            errorEWY[layerLoop] += track.truthPositions[layerLoop].Y() - track.energyWeightedXY[layerLoop].Y();       

            //Mean number of empty cells
            meanNumEmptyIn3x3[layerLoop] += 9 - track.hitsByLayer3x3[layerLoop].size();
            //Fraction no hits
            if (track.distsFromHitCentre[layerLoop].X() > 9000) {fracNoHits[layerLoop] += 1.0;} 
            //Fraction where all 9 in the 3x3 have energy
            if (track.hitsByLayer3x3[layerLoop].size() == 9) {fracAll3x3Energy[layerLoop] += 1.0;}       

        }

        errorEWX[layerLoop] /= (float)photonCount;
        errorEWY[layerLoop] /= (float)photonCount;
        meanEnergyDeposited[layerLoop] /= (float)photonCount;
        meanNumEmptyIn3x3[layerLoop] /= (float)photonCount;
        fracNoHits[layerLoop] /= (float)photonCount;
        fracAll3x3Energy[layerLoop] /= (float)photonCount;

    }

//Phi slice analysis
    //20 histos distributed in phi
    unsigned int nPhiSegments(20);
    float phiSegment = 3.14159/(float)nPhiSegments;
    float dPhi = 3.141459*10/180;
    unsigned layerStart(8), layerEnd(15);
    std::vector<TH2F*> phiSegmentEWXHistos(nPhiSegments);
    std::vector<TH2F*> phiSegmentEWYHistos(nPhiSegments);

    for (unsigned segment(0);segment<nPhiSegments;segment++) {

        float centralPhi = phiSegment*segment;

        std::ostringstream name;
        name << "EWX_Segment_" << centralPhi;
        TString histName = name.str(); 
        phiSegmentEWXHistos[segment] = new TH2F(histName,histName,50,-0.5,0.5,50,-3,3);   

        name.str(std::string());
        name << "EWY_Segment_" << centralPhi;
        histName = name.str(); 
        phiSegmentEWYHistos[segment] = new TH2F(histName,histName,50,-0.5,0.5,50,-3,3);   

        for (unsigned trackLoop(0);trackLoop<tracks.size();trackLoop++) {
            TrackTruth track = tracks[trackLoop];       
            if (track.particleInfo.phi() > centralPhi-dPhi && track.particleInfo.phi() < centralPhi+dPhi) { 
                for (unsigned layerLoop(layerStart);layerLoop<layerEnd+1;layerLoop++) {
                    phiSegmentEWXHistos[segment]->Fill(track.distsFromHitCentre[layerLoop].X(),track.truthPositions[layerLoop].X() - track.energyWeightedXY[layerLoop].X());
                    phiSegmentEWYHistos[segment]->Fill(track.distsFromHitCentre[layerLoop].Y(),track.truthPositions[layerLoop].Y() - track.energyWeightedXY[layerLoop].Y());
                } 
            }
        }

    }


    const Int_t numLayers = nLayers;
    TGraph * errorEWXGraph = new TGraph(numLayers,layerNumbers,errorEWX);
    errorEWXGraph->SetTitle("errorEWX");
    TGraph * errorEWYGraph = new TGraph(numLayers,layerNumbers,errorEWY);
    errorEWYGraph->SetTitle("errorEWY");
    TGraph * meanEnergyDepositedGraph = new TGraph(numLayers,layerNumbers,meanEnergyDeposited);
    meanEnergyDepositedGraph->SetTitle("meanEnergyDeposited");
    TGraph * meanNumEmptyIn3x3Graph = new TGraph(numLayers,layerNumbers,meanNumEmptyIn3x3);
    meanNumEmptyIn3x3Graph->SetTitle("meanNumEmptyIn3x3");
    TGraph * fracNoHitsGraph = new TGraph(numLayers,layerNumbers,fracNoHits);
    fracNoHitsGraph->SetTitle("fracNoHits");
    TGraph * fracAll3x3EnergyGraph = new TGraph(numLayers,layerNumbers,fracAll3x3Energy);
    fracAll3x3EnergyGraph->SetTitle("fracAll3x3Energy");

    errorEWXGraph->Write(); 
    errorEWYGraph->Write(); 
    meanEnergyDepositedGraph->Write(); 
    meanNumEmptyIn3x3Graph->Write(); 
    fracNoHitsGraph->Write(); 
    fracAll3x3EnergyGraph->Write(); 

    //Analyze x_t-x_ew vs y_t-y_ew averages
    TH2F *meanBias = new TH2F("meanBias","meanBias",50,-0.5,0.5,50,-0.5,0.5);
    Double_t meanBiasX[nLayers];
    Double_t meanBiasY[nLayers];
    for (unsigned int layerLoop(0);layerLoop<nLayers;layerLoop++) {
        std::cout << setw(12) << layerEWXHistos[layerLoop]->GetMean(2);
        std::cout << setw(12) << layerEWYHistos[layerLoop]->GetMean(2) << std::endl;   
        meanBiasX[layerLoop] = layerEWXHistos[layerLoop]->GetMean(2);
        meanBiasY[layerLoop] = layerEWYHistos[layerLoop]->GetMean(2);
    }
    TGraph * meanBiasXGraph = new TGraph(numLayers,layerNumbers,meanBiasX);
    TGraph * meanBiasYGraph = new TGraph(numLayers,layerNumbers,meanBiasY);

    TCanvas c1("c1");
    for (unsigned int layerLoop(1);layerLoop<nLayers;layerLoop++) {
        layerEWXHistos[layerLoop]->Draw();        
        c1.Print("gifs/EWX.gif+20");
    }
    c1.Print("gifs/EWX.gif++");

    c1.Clear();
    for (unsigned int layerLoop(1);layerLoop<nLayers;layerLoop++) {
        layerEWYHistos[layerLoop]->Draw();
        c1.Print("gifs/EWY.gif+20");
    }
    c1.Print("gifs/EWY.gif++");

    c1.Clear();
    for (unsigned int layerLoop(1);layerLoop<nLayers;layerLoop++) {
        layerEWrHistos[layerLoop]->Draw();
        c1.Print("gifs/EWr.gif+20");
    }
    c1.Print("gifs/EWr.gif++");

    c1.Clear();
    for (unsigned int layerLoop(1);layerLoop<nLayers;layerLoop++) {
        layerEWrPhiHistos[layerLoop]->Draw();
        c1.Print("gifs/EWrPhi.gif+20");
    }
    c1.Print("gifs/EWrPhi.gif++");

    c1.Clear();
    for (unsigned segment(0);segment<nPhiSegments;segment++) {
        phiSegmentEWXHistos[segment]->Draw();
        c1.Print("segmentEWX.gif+20");
    }
    c1.Print("gifs/segmentEWX.gif++");

    c1.Clear();
    for (unsigned segment(0);segment<nPhiSegments;segment++) {
        phiSegmentEWYHistos[segment]->Draw();
        c1.Print("gifs/segmentEWY.gif+20");
    }
    c1.Print("gifs/segmentEWY.gif++");

    phiSegmentEWXHistos[0]->Draw();
    c1.Print("EWX0.pdf");
    phiSegmentEWXHistos[10]->Draw();
    c1.Print("EWX10.pdf");

    c1.Clear();
    meanBiasXGraph->SetLineColor(2);
    meanBiasYGraph->Draw();
    meanBiasXGraph->SetLineColor(4);
    meanBiasXGraph->Draw("same");
    c1.Print("meanBiases.pdf");

*/

