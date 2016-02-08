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
        << " -- Number of runs " << nRuns << std::endl
        << " -- Sim filename "  << simFileName  << std::endl
        << " -- Reco filename " << recoFileName << std::endl
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
//    lstr << inputsim.str() << "_run" << i << ".root";
      lstr << filePath << "run_" << i << "/" << simFileName << "_run" << i << ".root";
      std::cout << "\n\n" << lstr.str() << "\n\n";

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
//      lstr << inputrec.str() << "_run" << i << ".root";
      lstr << filePath << "run_" << i << "/" << recoFileName << "_run" << i << ".root";

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

//TTree prep
    TrackInfo trackStruct;
    TString treeLeaves;
    TTree *trackTree = new TTree("tracks","position finding tree");
    treeLeaves  = TString("showerStart/I:energyWeightedX[28]/F:energyWeightedY[28]/F:truthX[28]/F:truthY[28]/F:distsFromHitCentreX[28]/F:");
    treeLeaves += TString("distsFromHitCentreY[28]/F:distsFromTileEdgesX[28]/F:distsFromTileEdgesY[28]/F:centralE[28]/F:totalE[28]/F:numHitsInLayer[28]/I:");
    treeLeaves += TString("energyWeightedU[28]/F:energyWeightedV[28]/F:energyWeightedW[28]/F:truthU[28]/F:truthV[28]/F:truthW[28]/F:");
    treeLeaves += TString("distsFromHitCentreU[28]/F:distsFromHitCentreV[28]/F:distsFromHitCentreW[28]/F:");
    treeLeaves += TString("distsFromHitCentreUPerp[28]/F:distsFromHitCentreVPerp[28]/F:distsFromHitCentreWPerp[28]/F");
    trackTree->Branch("truthInfo",&trackStruct.showerStart,treeLeaves);
    /*
    Float_t energyWeightedU[28];
    Float_t energyWeightedV[28];
    Float_t energyWeightedW[28];
    Float_t truthU[28];
    Float_t truthV[28];
    Float_t truthW[28];
    Float_t distsFromHitCentreU[28];
    Float_t distsFromHitCentreV[28];
    Float_t distsFromHitCentreW[28];
    Float_t distsFromHitCentreUPerp[28];
    Float_t distsFromHitCentreVPerp[28];
    Float_t distsFromHitCentreWPerp[28];
    */

//Calculating track truth
    std::vector<TrackTruth> tracks;
    bool ttpDebug = false;

    TrackTruthProducer trackTruthProducer(ttpDebug,nLayers,versionNumber);
    if (!trackTruthProducer.layersZsLoaded()) {
        std::cout << "Generating Z positions" << std::endl;
        trackTruthProducer.getZpositions(versionNumber, lSimTree, nEvts, nSiLayers);
    }


    unsigned photonCount(0);
    
    for (unsigned ievt(0); ievt<nEvts; ++ievt){

        std::cout << "\n... Processing entry: " << ievt << std::endl;
        std::cout << "Reading info from trees... ";

        lSimTree->GetEntry(ievt);
        lRecTree->GetEntry(ievt);

        std::cout << "Done!" << std::endl;
        std::cout << "Producing... " << std::endl;

        if ((*genvec).size() == 1) {
            photonCount++;
            trackTruthProducer.produce(genvec,rechitvec,geomConv,0.5,5,2);
            tracks.push_back(trackTruthProducer.getTrack(0));
            trackStruct = trackTruthProducer.trackStruct(0);
            trackTree->Fill();
        }else{
            std::cout << "Not a single photon -- Skipping event." << std::endl;
        }
        std::cout << "Done!" << std::endl;
    }
    std::cout << "There were " << photonCount << " photons" << std::endl;

    outputFile->cd();
    trackTree->Write();

    outputFile->Close();
  
    return 0;



}//main
