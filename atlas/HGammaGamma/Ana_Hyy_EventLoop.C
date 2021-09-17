// 
// Run this script by simply do
// > root -l -b -q Ana_Hyy_EventLoop.C
//
#include <TChain.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

using namespace std;

void Analyze(TString, TString);

void Ana_Hyy_EventLoop(){
  TString inDir = "/Users/nurfikri89/WorkspaceOpenData/atlas/samples/2020/GamGam/";  

  //
  // MC samples
  //
  Analyze(inDir+"MC/mc_343981.ggH125_gamgam.GamGam.root",       "MC.ggH125_gamgam");
  Analyze(inDir+"MC/mc_345041.VBFH125_gamgam.GamGam.root",      "MC.VBFH125_gamgam");
  Analyze(inDir+"MC/mc_345319.ZH125J_Zincl_gamgam.GamGam.root", "MC.ZH125J_Zincl_gamgam");
  Analyze(inDir+"MC/mc_345318.WpH125J_Wincl_gamgam.GamGam.root","MC.WpH125J_Wincl_gamgam");
  Analyze(inDir+"MC/mc_341081.ttH125_gamgam.GamGam.root",       "MC.ttH125_gamgam");

  //
  // Data samples
  //
  Analyze(inDir+"Data/data_A.GamGam.root",       "Data.A_GamGam");
  Analyze(inDir+"Data/data_B.GamGam.root",       "Data.B_GamGam");
  Analyze(inDir+"Data/data_C.GamGam.root",       "Data.C_GamGam");
  Analyze(inDir+"Data/data_D.GamGam.root",       "Data.D_GamGam");
}

void Analyze(TString inFile,TString sampleName){
  cout << "======================================================" << endl;
  cout << "Analyzing sample:  " << sampleName << endl;

  //
  // Setup TChain
  //
  TChain tree("mini");
  tree.Add(inFile);
  
  //
  // Setup TTreeReader
  //
  TTreeReader fReader;
  fReader.SetTree(&tree);

  TTreeReaderValue<Float_t> scaleFactor_PHOTON = {fReader, "scaleFactor_PHOTON"};
  TTreeReaderValue<Float_t> scaleFactor_PhotonTRIGGER = {fReader, "scaleFactor_PhotonTRIGGER"};
  TTreeReaderValue<Float_t> scaleFactor_PILEUP = {fReader, "scaleFactor_PILEUP"};
  TTreeReaderValue<Float_t> mcWeight = {fReader, "mcWeight"};
  TTreeReaderValue<Bool_t>  trigP = {fReader, "trigP"};
  
  TTreeReaderValue<UInt_t> photon_n  = {fReader, "photon_n"};
  TTreeReaderArray<float> photon_pt  = {fReader, "photon_pt"};
  TTreeReaderArray<float> photon_eta = {fReader, "photon_eta"};
  TTreeReaderArray<float> photon_phi = {fReader, "photon_phi"};
  TTreeReaderArray<float> photon_E   = {fReader, "photon_E"};
  TTreeReaderArray<float> photon_ptcone30 = {fReader, "photon_ptcone30"};
  TTreeReaderArray<float> photon_etcone20 = {fReader, "photon_etcone20"};
  TTreeReaderValue<vector<bool>> photon_isTightID = {fReader, "photon_isTightID"};
  
  uint nEventsInTree  = tree.GetEntries();
  int iEvt = 0;

  TH1F hist_mYY_bin1("hist_mYY_bin1","Diphoton invariant mass; m_{#gamma#gamma} [GeV];Events / bin", 30, 105, 160.);


  //
  // This is where we loop each event.
  //
  while (fReader.Next()) {
    if((iEvt)%10000 == 0){
      cout << iEvt << "/" << nEventsInTree << endl;
    }

    iEvt++;

    float scaleFactor = (*scaleFactor_PHOTON) * (*scaleFactor_PhotonTRIGGER) * (*scaleFactor_PILEUP);
    float weight = scaleFactor * (*mcWeight);

    if (sampleName.Contains("Data")){
      weight = 1.0f;
    }

    //
    // Skip event if it does not pass this Photon trigger
    //
    if (!(*trigP)) continue;
    
    //
    // Preselection of "Good" photons
    //

    std::vector<int> goodphoton_index; 
    int goodphoton_n = 0;
    
    for(unsigned int i=0; i < (*photon_n); i++){
      if(!(photon_isTightID->at(i))) continue; // photons are tight
      if( photon_pt[i] < 25000.) continue;
      if(!(TMath::Abs(photon_eta[i])<2.37 && (TMath::Abs(photon_eta[i]) < 1.37 || TMath::Abs(photon_eta[i]) > 1.52 ))) continue;
      goodphoton_n++;
      goodphoton_index.emplace_back(i);
    }

    //
    // If number of "Good" photons is not exactly 2, skip event.
    //
    if (goodphoton_n != 2) continue;

    int goodphoton1_index = goodphoton_index[0];
    int goodphoton2_index = goodphoton_index[1];
        
    // isolation requirement for leading photon
    if(!((photon_ptcone30[goodphoton1_index]/photon_pt[goodphoton1_index]) < 0.065)) continue;
    if(!((photon_etcone20[goodphoton1_index]/photon_pt[goodphoton1_index]) < 0.065)) continue; 
    
    // isolated requirement for sub-leading photon
    if(!((photon_ptcone30[goodphoton2_index]/photon_pt[goodphoton2_index]) < 0.065)) continue;
    if(!((photon_etcone20[goodphoton2_index]/photon_pt[goodphoton2_index]) < 0.065)) continue;

    // create 2 vectors          
    TLorentzVector Photon_1  = TLorentzVector();
    TLorentzVector Photon_2  = TLorentzVector();
    
    Photon_1.SetPtEtaPhiE(photon_pt[goodphoton1_index], photon_eta[goodphoton1_index], photon_phi[goodphoton1_index],photon_E[goodphoton1_index]);
    Photon_2.SetPtEtaPhiE(photon_pt[goodphoton2_index], photon_eta[goodphoton2_index], photon_phi[goodphoton2_index],photon_E[goodphoton2_index]);
    
    // calculate dPhi(photon-photon)
    float dPhi_yy = TMath::Abs(photon_phi[goodphoton1_index] - photon_phi[goodphoton2_index]);
    dPhi_yy       = dPhi_yy < TMath::Pi() ? dPhi_yy : 2*TMath::Pi() - dPhi_yy;
    
    // diphoton mass
    float m_yy  = sqrt( 2 * Photon_1.Pt()/1000. * Photon_2.Pt()/1000. * (cosh( Photon_1.Eta() - Photon_2.Eta()) - cos(dPhi_yy)));
    // kinematics
    float Photon_1_kin = Photon_1.Pt()/1000. / m_yy;
    float Photon_2_kin = Photon_2.Pt()/1000. / m_yy;

    if ( Photon_1_kin < 0.35) continue; 
    if ( Photon_2_kin < 0.25) continue;

    if(m_yy < 105) continue;
    if(m_yy > 160) continue;

    hist_mYY_bin1.Fill(m_yy, weight);
  }

  //
  // create output directory
  //
  system(("mkdir -p output")); 

  //
  // Save histogram in a root file
  //
  TFile outFile("output/Histo_"+sampleName+".root","recreate");
  outFile.cd();
  hist_mYY_bin1.Write();
  outFile.Close();
  cout << "" << endl;
}

