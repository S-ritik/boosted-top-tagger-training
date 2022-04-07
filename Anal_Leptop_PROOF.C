/*
  weight = event_weight; 
  should be in the same ntuple, event by event we are going to read it.
   
  Order lepton according to the matching of leading Ak8 jets. 
  Lepton tagging
  
*/

#define Anal_Leptop_PROOF_cxx
//#include "Anal_Leptop_PROOF.h"
#include "getobjects.h"

#include <TH2.h>
#include <TStyle.h>
#include <TVector3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <fstream>
#include <TProofOutputFile.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <TProofServ.h>


void Anal_Leptop_PROOF::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();

}
 
void Anal_Leptop_PROOF::SlaveBegin(TTree * /*tree*/)
{
  //The SlaveBegin() function is called after the Begin() function.
  //When running with PROOF SlaveBegin() is called on each slave server.
  //The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
  
  OutFile = new TProofOutputFile("test_output.root");
  
  fileOut = OutFile->OpenFile("RECREATE");
  if ( !(fileOut = OutFile->OpenFile("RECREATE")) )
    {
      Warning("SlaveBegin", "problems opening file: %s/%s",OutFile->GetDir(), OutFile->GetFileName());
    }
   
  isTT = true;
  isbQCD = false;

  
  Tout = new TTree("leptop_e","leptop_e");
  Tout->Branch("event_pt_weight",&event_pt_weight,"event_pt_weight/F");
  Tout->Branch("weight",&weight,"weight/F");

  Tout->Branch("selpfjetAK8NHadF",&selpfjetAK8NHadF,"selpfjetAK8NHadF/F");
  Tout->Branch("selpfjetAK8neunhadfrac",&selpfjetAK8neunhadfrac,"selpfjetAK8neunhadfrac/F");
  Tout->Branch("selpfjetAK8subhaddiff",&selpfjetAK8subhaddiff,"selpfjetAK8subhaddiff/F");
  Tout->Branch("selpfjetAK8tau21",&selpfjetAK8tau21,"selpfjetAK8tau21/F");
  Tout->Branch("selpfjetAK8chrad",&selpfjetAK8chrad,"selpfjetAK8chrad/F");
  Tout->Branch("selpfjetAK8sdmass",&selpfjetAK8sdmass,"selpfjetAK8sdmass/F");

  Tout->Branch("selpfjetAK8matchedeldxy_sv",&selpfjetAK8matchedeldxy_sv,"selpfjetAK8matchedeldxy_sv/F");
  Tout->Branch("selpfjetAK8matchedelpt",&selpfjetAK8matchedelpt,"selpfjetAK8matchedelpt/F");
  Tout->Branch("selpfjetAK8matchedelcleta",&selpfjetAK8matchedelcleta,"selpfjetAK8matchedelcleta/F");
  Tout->Branch("selpfjetAK8matchedelsigmaieta",&selpfjetAK8matchedelsigmaieta,"selpfjetAK8matchedelsigmaieta/F");
  Tout->Branch("selpfjetAK8matchedelsigmaiphi",&selpfjetAK8matchedelsigmaiphi,"selpfjetAK8matchedelsigmaiphi/F");

  Tout->Branch("selpfjetAK8matchedelr9full",&selpfjetAK8matchedelr9full,"selpfjetAK8matchedelr9full/F");
  Tout->Branch("selpfjetAK8matchedelsupcl_etaw",&selpfjetAK8matchedelsupcl_etaw,"selpfjetAK8matchedelsupcl_etaw/F");
  Tout->Branch("selpfjetAK8matchedelsupcl_phiw",&selpfjetAK8matchedelsupcl_phiw,"selpfjetAK8matchedelsupcl_phiw/F");
  Tout->Branch("selpfjetAK8matchedelhcaloverecal",&selpfjetAK8matchedelhcaloverecal,"selpfjetAK8matchedelhcaloverecal/F");
  Tout->Branch("selpfjetAK8matchedelcloctftrkn",&selpfjetAK8matchedelcloctftrkn,"selpfjetAK8matchedelcloctftrkn/F");
  Tout->Branch("selpfjetAK8matchedelcloctftrkchi2",&selpfjetAK8matchedelcloctftrkchi2,"selpfjetAK8matchedelcloctftrkchi2/F");
  Tout->Branch("selpfjetAK8matchedele1x5bye5x5",&selpfjetAK8matchedele1x5bye5x5,"selpfjetAK8matchedele1x5bye5x5/F");
  Tout->Branch("selpfjetAK8matchedelnormchi2",&selpfjetAK8matchedelnormchi2,"selpfjetAK8matchedelnormchi2/F");
  Tout->Branch("selpfjetAK8matchedelhitsmiss",&selpfjetAK8matchedelhitsmiss,"selpfjetAK8matchedelhitsmiss/F");
  Tout->Branch("selpfjetAK8matchedeltrkmeasure",&selpfjetAK8matchedeltrkmeasure,"selpfjetAK8matchedeltrkmeasure/F");
 
  Tout->Branch("selpfjetAK8matchedelecloverpout",&selpfjetAK8matchedelecloverpout,"selpfjetAK8matchedelecloverpout/F");
  Tout->Branch("selpfjetAK8matchedelecaletrkmomentum",&selpfjetAK8matchedelecaletrkmomentum,"selpfjetAK8matchedelecaletrkmomentum/F");
  Tout->Branch("selpfjetAK8matchedeldeltaetacltrkcalo",&selpfjetAK8matchedeldeltaetacltrkcalo,"selpfjetAK8matchedeldeltaetacltrkcalo/F");
  Tout->Branch("selpfjetAK8matchedelsupcl_preshvsrawe",&selpfjetAK8matchedelsupcl_preshvsrawe,"selpfjetAK8matchedelsupcl_preshvsrawe/F");
  Tout->Branch("selpfjetAK8matchedelpfisolsumphet",&selpfjetAK8matchedelpfisolsumphet,"selpfjetAK8matchedelpfisolsumphet/F");
  Tout->Branch("selpfjetAK8matchedelpfisolsumchhadpt",&selpfjetAK8matchedelpfisolsumchhadpt,"selpfjetAK8matchedelpfisolsumchhadpt/F");
  Tout->Branch("selpfjetAK8matchedelpfisolsumneuhadet",&selpfjetAK8matchedelpfisolsumneuhadet,"selpfjetAK8matchedelpfisolsumneuhadet/F");
  Tout->Branch("selpfjetAK8matchedeletain",&selpfjetAK8matchedeletain,"selpfjetAK8matchedeletain/F");
  Tout->Branch("selpfjetAK8matchedelphiin",&selpfjetAK8matchedelphiin,"selpfjetAK8matchedelphiin/F");
  Tout->Branch("selpfjetAK8matchedelfbrem",&selpfjetAK8matchedelfbrem,"selpfjetAK8matchedelfbrem/F");
  Tout->Branch("selpfjetAK8matchedeleoverp",&selpfjetAK8matchedeleoverp,"selpfjetAK8matchedeleoverp/F");
  Tout->Branch("selpfjetAK8matchedelhovere",&selpfjetAK8matchedelhovere,"selpfjetAK8matchedelhovere/F");
  Tout->Branch("selpfjetAK8matchedelRho", &selpfjetAK8matchedelRho,"selpfjetAK8matchedelRho/F");
  Tout->Branch("selpfjetAK8matchedelptrel", &selpfjetAK8matchedelptrel,"selpfjetAK8matchedelptrel/F");
  Tout->Branch("matched", &matched,"matched/O");
  Tout->Branch("matchel", &matchel,"matchel/I");
    
 
 
  
  char name[1000];

}

Bool_t Anal_Leptop_PROOF::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either Anal_Leptop_PROOF::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  GetEntry(entry);
  if(isMC){
    weight = event_weight;
  }else{
    weight = 1;
  }


  //Here you get electrons with your criteria
  vector <Electron> velectrons;
  getelectrons(velectrons,30,absetacut);    // electron_pt_cut & absetacut defined in Proof.h
  
  if((int)velectrons.size() < 1)
    return kFALSE;
  
  //  vector <Muon> vmuons;
  //  getmuons(vmuons,0,absetacut);
  //  getLeptons(vleptons,vmuons,velectrons,0);
  
  //Here you get AK4 jets with your criteria                                     
  //vector <AK4Jet> Jets;
  //getAK4jets(Jets,AK4jet_pt_cut,absetacut,isMC);
   
   
  //Here you get AK8 jets with your criteria 
  vector <AK8Jet> LJets;
  getAK8jets(LJets,300,absetacut,isMC);
  
  if((int)LJets.size() < 1)
    return kFALSE;
  // Get generator-level particles //
       
  vector<GenParton> genpartons;
  vector<GenParton> LHEtops;
  vector<TopQuark> gentops;
  
  vector<LHEparticle> lheparticles;
  getPartons(genpartons);
  getLHEParticles(lheparticles);

  if(isTT){
    // Get GEN-level top quarks                                                         
    getLHETops(LHEtops,genpartons); // before shower (get original top quarks which have decayed) --> will be usedto derive top pt reweighting                                   

    getGENTops(gentops,genpartons); // after shower (get top quarks from its daughters) --> will tell details about the signature of ttbar events at GEN level
    
    // top pt reweighting //
    float toppt_wt = 1;
    if(LHEtops.size()==2){
      toppt_wt = SF_TOP(0.0615,0.0005,TMath::Min(float(500),float(LHEtops[0].pt)),TMath::Min(float(500),float(LHEtops[1].pt)));
      weight *= toppt_wt;
    }
  }

  //Get RECO-level objects //
  //First sorted out electron and muon2
  
  
  //Get index of AK8 jet nearest to each lepton                                  
  //Get indices of nearest lepton, AK4 jet for each AK8 jet                        

  for(int ijet=0; ijet < int(LJets.size()); ijet++)
    {
      //Match lepton with AK8 jets
      float maxpt = -100;
      matchel = -100;
      for(int ie=0; ie < int(velectrons.size()); ie++)
	{ 
	  float dr = delta2R(velectrons[ie].p4,LJets[ijet].p4);
	  if(dr <0.7 && velectrons[ie].pt > maxpt)
	    {
	      matchel = ie;
	      maxpt = velectrons[ie].pt;
	    }
	}
      if(matchel < -1)
	continue;
      
      matched = false;
      if(isTT)
	{
	  for(int it=0; it< (int)gentops.size(); it++)
	    {
	      if(abs(gentops[it].daughter[0].pdgId) == 11)
		if(delta2R(LJets[ijet].p4,gentops[it].daughter[0].p4) < 0.6 && delta2R(LJets[ijet].p4,gentops[it].daughter[2].p4) < 0.6)
		  {
		    matched = true;
		    break;
		  }
	    }
	}
      
      else if(isbQCD)
	{
	  for(int ib=0; ib < (int)genpartons.size(); ib++)
	    {
	      if(abs(genpartons[ib].pdgId) == 5 && genpartons[ib].fromhard)
		{
		  if(delta2R(LJets[ijet].p4,genpartons[ib].p4) < 0.6)
		    {
		      matched = true;
		      break;
		    }
		}
	    }
	}
      
      if(!matched)
	continue;
	
      selpfjetAK8NHadF = LJets[ijet].NHadF;
      selpfjetAK8neunhadfrac = LJets[ijet].neunhadfrac;
      selpfjetAK8subhaddiff = LJets[ijet].subhaddiff;
      selpfjetAK8tau21 = LJets[ijet].tau21;
      selpfjetAK8chrad = LJets[ijet].chrad;
      selpfjetAK8sdmass = LJets[ijet].sdmass;
	
      int iel = matchel;

      selpfjetAK8matchedelcleta = velectrons[iel].supcl_eta; //elsupcl_eta[nearest];
      selpfjetAK8matchedelpt = fabs(velectrons[iel].pt); //elpt[nearest]);
      selpfjetAK8matchedelsigmaieta = velectrons[iel].sigmaieta; // elsigmaieta[nearest];
      selpfjetAK8matchedelsigmaiphi = velectrons[iel].sigmaiphi; // elsigmaiphi[nearest];
      selpfjetAK8matchedelr9full = velectrons[iel].r9full; // elr9full[nearest];
      selpfjetAK8matchedelsupcl_etaw = velectrons[iel].supcl_etaw; // elsupcl_etaw[nearest];
      selpfjetAK8matchedelsupcl_phiw = velectrons[iel].supcl_phiw; // elsupcl_phiw[nearest];
      selpfjetAK8matchedelhcaloverecal = velectrons[iel].hcaloverecal; // elhcaloverecal[nearest];
      selpfjetAK8matchedelcloctftrkn = velectrons[iel].cloctftrkn; // elcloctftrkn[nearest];
      selpfjetAK8matchedelcloctftrkchi2 = velectrons[iel].cloctftrkchi2; // elcloctftrkchi2[nearest];
      selpfjetAK8matchedele1x5bye5x5 = velectrons[iel].e1x5bye5x5; // ele1x5bye5x5[nearest];
      selpfjetAK8matchedelnormchi2 = velectrons[iel].normchi2; // elnormchi2[nearest];
      selpfjetAK8matchedelhitsmiss = velectrons[iel].hitsmiss; // elhitsmiss[iel];
      selpfjetAK8matchedeltrkmeasure = velectrons[iel].trkmeasure; // eltrkmeasure[nearest];
      selpfjetAK8matchedelecloverpout = velectrons[iel].ecloverpout; // elecloverpout[nearest];
      selpfjetAK8matchedelecaletrkmomentum = velectrons[iel].ecaletrkmomentum; //elecaletrkmomentum[nearest];
      selpfjetAK8matchedeldeltaetacltrkcalo = velectrons[iel].deltaetacltrkcalo; // eldeltaetacltrkcalo[nearest];
      selpfjetAK8matchedelsupcl_preshvsrawe = velectrons[iel].supcl_preshvsrawe; // elsupcl_preshvsrawe[nearest];
      selpfjetAK8matchedelpfisolsumphet = velectrons[iel].pfisolsumphet; // elpfisolsumphet[nearest];
      selpfjetAK8matchedelpfisolsumchhadpt = velectrons[iel].pfisolsumchhadpt; // elpfisolsumchhadpt[nearest];
      selpfjetAK8matchedelpfisolsumneuhadet = velectrons[iel].pfsiolsumneuhadet; //elpfsiolsumneuhadet[nearest];
      selpfjetAK8matchedeletain = velectrons[iel].etain; // eletain[nearest];
      selpfjetAK8matchedelphiin = velectrons[iel].phiin; // elphiin[nearest];
      selpfjetAK8matchedelfbrem = velectrons[iel].fbrem; //elfbrem[nearest];
      selpfjetAK8matchedeleoverp = velectrons[iel].eoverp; // eleoverp[nearest];
      selpfjetAK8matchedelhovere = velectrons[iel].hovere; //elhovere[nearest];
      selpfjetAK8matchedelRho = Rho;
      selpfjetAK8matchedeldxy_sv = velectrons[iel].eldxy_sv;
      //selpfjetAK8matchedelptrel

      if(selpfjetAK8sdmass <0 || isnan(selpfjetAK8NHadF) || isnan(selpfjetAK8neunhadfrac) || isnan(selpfjetAK8subhaddiff) || isnan(selpfjetAK8tau21) || isnan(selpfjetAK8chrad) || isnan(selpfjetAK8sdmass) || isnan(selpfjetAK8matchedelcleta) || isnan(selpfjetAK8matchedelpt) || isnan(selpfjetAK8matchedelsigmaieta) || isnan(selpfjetAK8matchedelsigmaiphi) || isnan(selpfjetAK8matchedelr9full) || isnan(selpfjetAK8matchedelsupcl_etaw) || isnan(selpfjetAK8matchedelhcaloverecal) || isnan(selpfjetAK8matchedelcloctftrkn) || isnan(selpfjetAK8matchedelcloctftrkchi2) || isnan(selpfjetAK8matchedele1x5bye5x5) || isnan(selpfjetAK8matchedelnormchi2) || isnan(selpfjetAK8matchedelnormchi2) || isnan(selpfjetAK8matchedelhitsmiss) || isnan(selpfjetAK8matchedeltrkmeasure) || isnan(selpfjetAK8matchedelecloverpout) || isnan(selpfjetAK8matchedelecaletrkmomentum) || isnan(selpfjetAK8matchedeldeltaetacltrkcalo) || isnan(selpfjetAK8matchedelsupcl_preshvsrawe) || isnan(selpfjetAK8matchedelpfisolsumphet) || isnan(selpfjetAK8matchedelpfisolsumchhadpt) || isnan(selpfjetAK8matchedelpfisolsumneuhadet) || isnan(selpfjetAK8matchedeletain) || isnan(selpfjetAK8matchedelphiin) || isnan(selpfjetAK8matchedelfbrem) || isnan(selpfjetAK8matchedeleoverp) || isnan(selpfjetAK8matchedelhovere) || isnan(selpfjetAK8matchedelRho) || isnan(selpfjetAK8matchedeldxy_sv) || isnan(weight) )
	continue;
      Tout->Fill(); 
    }
  return kTRUE;     
}

void Anal_Leptop_PROOF::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  fileOut->cd();
  fileOut->Write();
  fOutput->Add(OutFile);
  fileOut->Close();
}

void Anal_Leptop_PROOF::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
}
