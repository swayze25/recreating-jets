//JETS (and subjets)!

#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "TH1.h"
#include <TF2.h>
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TCanvas.h"
#include <TGraph2D.h>
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include <TLorentzVector.h>
#include <iostream> 
#include <sstream>
#include <vector> 
#include <cstdio>

using namespace fastjet;
using namespace std;

void jet_construction (const vector<fastjet::PseudoJet> &);
double EtMin=4.0, R=0.4, rapmin=-1.0, rapmax=2.0;
double ycut = 0.01, dcut;
fastjet::Strategy strategy = fastjet::Best;
fastjet::RecombinationScheme recomb_scheme=fastjet::E_scheme;
fastjet::JetAlgorithm jet_alg = fastjet::antikt_algorithm;
fastjet::JetDefinition jet_def(jet_alg, R, recomb_scheme, strategy);

void print_jets (const vector<fastjet::PseudoJet> &, int flag);

//Reading  the input events from the root file.
TFile *myFile = new TFile("eventdata.root","READ");
TTreeReader myReader("events", myFile);
TTreeReaderArray<Float_t> px(myReader, "px");
TTreeReaderArray<Float_t> py(myReader, "py");
TTreeReaderArray<Float_t> pz(myReader, "pz");
TTreeReaderArray<Float_t> e(myReader, "e");


//Writing an Output ROOT file:
TFile *jets = new TFile("jets.root","recreate");
std::vector<Float_t> pxj,pyj,pzj,ej,njets;
TTree *jet = new TTree("jet", "");
TH1F *nJets = new TH1F("nJets","#Jets; Number of Jets; Events",10,0,10);
std::vector<Float_t> pxsj,pysj,pzsj,esj,nsjets;
TTree *subjet = new TTree("subjet", "");
TH1F *nsJets = new TH1F("nsubJets","#subJets; Number of subJets; Events",10,0,10);

int ev=0, n=0, sn=0;

int main(){

	vector<fastjet::PseudoJet> input_particles;
	
	//Jets:
	jet->Branch("px_jets", "std::vector<Float_t>",  &pxj);
	jet->Branch("py_jets", "std::vector<Float_t>",  &pyj);
	jet->Branch("pz_jets", "std::vector<Float_t>",  &pzj);
	jet->Branch("E_jets", "std::vector<Float_t>",  &ej);
	jet->Branch("njets", "std::vector<Float_t>",  &njets);
	
	//SubJets:
	subjet->Branch("px_subjets", "std::vector<Float_t>",  &pxsj);
	subjet->Branch("py_subjets", "std::vector<Float_t>",  &pysj);
	subjet->Branch("pz_subjets", "std::vector<Float_t>",  &pzsj);
	subjet->Branch("E_subjets", "std::vector<Float_t>",  &esj);
	subjet->Branch("nsjets", "std::vector<Float_t>",  &nsjets);
	
	// a "header" for the output
  cout << "Ran " << jet_def.description() << endl;
	cout << "Minimum Transverse Energy: " << EtMin << " GeV" << endl;
	cout << "Rapidity Range: " <<  rapmin  << " to " << rapmax << endl;
	cout << "Showing the subjets with ycut: " << ycut << endl << endl;
	
	while (myReader.Next()){
		ev++;
		//------------------------------------INPUT PARTICLES------------------------------------
		for (int x = 0; x < px.GetSize(); x++){
			input_particles.push_back(fastjet::PseudoJet(px[x],py[x],pz[x],e[x]));
		}
		jet_construction(input_particles);			
		input_particles.clear();
  }
	
  cout << "\nNumber of Events Read : " << ev << endl;
  cout << "\n#Jets formed : " << 	n << endl;
  cout << "		Etmin : " << EtMin << " GeV" << endl;
  cout << "		Jet Radius: " << R << endl;
  cout << "\n#Sub-Jets formed : " << sn << endl;
  cout << "		ycut : " << ycut << endl;
  
  jets->Write();
	jets->Close();
	delete jets;
}

//-----------------------------------------JETS-----------------------------------------
void jet_construction (const vector<fastjet::PseudoJet> & input_particles){

  fastjet::JetDefinition jet_def(kt_algorithm,R);
	fastjet::ClusterSequence clust_seq(input_particles, jet_def);
	Selector jet_selector = SelectorEtMin(EtMin) && SelectorRapRange(rapmin, rapmax);  
	vector<PseudoJet> inclusive_jets = sorted_by_pt(jet_selector(clust_seq.inclusive_jets()));
	
	if(inclusive_jets.size()!=0){
		njets.push_back(inclusive_jets.size());
		nJets->Fill(inclusive_jets.size());
		n+=inclusive_jets.size();
		//Print Jets:
		cout << "\n";
		cout << "===========================================================================\n";
		cout << "								Printing Jets and Subjets for Event "<< ev << "\n";
		cout << "===========================================================================\n";
		print_jets(inclusive_jets,1);
	}

	for (int i = 0; i < inclusive_jets.size(); i++){
		pxj.push_back(inclusive_jets[i].px());
		pyj.push_back(inclusive_jets[i].py());
		pzj.push_back(inclusive_jets[i].pz());
		ej.push_back(inclusive_jets[i].e());
		int jconstituents = inclusive_jets[i].constituents().size();
		
//--------------------------------------SUBJETS--------------------------------------
		dcut =  ycut*pow(inclusive_jets[i].Et(),2);
		vector<PseudoJet> subJets = inclusive_jets[i].exclusive_subjets(dcut);
		
		if (subJets.size() != 0){
			nsjets.push_back(subJets.size());
			nsJets->Fill(subJets.size());
			sn+=subJets.size();
			//Print Subjets:
			cout << "\nPrinting subjets with dcut = " << dcut << " GeV^2 || " << "Jet Constituents: " << jconstituents << "\n";
			cout << "---------------------------------------------------------------------------\n";
			print_jets(subJets,2);
		}
		
		for (int j=0; j<subJets.size(); j++){
			pxsj.push_back(subJets[j].px());
			pysj.push_back(subJets[j].py());
			pzsj.push_back(subJets[j].pz());
			esj.push_back(subJets[j].e());
		}
	}
	
	jet->Fill(); subjet->Fill();
	pxj.clear(); pyj.clear(); pzj.clear(); ej.clear(); njets.clear();
	pxsj.clear(); pysj.clear(); pzsj.clear(); esj.clear(); nsjets.clear();
}


//------------------------------------PRINT JETS------------------------------------
void print_jets (const vector<fastjet::PseudoJet> & jets, int flag) {
  // sort jets into increasing pt
  vector<fastjet::PseudoJet> sorted_jets = sorted_by_pt(jets);  

  // label the columns
  if (flag == 1)
  printf("%5s %15s %15s %15s %15s\n","Jet", "rapidity", "phi", "pt", "constituents");
	if (flag == 2)
	printf("%5s %15s %15s %15s %15s\n","sJet", "rapidity", "phi", "pt", "constituents");
  
  // print out the details for each jet
  for (unsigned int i = 0; i < sorted_jets.size(); i++) {
    int n_constituents = sorted_jets[i].constituents().size();
    printf("%5u %15.8f %15.8f %15.8f %8u\n",
	   i, sorted_jets[i].rap(), sorted_jets[i].phi(),
	   sorted_jets[i].perp(), n_constituents);
  }
  cout << endl;
}

