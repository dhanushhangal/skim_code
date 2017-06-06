#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2D.h"
#include "TF1.h"

#include "TH2F.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TCut.h>

#include "assert.h"
#include <fstream>
#include "TMath.h"
#include <vector>

using namespace std;


enum enum_dataset_types {e_Data2015,e_Data_pp,e_HydJet15,e_HydJet30,e_HydJet50, e_HydJet80, e_HydJet120,e_HydJet170,e_HydJet220,e_HydJet280, e_HydJet370,e_Pythia15,e_Pythia30,e_Pythia50, e_Pythia80, e_Pythia120,e_Pythia170,e_Pythia220,e_Pythia280, e_Pythia370, e_n_dataset_types};
TString dataset_type_strs[e_n_dataset_types] = {"Data2015","Data_pp","HydJet","HydJet30","HydJet50","HydJet80", "HydJet120", "HydJet170","HydJet220","HydJet280","HydJet370","Pythia15","Pythia30","Pythia50","Pythia80", "Pythia120", "Pythia170","Pythia220","Pythia280","Pythia370"};

int dataset_pthats[e_n_dataset_types+1] = {0,0,15,30,50,80,120,170,220,280,370,15,30,50,80,120,170,220,280,370,999};

int dataset_type_code = -999;

//arg 1 = which data set, arg 2 = output file number
void jet_skim(bool doCrab=0, int jobID=0, int endfile = 20, int dataset_type_code = 2, int output_file_num = 1)
{
    double pthatbins[10] = {15,30,50,80,120,170,220,280,370,9999};
	double xsecs[10] = {5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05, 1.018E-05, 2.477E-06, 6.160E-07, 1.088E-07, 0};
	double pthatEntries[9] = {0, 444180, 485654, 484515, 479933, 448829, 258987, 234835, 50644}; //5 TeV PbPb MC

	bool is_data = false;

	if(dataset_type_code == 0 || dataset_type_code == 1) is_data = true;

	cout << "dataset code: " << dataset_type_code << endl;
	
	bool do_PbPb=1, do_pp_tracking=0;

	int radius = 4;
	
	if(dataset_type_code== 1 || dataset_type_code > 10){do_PbPb = 0;   do_pp_tracking = 1;}

	cout<<"do_PbPb = "<<do_PbPb<<endl;
	cout<<"do_pp_tracking = "<<do_pp_tracking<<endl;


	//////////###### centrality Bins ###########///////////////

	TTree *inp_tree;
	TTree *inp_tree2;
	TTree *inp_tree3;
	TTree *inp_tree4;
	TTree *inp_tree5;
	TTree *inp_tree6;
	TTree *inp_tree_CS;

	string in_file_name;

	if(doCrab){
		in_file_name = Form("job_input_file_list_%d.txt",jobID);
	}
	else if(is_data&&!do_PbPb){
		in_file_name = "ppdata_filelist.txt";
	}else if(is_data&&do_PbPb){
		in_file_name = "PbPbData_2015_lxplusskim_MartaAnalysis.txt";
	}else if(dataset_type_code > 10){
		in_file_name = "ppMC_pthat80_filelist.txt";
	}else if(dataset_type_code > 1&&dataset_type_code <11){
		in_file_name = "PbPb_MC_filelist_withCSfix.txt";
	}else{
		cerr<<"need to set up to run on that sample..."<<endl;
	}

	cout << "trying a filelist named "<< in_file_name << endl;

	//MC
	TString output_file_base = "./";

	output_file_base +=dataset_type_strs[dataset_type_code];

	TString output_file_extension = "";   
	//output_file_extension += output_file_num;   
	output_file_extension += ".root";
	TFile *output_file = new TFile((TString) (output_file_base+output_file_extension), "RECREATE");
	TTree *ftree = new TTree("unzipMixTree", "");
   
    double weight;    
	const int MAXPARTICLES = 100000;

	float calo_jteta, calo_jtphi, calo_jtpt;
	float calo_refpt;
	int calo_parton_flavor;

	Int_t HBHENoiseFilter, HBHENoiseFilterResultRun2Loose, eventSelection, pvFilter;
	Int_t HLT_Jet80, HLT_Jet100;
	Int_t pPAcollisionEventSelectionPA = -999;
	Int_t hiBin = -999;
	Float_t pthat = -999;
	Float_t vz = -999;
	Int_t nPFpart, nCSpart, calo_nref;
	Int_t nCSpart_perJet2_id1, nPFpart_perJet2_id1;

	ftree->Branch("vz", &vz);
	if(do_PbPb){
		ftree->Branch("hiBin", &hiBin);
	}

	ftree->Branch("jteta", &calo_jteta);
	ftree->Branch("jtphi", &calo_jtphi);
	ftree->Branch("jtpt", &calo_jtpt);
	
	if(!is_data){
	    ftree->Branch("refpt", &calo_refpt);
		ftree->Branch("pthat", &pthat);
     	ftree->Branch("refparton_flavor",&calo_parton_flavor);
	    ftree->Branch("weight", &weight);
	}

	ftree->Branch("nCScand",&nCSpart_perJet2_id1);
    ftree->Branch("nPFcand",&nPFpart_perJet2_id1);

    std::ifstream instr(in_file_name.c_str(), std::ifstream::in);
	if(!instr.is_open()) cout << "filelist not found!! Exiting..." << endl;
	std::string filename;
	int ifile=0;

	const int MAXJETS = 500;
	
	Float_t t_calo_jtpt[MAXJETS], t_calo_jteta[MAXJETS], t_calo_jtphi[MAXJETS];
	if(!is_data) Float_t t_calo_refpt[MAXJETS];
	if(!is_data) Int_t t_calo_parton_flavor[MAXJETS];	

	vector<float> *csCandEta=0, *csCandPhi=0, *csCandPt=0;
	vector<int> *csCandId=0;

	vector<float> *pfCandEta=0, *pfCandPhi=0, *pfCandPt=0;
    vector<int> *pfCandId=0;

/*
	while(instr>>filename && ifile<endfile){
		filename.erase(std::remove(filename.begin(), filename.end(), '"'), filename.end());
		filename.erase(std::remove(filename.begin(), filename.end(), ','), filename.end());
		filename.erase(std::remove(filename.begin(), filename.end(), '['), filename.end());
		filename.erase(std::remove(filename.begin(), filename.end(), ']'), filename.end());
		cout<<"File name is "<< filename <<endl;
		ifile++;

		//TFile *my_file = TFile::Open(filename.c_str());

		//if(!my_file){
			int pos = filename.find_first_of('s');
			string reducedfn = filename.substr(pos-1);
			string xrdPrefix = "root://cmsxrootd.fnal.gov//";
			cout << "local file not detected. Trying " << xrdPrefix+reducedfn << endl;
			TFile *my_file = TFile::Open((xrdPrefix+reducedfn).c_str());
			//TFile::Open((xrdPrefix+reducedfn).c_str());
		//}
*/
        TFile *my_file = TFile::Open("/data/forests/temp/HiForestAOD_80.root");

/*		
		if(!my_file){ cout << "File cannot be found!!" << endl; exit(1); }	

		if(my_file->IsZombie()) { 
			std::cout << "Is zombie" << std::endl;
		}    
*/
		if(do_PbPb){
			inp_tree = (TTree*)  my_file->Get(Form("akPu%dCaloJetAnalyzer/t",radius));
			inp_tree_CS = (TTree*) my_file->Get("pfcandAnalyzerCS/pfTree");
		}else{
			inp_tree = (TTree*)  my_file->Get(Form("ak%dCaloJetAnalyzer/t",radius));
			inp_tree_CS = (TTree*) my_file->Get("pfcandAnalyzerCS/pfTree");
		}

		inp_tree2 = (TTree*)  my_file->Get("pfcandAnalyzer/pfTree");
		if(!inp_tree2){ cout << "PFCand Tree not found!! Exiting..." << endl; exit(1); }

		inp_tree3 = (TTree*) my_file->Get("hiEvtAnalyzer/HiTree");
		if(!inp_tree3){ cout << "Evt Tree not found!! Exiting..." << endl; exit(1); }

		inp_tree4 = (TTree*) my_file->Get("skimanalysis/HltTree");
		if(!inp_tree4){ cout << "Skim Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree4);

		inp_tree6 = (TTree*) my_file->Get("hltanalysis/HltTree");
		if(is_data && !inp_tree6){ cout << "HLT Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree6);

		cout << "trees loaded" << endl;

		inp_tree3->SetBranchAddress("vz",&vz);
		if(do_PbPb){
			inp_tree3->SetBranchAddress("hiBin",&hiBin);
		}
		inp_tree->SetBranchAddress("nref",&calo_nref);
		inp_tree->SetBranchAddress("jtpt",t_calo_jtpt);
		inp_tree->SetBranchAddress("jteta",t_calo_jteta);
		inp_tree->SetBranchAddress("jtphi",t_calo_jtphi);
		
		if(!is_data){
    		inp_tree->SetBranchAddress("refpt",t_calo_refpt);
	    	inp_tree->SetBranchAddress("refparton_flavor",t_calo_parton_flavor);
    		inp_tree3->SetBranchAddress("pthat",&pthat);
        }
		//Run2
		inp_tree2->SetBranchAddress("nPFpart",&nPFpart);
		inp_tree2->SetBranchAddress("pfPt",&pfCandPt);
                inp_tree2->SetBranchAddress("pfEta",&pfCandEta);
                inp_tree2->SetBranchAddress("pfPhi",&pfCandPhi);
                inp_tree2->SetBranchAddress("pfId",&pfCandId);
                
		if(do_PbPb){
			inp_tree_CS->SetBranchAddress("nPFpart",&nCSpart);
			inp_tree_CS->SetBranchAddress("pfPt",&csCandPt);
			inp_tree_CS->SetBranchAddress("pfEta",&csCandEta);
			inp_tree_CS->SetBranchAddress("pfPhi",&csCandPhi);
			inp_tree_CS->SetBranchAddress("pfId",&csCandId);
		}


		inp_tree->SetBranchAddress("HBHENoiseFilterResult",&HBHENoiseFilter);
		inp_tree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&HBHENoiseFilterResultRun2Loose);

		if(!do_PbPb){
			inp_tree->SetBranchAddress("pprimaryVertexFilter",&pvFilter);
		}
		else{	
			inp_tree->SetBranchAddress("pcollisionEventSelection",&eventSelection);			
			inp_tree->SetBranchAddress("pprimaryVertexFilter",&pvFilter);
		}

		if(is_data){
			if(!do_PbPb){ 
				inp_tree->SetBranchAddress("HLT_AK4PFJet80_Eta5p1_v1", &HLT_Jet80);
				inp_tree->SetBranchAddress("HLT_AK4PFJet100_Eta5p1_v1", &HLT_Jet100);
			}
			else{
				inp_tree->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_v1", &HLT_Jet80);
				inp_tree->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_v1", &HLT_Jet100);
				inp_tree->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_v1_Prescl", &HLT_Jet80_ps);
				inp_tree->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_v1_Prescl", &HLT_Jet100_ps);
			}
		}


		int n_evt = inp_tree->GetEntriesFast();

		//n_evt=2;
		cout << "Entries: "<< n_evt << endl;
		for(int evi = 0; evi < n_evt; evi++) {

			inp_tree->GetEntry(evi);
			inp_tree2->GetEntry(evi);
			inp_tree3->GetEntry(evi);
			if(do_PbPb) inp_tree_CS->GetEntry(evi);

			if(is_data && !HLT_Jet80 && !HLT_Jet100) continue;
			if(!do_PbPb && (!pvFilter)) continue;

			if(do_PbPb && (!pvFilter || !HBHENoiseFilter || !eventSelection)) continue;

            nCSpart_perJet2_id1 = 0;
            nPFpart_perJet2_id1 = 0;

			//start calo jet loop
			for(int j4i = 0; j4i < calo_nref ; j4i++) {

				if( fabs(t_calo_jteta[j4i]) > 1.6 || t_calo_jtpt[j4i] < 50. || pthat<50) continue;

                calo_jtpt = t_calo_jtpt[j4i];
                calo_jteta = t_calo_jteta[j4i];
                calo_jtphi = t_calo_jtphi[j4i];
                calo_refpt = t_calo_refpt[j4i];
                calo_parton_flavor = t_calo_parton_flavor[j4i];
				
                //ncs cand for PbPb and P+H
				int ncs2_id1 = 0;
				if(do_PbPb){
					for(int icand=0; icand<nCSpart; icand++){
					
						if(abs(csCandEta->at(icand))>2.4) continue;
						double dr = sqrt(pow(csCandEta->at(icand) - calo_jteta, 2) + pow(acos(cos(calo_jtphi - csCandPhi->at(icand))),2));
						if(dr < 0.4 && csCandPt->at(icand)>2 && csCandId->at(icand)==1){
                          nCSpart_perJet2_id1++;
						} 
					}
				}

                //npf cand for PbPb and P+H
				int npf2_id1 = 0;
				if(do_PbPb){
					for(int icand=0; icand<nPFpart; icand++){

						if(abs(pfCandEta->at(icand))>2.4) continue;
						double dr = sqrt(pow(pfCandEta->at(icand) - calo_jteta, 2) + pow(acos(cos(calo_jtphi - pfCandPhi->at(icand))),2));
						if(dr < 0.4 && pfCandPt->at(icand)>2 && pfCandId->at(icand)==1){
						  nPFpart_perJet2_id1++; 
						}
					}
				}
				
     			int ibin=0;
                while(pthat>pthatbins[ibin+1]) ibin++;
                weight = (xsecs[ibin]-xsecs[ibin+1])/pthatEntries[ibin];
                     
				///// Fill it
		    	ftree->Fill();
			
			} /// calo jet loop

		}  ///event loop

		my_file->Close();
	//}

	cout<<"writing"<<endl;

	output_file->cd();
	ftree->Write();

	output_file->Close();

	cout<<"done"<<endl;

}