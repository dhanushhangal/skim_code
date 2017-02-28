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
TString dataset_type_strs[e_n_dataset_types] = {"Data2015","Data_pp","HydJet15","HydJet30","HydJet50","HydJet80", "HydJet120", "HydJet170","HydJet220","HydJet280","HydJet370","Pythia15","Pythia30","Pythia50","Pythia80", "Pythia120", "Pythia170","Pythia220","Pythia280","Pythia370"};

int dataset_pthats[e_n_dataset_types+1] = {0,0,15,30,50,80,120,170,220,280,370,15,30,50,80,120,170,220,280,370,999};

int dataset_type_code = -999;

void ReadFileList(std::vector<TString> &my_file_names, TString file_of_names, bool debug)
{
  ifstream file_stream(file_of_names);
  std::string line;
  my_file_names.clear();
  if( debug ) std::cout << "Open file " << file_of_names << " to extract files to run over" << std::endl;
  if( file_stream.is_open() ) {
    if( debug ) std::cout << "Opened " << file_of_names << " for reading" << std::endl;
    int line_num = 0;
    while( !file_stream.eof() ) {
      getline(file_stream, line);
      //if( debug ) std::cout << line_num << ": " << line << std::endl;
      TString tstring_line(line);
      if( tstring_line.CompareTo("", TString::kExact) != 0 ) my_file_names.push_back(tstring_line);
      line_num++;
    }
  } else {
    std::cout << "Error, could not open " << file_of_names << " for reading" << std::endl;
    assert(0);
  }
}

//arg 1 = which data set, arg 2 = output file number
void mini_skim(bool doCrab=1, int jobID=0, int endfile = 1, int dataset_type_code = 0, int output_file_num = 1)
{

	bool is_data = false;

	//if(dataset_type_code == 0 || dataset_type_code == 1) is_data = true;

	cout << "dataset code: " << dataset_type_code << endl;

	// assert(!is_data); //for now I'm interested in MC

	//-----------------------------
	// Set JFF-dependent corrections
	//-----------------------------

	float reco_eta, reco_phi, reco_pt;

	bool do_PbPb=1, do_pp_tracking=0;

	int radius = 4;

	if(dataset_type_code== 1 || dataset_type_code > 10){do_PbPb = 0;   do_pp_tracking = 1;}

	cout<<"do_PbPb = "<<do_PbPb<<endl;
	cout<<"do_pp_tracking = "<<do_pp_tracking<<endl;

	vector<float> corr_pt; 

	//--------------------------------

	//////////###### centrality Bins ###########///////////////

	TTree *inp_tree;
	TTree *inp_tree2;
	TTree *inp_tree3;
	TTree *inp_tree4;
	TTree *inp_tree5;
	TTree *inp_tree6;
	TTree *inp_tree7=0;
	TTree *pftree;
	TTree *inp_tree_CS;

	string in_file_name;

	if(doCrab){
		in_file_name = Form("job_input_file_list_%d.txt",jobID);
		//in_file_name = "Pythia6+Hydjet.txt";
	}
	else if(is_data&&!do_PbPb){
		in_file_name = "pp5TeV_HighPtPD_Apr2016.txt";
	}else if(is_data&&do_PbPb){
		in_file_name = "PbPbData_2015_lxplusskim_MartaAnalysis.txt";
	}else if(dataset_type_code > 10){
		in_file_name = "pp_5TeV_MC_NovSkimWithCS.txt";
		//in_file_name = "Pythia6_2p76TeV_MergedList.txt";
	}else if(dataset_type_code > 1&&dataset_type_code <11){
		//in_file_name = "pythiaHydjet_2p76TeV_forest.txt";
		in_file_name = "Pythia6_Hydjet_PbPbMix_PurdueList.txt";
	}else{
		cerr<<"need to set up to run on that sample..."<<endl;
	}

	cout << "trying a filelist named "<< in_file_name << endl;

	//MC
	TString output_file_base = "./unzippedSkim_new";

	//output_file_base +=dataset_type_strs[dataset_type_code];

	TString output_file_extension = "";   
	//output_file_extension += output_file_num;   
	output_file_extension += ".root";
	TFile *output_file = new TFile((TString) (output_file_base+output_file_extension), "RECREATE");
	TTree *mixing_tree = new TTree("mixing_tree", "");

	const int MAXPARTICLES = 100000;

	vector<float> calo_jteta, calo_jtphi, calo_jtpt, calo_corrpt;
	vector<float> pf_refpt, calo_refpt;
	vector<bool> highPurity;
	vector<int> calo_parton_flavor;

	Int_t HBHENoiseFilter, HBHENoiseFilterResultRun1, HBHENoiseFilterResultRun2Loose, HBHENoiseFilterResultRun2Tight, HBHEIsoNoiseFilterResult;
	Int_t phfCoincFilter1, phfCoincFilter2, phfCoincFilter3, phfCoincFilter4, phfCoincFilter5;
	Int_t eventSelection, pvFilter;
	Int_t HLT_Jet80, HLT_Jet100, HLT_Jet80_ps, HLT_Jet100_ps;
	Int_t pPAcollisionEventSelectionPA = -999;
	Int_t beamScrapeFilter;
	Int_t pClusterCompatibilityFilter, pVertexFilterCutGplus;
	Int_t hiBin = -999;
	Float_t evtPlane_HF2 = -999;
	Float_t pthat = -999;
	Float_t vz = -999;//, sumpt[15];
	Int_t nPFpart, nCSpart;
	vector<int> nCSpart_perJet, nCSpart_perJet0p5, nCSpart_perJet1, nCSpart_perJet2, nCSpart_perJet3, nCSpart_perJet0p5_id1, nCSpart_perJet1_id1, nCSpart_perJet2_id1, nCSpart_perJet3_id1;
	vector<float> calo_discr_ssvHighEff, calo_discr_ssvHighPur, calo_discr_csvV1, calo_discr_prob, calo_svtxm, calo_svtxpt, calo_svtxmcorr, calo_svtxdl, calo_svtxdls;

	/// higenparticles

	// Event skims commented out for Run1 data/mc
	mixing_tree->Branch("HLT_ak4CaloJet80", &HLT_Jet80);
	mixing_tree->Branch("HLT_ak4CaloJet100", &HLT_Jet100);
	mixing_tree->Branch("HLT_ak4CaloJet80_Prescale", &HLT_Jet80_ps);
	mixing_tree->Branch("HLT_ak4CaloJet100_Prescale", &HLT_Jet100_ps);

	mixing_tree->Branch("HBHENoiseFilterResult",&HBHENoiseFilter);
	mixing_tree->Branch("HBHENoiseFilterResultRun1",&HBHENoiseFilterResultRun1);
	mixing_tree->Branch("HBHENoiseFilterResultRun2Loose",&HBHENoiseFilterResultRun2Loose);
	mixing_tree->Branch("HBHENoiseFilterResultRun2Tight",&HBHENoiseFilterResultRun2Tight);
	mixing_tree->Branch("HBHEIsoNoiseFilterResult",&HBHEIsoNoiseFilterResult);

	if(!do_PbPb){
		mixing_tree->Branch("pPAprimaryVertexFilter",&pvFilter);
		//mixing_tree->Branch("pBeamScrapingFilter",&beamScrapeFilter);
		//mixing_tree->Branch("pVertexFilterCutGplus",&pVertexFilterCutGplus);
	}
	else{	
		mixing_tree->Branch("pcollisionEventSelection",&eventSelection);			
		mixing_tree->Branch("pprimaryVertexFilter",&pvFilter);
		mixing_tree->Branch("pclusterCompatibilityFilter",&pClusterCompatibilityFilter);
		mixing_tree->Branch("phfCoincFilter1",&phfCoincFilter1);
		mixing_tree->Branch("phfCoincFilter2",&phfCoincFilter2);
		mixing_tree->Branch("phfCoincFilter3",&phfCoincFilter3);
		mixing_tree->Branch("phfCoincFilter4",&phfCoincFilter4);
		mixing_tree->Branch("phfCoincFilter5",&phfCoincFilter5);
	}
	mixing_tree->Branch("highPurity", &highPurity);
	mixing_tree->Branch("vz", &vz);

	if(do_PbPb){
		mixing_tree->Branch("hiBin", &hiBin);
		mixing_tree->Branch("evtPlane_HF2",&evtPlane_HF2);
	}

	mixing_tree->Branch("calo_jteta", &calo_jteta);
	mixing_tree->Branch("calo_jtphi", &calo_jtphi);
	mixing_tree->Branch("calo_jtpt", &calo_jtpt);

	if(!is_data){
		mixing_tree->Branch("calo_refpt",&calo_refpt);
		mixing_tree->Branch("pthat", &pthat);
	}

	mixing_tree->Branch("nPFpart", &nPFpart, "nPFpart/I");
	mixing_tree->Branch("nCSpart", &nCSpart_perJet);
	mixing_tree->Branch("nCSpartGT0p5",&nCSpart_perJet0p5);
	mixing_tree->Branch("nCSpartGT1",&nCSpart_perJet1);
	mixing_tree->Branch("nCSpartGT2",&nCSpart_perJet2);
	mixing_tree->Branch("nCSpartGT3",&nCSpart_perJet3);
	mixing_tree->Branch("nCSpartGT0p5_id1",&nCSpart_perJet0p5_id1);
	mixing_tree->Branch("nCSpartGT1_id1",&nCSpart_perJet1_id1);
	mixing_tree->Branch("nCSpartGT2_id1",&nCSpart_perJet2_id1);
	mixing_tree->Branch("nCSpartGT3_id1",&nCSpart_perJet3_id1);
	
	if(!is_data){
		mixing_tree->Branch("calo_refparton_flavor",&calo_parton_flavor);
	}


	std::ifstream instr(in_file_name.c_str(), std::ifstream::in);
	if(!instr.is_open()) cout << "filelist not found!! Exiting..." << endl;
	std::string filename;
	int ifile=0;

	const int MAXJETS = 500;
	Float_t t_calo_jtpt[MAXJETS], t_calo_jteta[MAXJETS], t_calo_jtphi[MAXJETS], t_calo_discr_ssvHighEff[MAXJETS], t_calo_discr_ssvHighPur[MAXJETS], t_calo_discr_csvV1[MAXJETS], t_calo_discr_prob[MAXJETS], t_calo_svtxm[MAXJETS], t_calo_svtxpt[MAXJETS], t_calo_svtxmcorr[MAXJETS], t_calo_svtxdl[MAXJETS], t_calo_svtxdls[MAXJETS];
	Float_t t_calo_refpt[MAXJETS];

	Int_t t_calo_parton_flavor[MAXJETS];	
	Bool_t t_highPurity[MAXPARTICLES];

	Float_t t_hiEvtPlane[50];

	vector<float> *csCandEta=0, *csCandPhi=0, *csCandPt=0;
	vector<int> *csCandId=0;

	Int_t calo_nref;


    std::vector<TString> file_names;   file_names.clear();

    ReadFileList( file_names, in_file_name, true);

    

    for(int fi = 0; fi < (int) file_names.size(); fi++) {
    //for(int fi = 0; fi < 10; fi++) {
    
      TFile *my_file = TFile::Open(file_names.at(fi));
      //std::cout << "Current file: " << ", file_name: " << file_names.at(fi) << ", number " << fi << " of " << file_names.size() << std::endl;
      //if(my_file->IsZombie()) {
      //  std::cout << "Is zombie" << std::endl;
      //}
    cout<<"got file"<<endl;  
		if(do_PbPb){
			inp_tree = (TTree*)  my_file->Get(Form("akPu%dCaloJetAnalyzer/t",radius));
			inp_tree_CS = (TTree*) my_file->Get("pfcandAnalyzerCS/pfTree");
			pftree = (TTree*) my_file->Get(Form("akPu%dPFJetAnalyzer/t",radius));
		}else{
			inp_tree = (TTree*)  my_file->Get(Form("ak%dCaloJetAnalyzer/t",radius));
			inp_tree_CS = (TTree*) my_file->Get("pfcandAnalyzerCS/pfTree");
			pftree = (TTree*) my_file->Get(Form("ak%dPFJetAnalyzer/t",radius));
		}

		inp_tree2 = (TTree*)  my_file->Get("pfcandAnalyzer/pfTree");
		if(!inp_tree2){ cout << "PFCand Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree2);

		inp_tree3 = (TTree*) my_file->Get("hiEvtAnalyzer/HiTree");
		if(!inp_tree3){ cout << "Evt Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree3);

		inp_tree4 = (TTree*) my_file->Get("skimanalysis/HltTree");
		if(!inp_tree4){ cout << "Skim Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree4);

		if(do_PbPb){
			inp_tree5 = (TTree*) my_file->Get("anaTrack/trackTree");
		}else{
			inp_tree5 = (TTree*) my_file->Get("ppTrack/trackTree");
		}
		if(!inp_tree5){ cout << "track Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree5);

		inp_tree6 = (TTree*) my_file->Get("hltanalysis/HltTree");
		if(is_data && !inp_tree6){ cout << "HLT Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree6);

		if(!is_data){ 
			inp_tree7 = (TTree*) my_file->Get("HiGenParticleAna/hi");
			if(!inp_tree7){ cout << "GenPart Tree not found!! Exiting..." << endl; exit(1); }
			else inp_tree->AddFriend(inp_tree7);
		}

		cout << "trees loaded" << endl;

		inp_tree->SetBranchAddress("highPurity",t_highPurity);
		inp_tree->SetBranchAddress("vz",&vz);
		if(do_PbPb){
			inp_tree->SetBranchAddress("hiBin",&hiBin);
			inp_tree->SetBranchAddress("hiEvtPlanes",t_hiEvtPlane);
		}
		inp_tree->SetBranchAddress("nref",&calo_nref);
		inp_tree->SetBranchAddress("jtpt",t_calo_jtpt);
		inp_tree->SetBranchAddress("jteta",t_calo_jteta);
		inp_tree->SetBranchAddress("jtphi",t_calo_jtphi);
		if(!is_data) inp_tree->SetBranchAddress("refpt",t_calo_refpt);

		if(do_PbPb){
			inp_tree_CS->SetBranchAddress("nPFpart",&nCSpart);
			inp_tree_CS->SetBranchAddress("pfPt",&csCandPt);
			inp_tree_CS->SetBranchAddress("pfEta",&csCandEta);
			inp_tree_CS->SetBranchAddress("pfPhi",&csCandPhi);
			inp_tree_CS->SetBranchAddress("pfId",&csCandId);
		}

		inp_tree->SetBranchAddress("discr_ssvHighEff", t_calo_discr_ssvHighEff);
		inp_tree->SetBranchAddress("discr_ssvHighPur", t_calo_discr_ssvHighPur);
		if(do_PbPb || !is_data) inp_tree->SetBranchAddress("discr_csvV1", t_calo_discr_csvV1);
		else{
			inp_tree->SetBranchAddress("discr_csvSimple", t_calo_discr_csvV1);
		}
		if(!is_data){
			inp_tree->SetBranchAddress("refparton_flavor",t_calo_parton_flavor);
		}
		inp_tree->SetBranchAddress("discr_prob", t_calo_discr_prob);
		inp_tree->SetBranchAddress("svtxdl", t_calo_svtxdl);
		inp_tree->SetBranchAddress("svtxdls", t_calo_svtxdls);
		inp_tree->SetBranchAddress("svtxm", t_calo_svtxm);
		inp_tree->SetBranchAddress("svtxpt", t_calo_svtxpt);
		inp_tree->SetBranchAddress("svtxmcorr", t_calo_svtxmcorr);

		inp_tree->SetBranchAddress("HBHENoiseFilterResult",&HBHENoiseFilter);
		inp_tree->SetBranchAddress("HBHENoiseFilterResultRun1",&HBHENoiseFilterResultRun1);
		inp_tree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&HBHENoiseFilterResultRun2Loose);
		inp_tree->SetBranchAddress("HBHENoiseFilterResultRun2Tight",&HBHENoiseFilterResultRun2Tight);
		inp_tree->SetBranchAddress("HBHEIsoNoiseFilterResult",&HBHEIsoNoiseFilterResult);

		if(!do_PbPb){
			inp_tree->SetBranchAddress("pPAprimaryVertexFilter",&pvFilter); //Run2
			inp_tree->SetBranchAddress("pBeamScrapingFilter",&beamScrapeFilter); //Run2
			inp_tree->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus); //Run2
			//inp_tree->SetBranchAddress("pPAcollisionEventSelectionPA",&pvFilter); //Run1
		}
		else{	
			inp_tree->SetBranchAddress("pcollisionEventSelection",&eventSelection);			
			inp_tree->SetBranchAddress("pprimaryVertexFilter",&pvFilter);
			inp_tree->SetBranchAddress("pclusterCompatibilityFilter",&pClusterCompatibilityFilter);
			inp_tree->SetBranchAddress("phfCoincFilter1",&phfCoincFilter1);
			inp_tree->SetBranchAddress("phfCoincFilter2",&phfCoincFilter2);
			inp_tree->SetBranchAddress("phfCoincFilter3",&phfCoincFilter3);
			inp_tree->SetBranchAddress("phfCoincFilter4",&phfCoincFilter4);
			inp_tree->SetBranchAddress("phfCoincFilter5",&phfCoincFilter5);
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

		//n_evt=10000;
		cout << "Entries: "<< n_evt << endl;
		for(int evi = 0; evi < n_evt; evi++) {

			inp_tree->GetEntry(evi);
			if(do_PbPb) inp_tree_CS->GetEntry(evi);
			pftree->GetEntry(evi);

			if(is_data && !HLT_Jet80 && !HLT_Jet100) continue;
			if(!do_PbPb && (!pvFilter)) continue;

			if(do_PbPb) evtPlane_HF2 = t_hiEvtPlane[8]; 
			else evtPlane_HF2 = -999;

			// Removing event selection for PbPb skims so we can test various selection criteria
			if(do_PbPb && (!pvFilter || !HBHENoiseFilterResultRun2Loose || !eventSelection)) continue;

			//start calo jet loop
			for(int j4i = 0; j4i < calo_nref ; j4i++) {

				if( fabs(t_calo_jteta[j4i]) > 2. ) continue;

				//-----------------------------------------------------------------------------------------
				// Jet Energy Corrections (JFF-dependent, store final corrected values in vector corr_pt)
				//----------------------------------------------------------------------------------------


				reco_pt = t_calo_jtpt[j4i];
				reco_phi = t_calo_jtphi[j4i];
				reco_eta = t_calo_jteta[j4i];

				int ncs=0, ncs0p5=0, ncs1=0, ncs2=0, ncs3=0;
				int ncs0p5_id1=0, ncs1_id1=0, ncs2_id1=0, ncs3_id1=0;
				if(do_PbPb){
					for(int icand=0; icand<nCSpart; icand++){
						if(abs(csCandEta->at(icand))>2.4) continue;
						double dr = sqrt(pow(csCandEta->at(icand) - reco_eta, 2) + pow(acos(cos(reco_phi - csCandPhi->at(icand))),2));
						if(dr < 0.4){
							ncs++;
							if(csCandPt->at(icand)>0.5) ncs0p5++;
							if(csCandPt->at(icand)>1) ncs1++;
							if(csCandPt->at(icand)>2) ncs2++;
							if(csCandPt->at(icand)>3) ncs3++;
							if(csCandPt->at(icand)>0.5 && csCandId->at(icand)==1) ncs0p5_id1++;
							if(csCandPt->at(icand)>1 && csCandId->at(icand)==1) ncs1_id1++;
							if(csCandPt->at(icand)>2 && csCandId->at(icand)==1) ncs2_id1++;
							if(csCandPt->at(icand)>3 && csCandId->at(icand)==1) ncs3_id1++;
						}
					}
					nCSpart_perJet.push_back(ncs);
					nCSpart_perJet0p5.push_back(ncs0p5);
					nCSpart_perJet1.push_back(ncs1);
					nCSpart_perJet2.push_back(ncs2);
					nCSpart_perJet3.push_back(ncs3);
					nCSpart_perJet0p5_id1.push_back(ncs0p5_id1);
					nCSpart_perJet1_id1.push_back(ncs1_id1);
					nCSpart_perJet2_id1.push_back(ncs2_id1);
					nCSpart_perJet3_id1.push_back(ncs3_id1);
				}				
		
				if( t_calo_jtpt[j4i] < 25) continue;

				calo_jteta.push_back(reco_eta);
				calo_jtphi.push_back(reco_phi);
				calo_jtpt.push_back(reco_pt);
				if(!is_data){ 
					calo_refpt.push_back(t_calo_refpt[j4i]);
					calo_parton_flavor.push_back(t_calo_parton_flavor[j4i]);
				}
				calo_discr_ssvHighEff.push_back(t_calo_discr_ssvHighEff[j4i]);
				calo_discr_ssvHighPur.push_back(t_calo_discr_ssvHighPur[j4i]);
				calo_discr_csvV1.push_back(t_calo_discr_csvV1[j4i]);
				calo_discr_prob.push_back(t_calo_discr_prob[j4i]);
				calo_svtxm.push_back(t_calo_svtxm[j4i]);
				calo_svtxpt.push_back(t_calo_svtxpt[j4i]);
				calo_svtxmcorr.push_back(t_calo_svtxmcorr[j4i]);
				calo_svtxdl.push_back(t_calo_svtxdl[j4i]);
				calo_svtxdls.push_back(t_calo_svtxdls[j4i]);


			} /// calo jet loop


			///// Fill it
			mixing_tree->Fill();

			///// Reset
			highPurity.clear();

			nCSpart_perJet.clear();
			nCSpart_perJet0p5.clear();
			nCSpart_perJet1.clear();
			nCSpart_perJet2.clear();
			nCSpart_perJet3.clear();
			nCSpart_perJet0p5_id1.clear();
			nCSpart_perJet1_id1.clear();
			nCSpart_perJet2_id1.clear();
			nCSpart_perJet3_id1.clear();
	
			calo_jteta.clear();
			calo_jtphi.clear();
			calo_jtpt.clear();

			calo_parton_flavor.clear();
			calo_discr_ssvHighEff.clear();
			calo_discr_ssvHighPur.clear();
			calo_discr_csvV1.clear();
			calo_discr_prob.clear();
			calo_svtxm.clear();
			calo_svtxpt.clear();
			calo_svtxmcorr.clear();
			calo_svtxdl.clear();
			calo_svtxdls.clear();

			if(!is_data){
				calo_refpt.clear();
			}

		}  ///event loop

		my_file->Close();
	}


	cout<<"writing"<<endl;

	//mixing_tree->Write();
	output_file->cd();
	mixing_tree->Write();

	output_file->Close();

	cout<<"done"<<endl;

}
