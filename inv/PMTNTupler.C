#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include <time.h>
#include <math.h>
#include <string.h>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TString.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TSpectrum.h"

#include "TVirtualFFT.h"
#include "Math/Functor.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"
#include "Math/Vector3D.h"
#include "Fit/Fitter.h"

#include <netinet/in.h>

#define NPMTs 8

using namespace std;
using namespace ROOT::Math;

struct event
{
	int EventNo;
	vector <float> signal[256];
};
event e;

struct edata
{
	vector <vector<int>>* TPCWF;
	vector <vector<int>>* PMTWF;
};
edata ed;

int RunNo=0;
int EventNo=0;
int opt=0;
ifstream::pos_type size=2097416;
char hname[1000];
char hname2[1000];
char TPCinFilePath[1000];
char PMTinFilePath[1000];
char AnalysisFilePath[1000];
char NTupleFilePath[1000];
ifstream file;
int data0loc=0;
int TPCBaselineRun=0;
int PMTBaselineRun=0;

int ThInd=7;
int ThCol=7;
float ThPMT=2.5;

float TPCBaselines[2][128][2]={{{0}}};
float PMTBaselines[3][2]={{0}};

int IC[256]={0};
int CH[256]={0};

double pStart[5] = {0};
double fitParams[5]={0};

int peakSearchWidth=1;

string wt[2]={"Collection","Induction"};

vector <int> B;

void WriteNTuple()
{
	ifstream pmtinfile[NPMTs];
	for(int i1=0;i1<NPMTs;i1++)
	{
		sprintf(hname,"cp %s/waveforms_ch%d_run%d.txt .;wait;",PMTinFilePath,i1,RunNo);
		system(hname);
		sprintf(hname,"waveforms_ch%d_run%d.txt",i1,RunNo);
		pmtinfile[i1].open(hname);
	}
	
	sprintf(hname,"Run_%d.root",RunNo);
	TFile* outroot=new TFile(hname,"recreate");
	TTree* T =  new TTree("T","T");
	T->Branch("TPCWF",&ed.TPCWF);
	T->Branch("PMTWF",&ed.PMTWF);
	
	vector <int> bb;
	string strarr;
	int Ev=0;
	int ii=0;
	
	while(getline(pmtinfile[0],strarr,'\n'))
	{
		ed.PMTWF->clear();
		istringstream iss(strarr);
		ii=0;
		bb.clear();
		for(std::string s; iss >> s; )
		{
			if(ii>0)
			{
				bb.push_back(atoi(s.c_str()));
			}
			ii++;
		}
		ed.PMTWF->push_back(bb);
		
		for(int i1=1;i1<NPMTs;i1++)
		{
			getline(pmtinfile[i1],strarr,'\n');
			istringstream iss(strarr);
			ii=0;
			bb.clear();
			for(std::string s; iss >> s; )
			{
				if(ii>0)
				{
					bb.push_back(atoi(s.c_str()));
				}
				ii++;
			}
			ed.PMTWF->push_back(bb);
		}
		T->Fill();
		Ev++;
	}
	for(int i1=0;i1<NPMTs;i1++)
	{
		pmtinfile[i1].close();
	}
	outroot->cd();
	T->Write();
	outroot->Close();
	file.close();
	
	sprintf(hname,"mv Run_%d.root %s/Run_%d.root;wait;",RunNo,NTupleFilePath,RunNo);
	system(hname);
}

int main( int argc, const char* argv[] )
{
	RunNo=atoi(argv[1]);
	
	//local
//	sprintf(TPCinFilePath,"/home/bbilki/D/Data/Flic/TPC");//where the Run00X.dat files reside
//	sprintf(PMTinFilePath,"/home/bbilki/D/Data/Flic/PMT");//where the waveforms_chX.txt files reside
// 	sprintf(AnalysisFilePath,".");//where Files and Histos will be placed
	
// 	//pc04-warp
//  	sprintf(TPCinFilePath,"/Data1/completed");//where the Run00X.dat files reside
//  	sprintf(PMTinFilePath,"/home/icadaq/CAENV1751Daq/PMTDAQCode_updated/mcrbee10");//where the waveforms_chX.txt files reside
// 	sprintf(AnalysisFilePath,".");//where Files and Histos will be placed
	
// 	//eos
 	sprintf(TPCinFilePath,"/eos/project/f/flic2019/Data/TPC/Runs");//where the Run00X.dat files reside
 	sprintf(PMTinFilePath,"/eos/project/f/flic2019/Data/PMT/WaveForms");//where the waveforms_chX.txt files reside
	sprintf(AnalysisFilePath,"/eos/project/f/flic2019/Analysis");//where Files and Histos will be placed
	sprintf(NTupleFilePath,"/eos/project/f/flic2019/Data/NTuple");//where the Run_XX.root files reside
	
	WriteNTuple();
}

