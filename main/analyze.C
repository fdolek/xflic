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
#include "TLegend.h"
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

using namespace std;
using namespace ROOT::Math;

struct edata
{
	vector <vector<int>>* TPCWF;
	vector <vector<int>>* PMTWF;
};
edata ed;

struct hitdata
{
	int E;//event
	vector <vector<int>>* ColIndT;
	vector <vector<float>>* Int;
	float QColTot;
	float QColTotZS;
	float QHitTot;
	vector <float>* PMTIntegral;
	vector <int>* ColID;
	vector <float>* ColT;
	vector <float>* ColA;
	vector <int>* ColInt;
	vector <int>* Colw;
	vector <int>* IndID;
	vector <float>* IndT;
	vector <float>* IndA;
	vector <int>* Indw;
	int EventType;
	//0 - single track traversing the entire drift region; 
	//1 - single track travesing partially;
	//2 - showers (includes partial showers and multiple tracks)
	//3 - accidentals (scattered shower fragments triggering)
};
hitdata hd;

struct trackdata
{
	int E;//event
	vector <int>* StartEndColIndT;//start col, start ind, start t, end col, end ind, end t
	vector <float>* FitParams;
	float FitNormChi2;
	float QColTot;
	int NHits;
	int Nexcl;
	vector <float>* PMTIntegral;
	vector <float>* ColTStartEnd;
	vector <float>* ColHitTStartEnd;
};
trackdata td;

struct peaks
{
	int IC;
	int CH;
	vector <int> t;
};
peaks p;

struct hits
{
	int col;
	int ind;
	int t;
	int Int;
	int Int2;
	int Int3;
	
	float colT;
	float indT;
	int colWidth;
	int indWidth;
	float colAmp;
	float indAmp;
};
hits hit;

struct hitcandidates
{
	vector <int> col;
	vector <int> ind;
	vector <int> colstartw;
	vector <int> colendw;
	int t;
	vector <int> n;
	vector <int> Int2;//induction
	vector <int> Int3;//collection
};
hitcandidates hc;

struct hitcandidateshort
{
	int ci;//collection 0 induction 1
	int wid;//wire id
	float t;//positive peak time
	float t2;//zero crossing time
	float t3;//negative peak time
	int Int;//default collection integral
	int Int2;//undershoot corrected collection integral
	int w;//width
	float a;//positive peak amplitude
	float a2;//pos peak to neg peak
	float a3;//negative peak amplitude
	float us;//undershoot
};
hitcandidateshort hcs;

struct clusters
{
	vector <int> hi;//hit id
};
clusters cl;

struct tlist
{
	int t;
	vector <int> ch;//0-127 col; 128-255 ind
};
tlist tl;

struct fullevent
{
	int E;
	float fitParams[5];
	float fitnormchi2;
	int clsize;
	float PMTint[3];
	vector <hits> H;
};
fullevent fe;

Double_t fitf2f(Double_t *x,Double_t *par)
{
	Double_t fitval=0.;
	Double_t sigma=0.;
	Double_t tt=(x[0]-par[1]);
	if(x[0]<par[1])
	{
		sigma=exp(par[2]*x[0]);
	}
	else
	{
		sigma=exp(-x[0]/par[3])+exp(-x[0]/par[4])+exp(-x[0]/par[5]);
	}
	fitval=par[0]*sigma;
	return fitval;
}

Double_t fitf3f(Double_t *x,Double_t *par)
{
	Double_t fitval=0.;
	Double_t sigma=0.;
	Double_t tt=(x[0]-par[1]);
	if(x[0]<par[1])
	{
		sigma=par[2]+par[3]*pow(par[1]-x[0],par[4]);
	}
	else
	{
		sigma=par[2]+par[5]*pow(x[0]-par[1],par[6]);
	}
	fitval=par[0]*exp(-0.5*pow(tt/sigma,2));
	return fitval;
}

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
int TimeRange[2]={750,2500};
float HitColIndDist=5.;
int NPMTs=3;//number of PMTs

// int ThInd=7;
// int ThCol=7;
// int ThInd=3;
// int ThCol=4;
// float ThInd=3;
// float ThCol=4;
float ThColInd[2]={0};//Col ind
float ThColIndW[2][128]={{0}};//Col ind
float ThPMT[3]={0.};

float TPCBaselines[2][128][2]={{{0}}};
float PMTBaselines[3][2]={{0}};
float ThSigma[2]={4.5,4.5};//threshold for collection and induction signal sigmas above baseline for hit reconstruction

int IC[256]={0};
int CH[256]={0};

double pStart[5] = {0};
double fitParams[5]={0};

int peakSearchWidth=1;

string wt[2]={"Collection","Induction"};
// string PMTNames[3]={"LAr w/TPB","LAr no TPB","LXe"};
string PMTNames[6]={"LAr_wTPB","LAr_noTPB","LXe","LAr_wTPB_new","SC1","SC2"};

vector <int> B;

#include "core.cc"

void ColWireLines()
{
	sprintf(hname,"ColWireLines_%d.root",RunNo);
	TFile* outroot=new TFile(hname,"recreate");
	outroot->mkdir("Trends");
	
	sprintf(hname,"%s/Run_%d.root",NTupleFilePath,RunNo);
	TFile* inroot=new TFile(hname);
	TTree* T =  (TTree*) inroot->Get("T");
	T->SetBranchAddress("TPCWF",&ed.TPCWF);
	T->SetBranchAddress("PMTWF",&ed.PMTWF);
	
	TH1F* HitCol=new TH1F("HitCol","HitCol",128,-0.5,127.5);
	TH1F* HitTime=new TH1F("HitTime","HitTime",4096,-0.5,4095.5);
	TH1F* HitIntegral1=new TH1F("HitIntegral1","HitIntegral1",1024,0,1023.5);
	TH1F* HitIntegral2=new TH1F("HitIntegral2","HitIntegral2",1024,0,1023.5);
	TH1F* HitIntegral3=new TH1F("HitIntegral3","HitIntegral3",500,0,3000);
	TH1F* NHits=new TH1F("NHits","NHits",2000,-0.5,1999.5);
	
	TH1F* StartTime=new TH1F("StartTime","StartTime",512,-0.5,4095.5);
	TH1F* EndTime=new TH1F("EndTime","EndTime",512,-0.5,4095.5);
	TH1F* StartCollection=new TH1F("StartCollection","StartCollection",128,-0.5,127.5);
	TH1F* EndCollection=new TH1F("EndCollection","EndCollection",128,-0.5,127.5);
	TH1F* DeltaT=new TH1F("DeltaT","DeltaT",1024,-0.5,4095.5);
	TH1F* DeltaCol=new TH1F("DeltaCol","DeltaCol",128,-0.5,127.5);
	TH1F* FitNormChi2=new TH1F("FitNormChi2","FitNormChi2",500,0.,50.);
	TH1F* LinePar0=new TH1F("LinePar0","LinePar0",800,-200.,200.);
	TH1F* LinePar1=new TH1F("LinePar1","LinePar1",1000,-1000.,1000.);
	
	TH2F* Col_vs_t[2];
	Col_vs_t[0]=new TH2F("Col_vs_t","Col_vs_t",4096,-0.5,4095.5,128,-0.5,127.5);
	Col_vs_t[1]=new TH2F("Col_vs_t_norm","Col_vs_t_norm",4096,-0.5,4095.5,128,-0.5,127.5);
	
	TProfile* Qvst=new TProfile("Qvst","Qvst",500,0,4000,0,2000);
	TGraph* tg1;
	
	TF1* lin1=new TF1("lin1","[0]+[1]*x",0.,200.);
	
	ReadTPCBaselines();
	
	TSpectrum *s;
	Double_t *xpeaks;
	TH1F* hh[2][128];TH1F* hh2[2][128];
	TH1F* ht[4096];TH1F* ht2[4096];
	
	for(int i1=0;i1<2;i1++)
	{
		for(int i2=0;i2<128;i2++)
		{
			sprintf(hname,"WF_%s_%d",wt[i1].c_str(),i2);
			hh[i1][i2]=new TH1F(hname,hname,4096,-0.5,4095.5);
			sprintf(hname,"WFZS_%s_%d",wt[i1].c_str(),i2);
			hh2[i1][i2]=new TH1F(hname,hname,4096,-0.5,4095.5);
		}
	}
	for(int i1=0;i1<4096;i1++)
	{	
		sprintf(hname,"WFt_%d",i1);
		ht[i1]=new TH1F(hname,hname,256,-0.5,255.5);
		sprintf(hname,"WFtZS_%d",i1);
		ht2[i1]=new TH1F(hname,hname,256,-0.5,255.5);
	}
	
	vector <hits> Hits;
	int nh=0;
	vector <hits> H;
	
	vector <int> bb;
	vector <float> bc;
	int Ev=0;
	int integral=0;int integral2=0;int integral3=0;
	int peakbin=0;
	int ic=0;int gch=0;
	float QColTot=0;
	float QColTotZS=0;
	float nchi2=0.;
	for(int I=0;I<T->GetEntries();I++)
	{
		T->GetEntry(I);
		int npeaks = 200;
		int nfound=0;
		float ss=0;
		bool found=true;
		bool HitOK=true;
		
		QColTot=0;QColTotZS=0;
		for(int i1=0;i1<128;i1++)
		{
			for(int i2=0;i2<4096;i2++)
			{
				ss=((float)ed.TPCWF->at(i1)[i2])-TPCBaselines[0][i1][0];
				ht[i2]->SetBinContent(ht[i2]->FindBin(i1),ss);
				if(ss>ThColInd[0])
				{
					ht2[i2]->SetBinContent(ht2[i2]->FindBin(i1),ss);
					hh2[0][i1]->SetBinContent(hh2[0][i1]->FindBin(i2),ss);
				}
				QColTot+=ss;
				if(ss>0) QColTotZS+=ss;
			}
		}
		int nh=0;
		for(int i1=0;i1<128;i1++)
		{
			s = new TSpectrum(2*npeaks);
			nfound = s->Search(hh2[0][i1],3.,"",0.1);
			xpeaks = s->GetPositionX();
			for(int is1=0;is1<nfound;is1++)
			{
				peakbin=hh2[0][i1]->FindBin(xpeaks[is1]);
				integral3=0;
				for(int ip1=peakbin;ip1<=hh2[0][i1]->GetNbinsX();ip1++)
				{
					if(hh2[0][i1]->GetBinContent(ip1)>ThColInd[0])
					{
						integral3+=hh2[0][i1]->GetBinContent(ip1);
					}
					else break;
				}
				for(int ip1=peakbin-1;ip1>=1;ip1--)
				{
					if(hh2[0][i1]->GetBinContent(ip1)>ThColInd[0])
					{
						integral3+=hh2[0][i1]->GetBinContent(ip1);
					}
					else break;
				}
				Col_vs_t[0]->Fill(xpeaks[is1],i1,integral3);
				Col_vs_t[1]->Fill(xpeaks[is1],i1);
				nh++;
			}
			delete s;
		}
		
// 		for(int i1=0;i1<4096;i1++)
// 		{
// 			s = new TSpectrum(2*npeaks);
// 			nfound = s->Search(ht2[i1],3.,"",0.1);
// 			xpeaks = s->GetPositionX();
// 			
// 			for(int is1=0;is1<nfound;is1++)
// 			{
// 				peakbin=ht[i1]->FindBin(xpeaks[is1]);
// 				integral=0;integral2=0;integral3=0;
// 				for(int ip1=peakbin;ip1<=ht[i1]->GetNbinsX();ip1++)
// 				{
// 					if(ht[i1]->GetBinContent(ip1)>ThColInd[0])
// 					{
// 						integral3+=ht[i1]->GetBinContent(ip1);
// 					}
// 					else break;
// 				}
// 				for(int ip1=peakbin-1;ip1>=1;ip1--)
// 				{
// 					if(ht[i1]->GetBinContent(ip1)>ThColInd[0])
// 					{
// 						integral3+=ht[i1]->GetBinContent(ip1);
// 					}
// 					else break;
// 				}
// 				
// 				hit.col=((int)xpeaks[is1]);
// 				hit.ind=0;
// 				hit.t=i1;
// 				hit.Int=0;//induction amplitude
// 				hit.Int2=ht2[i1]->GetBinContent(ht2[i1]->FindBin(xpeaks[is1]));//collection amplitude
// 				hit.Int3=integral3;
// 				
// 				Hits.push_back(hit);
// 				nh++;
// 			}
// 			delete s;
// 		}
// 		tg1=new TGraph();
// 		
// 		int mins[3]={5000,5000,5000};int maxs[3]={0,0,0};//Col, ind, t
// 		for(int i1=0;i1<Hits.size();i1++)
// 		{
// 			tg1->SetPoint(tg1->GetN(),Hits[i1].col,Hits[i1].t);
// 			if(Hits[i1].col<mins[0]) mins[0]=Hits[i1].col;if(Hits[i1].col>maxs[0]) maxs[0]=Hits[i1].col;
// 			if(Hits[i1].t<mins[2]) mins[2]=Hits[i1].t;if(Hits[i1].t>maxs[2]) maxs[2]=Hits[i1].t;
// 		}
// 		tg1->Fit(lin1,"q","q",0.,200.);
// 		nchi2=lin1->GetChisquare()/lin1->GetNDF();
// 		LinePar0->Fill(lin1->GetParameter(0));
// 		LinePar1->Fill(lin1->GetParameter(1));
// 		FitNormChi2->Fill(nchi2);
// 		outroot->cd("Trends");
// 		tg1->SetMarkerStyle(20);
// 		tg1->SetMarkerColor(1);
// 		sprintf(hname,"Ev_%d",I);
// 		tg1->SetName(hname);
// 		tg1->SetTitle(hname);
// 		tg1->Write();
// // 		if(nchi2<1.)
// 		{
// 			NHits->Fill(nh);
// 			StartTime->Fill(mins[2]);
// 			EndTime->Fill(maxs[2]);
// 			StartCollection->Fill(mins[0]);
// 			EndCollection->Fill(maxs[0]);
// 			DeltaT->Fill(maxs[2]-mins[2]);
// 			DeltaCol->Fill(maxs[0]-mins[0]);
// 			for(int i1=0;i1<Hits.size();i1++)
// 			{
// 				Qvst->Fill(Hits[i1].t,Hits[i1].Int3);
// 				HitCol->Fill(Hits[i1].col);
// 				HitTime->Fill(Hits[i1].t);
// 				HitIntegral3->Fill(Hits[i1].Int3);
// 			}
// 		}
		for(int i1=0;i1<4096;i1++)
		{
			ht[i1]->Reset();
			ht2[i1]->Reset();
		}
		for(int i1=0;i1<128;i1++)
		{
			hh2[0][i1]->Reset();
		}
		Hits.clear();
		if(I%10==0) cout<<"Event: "<<I<<endl;
// 		if(I==1000) break;
	}
	outroot->cd();
	HitCol->Write();
	HitTime->Write();
	HitIntegral1->Write();
	HitIntegral2->Write();
	HitIntegral3->Write();
	NHits->Write();
	StartTime->Write();
	EndTime->Write();
	StartCollection->Write();
	EndCollection->Write();
	DeltaT->Write();
	DeltaCol->Write();
	FitNormChi2->Write();
	LinePar0->Write();
	LinePar1->Write();
	Col_vs_t[0]->Divide(Col_vs_t[1]);
	Col_vs_t[0]->Write();
	Col_vs_t[1]->Write();
	
	outroot->Close();
	inroot->Close();
	
	sprintf(hname,"cp ColWireLines_%d.root %s/Histos/ColWireLines_%d.root;wait;",RunNo,AnalysisFilePath,RunNo);system(hname);
}

void CompareTPCBaselines()
{
// 	int nr=7;
// 	int runnos[7]={7075,7077,7080,7083,7086,7088,7090};
// 	
// 	float TPCBL[2][128][2][20]={{{{0.}}}};
// 	
// 	for(int i1=0;i1<nr;i1++)
// 	{
// 		TPCBaselineRun=runnos[i1];
// 		ReadTPCBaselines();
// 		for(int i2=0;i2<2;i2++)
// 		{
// 			for(int i3=0;i3<128;i3++)
// 			{
// 				for(int i4=0;i4<2;i4++)
// 				{
// 					TPCBL[i2][i3][i4][i1]=TPCBaselines[i2][i3][i4];
// 				}
// 			}
// 		}
// 	}
// 	sprintf(hname,"%s/Histos/CompareTPCBaselines.root",AnalysisFilePath);
// 	TFile* outroot=new TFile(hname,"recreate");
// 	
// 	TH1F* Diff[2];
// 	Diff[0]=new TH1F("Diff1","Diff1",11,-5.5,5.5);
// 	Diff[1]=new TH1F("Diff2","Diff2",22,-5.5,5.5);
// 	for(int i1=1;i1<nr;i1++)
// 	{
// 		for(int i2=0;i2<2;i2++)
// 		{
// 			for(int i3=0;i3<128;i3++)
// 			{
// 				for(int i4=0;i4<2;i4++)
// 				{
// 					Diff[i4]->Fill(TPCBL[i2][i3][i4][i1]-TPCBL[i2][i3][i4][0]);
// 				}
// 			}
// 		}
// 	}
// 	outroot->Write();
// 	outroot->Close();
}

void CombinedTPCPMTPlot()
{
// 	vector <fullevent> FE;
// // 	sprintf(hname,"%s/Files/TrackParameters_%d.txt",AnalysisFilePath,RunNo);
// 	sprintf(hname,"cp %s/Files/TrackParameters_%d.txt .;wait;",AnalysisFilePath,RunNo);system(hname);
// 	sprintf(hname,"TrackParameters_%d.txt",RunNo);
// 	ifstream inTracks(hname);
// 	while(!inTracks.eof())
// 	{
// 		inTracks>>fe.E>>fe.fitParams[0]>>fe.fitParams[1]>>fe.fitParams[2]>>fe.fitParams[3]>>fe.fitParams[4]>>fe.fitnormchi2>>fe.clsize;
// 		fe.PMTint[0]=0;fe.PMTint[1]=0;fe.PMTint[2]=0;
// 		FE.push_back(fe);
// 	}
// 	inTracks.close();
// 	if(FE[FE.size()-1].E==FE[FE.size()-2].E) FE.pop_back();
// 	int nEvPMT=0;
// // 	sprintf(hname,"%s/Files/PMTIntegrals_%d.txt",AnalysisFilePath,RunNo);
// 	sprintf(hname,"cp %s/Files/PMTIntegrals_%d.txt .;wait;",AnalysisFilePath,RunNo);system(hname);
// 	sprintf(hname,"PMTIntegrals_%d.txt",RunNo);
// 	ifstream inPMTs(hname);
// 	while(!inPMTs.eof())
// 	{
// 		inPMTs>>fe.E>>fe.PMTint[0]>>fe.PMTint[1]>>fe.PMTint[2];
// 		for(int i1=(fe.E>1?fe.E-2:0);i1<FE.size();i1++)
// 		{
// 			if(FE[i1].E==fe.E)
// 			{
// 				FE[i1].PMTint[0]=fe.PMTint[0];
// 				FE[i1].PMTint[1]=fe.PMTint[1];
// 				FE[i1].PMTint[2]=fe.PMTint[2];
// 				break;
// 			}
// 		}
// 		nEvPMT++;
// 	}
// 	inPMTs.close();
// 	
// 	vector <int> eventlist;
// 	for(int i1=0;i1<nEvPMT;i1++)
// 	{
// 		eventlist.push_back(-1);
// 		for(int i2=0;i2<FE.size();i2++)
// 		{
// 			if(i1==FE[i2].E)
// 			{
// 				eventlist[eventlist.size()-1]=i2;
// 				break;
// 			}
// 		}
// 	}
// 	
// // 	sprintf(hname,"%s/Files/Hits_%d.txt",AnalysisFilePath,RunNo);
// 	sprintf(hname,"cp %s/Files/Hits_%d.txt .;wait;",AnalysisFilePath,RunNo);system(hname);
// 	sprintf(hname,"Hits_%d.txt",RunNo);
// 	ifstream inHits(hname);
// 	int b[7]={0};
// 	int Ev=0;
// 	bool found=false;
// 	bool skipEvent=false;
// 	int curInd=0;
// 	
// 	while(!inHits.eof())
// 	{
// 		inHits>>b[0]>>b[1]>>b[2]>>b[3]>>b[4]>>b[5]>>b[6];
// 		if(eventlist[b[0]]>=0)
// 		{
// 			hit.col=b[1];
// 			hit.ind=b[2];
// 			hit.t=b[3];
// 			hit.Int=b[4];
// 			hit.Int2=b[5];
// 			hit.Int3=b[6];
// 			FE[eventlist[b[0]]].H.push_back(hit);
// 		}
// 	}
// 	inHits.close();
// 	
// // 	sprintf(hname,"%s/Histos/PMTWaveforms_%d.root",AnalysisFilePath,RunNo);
// 	sprintf(hname,"cp %s/Histos/PMTWaveforms_%d.root .;wait;",AnalysisFilePath,RunNo);system(hname);
// 	sprintf(hname,"PMTWaveforms_%d.root",RunNo);
// 	TFile* inrootPMTWF=new TFile(hname);
// 	TH1F* hPMT[3];
// 	string PMTNames[3]={"LAr w/TPB","LAr no TPB","LXe"};
// 	int PMTlc[3]={1,2,4};
// 	
// 	sprintf(hname,"CombinedPlots_%d.root",RunNo);
// 	TFile* outroot=new TFile(hname,"recreate");
// 	
// 	TH2F* EntryXY=new TH2F("EntryXY","EntryXY",128,-0.5,127.5,128,-0.5,127.5);
// 	TH1F* EntryT=new TH1F("EntryT","EntryT",256,-0.5,4095.5);
// 	TH2F* ExitXY=new TH2F("ExitXY","ExitXY",128,-0.5,127.5,128,-0.5,127.5);
// 	TH1F* ExitT=new TH1F("ExitT","ExitT",256,-0.5,4095.5);
// 	TH1F* DeltaT=new TH1F("DeltaT","DeltaT",256,-0.5,4095.5);
// 	TH1F* AllHitX=new TH1F("AllHitX","AllHitX",128,-0.5,127.5);
// 	TH1F* AllHitY=new TH1F("AllHitY","AllHitY",128,-0.5,127.5);
// 	TH1F* AllHitT=new TH1F("AllHitT","AllHitT",4096,-0.5,4095.5);
// 	TH1F* AllHitQ=new TH1F("AllHitQ","AllHitQ",10000,-0.5,9999.5);
// 	TH2F* Qvst=new TH2F("Qvst","Qvst",1400,1000,2400,1000,0,1000);
// 	
// 	TH2F* TrackHits[7][2];
// 	for(int i1=0;i1<7;i1++)
// 	{
// 		sprintf(hname,"MeanCharge_%d-%d",(i1*200)+1000,((i1+1)*200)+1000);
// 		TrackHits[i1][0]=new TH2F(hname,hname,128,-0.5,127.5,128,-0.5,127.5);
// 		sprintf(hname,"NHits_%d-%d",(i1*200)+1000,((i1+1)*200)+1000);
// 		TrackHits[i1][1]=new TH2F(hname,hname,128,-0.5,127.5,128,-0.5,127.5);
// 	}	
// 	
// 	TPolyLine3D *l;
// 	TH3F *frame3d=new TH3F("frame3d","frame3d",1,-0.5,127.5,1,-0.5,127.5,1,-0.5,4095.5);
// 	TPolyMarker3D *pm3d1;
// 	TPolyMarker3D *pm3d2;
// 	double x(0),y(0),z(0),t(0);
// 	double tlast(0);
// 	TCanvas* cc=new TCanvas("cc","cc",200,10,900,600);
// 	TPad *padTPC = new TPad("padTPC","TPC",0.01,0.01,0.65,0.99,0,1,0);
// 	TPad *padPMT[3];
// 	padPMT[0] = new TPad("pad0","LAr w/TPB",0.67,0.67,0.99,0.99,0,1,0);
// 	padPMT[1] = new TPad("pad1","LAr no TPB",0.67,0.35,0.99,0.65,0,1,0);
// 	padPMT[2] = new TPad("pad2","LXe",0.67,0.01,0.99,0.33,0,1,0);
// 	cc->cd();
// 	gStyle->SetOptStat(0);
// 	padTPC->Draw();
// 	gStyle->SetOptStat(0);
// 	padPMT[0]->Draw();padPMT[1]->Draw();padPMT[2]->Draw();
// 	bool foundfirst=true;bool foundlast=true;
// 	int tmin=5000;int indtmin=0;
// 	int tmax=0;int indtmax=0;
// 	
// 	for(int i1=0;i1<FE.size();i1++)
// 	{
// // 		cc->cd();
// 		padTPC->cd();
// // 		gStyle->SetTitleFontSize(0);
// 		frame3d->Draw();
// 		
// 		pm3d1 = new TPolyMarker3D(FE[i1].H.size());
// 		
// 		tmin=5000;tmax=0;
// 		indtmin=0;indtmax=0;
// 		for(int i2=0;i2<FE[i1].H.size();i2++)
// 		{
// 			pm3d1->SetPoint(i2,FE[i1].H[i2].col,FE[i1].H[i2].ind,FE[i1].H[i2].t);
// 			AllHitX->Fill(FE[i1].H[i2].col);
// 			AllHitY->Fill(FE[i1].H[i2].ind);
// 			AllHitT->Fill(FE[i1].H[i2].t);
// 			AllHitQ->Fill(FE[i1].H[i2].Int3);
// 			Qvst->Fill(FE[i1].H[i2].t,FE[i1].H[i2].Int3);
// 			if(FE[i1].H[i2].t<tmin){tmin=FE[i1].H[i2].t;indtmin=i2;}
// 			if(FE[i1].H[i2].t>tmax){tmax=FE[i1].H[i2].t;indtmax=i2;}
// 		}
// 		EntryXY->Fill(FE[i1].H[indtmin].col,FE[i1].H[indtmin].ind);
// 		ExitXY->Fill(FE[i1].H[indtmax].col,FE[i1].H[indtmax].ind);
// 		DeltaT->Fill(FE[i1].H[indtmax].t-FE[i1].H[indtmin].t);
// 		EntryT->Fill(FE[i1].H[indtmin].t);
// 		ExitT->Fill(FE[i1].H[indtmax].t);
// 		
// 		pm3d1->SetMarkerColor(kRed);
// 		pm3d1->SetMarkerStyle(24);   
// 		pm3d1->Draw();
// 		l = new TPolyLine3D(2);
// 		
// 		t=((FE[i1].fitParams[4]-(FE[i1].fitParams[0]/FE[i1].fitParams[1]))>(FE[i1].fitParams[4]-(FE[i1].fitParams[2]/FE[i1].fitParams[3]))?(FE[i1].fitParams[4]-(FE[i1].fitParams[0]/FE[i1].fitParams[1])):(FE[i1].fitParams[4]-(FE[i1].fitParams[2]/FE[i1].fitParams[3])));
// 		x=FE[i1].fitParams[0]+FE[i1].fitParams[1]*(t-FE[i1].fitParams[4]);
// 		y=FE[i1].fitParams[2]+FE[i1].fitParams[3]*(t-FE[i1].fitParams[4]);
// 		l->SetPoint(0,x,y,t);
// 		t=((FE[i1].fitParams[4]+((127-FE[i1].fitParams[0])/FE[i1].fitParams[1]))>(FE[i1].fitParams[4]+((127-FE[i1].fitParams[2])/FE[i1].fitParams[3]))?(FE[i1].fitParams[4]+((127-FE[i1].fitParams[2])/FE[i1].fitParams[3])):(FE[i1].fitParams[4]+((127-FE[i1].fitParams[0])/FE[i1].fitParams[1])));
// 		x=FE[i1].fitParams[0]+FE[i1].fitParams[1]*(t-FE[i1].fitParams[4]);
// 		y=FE[i1].fitParams[2]+FE[i1].fitParams[3]*(t-FE[i1].fitParams[4]);
// 		l->SetPoint(1,x,y,t);
// 		
// 		l->SetLineColor(kBlack);
// 		l->SetLineWidth(2);
// 		l->Draw("same");
// 		
// 		for(int is1=0;is1<3;is1++)
// 		{
// 			sprintf(hname,"Ch%d/WF_PMT_%d_Ev_%d",is1,is1,FE[i1].E);
// 			inrootPMTWF->GetObject(hname,hPMT[is1]);
// 			sprintf(hname2,"%s %s",PMTNames[is1].c_str(),hPMT[is1]->GetTitle());
// 			hPMT[is1]->SetTitle(hname2);
// 			hPMT[is1]->SetLineColor(PMTlc[is1]);
// 			padPMT[is1]->cd();
// 			hPMT[is1]->Draw();
// 		}
// 		
// 		sprintf(hname,"Ev_%d",FE[i1].E);
// 		cc->SetName(hname);cc->SetTitle(hname);
// 		outroot->cd();
// 		cc->Write();
// 		
// 		for(int i2=0;i2<FE[i1].H.size();i2++)
// 		{
// 			if(FE[i1].H[i2].t<1000 || FE[i1].H[i2].t>=2400) continue;
// 			TrackHits[(FE[i1].H[i2].t-1000)/200][0]->Fill(FE[i1].H[i2].col,FE[i1].H[i2].ind,FE[i1].H[i2].Int3);
// 			TrackHits[(FE[i1].H[i2].t-1000)/200][1]->Fill(FE[i1].H[i2].col,FE[i1].H[i2].ind);
// 		}
// 	}
// 	outroot->cd();
// 	EntryXY->Write();
// 	EntryT->Write();
// 	ExitXY->Write();
// 	ExitT->Write();
// 	DeltaT->Write();
// 	AllHitX->Write();
// 	AllHitY->Write();
// 	AllHitT->Write();
// 	AllHitQ->Write();
// 	Qvst->Write();
// 	
// 	for(int i1=0;i1<7;i1++)
// 	{
// 		TrackHits[i1][0]->Divide(TrackHits[i1][1]);
// 		TrackHits[i1][0]->Write();
// 		TrackHits[i1][1]->Write();
// 	}
// 	
// 	outroot->Close();
// 	inrootPMTWF->Close();
// 	
// 	sprintf(hname,"rm TrackParameters_%d.txt",RunNo);system(hname);
// 	sprintf(hname,"rm PMTIntegrals_%d.txt",RunNo);system(hname);
// 	sprintf(hname,"rm Hits_%d.txt",RunNo);system(hname);
// 	sprintf(hname,"mv CombinedPlots_%d.root %s/Histos/CombinedPlots_%d.root",RunNo,AnalysisFilePath,RunNo);system(hname);
}

void ShowerAnalysis()
{
	int tmin=1000;
	if(RunNo>7300) tmin=1500;
	sprintf(hname,"ShowerPlots_%d.root",RunNo);
	TFile* outroot=new TFile(hname,"recreate");
	TH1F* TotalCollectionCharge=new TH1F("TotalCollectionCharge","TotalCollectionCharge",5000,0,1000000);
	TH2F* TCCvsNHits=new TH2F("TCCvsNHits","TCCvsNHits",1000,0,1000,1000,0,500000);
	TH1F* AllHitX=new TH1F("AllHitX","AllHitX",128,-0.5,127.5);
	TH1F* AllHitY=new TH1F("AllHitY","AllHitY",128,-0.5,127.5);
	TH1F* AllHitT=new TH1F("AllHitT","AllHitT",4096,-0.5,4095.5);
	TH1F* AllHitQ=new TH1F("AllHitQ","AllHitQ",10000,-0.5,9999.5);
	TH2F* Qvst=new TH2F("Qvst","Qvst",1400,1000,2400,1000,0,1000);
	TH3F* NHitsvst3d=new TH3F("NHitsvst3d","NHitsvst3d",128,-0.5,127.5,128,-0.5,127.5,14,1000,2400);
// 	int hc[128][128][140]={{{0}}};
	TGraph2D* tg2d=new TGraph2D();
	NHitsvst3d->SetFillColor(4);
	TH2F* TrackHits[7][3];
	for(int i1=0;i1<7;i1++)
	{
		sprintf(hname,"MeanCharge_%d-%d",(i1*200)+tmin,((i1+1)*200)+tmin);
		TrackHits[i1][0]=new TH2F(hname,hname,128,-0.5,127.5,128,-0.5,127.5);
// 		sprintf(hname,"NHits_%d-%d",(i1*200)+1000,((i1+1)*200)+1000);
// 		TrackHits[i1][1]=new TH2F(hname,hname,128,-0.5,127.5,128,-0.5,127.5);
		
		sprintf(hname,"NHits_%d-%d",(i1*200)+tmin,((i1+1)*200)+tmin);
		TrackHits[i1][1]=new TH2F(hname,hname,128,-0.5,127.5,128,-0.5,127.5);
		
		sprintf(hname,"NHitsN_%d-%d",(i1*200)+tmin,((i1+1)*200)+tmin);
		TrackHits[i1][2]=new TH2F(hname,hname,128,-0.5,127.5,128,-0.5,127.5);
	}
	
	int b[7]={0};
	int Ev=0;
	bool found=false;
	bool skipEvent=false;
	int curInd=0;
	float qtot=0;
	int nh=0;
	
// 	TH3D* H3d=new TH3D("H3d","H3d",128,-0.5,127.5,128,-0.5,127.5,100,1000,2000);
	TH3D* H3d=new TH3D("H3d","H3d",128,-0.5,127.5,128,-0.5,127.5,100,tmin,tmin+1000);
	
	vector <hits> H;
	
	sprintf(hname,"%s/Histos/Hits_%d.root",AnalysisFilePath,RunNo);
	TFile* inroot2=new TFile(hname);
	TTree* THit =  (TTree*) inroot2->Get("T");
	THit->SetBranchAddress("E",&hd.E);
	THit->SetBranchAddress("ColIndT",&hd.ColIndT);
	THit->SetBranchAddress("Int",&hd.Int);
	THit->SetBranchAddress("QColTot",&hd.QColTot);
	THit->SetBranchAddress("QColTotZS",&hd.QColTotZS);
	THit->SetBranchAddress("PMTIntegral",&hd.PMTIntegral);
	
	for(int I=0;I<THit->GetEntries();I++)
	{
		THit->GetEntry(I);
		TotalCollectionCharge->Fill(hd.QColTot);
		TCCvsNHits->Fill(hd.ColIndT->size(),hd.QColTot);
// 		if(hd.QColTot>30000)
		{
			for(int i2=0;i2<hd.ColIndT->size();i2++)
			{
				AllHitX->Fill(hd.ColIndT->at(i2)[0]);
				AllHitY->Fill(hd.ColIndT->at(i2)[1]);
				AllHitT->Fill(hd.ColIndT->at(i2)[2]);
				AllHitQ->Fill(hd.Int->at(i2)[2]);
				Qvst->Fill(hd.ColIndT->at(i2)[2],hd.Int->at(i2)[2]);
				NHitsvst3d->Fill(hd.ColIndT->at(i2)[0],hd.ColIndT->at(i2)[1],hd.ColIndT->at(i2)[2]);
// 				if(!(hd.ColIndT->at(i2)[2]<1000 || hd.ColIndT->at(i2)[2]>=2000))
// 				if(!(hd.ColIndT->at(i2)[2]<1500 || hd.ColIndT->at(i2)[2]>=2500))
				if(!(hd.ColIndT->at(i2)[2]<tmin || hd.ColIndT->at(i2)[2]>=(tmin+1000)))
				{
// 					TrackHits[(hd.ColIndT->at(i2)[2]-1000)/200][0]->Fill(hd.ColIndT->at(i2)[0],hd.ColIndT->at(i2)[1],hd.Int->at(i2)[2]);
// 					TrackHits[(hd.ColIndT->at(i2)[2]-1000)/200][1]->Fill(hd.ColIndT->at(i2)[0],hd.ColIndT->at(i2)[1]);
// 					TrackHits[(hd.ColIndT->at(i2)[2]-1500)/200][0]->Fill(hd.ColIndT->at(i2)[0],hd.ColIndT->at(i2)[1],hd.Int->at(i2)[2]);
// 					TrackHits[(hd.ColIndT->at(i2)[2]-1500)/200][1]->Fill(hd.ColIndT->at(i2)[0],hd.ColIndT->at(i2)[1]);
					TrackHits[(hd.ColIndT->at(i2)[2]-tmin)/200][0]->Fill(hd.ColIndT->at(i2)[0],hd.ColIndT->at(i2)[1],hd.Int->at(i2)[2]);
					TrackHits[(hd.ColIndT->at(i2)[2]-tmin)/200][1]->Fill(hd.ColIndT->at(i2)[0],hd.ColIndT->at(i2)[1]);
					tg2d->SetPoint(tg2d->GetN(),hd.ColIndT->at(i2)[0],hd.ColIndT->at(i2)[1],hd.ColIndT->at(i2)[2]);
// 						hc[H[i2].col][H[i2].ind][(H[i2].t-1000)/100]++;
					H3d->AddBinContent(H3d->FindBin(hd.ColIndT->at(i2)[0],hd.ColIndT->at(i2)[1],hd.ColIndT->at(i2)[2]));
				}
			}
		}
	}
	
	for(int i1=0;i1<7;i1++)
	{
		float nmax=0.;
		for(int is1=1;is1<=TrackHits[i1][1]->GetNbinsX();is1++)
		{
			for(int is2=1;is2<=TrackHits[i1][1]->GetNbinsY();is2++)
			{
				if(TrackHits[i1][1]->GetBinContent(is1,is2)>nmax) nmax=TrackHits[i1][1]->GetBinContent(is1,is2);
			}
		}
		for(int is1=1;is1<=TrackHits[i1][2]->GetNbinsX();is1++)
		{
			for(int is2=1;is2<=TrackHits[i1][2]->GetNbinsY();is2++)
			{
				TrackHits[i1][2]->SetBinContent(is1,is2,TrackHits[i1][1]->GetBinContent(is1,is2)/nmax);
			}
		}
	}
	
	inroot2->Close();
	
	TotalCollectionCharge->Write();
	TCCvsNHits->Write();
	outroot->cd();
//	EntryXY->Write();
//	EntryT->Write();
//	ExitXY->Write();
//	ExitT->Write();
	AllHitX->Write();
	AllHitY->Write();
	AllHitT->Write();
	AllHitQ->Write();
	Qvst->Write();
	H3d->Write();
	
	TCanvas*cc=new TCanvas("cc","cc",600,600);
//	tg2d->Draw("surf4");
	H3d->Draw("iso");
	cc->SetName("H3d_iso");
	cc->Write();
	H3d->Draw("box");
	cc->SetName("H3d_box");
	cc->Write();
	H3d->Draw("box1");
	cc->SetName("H3d_box1");
	cc->Write();
	H3d->Draw("box2");
	cc->SetName("H3d_box2");
	cc->Write();
	H3d->Draw("box2z");
	cc->SetName("H3d_box2z");
	cc->Write();
	H3d->Draw("box3");
	cc->SetName("H3d_box3");
	cc->Write();
	H3d->Draw("lego");
	cc->SetName("H3d_lego");
	cc->Write();
	
	NHitsvst3d->Write();
	tg2d->Write();
	for(int i1=0;i1<7;i1++)
	{
		TrackHits[i1][0]->Divide(TrackHits[i1][1]);
		TrackHits[i1][0]->Write();
		TrackHits[i1][1]->Write();
		TrackHits[i1][2]->Write();
	}
	
	outroot->Close();
//	inrootPMTWF->Close();
}

void TrackAnalysis()
{
	sprintf(hname,"%s/Histos/Tracks_%d.root",AnalysisFilePath,RunNo);
	TFile* inroot1=new TFile(hname);
	TTree* TTrack =  (TTree*) inroot1->Get("T");
	TTrack->SetBranchAddress("E",&td.E);
	TTrack->SetBranchAddress("StartEndColIndT",&td.StartEndColIndT);
	TTrack->SetBranchAddress("FitParams",&td.FitParams);
	TTrack->SetBranchAddress("FitNormChi2",&td.FitNormChi2);
	TTrack->SetBranchAddress("QColTot",&td.QColTot);
	TTrack->SetBranchAddress("NHits",&td.NHits);
	TTrack->SetBranchAddress("Nexcl",&td.Nexcl);
	TTrack->SetBranchAddress("PMTIntegral",&td.PMTIntegral);
	
	sprintf(hname,"%s/Histos/Hits_%d.root",AnalysisFilePath,RunNo);
	TFile* inroot2=new TFile(hname);
	TTree* THit =  (TTree*) inroot2->Get("T");
	THit->SetBranchAddress("E",&hd.E);
	THit->SetBranchAddress("ColIndT",&hd.ColIndT);
	THit->SetBranchAddress("Int",&hd.Int);
	THit->SetBranchAddress("QColTot",&hd.QColTot);
	THit->SetBranchAddress("QColTotZS",&hd.QColTotZS);
	THit->SetBranchAddress("PMTIntegral",&hd.PMTIntegral);
	
	cout<<TTrack->GetEntries()<<" tracks "<<THit->GetEntries()<<" hits "<<endl;
	
	TH1F* hPMT[3];
// 	string PMTNames[3]={"LAr w/TPB","LAr no TPB","LXe"};
	int PMTlc[3]={1,2,4};
	
	sprintf(hname,"TrackAnalysis2_%d.root",RunNo);
	TFile* outroot=new TFile(hname,"recreate");
	
	TH2F* EntryXY=new TH2F("EntryXY","EntryXY",128,-0.5,127.5,128,-0.5,127.5);
	TH1F* EntryT=new TH1F("EntryT","EntryT",256,-0.5,4095.5);
	TH2F* ExitXY=new TH2F("ExitXY","ExitXY",128,-0.5,127.5,128,-0.5,127.5);
	TH1F* ExitT=new TH1F("ExitT","ExitT",256,-0.5,4095.5);
	TH1F* Trackt0=new TH1F("Trackt0","Trackt0",400,0,4000);
	TH1F* DeltaT=new TH1F("DeltaT","DeltaT",256,-0.5,4095.5);
	TH1F* TrackDeltaT=new TH1F("TrackDeltaT","TrackDeltaT",400,0,4000);
	TH1F* AllHitX=new TH1F("AllHitX","AllHitX",128,-0.5,127.5);
	TH1F* AllHitY=new TH1F("AllHitY","AllHitY",128,-0.5,127.5);
	TH1F* AllHitT=new TH1F("AllHitT","AllHitT",4096,-0.5,4095.5);
	TH1F* AllHitQ=new TH1F("AllHitQ","AllHitQ",500,0,1000);
	TH1F* TrackQ=new TH1F("TrackQ","TrackQ",100,0,100000);
	TH1F* TrackDeltaR=new TH1F("TrackDeltaR","TrackDeltaR",1000,0,2500);
	TH1F* TrackDeltaXY=new TH1F("TrackDeltaXY","TrackDeltaXY",200,0,200);
	TH2F* Qvst=new TH2F("Qvst","Qvst",1400,1000,2400,1000,0,1000);
	TH2F* Dt1vsDt2=new TH2F("Dt1vsDt2","Dt1vsDt2",500,-1000,1000,500,-1000,1000);
	TProfile* QvstPr=new TProfile("QvstPr","QvstPr",200,1000,2400,0,1000);
	TProfile* DeltaTvsDeltaR=new TProfile("DeltaTvsDeltaR","DeltaTvsDeltaR",200,0,180,0,1000);
	TProfile* DeltaTvsTrBin=new TProfile("DeltaTvsTrBin","DeltaTvsTrBin",200,0,1000,0,1000);
	TGraph* DeltaTvsTr=new TGraph();
	DeltaTvsTr->SetName("DeltaTvsTr");DeltaTvsTr->SetTitle("DeltaTvsTr");
	DeltaTvsTr->SetMarkerStyle(20);DeltaTvsTr->SetMarkerColor(1);
	TH1F* TrackHitCharge=new TH1F("TrackHitCharge","TrackHitCharge",1000,0,10000);
	TH1F* OutOfTrackHitCharge=new TH1F("OutOfTrackHitCharge","OutOfTrackHitCharge",1000,0,10000);
	TH2F* NHitsinvsoutT=new TH2F("NHitsinvsoutT","NHitsinvsoutT",500,-0.5,499.5,500,-0.5,499.5);
	
	TH2F* EntryXYPMT[3];TH2F* EntryXYPMTSat[3];TH2F* MeanPMTChargeEntryXY[3];
	TH1F* PMTInt[3];TH1F* PMTIntPE[3];TH1F* PMTIntAll[3];
	TProfile* PMTProf[3];TProfile* PMTProfPE[3];
	int PMTcols[3]={1,2,4};
	float PMTPEs[3]={34,1,23};
	for(int i1=0;i1<3;i1++)
	{
		sprintf(hname,"PMTIntAll_%d",i1);
		PMTIntAll[i1]=new TH1F(hname,hname,2000,-1000000,1000000);
		PMTIntAll[i1]->SetLineColor(PMTcols[i1]);
		sprintf(hname,"PMTInt_%d",i1);
		PMTInt[i1]=new TH1F(hname,hname,200,0,100000);
		PMTInt[i1]->SetLineColor(PMTcols[i1]);
		sprintf(hname,"PMTIntPE_%d",i1);
		PMTIntPE[i1]=new TH1F(hname,hname,300,0,3000);
		PMTIntPE[i1]->SetLineColor(PMTcols[i1]);
		sprintf(hname,"EntryXYPMT_%d",i1);
		EntryXYPMT[i1]=new TH2F(hname,hname,128,-0.5,127.5,128,-0.5,127.5);
		sprintf(hname,"EntryXYPMTSat_%d",i1);
		EntryXYPMTSat[i1]=new TH2F(hname,hname,128,-0.5,127.5,128,-0.5,127.5);
		sprintf(hname,"MeanPMTChargeEntryXY_%d",i1);
		MeanPMTChargeEntryXY[i1]=new TH2F(hname,hname,128,-0.5,127.5,128,-0.5,127.5);
		sprintf(hname,"PMTChageProf_%d",i1);
		PMTProf[i1]=new TProfile(hname,hname,50,0,150,0,100000);
		PMTProf[i1]->SetLineColor(PMTcols[i1]);
		sprintf(hname,"PMTPEProf_%d",i1);
		PMTProfPE[i1]=new TProfile(hname,hname,50,0,150,0,3000);
		PMTProfPE[i1]->SetLineColor(PMTcols[i1]);
	}
	
	TH2F* TrackHits[9][2];
	for(int i1=0;i1<9;i1++)
	{
		sprintf(hname,"MeanCharge_%d-%d",(i1*200)+800,((i1+1)*200)+800);
		TrackHits[i1][0]=new TH2F(hname,hname,128,-0.5,127.5,128,-0.5,127.5);
		sprintf(hname,"NHits_%d-%d",(i1*200)+800,((i1+1)*200)+800);
		TrackHits[i1][1]=new TH2F(hname,hname,128,-0.5,127.5,128,-0.5,127.5);
	}
	
	TH1F* QperHit[32];
	for(int i1=0;i1<32;i1++)
	{
		sprintf(hname,"QperHit_%d",i1);
		QperHit[i1]=new TH1F(hname,hname,500,0,1000);
	}
	
	TF1 *func = new TF1("fit",fitf2,0.,1000.,5);
	TF1 *exp1 = new TF1("exp1","[0]*exp(-x*[1])",0.,2400.);
	TF1 *land1 = new TF1("land1","landau",0.,100000.);
	TF1 *gplusl = new TF1("gplusl","gaus(0)+landau(3)",0.,1000.);
	TF1 *gaus1 = new TF1("gaus1","gaus",0.,1000.);
	
	bool foundfirst=true;bool foundlast=true;
	int tmin=5000;int indtmin=0;
	int tmax=0;int indtmax=0;
	double trackQ=0.;
	vector <float> tts;
	
	float x=0;float y=0;float ttmin=0;float ttmax=0;
	bool TrackOK=true;
	int trbin=0;
	float dr=0.;
	int NHitinTrack=0;int NHitoutsideTrack=0;
	vector <int> trhitinds;
	
	float PMTXY[3][2]={{96.5,63.5},{63.5,30.5},{63.5,96.5}};
	int I=0;int J=0;
	
	for(I=0;I<TTrack->GetEntries();I++)
	{
		TTrack->GetEntry(I);
// 		for(J=((I+5)<THit->GetEntries()?(I+5):THit->GetEntries());J>=0;J--)
		for(J=I;J<THit->GetEntries();J++)
		{
			THit->GetEntry(J);
			if(td.E==hd.E) break;
		}
// 		if(J==-1){cout<<I<<" / "<<TTrack->GetEntries()<<" : "<<I<<" "<<J<<endl;continue;}
// 		cout<<I<<" / "<<TTrack->GetEntries()<<" : "<<I<<" "<<J<<endl;
		
		tmin=5000;tmax=0;
		indtmin=0;indtmax=0;
		trackQ=0.;
		NHitinTrack=0;NHitoutsideTrack=0;
		trhitinds.clear();
		for(int i1=0;i1<hd.ColIndT->size();i1++)
		{
			AllHitQ->Fill(hd.Int->at(i1)[2]);
			AllHitX->Fill(hd.ColIndT->at(i1)[0]);
			AllHitY->Fill(hd.ColIndT->at(i1)[1]);
			AllHitT->Fill(hd.ColIndT->at(i1)[2]);
		
			x=td.FitParams->at(0)+td.FitParams->at(1)*(hd.ColIndT->at(i1)[2]-td.FitParams->at(4));
			y=td.FitParams->at(2)+td.FitParams->at(3)*(hd.ColIndT->at(i1)[2]-td.FitParams->at(4));
			if(sqrt(pow(hd.ColIndT->at(i1)[0]-x,2.)+pow(hd.ColIndT->at(i1)[1]-y,2.))<10)
			{
				if(hd.ColIndT->at(i1)[2]>tmax) {tmax=hd.ColIndT->at(i1)[2];indtmax=i1;}
				if(hd.ColIndT->at(i1)[2]<tmin) {tmin=hd.ColIndT->at(i1)[2];indtmin=i1;}
				TrackHitCharge->Fill(hd.Int->at(i1)[2]);
				NHitinTrack++;
				trhitinds.push_back(i1);
				trackQ+=hd.Int->at(i1)[2];
			}
			else
			{
				OutOfTrackHitCharge->Fill(hd.Int->at(i1)[2]);
				NHitoutsideTrack++;
			}
		}
		NHitsinvsoutT->Fill(NHitoutsideTrack,NHitinTrack);
		TrackOK=false;
		
		if(hd.ColIndT->at(indtmin)[0]>10 && hd.ColIndT->at(indtmin)[0]<117 && hd.ColIndT->at(indtmin)[1]>10 && hd.ColIndT->at(indtmin)[1]<117 && hd.ColIndT->at(indtmax)[0]>10 && hd.ColIndT->at(indtmax)[0]<117 && hd.ColIndT->at(indtmax)[1]>10 && hd.ColIndT->at(indtmax)[1]<117) TrackOK=true;
		if(NHitinTrack<50) TrackOK=false;
		
		dr=sqrt(pow(hd.ColIndT->at(indtmin)[0]-hd.ColIndT->at(indtmax)[0],2.)+pow(hd.ColIndT->at(indtmin)[1]-hd.ColIndT->at(indtmax)[1],2.)+pow(hd.ColIndT->at(indtmin)[2]-hd.ColIndT->at(indtmax)[2],2.));
// 		if(dr<700) TrackOK=false;
		
		EntryXY->Fill(hd.ColIndT->at(indtmin)[0],hd.ColIndT->at(indtmin)[1]);
		ExitXY->Fill(hd.ColIndT->at(indtmax)[0],hd.ColIndT->at(indtmax)[1]);
		
		if(TrackOK)
		{
			for(int i1=0;i1<trhitinds.size();i1++)
			{
				Qvst->Fill(hd.ColIndT->at(trhitinds[i1])[2],hd.Int->at(trhitinds[i1])[2]);
				QvstPr->Fill(hd.ColIndT->at(trhitinds[i1])[2],hd.Int->at(trhitinds[i1])[2]);
				if(hd.ColIndT->at(trhitinds[i1])[2]>=800 && hd.ColIndT->at(trhitinds[i1])[2]<2400)
				{
					QperHit[(hd.ColIndT->at(trhitinds[i1])[2]-800)/50]->Fill(hd.Int->at(trhitinds[i1])[2]);
				}
			}
			
			DeltaT->Fill(hd.ColIndT->at(indtmax)[2]-hd.ColIndT->at(indtmin)[2]);
// 			if((FE[i1].H[indtmax].t-FE[i1].H[indtmin].t)<500) cout<<FE[i1].E<<" "<<(FE[i1].H[indtmax].t-FE[i1].H[indtmin].t)<<" : "<<FE[i1].H[indtmin].t<<" "<<FE[i1].H[indtmax].t<<" "<<FE[i1].H.size()<<" "<<FE[i1].fitnormchi2<<endl;
			EntryT->Fill(hd.ColIndT->at(indtmin)[2]);
			ExitT->Fill(hd.ColIndT->at(indtmax)[2]);
			Trackt0->Fill(td.FitParams->at(4));
			TrackQ->Fill(trackQ);
			TrackDeltaR->Fill(dr);
			TrackDeltaXY->Fill(sqrt(pow(hd.ColIndT->at(indtmin)[0]-hd.ColIndT->at(indtmax)[0],2.)+pow(hd.ColIndT->at(indtmin)[1]-hd.ColIndT->at(indtmax)[1],2.)));
			Dt1vsDt2->Fill(td.FitParams->at(4)-hd.ColIndT->at(indtmin)[2],hd.ColIndT->at(indtmax)[2]-td.FitParams->at(4));
		
// 			tts.clear();
// 			tts.push_back(FE[i1].fitParams[4]-(FE[i1].fitParams[0]/FE[i1].fitParams[1]));
// 			tts.push_back(FE[i1].fitParams[4]-(FE[i1].fitParams[2]/FE[i1].fitParams[3]));
// 			tts.push_back(FE[i1].fitParams[4]+((127-FE[i1].fitParams[0])/FE[i1].fitParams[1]));
// 			tts.push_back(FE[i1].fitParams[4]+((127-FE[i1].fitParams[2])/FE[i1].fitParams[3]));
// 			sort(tts.begin(), tts.end());
// 			
// 			DeltaTvsTrBin->Fill(trbin/20,FE[i1].H[indtmax].t-FE[i1].H[indtmin].t);
// 			trbin++;
// 			DeltaTvsTr->SetPoint(DeltaTvsTr->GetN(),DeltaTvsTr->GetN(),FE[i1].H[indtmax].t-FE[i1].H[indtmin].t);
		
// 		ttmin=((FE[i1].fitParams[4]-(FE[i1].fitParams[0]/FE[i1].fitParams[1]))>(FE[i1].fitParams[4]-(FE[i1].fitParams[2]/FE[i1].fitParams[3]))?(FE[i1].fitParams[4]-(FE[i1].fitParams[0]/FE[i1].fitParams[1])):(FE[i1].fitParams[4]-(FE[i1].fitParams[2]/FE[i1].fitParams[3])));
// 		ttmax=((FE[i1].fitParams[4]+((127-FE[i1].fitParams[0])/FE[i1].fitParams[1]))>(FE[i1].fitParams[4]+((127-FE[i1].fitParams[2])/FE[i1].fitParams[3]))?(FE[i1].fitParams[4]+((127-FE[i1].fitParams[2])/FE[i1].fitParams[3])):(FE[i1].fitParams[4]+((127-FE[i1].fitParams[0])/FE[i1].fitParams[1])));
// 		TrackDeltaT->Fill(ttmax-ttmin);
// 			TrackDeltaT->Fill(tts[2]-tts[1]);
			
			DeltaTvsDeltaR->Fill(sqrt(pow(hd.ColIndT->at(indtmin)[0]-hd.ColIndT->at(indtmax)[0],2.)+pow(hd.ColIndT->at(indtmin)[1]-hd.ColIndT->at(indtmax)[1],2.)),hd.ColIndT->at(indtmax)[2]-hd.ColIndT->at(indtmin)[2]);
		
// 		cout<<FE[i1].E<<" "<<FE[i1].PMTint[0]<<" "<<FE[i1].PMTint[1]<<" "<<FE[i1].PMTint[2]<<endl;
			for(int i2=0;i2<3;i2++)
			{
				PMTIntAll[i2]->Fill(td.PMTIntegral->at(i2));
				if(td.PMTIntegral->at(i2)>=0)
				{
					PMTInt[i2]->Fill(td.PMTIntegral->at(i2));
					PMTIntPE[i2]->Fill(td.PMTIntegral->at(i2)/PMTPEs[i2]);
					EntryXYPMT[i2]->Fill(hd.ColIndT->at(indtmin)[0],hd.ColIndT->at(indtmin)[1]);
					MeanPMTChargeEntryXY[i2]->Fill(hd.ColIndT->at(indtmin)[0],hd.ColIndT->at(indtmin)[1],td.PMTIntegral->at(i2));
					PMTProf[i2]->Fill(sqrt(pow(hd.ColIndT->at(indtmin)[0]-PMTXY[i2][0],2.)+pow(hd.ColIndT->at(indtmin)[1]-PMTXY[i2][1],2.)),td.PMTIntegral->at(i2));
					PMTProfPE[i2]->Fill(sqrt(pow(hd.ColIndT->at(indtmin)[0]-PMTXY[i2][0],2.)+pow(hd.ColIndT->at(indtmin)[1]-PMTXY[i2][1],2.)),td.PMTIntegral->at(i2)/PMTPEs[i2]);
				}
				else
				{
					EntryXYPMTSat[i2]->Fill(hd.ColIndT->at(indtmin)[0],hd.ColIndT->at(indtmin)[1]);
				}
			}
		}
// 		for(int i2=0;i2<hd.ColIndT->size();i2++)
		for(int i2=0;i2<trhitinds.size();i2++)
		{
			
			if(hd.ColIndT->at(trhitinds[i2])[2]<800 || hd.ColIndT->at(trhitinds[i2])[2]>=2400) continue;
			TrackHits[(hd.ColIndT->at(trhitinds[i2])[2]-800)/200][0]->Fill(hd.ColIndT->at(trhitinds[i2])[0],hd.ColIndT->at(trhitinds[i2])[1],hd.Int->at(trhitinds[i2])[2]);
			TrackHits[(hd.ColIndT->at(trhitinds[i2])[2]-800)/200][1]->Fill(hd.ColIndT->at(trhitinds[i2])[0],hd.ColIndT->at(trhitinds[i2])[1]);
		}
	}
	outroot->cd();
	
	func->SetParLimits(0,0,1000.);
	func->SetParLimits(1,0,1000.);
	func->SetParLimits(2,0,100.);
	func->SetParLimits(3,-500,500);
	func->SetParLimits(4,-100,100);
	
	func->SetParameter(0,DeltaT->GetBinContent(DeltaT->GetMaximumBin()));
	func->SetParameter(1,DeltaT->GetBinCenter(DeltaT->GetMaximumBin()));
	func->SetParameter(2,10);
	func->SetParameter(3,0.1);
	func->SetParameter(4,0.1);
	DeltaT->Fit(func,"q","q",0.,1000.);
	
	QvstPr->Fit(exp1,"q","q",1000.,2400.);
	
// 	gplusl->SetParameter(0,AllHitQ->GetBinContent(AllHitQ->GetMaximumBin()));
// 	gplusl->SetParameter(1,AllHitQ->GetBinCenter(AllHitQ->GetMaximumBin()));
// 	gplusl->SetParameter(2,10);
// 	gplusl->SetParameter(3,AllHitQ->GetBinContent(AllHitQ->GetMaximumBin()));
// 	gplusl->SetParameter(4,AllHitQ->GetBinCenter(AllHitQ->GetMaximumBin()));
// 	gplusl->SetParameter(5,10);
// 	AllHitQ->Fit(gplusl,"q","q",0.,1000.);
	AllHitQ->Fit(land1,"q","q",0.,1000.);
	
	TrackQ->Fit(land1,"q","q",0.,100000.);
	
	EntryXY->Write();
	EntryT->Write();
	ExitXY->Write();
	ExitT->Write();
	DeltaT->Write();
	Trackt0->Write();
	TrackDeltaT->Write();
	AllHitX->Write();
	AllHitY->Write();
	AllHitT->Write();
	AllHitQ->Write();
	TrackQ->Write();
	TrackDeltaR->Write();
	TrackDeltaXY->Write();
	Dt1vsDt2->Write();
	Qvst->Write();
	QvstPr->Write();
	DeltaTvsDeltaR->Write();
	DeltaTvsTrBin->Write();
	DeltaTvsTr->Write();
	                                                                
	TrackHitCharge->Write();
	OutOfTrackHitCharge->Write();
	NHitsinvsoutT->Write();
	
	for(int i1=0;i1<3;i1++)
	{
		PMTInt[i1]->Scale(1./PMTInt[i1]->Integral());
		PMTInt[i1]->Write();
		PMTIntAll[i1]->Write();
		PMTIntPE[i1]->Write();
		EntryXYPMT[i1]->Write();
		EntryXYPMTSat[i1]->Write();
		MeanPMTChargeEntryXY[i1]->Divide(EntryXYPMT[i1]);
		MeanPMTChargeEntryXY[i1]->Write();
		PMTProf[i1]->Write();
		PMTProfPE[i1]->Write();
	}
	
	for(int i1=0;i1<9;i1++)
	{
		TrackHits[i1][0]->Divide(TrackHits[i1][1]);
		TrackHits[i1][0]->Write();
		TrackHits[i1][1]->Write();
	}
	
	TGraphErrors* tge=new TGraphErrors();
	tge->SetName("MeanHitQvsTime");tge->SetTitle("MeanHitQvsTime");
	tge->SetMarkerStyle(20);tge->SetMarkerColor(1);
	
	for(int i1=0;i1<32;i1++)
	{
// 		gplusl->SetParLimits(0,0,10000);
// 		gplusl->SetParLimits(1,0,100);
// 		gplusl->SetParLimits(2,0,100);
// 		gplusl->SetParLimits(3,0,1000);
// 		gplusl->SetParLimits(4,0,100);
// 		gplusl->SetParLimits(5,0,100);
// 		
// 		gplusl->SetParameter(0,QperHit[i1]->GetBinContent(QperHit[i1]->GetMaximumBin()));
// 		gplusl->SetParameter(1,QperHit[i1]->GetBinCenter(QperHit[i1]->GetMaximumBin()));
// 		gplusl->SetParameter(2,100);
// 		gplusl->SetParameter(3,QperHit[i1]->GetBinContent(QperHit[i1]->GetMaximumBin()));
// 		gplusl->SetParameter(4,QperHit[i1]->GetBinCenter(QperHit[i1]->GetMaximumBin()));
// 		gplusl->SetParameter(5,10);
		if(QperHit[i1]->GetEntries()>0)
		{
			if(QperHit[i1]->GetEntries()<2000) QperHit[i1]->Rebin(4);
			QperHit[i1]->Fit(gaus1,"q","q",QperHit[i1]->GetMean()-QperHit[i1]->GetRMS(),QperHit[i1]->GetMean()+QperHit[i1]->GetRMS());
			QperHit[i1]->Fit(gaus1,"q","q",gaus1->GetParameter(1)-gaus1->GetParameter(2),gaus1->GetParameter(1)+gaus1->GetParameter(2));
			tge->SetPoint(tge->GetN(),i1,gaus1->GetParameter(1));
			tge->SetPointError(tge->GetN()-1,0,gaus1->GetParError(1));
		}
		
// 		QperHit[i1]->Fit(gplusl,"q","q",0.,1000.);
// 		QperHit[i1]->Fit(land1,"q","q",0.,1000.);
		QperHit[i1]->Write();
	}
	
	tge->Fit(exp1,"q","q",1.,8.);
	tge->Write();
	
	outroot->Close();
// 	outfile.close();
	
	inroot1->Close();
	inroot2->Close();
	
// 	sprintf(hname,"rm TrackParameters_%d.txt",RunNo);system(hname);
// 	sprintf(hname,"rm PMTIntegrals_%d.txt",RunNo);system(hname);
// 	sprintf(hname,"rm Hits_%d.txt",RunNo);system(hname);
	sprintf(hname,"cp TrackAnalysis2_%d.root %s/Histos/TrackAnalysis2_%d.root;wait;",RunNo,AnalysisFilePath,RunNo);system(hname);
// 	sprintf(hname,"mv ValidTrackParameters_%d.txt %s/Files/ValidTrackParameters_%d.txt;wait;",RunNo,AnalysisFilePath,RunNo);system(hname);
}

void MultiViewPlot()
{
// 	vector <fullevent> FE;
// 	sprintf(hname,"cp %s/Files/NHitsperEvent_%d.txt .;wait;",AnalysisFilePath,RunNo);system(hname);
// 	sprintf(hname,"NHitsperEvent_%d.txt",RunNo);
// 	ifstream inNHits(hname);
// 	while(!inNHits.eof())
// 	{
// 		inNHits>>fe.E>>fe.clsize;
// 		fe.fitParams[0]=0;fe.fitParams[1]=0;fe.fitParams[2]=0;fe.fitParams[3]=0;fe.fitParams[4]=0;
// 		fe.fitnormchi2=-1;
// 		fe.PMTint[0]=0;fe.PMTint[1]=0;fe.PMTint[2]=0;
// 		FE.push_back(fe);
// 	}
// 	inNHits.close();
// 	if(FE[FE.size()-1].E==FE[FE.size()-2].E) FE.pop_back();
// 	
// 	struct fullevent
// 	{
// 		int E;
// 		float fitParams[5];
// 		float fitnormchi2;
// 		int clsize;
// 		float PMTint[3];
// 		vector <hits> H;
// 	};
// 	
// 	sprintf(hname,"cp %s/Files/TrackParameters_%d.txt .;wait;",AnalysisFilePath,RunNo);system(hname);
// 	sprintf(hname,"TrackParameters_%d.txt",RunNo);
// 	ifstream inTracks(hname);
// 	while(!inTracks.eof())
// 	{
// 		inTracks>>fe.E>>fe.fitParams[0]>>fe.fitParams[1]>>fe.fitParams[2]>>fe.fitParams[3]>>fe.fitParams[4]>>fe.fitnormchi2>>fe.clsize;
// 		for(int i1=0;i1<FE.size();i1++)
// 		{
// 			if(FE[i1].E==fe.E)
// 			{
// 				for(int i2=0;i2<5;i2++){FE[i1].fitParams[i2]=fe.fitParams[i2];}
// 				FE[i1].fitnormchi2=fe.fitnormchi2;
// 				FE[i1].clsize=fe.clsize;//if there is a track, clsize is the # of hits in the fit
// 				break;
// 			}
// 		}
// 	}
// 	inTracks.close();
// 	
// 	int nEvPMT=0;
// 	sprintf(hname,"cp %s/Files/PMTIntegrals_%d.txt .;wait;",AnalysisFilePath,RunNo);system(hname);
// 	sprintf(hname,"PMTIntegrals_%d.txt",RunNo);
// 	ifstream inPMTs(hname);
// 	while(!inPMTs.eof())
// 	{
// 		inPMTs>>fe.E>>fe.PMTint[0]>>fe.PMTint[1]>>fe.PMTint[2];
// 		for(int i1=0;i1<FE.size();i1++)
// 		{
// 			if(FE[i1].E==fe.E)
// 			{
// 				FE[i1].PMTint[0]=fe.PMTint[0];
// 				FE[i1].PMTint[1]=fe.PMTint[1];
// 				FE[i1].PMTint[2]=fe.PMTint[2];
// 				break;
// 			}
// 		}
// 		nEvPMT++;
// 	}
// 	inPMTs.close();
// 	
// 	sprintf(hname,"cp %s/Files/Hits_%d.txt .;wait;",AnalysisFilePath,RunNo);system(hname);
// 	sprintf(hname,"Hits_%d.txt",RunNo);
// 	ifstream inHits(hname);
// 	int b[7]={0};
// 	int Ev=0;
// 	bool found=false;
// 	bool skipEvent=false;
// 	int curInd=0;
// 	
// 	while(!inHits.eof())
// 	{
// 		inHits>>b[0]>>b[1]>>b[2]>>b[3]>>b[4]>>b[5]>>b[6];
// 		hit.col=b[1];
// 		hit.ind=b[2];
// 		hit.t=b[3];
// 		hit.Int=b[4];
// 		hit.Int2=b[5];
// 		hit.Int3=b[6];
// 		for(int i1=0;i1<FE.size();i1++)
// 		{
// 			if(FE[i1].E==b[0])
// 			{
// 				FE[i1].H.push_back(hit);
// 				break;
// 			}
// 		}
// 	}
// 	inHits.close();
// 	
// 	sprintf(hname,"cp %s/Histos/PMTWaveforms_%d.root .;wait;",AnalysisFilePath,RunNo);system(hname);
// 	sprintf(hname,"PMTWaveforms_%d.root",RunNo);
// 	TFile* inrootPMTWF=new TFile(hname);
// 	TH1F* hPMT[3];
// 	string PMTNames[3]={"LAr w/TPB","LAr no TPB","LXe"};
// 	int PMTlc[3]={1,2,4};
// 	
// 	sprintf(hname,"cp %s/Histos/TPCWaveforms2D_%d.root .;wait;",AnalysisFilePath,RunNo);system(hname);
// 	sprintf(hname,"TPCWaveforms2D_%d.root",RunNo);
// 	TFile* inrootTPCWF2D=new TFile(hname);
// 	TH2F* hTPCWF;
// 	
// // 	TCanvas * cca1=new TCanvas("cca1","cca1",300,300);
// // 	TCanvas * cca2=new TCanvas("cca2","cca2",300,300);
// // 	TCanvas * cca3=new TCanvas("cca3","cca3",300,300);
// 	TCanvas * cca1;
// 	TCanvas * cca2;
// 	TCanvas * cca3;
// 	gROOT->GetStyle("Default")->SetPalette(kRainBow);
// 	
// 	TPolyLine3D *l;
// 	TH3F *frame3d=new TH3F("Hits","Hits",1,-0.5,127.5,1,-0.5,127.5,1,-0.5,4095.5);
// 	TPolyMarker3D *pm3d1;
// 	TPolyMarker3D *pm3d2;
// 	double x(0),y(0),z(0),t(0);
// 	double tlast(0);
// 	
// 	bool foundfirst=true;bool foundlast=true;
// 	int tmin=5000;int indtmin=0;
// 	int tmax=0;int indtmax=0;
// 	
// // 	for(int i1=0;i1<FE.size();i1++)
// 	for(int i1=0;i1<FE.size()-2;i1++)
// // 	for(int i1=0;i1<10;i1++)
// 	{
// 		cca1=new TCanvas("cca1","cca1",300,300);
// 		cca2=new TCanvas("cca2","cca2",300,300);
// 		cca3=new TCanvas("cca3","cca3",300,300);
// 		
// 		sprintf(hname,"Ev2D_%d",FE[i1].E);
// 		inrootTPCWF2D->GetObject(hname,hTPCWF);
// 		cca1->cd();
// 		gStyle->SetOptStat(0);
// 		hTPCWF->Draw("colz");
// 		
// 		pm3d1 = new TPolyMarker3D(FE[i1].H.size());
// 		for(int i2=0;i2<FE[i1].H.size();i2++)
// 		{
// 			pm3d1->SetPoint(i2,FE[i1].H[i2].col,FE[i1].H[i2].ind,FE[i1].H[i2].t);
// 		}
// 		
// 		cca2->cd();
// 		gStyle->SetOptStat(0);
// 		frame3d->Draw();
// 		pm3d1->SetMarkerColor(kRed);
// 		pm3d1->SetMarkerStyle(24);   
// 		pm3d1->Draw();
// 		if(FE[i1].fitnormchi2>=0)
// 		{
// 			l = new TPolyLine3D(2);
// 			
// 			t=((FE[i1].fitParams[4]-(FE[i1].fitParams[0]/FE[i1].fitParams[1]))>(FE[i1].fitParams[4]-(FE[i1].fitParams[2]/FE[i1].fitParams[3]))?(FE[i1].fitParams[4]-(FE[i1].fitParams[0]/FE[i1].fitParams[1])):(FE[i1].fitParams[4]-(FE[i1].fitParams[2]/FE[i1].fitParams[3])));
// 			x=FE[i1].fitParams[0]+FE[i1].fitParams[1]*(t-FE[i1].fitParams[4]);
// 			y=FE[i1].fitParams[2]+FE[i1].fitParams[3]*(t-FE[i1].fitParams[4]);
// 			l->SetPoint(0,x,y,t);
// 			t=((FE[i1].fitParams[4]+((127-FE[i1].fitParams[0])/FE[i1].fitParams[1]))>(FE[i1].fitParams[4]+((127-FE[i1].fitParams[2])/FE[i1].fitParams[3]))?(FE[i1].fitParams[4]+((127-FE[i1].fitParams[2])/FE[i1].fitParams[3])):(FE[i1].fitParams[4]+((127-FE[i1].fitParams[0])/FE[i1].fitParams[1])));
// 			x=FE[i1].fitParams[0]+FE[i1].fitParams[1]*(t-FE[i1].fitParams[4]);
// 			y=FE[i1].fitParams[2]+FE[i1].fitParams[3]*(t-FE[i1].fitParams[4]);
// 			l->SetPoint(1,x,y,t);
// 			
// 			l->SetLineColor(kBlack);
// 			l->SetLineWidth(2);
// 			l->Draw("same");
// 		}
// 		
// 		double PMTmaxAmp=0.;
// 		for(int is1=0;is1<3;is1++)
// 		{
// 			sprintf(hname,"Ch%d/WF_PMT_%d_Ev_%d",is1,is1,FE[i1].E);
// 			inrootPMTWF->GetObject(hname,hPMT[is1]);
// 			sprintf(hname2,"%s %s",PMTNames[is1].c_str(),hPMT[is1]->GetTitle());
// 			hPMT[is1]->SetTitle(hname2);
// 			hPMT[is1]->SetLineColor(PMTlc[is1]);
// 			if(hPMT[is1]->GetBinContent(hPMT[is1]->GetMaximumBin())>PMTmaxAmp) PMTmaxAmp=hPMT[is1]->GetBinContent(hPMT[is1]->GetMaximumBin());
// 		}
// 		PMTmaxAmp*=1.2;
// 		cca3->cd();
// 		gStyle->SetOptStat(0);
// 		
// 		sprintf(hname,"%5.2f %5.2f %5.2f",FE[i1].PMTint[0],FE[i1].PMTint[1],FE[i1].PMTint[2]);
// 		hPMT[0]->SetTitle(hname);
// 		hPMT[0]->Draw();
// 		hPMT[0]->GetYaxis()->SetRangeUser(1,PMTmaxAmp);
// 		gPad->SetLogy(1);
// 		hPMT[1]->Draw("same");
// 		hPMT[2]->Draw("same");
// 		
// 		sprintf(hname,"F_1_E_%d.png",FE[i1].E);cca1->SaveAs(hname);
// 		sprintf(hname,"F_2_E_%d.png",FE[i1].E);cca2->SaveAs(hname);
// 		sprintf(hname,"F_3_E_%d.png",FE[i1].E);cca3->SaveAs(hname);
// 		
// 		delete hTPCWF;delete hPMT[0];delete hPMT[1];delete hPMT[2];
// 		delete pm3d1;
// 		delete cca1;delete cca2;delete cca3;
// 	}
// 	inrootPMTWF->Close();
// 	inrootTPCWF2D->Close();
// 	
// 	sprintf(hname,"rm TrackParameters_%d.txt",RunNo);system(hname);
// 	sprintf(hname,"rm PMTIntegrals_%d.txt",RunNo);system(hname);
// 	sprintf(hname,"rm Hits_%d.txt",RunNo);system(hname);
// 	sprintf(hname,"rm NHitsperEvent_%d.txt",RunNo);system(hname);
// // 	sprintf(hname,"if [ ! -d /eos/project/f/flic2019/www/Run1/%d ] then mkdir /eos/project/f/flic2019/www/Run1/%d;fi;wait;",RunNo,RunNo);system(hname);
// 	sprintf(hname,"mkdir /eos/project/f/flic2019/www/Run1/%d;wait;",RunNo);system(hname);
// 	sprintf(hname,"mv *.png /eos/project/f/flic2019/www/Run1/%d;wait;",RunNo);system(hname);
// 	
// 	FE.clear();
}

void PrintNevt()
{
	sprintf(hname,"%s/Run_%d.root",NTupleFilePath,RunNo);
	TFile* inroot=new TFile(hname);
	TTree* T =  (TTree*) inroot->Get("T");
	cout<<"Run : "<<RunNo<<" NEvt : "<<T->GetEntries()<<endl;
	inroot->Close();
}

void TrackAnalysis3()
{
	vector <int> RunNos;
	if(RunNo==1000)
	{
		RunNos.push_back(7074);
		RunNos.push_back(7076);
		RunNos.push_back(7078);
		RunNos.push_back(7082);
		RunNos.push_back(7084);
		RunNos.push_back(7092);
		RunNos.push_back(7094);
		RunNos.push_back(7096);
		RunNos.push_back(7104);
		RunNos.push_back(7106);
		RunNos.push_back(7107);
		RunNos.push_back(7110);
		RunNos.push_back(7121);
	}
	else if(RunNo==2000)
	{
		RunNos.push_back(7252);
		RunNos.push_back(7254);
		RunNos.push_back(7256);
		RunNos.push_back(7259);
		RunNos.push_back(7268);
		RunNos.push_back(7271);
		RunNos.push_back(7273);
	}
	else if(RunNo==3000)
	{
		RunNos.push_back(7308);
		RunNos.push_back(7312);
		RunNos.push_back(7314);
		RunNos.push_back(7316);
		RunNos.push_back(7318);
		RunNos.push_back(7320);
		RunNos.push_back(7322);
		RunNos.push_back(7325);
		RunNos.push_back(7327);
		RunNos.push_back(7385);
		RunNos.push_back(7388);
	}
	
	TH1F* hPMT[3];
// 	string PMTNames[3]={"LAr w/TPB","LAr no TPB","LXe"};
	int PMTlc[3]={1,2,4};
	
	sprintf(hname,"TrackAnalysis3_%d.root",RunNo);
	TFile* outroot=new TFile(hname,"recreate");
	
	TH2F* EntryXY=new TH2F("EntryXY","EntryXY",128,-0.5,127.5,128,-0.5,127.5);
	TH1F* EntryT=new TH1F("EntryT","EntryT",256,-0.5,4095.5);
	TH2F* ExitXY=new TH2F("ExitXY","ExitXY",128,-0.5,127.5,128,-0.5,127.5);
	TH1F* ExitT=new TH1F("ExitT","ExitT",256,-0.5,4095.5);
	TH1F* Trackt0=new TH1F("Trackt0","Trackt0",400,0,4000);
	TH1F* DeltaT=new TH1F("DeltaT","DeltaT",256,-0.5,4095.5);
	TH1F* TrackDeltaT=new TH1F("TrackDeltaT","TrackDeltaT",400,0,4000);
	TH1F* AllHitX=new TH1F("AllHitX","AllHitX",128,-0.5,127.5);
	TH1F* AllHitY=new TH1F("AllHitY","AllHitY",128,-0.5,127.5);
	TH1F* AllHitT=new TH1F("AllHitT","AllHitT",4096,-0.5,4095.5);
	TH1F* AllHitQ=new TH1F("AllHitQ","AllHitQ",500,0,1000);
	TH1F* TrackQ=new TH1F("TrackQ","TrackQ",100,0,100000);
	TH1F* TrackDeltaR=new TH1F("TrackDeltaR","TrackDeltaR",1000,0,2500);
	TH1F* TrackDeltaXY=new TH1F("TrackDeltaXY","TrackDeltaXY",200,0,200);
	TH2F* Qvst=new TH2F("Qvst","Qvst",1400,1000,2400,1000,0,1000);
	TH2F* Dt1vsDt2=new TH2F("Dt1vsDt2","Dt1vsDt2",500,-1000,1000,500,-1000,1000);
	TProfile* QvstPr=new TProfile("QvstPr","QvstPr",200,1000,2400,0,1000);
	TProfile* DeltaTvsDeltaR=new TProfile("DeltaTvsDeltaR","DeltaTvsDeltaR",200,0,180,0,1000);
	TProfile* DeltaTvsTrBin=new TProfile("DeltaTvsTrBin","DeltaTvsTrBin",200,0,1000,0,1000);
	TGraph* DeltaTvsTr=new TGraph();
	DeltaTvsTr->SetName("DeltaTvsTr");DeltaTvsTr->SetTitle("DeltaTvsTr");
	DeltaTvsTr->SetMarkerStyle(20);DeltaTvsTr->SetMarkerColor(1);
	TH1F* TrackHitCharge=new TH1F("TrackHitCharge","TrackHitCharge",1000,0,10000);
	TH1F* OutOfTrackHitCharge=new TH1F("OutOfTrackHitCharge","OutOfTrackHitCharge",1000,0,10000);
	TH2F* NHitsinvsoutT=new TH2F("NHitsinvsoutT","NHitsinvsoutT",500,-0.5,499.5,500,-0.5,499.5);
	
	TH2F* EntryXYPMT[3];TH2F* EntryXYPMTSat[3];TH2F* MeanPMTChargeEntryXY[3];
	TH1F* PMTInt[3];TH1F* PMTIntPE[3];TH1F* PMTIntAll[3];
	TProfile* PMTProf[3];TProfile* PMTProfPE[3];
	int PMTcols[3]={1,2,4};
	float PMTPEs[3]={34,1,23};
	for(int i1=0;i1<3;i1++)
	{
		sprintf(hname,"PMTIntAll_%d",i1);
		PMTIntAll[i1]=new TH1F(hname,hname,2000,-1000000,1000000);
		PMTIntAll[i1]->SetLineColor(PMTcols[i1]);
		sprintf(hname,"PMTInt_%d",i1);
		PMTInt[i1]=new TH1F(hname,hname,200,0,100000);
		PMTInt[i1]->SetLineColor(PMTcols[i1]);
		sprintf(hname,"PMTIntPE_%d",i1);
		PMTIntPE[i1]=new TH1F(hname,hname,300,0,3000);
		PMTIntPE[i1]->SetLineColor(PMTcols[i1]);
		sprintf(hname,"EntryXYPMT_%d",i1);
		EntryXYPMT[i1]=new TH2F(hname,hname,128,-0.5,127.5,128,-0.5,127.5);
		sprintf(hname,"EntryXYPMTSat_%d",i1);
		EntryXYPMTSat[i1]=new TH2F(hname,hname,128,-0.5,127.5,128,-0.5,127.5);
		sprintf(hname,"MeanPMTChargeEntryXY_%d",i1);
		MeanPMTChargeEntryXY[i1]=new TH2F(hname,hname,128,-0.5,127.5,128,-0.5,127.5);
		sprintf(hname,"PMTChageProf_%d",i1);
		PMTProf[i1]=new TProfile(hname,hname,50,0,150,0,100000);
		PMTProf[i1]->SetLineColor(PMTcols[i1]);
		sprintf(hname,"PMTPEProf_%d",i1);
		PMTProfPE[i1]=new TProfile(hname,hname,50,0,150,0,3000);
		PMTProfPE[i1]->SetLineColor(PMTcols[i1]);
	}
	
	TH2F* TrackHits[9][2];
	for(int i1=0;i1<9;i1++)
	{
		sprintf(hname,"MeanCharge_%d-%d",(i1*200)+800,((i1+1)*200)+800);
		TrackHits[i1][0]=new TH2F(hname,hname,128,-0.5,127.5,128,-0.5,127.5);
		sprintf(hname,"NHits_%d-%d",(i1*200)+800,((i1+1)*200)+800);
		TrackHits[i1][1]=new TH2F(hname,hname,128,-0.5,127.5,128,-0.5,127.5);
	}
	
	TH1F* QperHit[32];
	for(int i1=0;i1<32;i1++)
	{
		sprintf(hname,"QperHit_%d",i1);
		QperHit[i1]=new TH1F(hname,hname,500,0,1000);
	}
	
	TF1 *func = new TF1("fit",fitf2,0.,1000.,5);
	TF1 *exp1 = new TF1("exp1","[0]*exp(-x*[1])",0.,2400.);
	TF1 *land1 = new TF1("land1","landau",0.,100000.);
	TF1 *gplusl = new TF1("gplusl","gaus(0)+landau(3)",0.,1000.);
	TF1 *gaus1 = new TF1("gaus1","gaus",0.,1000.);
	
	bool foundfirst=true;bool foundlast=true;
	int tmin=5000;int indtmin=0;
	int tmax=0;int indtmax=0;
	double trackQ=0.;
	vector <float> tts;
	
	float x=0;float y=0;float ttmin=0;float ttmax=0;
	bool TrackOK=true;
	int trbin=0;
	float dr=0.;
	int NHitinTrack=0;int NHitoutsideTrack=0;
	vector <int> trhitinds;
	
// 	float PMTXY[3][2]={{96.5,63.5},{63.5,30.5},{63.5,96.5}};
	float PMTXY[3][2]={{95.,63.5},{63.5,32.},{63.5,95.}};
	int I=0;int J=0;int II=0;
	
	TFile* inroot;
	TTree* TH;
	TTree* TT;
	bool PMTsat=false;
	for(int J=0;J<RunNos.size();J++)
	{
		RunNo=RunNos[J];
		cout<<"Run "<<RunNo<<endl;
		sprintf(hname,"%s/Histos/HitsAndTracks_%d.root",AnalysisFilePath,RunNo);
		inroot=new TFile(hname);
		TH =  (TTree*) inroot->Get("Hits");
		TH->SetBranchAddress("E",&hd.E);
		TH->SetBranchAddress("ColIndT",&hd.ColIndT);
		TH->SetBranchAddress("Int",&hd.Int);
		TH->SetBranchAddress("QColTot",&hd.QColTot);
		TH->SetBranchAddress("QColTotZS",&hd.QColTotZS);
		TH->SetBranchAddress("QHitTot",&hd.QHitTot);
		TH->SetBranchAddress("PMTIntegral",&hd.PMTIntegral);
		TH->SetBranchAddress("ColID",&hd.ColID);
		TH->SetBranchAddress("ColT",&hd.ColT);
		TH->SetBranchAddress("ColA",&hd.ColA);
		TH->SetBranchAddress("ColInt",&hd.ColInt);
		TH->SetBranchAddress("Colw",&hd.Colw);
		TH->SetBranchAddress("IndID",&hd.IndID);
		TH->SetBranchAddress("IndT",&hd.IndT);
		TH->SetBranchAddress("IndA",&hd.IndA);
		TH->SetBranchAddress("Indw",&hd.Indw);
		TH->SetBranchAddress("EventType",&hd.EventType);
		TT =  (TTree*) inroot->Get("Tracks");
		TT->SetBranchAddress("E",&td.E);
		TT->SetBranchAddress("StartEndColIndT",&td.StartEndColIndT);
		TT->SetBranchAddress("FitParams",&td.FitParams);
		TT->SetBranchAddress("FitNormChi2",&td.FitNormChi2);
		TT->SetBranchAddress("NHits",&td.NHits);
		TT->SetBranchAddress("Nexcl",&td.Nexcl);
		TT->SetBranchAddress("PMTIntegral",&td.PMTIntegral);
		TT->SetBranchAddress("ColTStartEnd",&td.ColTStartEnd);
		TT->SetBranchAddress("ColHitTStartEnd",&td.ColHitTStartEnd);
		
		for(int I=0;I<TT->GetEntries();I++)
		{
			cout<<I<<endl;
			TT->GetEntry(I);
			
			for(int II=td.E;II<TH->GetEntries();II++)
			{
				TH->GetEntry(II);
				if(hd.E==td.E) break;
			}
			TrackOK=true;
			if(hd.EventType!=0) TrackOK=false;//choose tracks traversing the entire drift distance
			if(td.Nexcl>=10) TrackOK=false;
			if(!TrackOK) continue;
			
			EntryXY->Fill(td.StartEndColIndT->at(0),td.StartEndColIndT->at(1));
			ExitXY->Fill(td.StartEndColIndT->at(3),td.StartEndColIndT->at(4));
			
// 			for(int i1=0;i1<trhitinds.size();i1++)
// 			{
// 				Qvst->Fill(hd.ColIndT->at(trhitinds[i1])[2],hd.Int->at(trhitinds[i1])[2]);
// 				QvstPr->Fill(hd.ColIndT->at(trhitinds[i1])[2],hd.Int->at(trhitinds[i1])[2]);
// 				if(hd.ColIndT->at(trhitinds[i1])[2]>=800 && hd.ColIndT->at(trhitinds[i1])[2]<2400)
// 				{
// 					QperHit[(hd.ColIndT->at(trhitinds[i1])[2]-800)/50]->Fill(hd.Int->at(trhitinds[i1])[2]);
// 				}
// 			}
// 			
// 			DeltaT->Fill(hd.ColIndT->at(indtmax)[2]-hd.ColIndT->at(indtmin)[2]);
// // 			if((FE[i1].H[indtmax].t-FE[i1].H[indtmin].t)<500) cout<<FE[i1].E<<" "<<(FE[i1].H[indtmax].t-FE[i1].H[indtmin].t)<<" : "<<FE[i1].H[indtmin].t<<" "<<FE[i1].H[indtmax].t<<" "<<FE[i1].H.size()<<" "<<FE[i1].fitnormchi2<<endl;
// 			EntryT->Fill(hd.ColIndT->at(indtmin)[2]);
// 			ExitT->Fill(hd.ColIndT->at(indtmax)[2]);
// 			Trackt0->Fill(td.FitParams->at(4));
// 			TrackQ->Fill(trackQ);
// 			TrackDeltaR->Fill(dr);
// 			TrackDeltaXY->Fill(sqrt(pow(hd.ColIndT->at(indtmin)[0]-hd.ColIndT->at(indtmax)[0],2.)+pow(hd.ColIndT->at(indtmin)[1]-hd.ColIndT->at(indtmax)[1],2.)));
// 			Dt1vsDt2->Fill(td.FitParams->at(4)-hd.ColIndT->at(indtmin)[2],hd.ColIndT->at(indtmax)[2]-td.FitParams->at(4));
			
// 			DeltaTvsDeltaR->Fill(sqrt(pow(hd.ColIndT->at(indtmin)[0]-hd.ColIndT->at(indtmax)[0],2.)+pow(hd.ColIndT->at(indtmin)[1]-hd.ColIndT->at(indtmax)[1],2.)),hd.ColIndT->at(indtmax)[2]-hd.ColIndT->at(indtmin)[2]);
		
// 		cout<<FE[i1].E<<" "<<FE[i1].PMTint[0]<<" "<<FE[i1].PMTint[1]<<" "<<FE[i1].PMTint[2]<<endl;
			for(int i2=0;i2<3;i2++)
			{
				PMTIntAll[i2]->Fill(td.PMTIntegral->at(i2));
				if(td.PMTIntegral->at(i2)>=0)
				{
					PMTInt[i2]->Fill(td.PMTIntegral->at(i2));
					PMTIntPE[i2]->Fill(td.PMTIntegral->at(i2)/PMTPEs[i2]);
					EntryXYPMT[i2]->Fill(td.StartEndColIndT->at(0),td.StartEndColIndT->at(1));
					MeanPMTChargeEntryXY[i2]->Fill(td.StartEndColIndT->at(0),td.StartEndColIndT->at(1),td.PMTIntegral->at(i2));
					PMTProf[i2]->Fill(sqrt(pow(td.StartEndColIndT->at(0)-PMTXY[i2][0],2.)+pow(td.StartEndColIndT->at(1)-PMTXY[i2][1],2.)),td.PMTIntegral->at(i2));
					PMTProfPE[i2]->Fill(sqrt(pow(td.StartEndColIndT->at(0)-PMTXY[i2][0],2.)+pow(td.StartEndColIndT->at(1)-PMTXY[i2][1],2.)),td.PMTIntegral->at(i2)/PMTPEs[i2]);
				}
				else
				{
					EntryXYPMTSat[i2]->Fill(td.StartEndColIndT->at(0),td.StartEndColIndT->at(1));
				}
			}
		}
		
		
		
	
	
	
	
	
	
	
// 	for(I=0;I<TTrack->GetEntries();I++)
// 	{
// 		TTrack->GetEntry(I);
// // 		for(J=((I+5)<THit->GetEntries()?(I+5):THit->GetEntries());J>=0;J--)
// 		for(J=I;J<THit->GetEntries();J++)
// 		{
// 			THit->GetEntry(J);
// 			if(td.E==hd.E) break;
// 		}
// // 		if(J==-1){cout<<I<<" / "<<TTrack->GetEntries()<<" : "<<I<<" "<<J<<endl;continue;}
// // 		cout<<I<<" / "<<TTrack->GetEntries()<<" : "<<I<<" "<<J<<endl;
// 		
// 		tmin=5000;tmax=0;
// 		indtmin=0;indtmax=0;
// 		trackQ=0.;
// 		NHitinTrack=0;NHitoutsideTrack=0;
// 		trhitinds.clear();
// 		for(int i1=0;i1<hd.ColIndT->size();i1++)
// 		{
// 			AllHitQ->Fill(hd.Int->at(i1)[2]);
// 			AllHitX->Fill(hd.ColIndT->at(i1)[0]);
// 			AllHitY->Fill(hd.ColIndT->at(i1)[1]);
// 			AllHitT->Fill(hd.ColIndT->at(i1)[2]);
// 		
// 			x=td.FitParams->at(0)+td.FitParams->at(1)*(hd.ColIndT->at(i1)[2]-td.FitParams->at(4));
// 			y=td.FitParams->at(2)+td.FitParams->at(3)*(hd.ColIndT->at(i1)[2]-td.FitParams->at(4));
// 			if(sqrt(pow(hd.ColIndT->at(i1)[0]-x,2.)+pow(hd.ColIndT->at(i1)[1]-y,2.))<10)
// 			{
// 				if(hd.ColIndT->at(i1)[2]>tmax) {tmax=hd.ColIndT->at(i1)[2];indtmax=i1;}
// 				if(hd.ColIndT->at(i1)[2]<tmin) {tmin=hd.ColIndT->at(i1)[2];indtmin=i1;}
// 				TrackHitCharge->Fill(hd.Int->at(i1)[2]);
// 				NHitinTrack++;
// 				trhitinds.push_back(i1);
// 				trackQ+=hd.Int->at(i1)[2];
// 			}
// 			else
// 			{
// 				OutOfTrackHitCharge->Fill(hd.Int->at(i1)[2]);
// 				NHitoutsideTrack++;
// 			}
// 		}
// 		NHitsinvsoutT->Fill(NHitoutsideTrack,NHitinTrack);
// 		TrackOK=false;
// 		
// 		if(hd.ColIndT->at(indtmin)[0]>10 && hd.ColIndT->at(indtmin)[0]<117 && hd.ColIndT->at(indtmin)[1]>10 && hd.ColIndT->at(indtmin)[1]<117 && hd.ColIndT->at(indtmax)[0]>10 && hd.ColIndT->at(indtmax)[0]<117 && hd.ColIndT->at(indtmax)[1]>10 && hd.ColIndT->at(indtmax)[1]<117) TrackOK=true;
// 		if(NHitinTrack<50) TrackOK=false;
// 		
// 		dr=sqrt(pow(hd.ColIndT->at(indtmin)[0]-hd.ColIndT->at(indtmax)[0],2.)+pow(hd.ColIndT->at(indtmin)[1]-hd.ColIndT->at(indtmax)[1],2.)+pow(hd.ColIndT->at(indtmin)[2]-hd.ColIndT->at(indtmax)[2],2.));
// // 		if(dr<700) TrackOK=false;
// 		
// 		EntryXY->Fill(hd.ColIndT->at(indtmin)[0],hd.ColIndT->at(indtmin)[1]);
// 		ExitXY->Fill(hd.ColIndT->at(indtmax)[0],hd.ColIndT->at(indtmax)[1]);
// 		
// 		if(TrackOK)
// 		{
// 			for(int i1=0;i1<trhitinds.size();i1++)
// 			{
// 				Qvst->Fill(hd.ColIndT->at(trhitinds[i1])[2],hd.Int->at(trhitinds[i1])[2]);
// 				QvstPr->Fill(hd.ColIndT->at(trhitinds[i1])[2],hd.Int->at(trhitinds[i1])[2]);
// 				if(hd.ColIndT->at(trhitinds[i1])[2]>=800 && hd.ColIndT->at(trhitinds[i1])[2]<2400)
// 				{
// 					QperHit[(hd.ColIndT->at(trhitinds[i1])[2]-800)/50]->Fill(hd.Int->at(trhitinds[i1])[2]);
// 				}
// 			}
// 			
// 			DeltaT->Fill(hd.ColIndT->at(indtmax)[2]-hd.ColIndT->at(indtmin)[2]);
// // 			if((FE[i1].H[indtmax].t-FE[i1].H[indtmin].t)<500) cout<<FE[i1].E<<" "<<(FE[i1].H[indtmax].t-FE[i1].H[indtmin].t)<<" : "<<FE[i1].H[indtmin].t<<" "<<FE[i1].H[indtmax].t<<" "<<FE[i1].H.size()<<" "<<FE[i1].fitnormchi2<<endl;
// 			EntryT->Fill(hd.ColIndT->at(indtmin)[2]);
// 			ExitT->Fill(hd.ColIndT->at(indtmax)[2]);
// 			Trackt0->Fill(td.FitParams->at(4));
// 			TrackQ->Fill(trackQ);
// 			TrackDeltaR->Fill(dr);
// 			TrackDeltaXY->Fill(sqrt(pow(hd.ColIndT->at(indtmin)[0]-hd.ColIndT->at(indtmax)[0],2.)+pow(hd.ColIndT->at(indtmin)[1]-hd.ColIndT->at(indtmax)[1],2.)));
// 			Dt1vsDt2->Fill(td.FitParams->at(4)-hd.ColIndT->at(indtmin)[2],hd.ColIndT->at(indtmax)[2]-td.FitParams->at(4));
// 		
// // 			tts.clear();
// // 			tts.push_back(FE[i1].fitParams[4]-(FE[i1].fitParams[0]/FE[i1].fitParams[1]));
// // 			tts.push_back(FE[i1].fitParams[4]-(FE[i1].fitParams[2]/FE[i1].fitParams[3]));
// // 			tts.push_back(FE[i1].fitParams[4]+((127-FE[i1].fitParams[0])/FE[i1].fitParams[1]));
// // 			tts.push_back(FE[i1].fitParams[4]+((127-FE[i1].fitParams[2])/FE[i1].fitParams[3]));
// // 			sort(tts.begin(), tts.end());
// // 			
// // 			DeltaTvsTrBin->Fill(trbin/20,FE[i1].H[indtmax].t-FE[i1].H[indtmin].t);
// // 			trbin++;
// // 			DeltaTvsTr->SetPoint(DeltaTvsTr->GetN(),DeltaTvsTr->GetN(),FE[i1].H[indtmax].t-FE[i1].H[indtmin].t);
// 		
// // 		ttmin=((FE[i1].fitParams[4]-(FE[i1].fitParams[0]/FE[i1].fitParams[1]))>(FE[i1].fitParams[4]-(FE[i1].fitParams[2]/FE[i1].fitParams[3]))?(FE[i1].fitParams[4]-(FE[i1].fitParams[0]/FE[i1].fitParams[1])):(FE[i1].fitParams[4]-(FE[i1].fitParams[2]/FE[i1].fitParams[3])));
// // 		ttmax=((FE[i1].fitParams[4]+((127-FE[i1].fitParams[0])/FE[i1].fitParams[1]))>(FE[i1].fitParams[4]+((127-FE[i1].fitParams[2])/FE[i1].fitParams[3]))?(FE[i1].fitParams[4]+((127-FE[i1].fitParams[2])/FE[i1].fitParams[3])):(FE[i1].fitParams[4]+((127-FE[i1].fitParams[0])/FE[i1].fitParams[1])));
// // 		TrackDeltaT->Fill(ttmax-ttmin);
// // 			TrackDeltaT->Fill(tts[2]-tts[1]);
// 			
// 			DeltaTvsDeltaR->Fill(sqrt(pow(hd.ColIndT->at(indtmin)[0]-hd.ColIndT->at(indtmax)[0],2.)+pow(hd.ColIndT->at(indtmin)[1]-hd.ColIndT->at(indtmax)[1],2.)),hd.ColIndT->at(indtmax)[2]-hd.ColIndT->at(indtmin)[2]);
// 		
// // 		cout<<FE[i1].E<<" "<<FE[i1].PMTint[0]<<" "<<FE[i1].PMTint[1]<<" "<<FE[i1].PMTint[2]<<endl;
// 			for(int i2=0;i2<3;i2++)
// 			{
// 				PMTIntAll[i2]->Fill(td.PMTIntegral->at(i2));
// 				if(td.PMTIntegral->at(i2)>=0)
// 				{
// 					PMTInt[i2]->Fill(td.PMTIntegral->at(i2));
// 					PMTIntPE[i2]->Fill(td.PMTIntegral->at(i2)/PMTPEs[i2]);
// 					EntryXYPMT[i2]->Fill(hd.ColIndT->at(indtmin)[0],hd.ColIndT->at(indtmin)[1]);
// 					MeanPMTChargeEntryXY[i2]->Fill(hd.ColIndT->at(indtmin)[0],hd.ColIndT->at(indtmin)[1],td.PMTIntegral->at(i2));
// 					PMTProf[i2]->Fill(sqrt(pow(hd.ColIndT->at(indtmin)[0]-PMTXY[i2][0],2.)+pow(hd.ColIndT->at(indtmin)[1]-PMTXY[i2][1],2.)),td.PMTIntegral->at(i2));
// 					PMTProfPE[i2]->Fill(sqrt(pow(hd.ColIndT->at(indtmin)[0]-PMTXY[i2][0],2.)+pow(hd.ColIndT->at(indtmin)[1]-PMTXY[i2][1],2.)),td.PMTIntegral->at(i2)/PMTPEs[i2]);
// 				}
// 				else
// 				{
// 					EntryXYPMTSat[i2]->Fill(hd.ColIndT->at(indtmin)[0],hd.ColIndT->at(indtmin)[1]);
// 				}
// 			}
// 		}
// // 		for(int i2=0;i2<hd.ColIndT->size();i2++)
// 		for(int i2=0;i2<trhitinds.size();i2++)
// 		{
// 			
// 			if(hd.ColIndT->at(trhitinds[i2])[2]<800 || hd.ColIndT->at(trhitinds[i2])[2]>=2400) continue;
// 			TrackHits[(hd.ColIndT->at(trhitinds[i2])[2]-800)/200][0]->Fill(hd.ColIndT->at(trhitinds[i2])[0],hd.ColIndT->at(trhitinds[i2])[1],hd.Int->at(trhitinds[i2])[2]);
// 			TrackHits[(hd.ColIndT->at(trhitinds[i2])[2]-800)/200][1]->Fill(hd.ColIndT->at(trhitinds[i2])[0],hd.ColIndT->at(trhitinds[i2])[1]);
// 		}
	
		inroot->Close();

	}
	outroot->cd();
	
// 	func->SetParLimits(0,0,1000.);
// 	func->SetParLimits(1,0,1000.);
// 	func->SetParLimits(2,0,100.);
// 	func->SetParLimits(3,-500,500);
// 	func->SetParLimits(4,-100,100);
// 	
// 	func->SetParameter(0,DeltaT->GetBinContent(DeltaT->GetMaximumBin()));
// 	func->SetParameter(1,DeltaT->GetBinCenter(DeltaT->GetMaximumBin()));
// 	func->SetParameter(2,10);
// 	func->SetParameter(3,0.1);
// 	func->SetParameter(4,0.1);
// 	DeltaT->Fit(func,"q","q",0.,1000.);
// 	
// 	QvstPr->Fit(exp1,"q","q",1000.,2400.);
// 	
// // 	gplusl->SetParameter(0,AllHitQ->GetBinContent(AllHitQ->GetMaximumBin()));
// // 	gplusl->SetParameter(1,AllHitQ->GetBinCenter(AllHitQ->GetMaximumBin()));
// // 	gplusl->SetParameter(2,10);
// // 	gplusl->SetParameter(3,AllHitQ->GetBinContent(AllHitQ->GetMaximumBin()));
// // 	gplusl->SetParameter(4,AllHitQ->GetBinCenter(AllHitQ->GetMaximumBin()));
// // 	gplusl->SetParameter(5,10);
// // 	AllHitQ->Fit(gplusl,"q","q",0.,1000.);
// 	AllHitQ->Fit(land1,"q","q",0.,1000.);
// 	
// 	TrackQ->Fit(land1,"q","q",0.,100000.);
	
	EntryXY->Write();
	EntryT->Write();
	ExitXY->Write();
	ExitT->Write();
	DeltaT->Write();
	Trackt0->Write();
	TrackDeltaT->Write();
	AllHitX->Write();
	AllHitY->Write();
	AllHitT->Write();
	AllHitQ->Write();
	TrackQ->Write();
	TrackDeltaR->Write();
	TrackDeltaXY->Write();
	Dt1vsDt2->Write();
	Qvst->Write();
	QvstPr->Write();
	DeltaTvsDeltaR->Write();
	DeltaTvsTrBin->Write();
	DeltaTvsTr->Write();
	                                                                
	TrackHitCharge->Write();
	OutOfTrackHitCharge->Write();
	NHitsinvsoutT->Write();
	
	for(int i1=0;i1<3;i1++)
	{
		PMTInt[i1]->Scale(1./PMTInt[i1]->Integral());
		PMTInt[i1]->Write();
		PMTIntAll[i1]->Write();
		PMTIntPE[i1]->Write();
		EntryXYPMT[i1]->Write();
		EntryXYPMTSat[i1]->Write();
		MeanPMTChargeEntryXY[i1]->Divide(EntryXYPMT[i1]);
		MeanPMTChargeEntryXY[i1]->Write();
		PMTProf[i1]->Write();
		PMTProfPE[i1]->Write();
	}
	
	for(int i1=0;i1<9;i1++)
	{
		TrackHits[i1][0]->Divide(TrackHits[i1][1]);
		TrackHits[i1][0]->Write();
		TrackHits[i1][1]->Write();
	}
	
	TGraphErrors* tge=new TGraphErrors();
	tge->SetName("MeanHitQvsTime");tge->SetTitle("MeanHitQvsTime");
	tge->SetMarkerStyle(20);tge->SetMarkerColor(1);
	
	for(int i1=0;i1<32;i1++)
	{
// 		gplusl->SetParLimits(0,0,10000);
// 		gplusl->SetParLimits(1,0,100);
// 		gplusl->SetParLimits(2,0,100);
// 		gplusl->SetParLimits(3,0,1000);
// 		gplusl->SetParLimits(4,0,100);
// 		gplusl->SetParLimits(5,0,100);
// 		
// 		gplusl->SetParameter(0,QperHit[i1]->GetBinContent(QperHit[i1]->GetMaximumBin()));
// 		gplusl->SetParameter(1,QperHit[i1]->GetBinCenter(QperHit[i1]->GetMaximumBin()));
// 		gplusl->SetParameter(2,100);
// 		gplusl->SetParameter(3,QperHit[i1]->GetBinContent(QperHit[i1]->GetMaximumBin()));
// 		gplusl->SetParameter(4,QperHit[i1]->GetBinCenter(QperHit[i1]->GetMaximumBin()));
// 		gplusl->SetParameter(5,10);
		if(QperHit[i1]->GetEntries()>0)
		{
			if(QperHit[i1]->GetEntries()<2000) QperHit[i1]->Rebin(4);
			QperHit[i1]->Fit(gaus1,"q","q",QperHit[i1]->GetMean()-QperHit[i1]->GetRMS(),QperHit[i1]->GetMean()+QperHit[i1]->GetRMS());
			QperHit[i1]->Fit(gaus1,"q","q",gaus1->GetParameter(1)-gaus1->GetParameter(2),gaus1->GetParameter(1)+gaus1->GetParameter(2));
			tge->SetPoint(tge->GetN(),i1,gaus1->GetParameter(1));
			tge->SetPointError(tge->GetN()-1,0,gaus1->GetParError(1));
		}
		
// 		QperHit[i1]->Fit(gplusl,"q","q",0.,1000.);
// 		QperHit[i1]->Fit(land1,"q","q",0.,1000.);
		QperHit[i1]->Write();
	}
	
	tge->Fit(exp1,"q","q",1.,8.);
	tge->Write();
	
	outroot->Close();
// 	outfile.close();
	
// 	inroot1->Close();
// 	inroot2->Close();
	
// 	sprintf(hname,"rm TrackParameters_%d.txt",RunNo);system(hname);
// 	sprintf(hname,"rm PMTIntegrals_%d.txt",RunNo);system(hname);
// 	sprintf(hname,"rm Hits_%d.txt",RunNo);system(hname);
	sprintf(hname,"cp TrackAnalysis3_%d.root %s/Histos/TrackAnalysis3_%d.root;wait;",RunNo,AnalysisFilePath,RunNo);system(hname);
// 	sprintf(hname,"mv ValidTrackParameters_%d.txt %s/Files/ValidTrackParameters_%d.txt;wait;",RunNo,AnalysisFilePath,RunNo);system(hname);
}

void PMTAnalysis()
{
	int PRunNo=RunNo;
	vector <int> RunNos;
	if(PRunNo==1000)
	{
		//up to 7095 different voltages
// 		RunNos.push_back(7074);
// 		RunNos.push_back(7076);
// 		RunNos.push_back(7078);
// 		RunNos.push_back(7082);
// 		RunNos.push_back(7084);
// 		RunNos.push_back(7092);
// 		RunNos.push_back(7094);
		
// 		RunNos.push_back(7096);
		RunNos.push_back(7104);
		RunNos.push_back(7106);
		RunNos.push_back(7107);
		RunNos.push_back(7110);
		RunNos.push_back(7121);
	}
	else if(PRunNo==2000)
	{
		RunNos.push_back(7252);
		RunNos.push_back(7254);
		RunNos.push_back(7256);
		RunNos.push_back(7259);
		RunNos.push_back(7268);
		RunNos.push_back(7271);
		RunNos.push_back(7273);
	}
	else if(PRunNo==3000)
	{
		RunNos.push_back(7308);
		RunNos.push_back(7312);
		RunNos.push_back(7314);
		RunNos.push_back(7316);
		RunNos.push_back(7318);
		RunNos.push_back(7320);
		RunNos.push_back(7322);
		RunNos.push_back(7325);
		RunNos.push_back(7327);
		RunNos.push_back(7385);
		RunNos.push_back(7388);
	}
	
	TH1F* hPMT[3];
// 	string PMTNames[3]={"LAr w/TPB","LAr no TPB","LXe"};
	int PMTlc[3]={1,2,4};
	float TOTThresholds[20]={10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};
// 	float ThresholdOffset[3]={70,0,10};
	float ThresholdOffset[3]={0,0,0};
	float threshold=0;
	
	TF1* tf1=new TF1("gaus","gaus",0.,300.);
	TF1* tf2=new TF1("lin","[0]*x-[1]",0.,15000.);
	
	sprintf(hname,"PMTAnalysis_%d.root",RunNo);
	TFile* outroot=new TFile(hname,"recreate");
	
	TH1D* Amplitude[3];TH1D* Integral[3];TH1D* MaxAmplitude[3];
	TH1D* IntegralforThreshold[3][20];TH1D* DeltaTforThreshold[3][20];TH1D* DeltaTforThresholdSat[3][20];
	TProfile* TOT[3][20];TProfile* TOT2[3][20];TProfile* TOT3[3];
	
// 	TH1F* WFt0=new TH1F("WFt0","WFt0",15000,-0.5,14999.5);
	TH1F* AvgWFSat[3][2][15];
	TH1F* Ht0[3];TH1F* Ht02[3][15];
	for(int i1=0;i1<3;i1++)
	{
		for(int i2=0;i2<2;i2++)
		{
			for(int i3=0;i3<RunNos.size();i3++)
			{
				sprintf(hname,"AvgWFSat_%d_%d_R%d",i1,i2,RunNos[i3]);
				AvgWFSat[i1][i2][i3]=new TH1F(hname,hname,15000,-0.5,14999.5);
			}
		}
		sprintf(hname,"Amplitude_%d",i1);
		Amplitude[i1]=new TH1D(hname,hname,1000,-0.5,999.5);
		sprintf(hname,"Ht0_%d",i1);
		Ht0[i1]=new TH1F(hname,hname,15000,0,15000);
		for(int i2=0;i2<RunNos.size();i2++)
		{
			sprintf(hname,"Ht02_%d_R%d",i1,RunNos[i2]);
			Ht02[i1][i2]=new TH1F(hname,hname,15000,0,15000);
		}
		
		sprintf(hname,"MaxAmplitude_%d",i1);
		MaxAmplitude[i1]=new TH1D(hname,hname,1000,-0.5,999.5);
		sprintf(hname,"Integral_%d",i1);
		Integral[i1]=new TH1D(hname,hname,1000,0,1e5);
		for(int i2=0;i2<20;i2++)
		{
			threshold=TOTThresholds[i2]+ThresholdOffset[i1];
			sprintf(hname,"IntegralforThreshold_%d_%5.2f",i1,threshold);
			IntegralforThreshold[i1][i2]=new TH1D(hname,hname,1000,0,50000);
			sprintf(hname,"TOT_%d_%5.2f",i1,threshold);
// 			TOT[i1][i2]=new TProfile(hname,hname,300,0.,300.,0.,30000);
			TOT[i1][i2]=new TProfile(hname,hname,400,0.,100.,0.,20000);
			sprintf(hname,"TOT2_%d_%5.2f",i1,threshold);
			TOT2[i1][i2]=new TProfile(hname,hname,400,0.,100.,0.,20000);
			sprintf(hname,"DeltaTforThreshold_%d_%5.2f",i1,threshold);
			DeltaTforThreshold[i1][i2]=new TH1D(hname,hname,300,0,300);
			sprintf(hname,"DeltaTforThresholdSat_%d_%5.2f",i1,threshold);
			DeltaTforThresholdSat[i1][i2]=new TH1D(hname,hname,300,0,300);
		}
		sprintf(hname,"TOT_%d_400",i1,threshold);
		TOT3[i1]=new TProfile(hname,hname,400,0.,100.,0.,20000);
	}
	
// 	vector <TH1F*> HDeltaT[3];
// 	TH1F* HDeltaT[3][40];
// 	vector <float> bmin[3];
// 	vector <float> bmax[3];
// 	TH1F* hdt=new TH1F("hd","hd",300,0,300);
// 	for(int i1=0;i1<40;i1++)
// 	{
// 		for(int i2=0;i2<3;i2++)
// 		{
// 			sprintf(hname,"Th_%d_%d",5+(i1*5),i2);
// 			HDeltaT[i2][i1]=new TH1F(hname,hname,300,0,300);
// // 			hdt->SetName(hname);hdt->SetTitle(hname);
// // 			HDeltaT[i2].push_back(hdt);
// 			bmin[i2].push_back(0);
// 			bmax[i2].push_back(0);
// 		}
// 	}
	
	float maxamp=0;float amp10p=0;float amp90p=0;float amp25p=0;float amp75p=0;
// 	float PMTXY[3][2]={{96.5,63.5},{63.5,30.5},{63.5,96.5}};
	float PMTXY[3][2]={{95.,63.5},{63.5,32.},{63.5,95.}};
	int I=0;int J=0;int II=0;
	float qint=0;float sg=0;
	float qintTOT[20]={0};float qintTOT400=0;int binmin400=0;int binmax400=0;
	int binmin[20]={0};int binmax[20]={0};
	
	TFile* inroot;
	TFile* inroot2;
	TTree* TH;
	TTree* TT;
	TTree* T;
	bool PMTsat=false;
	bool foundt0=false;int tt0=0;
	for(int J=0;J<RunNos.size();J++)
	{
		RunNo=RunNos[J];
		PMTBaselineRun=FindBaselineRun(RunNo);
		cout<<"Run: "<<RunNo<<" PMTBaselineRun: "<<PMTBaselineRun<<endl;
		ReadPMTBaselines();
		sprintf(hname,"%s/Histos/HitsAndTracks_%d.root",AnalysisFilePath,RunNo);
		inroot=new TFile(hname);
		TH =  (TTree*) inroot->Get("Hits");
		TH->SetBranchAddress("E",&hd.E);
		TH->SetBranchAddress("ColIndT",&hd.ColIndT);
		TH->SetBranchAddress("Int",&hd.Int);
		TH->SetBranchAddress("QColTot",&hd.QColTot);
		TH->SetBranchAddress("QColTotZS",&hd.QColTotZS);
		TH->SetBranchAddress("QHitTot",&hd.QHitTot);
		TH->SetBranchAddress("PMTIntegral",&hd.PMTIntegral);
		TH->SetBranchAddress("ColID",&hd.ColID);
		TH->SetBranchAddress("ColT",&hd.ColT);
		TH->SetBranchAddress("ColA",&hd.ColA);
		TH->SetBranchAddress("ColInt",&hd.ColInt);
		TH->SetBranchAddress("Colw",&hd.Colw);
		TH->SetBranchAddress("IndID",&hd.IndID);
		TH->SetBranchAddress("IndT",&hd.IndT);
		TH->SetBranchAddress("IndA",&hd.IndA);
		TH->SetBranchAddress("Indw",&hd.Indw);
		TH->SetBranchAddress("EventType",&hd.EventType);
		TT =  (TTree*) inroot->Get("Tracks");
		TT->SetBranchAddress("E",&td.E);
		TT->SetBranchAddress("StartEndColIndT",&td.StartEndColIndT);
		TT->SetBranchAddress("FitParams",&td.FitParams);
		TT->SetBranchAddress("FitNormChi2",&td.FitNormChi2);
		TT->SetBranchAddress("NHits",&td.NHits);
		TT->SetBranchAddress("Nexcl",&td.Nexcl);
		TT->SetBranchAddress("PMTIntegral",&td.PMTIntegral);
		TT->SetBranchAddress("ColTStartEnd",&td.ColTStartEnd);
		TT->SetBranchAddress("ColHitTStartEnd",&td.ColHitTStartEnd);
		
		sprintf(hname,"%s/Run_%d.root",NTupleFilePath,RunNo);
// 		cout<<hname<<endl;
		inroot2=new TFile(hname);
		T =  (TTree*) inroot2->Get("T");
		T->SetBranchAddress("TPCWF",&ed.TPCWF);
		T->SetBranchAddress("PMTWF",&ed.PMTWF);
		
		for(int I=0;I<TH->GetEntries();I++)
		{
			TH->GetEntry(I);
			if(!(hd.EventType==0 || hd.EventType==1)) continue;
			T->GetEntry(hd.E);
			for(int i1=0;i1<3;i1++)
			{
				if(i1==1) continue;
				
				if(hd.PMTIntegral->at(i1)<0)//saturated
				{
					for(int i2=0;i2<ed.PMTWF->at(i1).size();i2++)
					{
						sg=(-1.)*(((float)ed.PMTWF->at(i1)[i2])-PMTBaselines[i1][0]);
						AvgWFSat[i1][0][J]->Fill(i2,sg);
						AvgWFSat[i1][1][J]->Fill(i2);
					}
				}
				
				qint=0.;
				maxamp=0;
				foundt0=false;
				if(hd.PMTIntegral->at(i1)>0)//unsaturated
				{
					for(int i2=1900;i2<2100;i2++)
					{
						sg=(-1.)*(((float)ed.PMTWF->at(i1)[i2])-PMTBaselines[i1][0]);
						if(sg>400)
						{
							if(binmin400==0) binmin400=i2;
							binmax400=i2;
							qintTOT400+=sg;
						}
						else
						{
							if(binmin400!=0) break;
						}
					}
				}
// 				for(int i2=0;i2<20;i2++){qintTOT[i2]=0;binmin[i2]=0;binmax[i2]=0;}
				for(int i2=1900;i2<2100;i2++)
				{
					sg=(-1.)*(((float)ed.PMTWF->at(i1)[i2])-PMTBaselines[i1][0]);
					if(hd.PMTIntegral->at(i1)>=0) Amplitude[i1]->Fill(sg);//unsaturated
					if(sg>maxamp) maxamp=sg;
					if(sg>10 && !foundt0){tt0=i2;foundt0=true;}
					qint+=sg;
					
// 					for(int i3=0;i3<20;i3++)
// 					{
// 						threshold=TOTThresholds[i3]+ThresholdOffset[i1];
// 						if(sg>threshold)
// 						{
// 							qintTOT[i3]+=sg;
// 							if(binmin[i3]==0) binmin[i3]=i2;
// 							binmax[i3]=i2;
// 						}
// 						else
// 						{
// 							if(binmin[i3]!=0) continue;
// 						}
// 					}
// 					for(int i3=0;i3<40;i3++)
// 					{
// 						if(sg>((i3*5)+5))
// 						{
// 							if(bmin[i1][i3]==0) bmin[i1][i3]=i2;
// 							bmax[i1][i3]=i2;
// 						}
// 						else
// 						{
// 							if(bmin[i1][i3]!=0) continue;
// 						}
// 					}
				}
				if(hd.PMTIntegral->at(i1)>=0)
				{
					Integral[i1]->Fill(qint);//unsaturated
					MaxAmplitude[i1]->Fill(maxamp);//unsaturated
// 					for(int i3=0;i3<40;i3++)
// 					{
// 						if((bmax[i1][i3]-bmin[i1][i3])>0)
// 						{
// 							HDeltaT[i1][i3]->Fill(bmax[i1][i3]-bmin[i1][i3]);
// 						}
// 						bmax[i1][i3]=0;bmin[i1][i3]=0;
// 					}
					
// 					amp25p=0.25*maxamp;
// 					amp75p=0.75*maxamp;
// 					WFt0->Reset();
// 					for(int i2=1900;i2<2100;i2++)
// 					{
// 						sg=(-1.)*(((float)ed.PMTWF->at(i1)[i2])-PMTBaselines[i1][0]);
// 						WFt0->SetBinContent(i2,sg);
// 					}
// 					WFt0->Fit(tf2,"q","q",WFt0->FindBin(amp25p),WFt0->FindBin(amp75p));
// 					Ht0[i1]->Fill(tf2->GetParameter(1)/tf2->GetParameter(0));
					Ht0[i1]->Fill(tt0);
					Ht02[i1][J]->Fill(tt0);
				}
				if(maxamp>600 && hd.PMTIntegral->at(i1)>=0)
				{
					TOT3[i1]->Fill(binmax400-binmin400,qintTOT400);
				}
// 				for(int i2=0;i2<20;i2++)
// 				{
// 					threshold=TOTThresholds[i2]+ThresholdOffset[i1];
// 					if(hd.PMTIntegral->at(i1)>=0) IntegralforThreshold[i1][i2]->Fill(qintTOT[i2]);
// 					if((binmax[i2]-binmin[i2])>2)
// 					{
// 						if(hd.PMTIntegral->at(i1)>=0)//unsaturated
// 						{
// 							DeltaTforThreshold[i1][i2]->Fill(binmax[i2]-binmin[i2]);
// 							TOT[i1][i2]->Fill(binmax[i2]-binmin[i2],qintTOT[i2]);
// 						}
// 						else//saturated
// 						{
// 							DeltaTforThresholdSat[i1][i2]->Fill(binmax[i2]-binmin[i2]);
// 						}
// 					}
// 				}
			}
		}
		inroot->Close();
		inroot2->Close();
	}
// 	TGraph* tg1[3];
// 	tg1[0]=new TGraph();tg1[2]=new TGraph();
// 	tg1[0]->SetName("chi2vstot_0");tg1[0]->SetTitle("chi2vstot_0");
// 	tg1[0]->SetMarkerStyle(20);tg1[0]->SetMarkerColor(1);
// 	tg1[0]->GetXaxis()->SetTitle("Time Over Threshold");tg1[0]->GetXaxis()->CenterTitle();
// 	tg1[0]->GetYaxis()->SetTitle("#chi^{2}/ndf ");tg1[0]->GetYaxis()->CenterTitle();
// 	tg1[2]->SetName("chi2vstot_2");tg1[2]->SetTitle("chi2vstot_2");
// 	tg1[2]->SetMarkerStyle(20);tg1[2]->SetMarkerColor(1);
// 	tg1[2]->GetXaxis()->SetTitle("Time Over Threshold");tg1[2]->GetXaxis()->CenterTitle();
// 	tg1[2]->GetYaxis()->SetTitle("#chi^{2}/ndf ");tg1[2]->GetYaxis()->CenterTitle();
// 	outroot->cd();
// 	for(int i1=0;i1<3;i1++)
// 	{
// 		if(i1==1) continue;
// 		for(int i3=0;i3<40;i3++)
// 		{
// 			HDeltaT[i1][i3]->Fit(tf1,"q","q",0.,300.);
// 			HDeltaT[i1][i3]->Fit(tf1,"q","q",tf1->GetParameter(1)-2*tf1->GetParameter(2),tf1->GetParameter(1)+2*tf1->GetParameter(2));
// 			tg1[i1]->SetPoint(tg1[i1]->GetN(),(i3*5)+5,tf1->GetChisquare()/tf1->GetNDF());
// 			HDeltaT[i1][i3]->Write();
// 		}
// 	}
	
// 	float TOTRange[3][20][2]={{{0}}};
// 	for(int i1=0;i1<3;i1++)
// 	{
// 		if(i1==1) continue;
// 		for(int i2=0;i2<20;i2++)
// 		{
// 			DeltaTforThreshold[i1][i2]->Fit(tf1,"q","q",0.,50.);
// 			DeltaTforThreshold[i1][i2]->Fit(tf1,"q","q",tf1->GetParameter(1)-2*tf1->GetParameter(2),tf1->GetParameter(1)+2*tf1->GetParameter(2));
// 			TOTRange[i1][i2][0]=tf1->GetParameter(1)-2.5*tf1->GetParameter(2);
// 			TOTRange[i1][i2][1]=tf1->GetParameter(1)+2.5*tf1->GetParameter(2);
// 		}
// 	}
// 	for(int J=0;J<RunNos.size();J++)
// 	{
// 		RunNo=RunNos[J];
// 		PMTBaselineRun=FindBaselineRun(RunNo);
// 		cout<<"Run: "<<RunNo<<" PMTBaselineRun: "<<PMTBaselineRun<<endl;
// 		ReadPMTBaselines();
// 		sprintf(hname,"%s/Histos/HitsAndTracks_%d.root",AnalysisFilePath,RunNo);
// 		inroot=new TFile(hname);
// 		TH =  (TTree*) inroot->Get("Hits");
// 		TH->SetBranchAddress("E",&hd.E);
// 		TH->SetBranchAddress("ColIndT",&hd.ColIndT);
// 		TH->SetBranchAddress("Int",&hd.Int);
// 		TH->SetBranchAddress("QColTot",&hd.QColTot);
// 		TH->SetBranchAddress("QColTotZS",&hd.QColTotZS);
// 		TH->SetBranchAddress("QHitTot",&hd.QHitTot);
// 		TH->SetBranchAddress("PMTIntegral",&hd.PMTIntegral);
// 		TH->SetBranchAddress("ColID",&hd.ColID);
// 		TH->SetBranchAddress("ColT",&hd.ColT);
// 		TH->SetBranchAddress("ColA",&hd.ColA);
// 		TH->SetBranchAddress("ColInt",&hd.ColInt);
// 		TH->SetBranchAddress("Colw",&hd.Colw);
// 		TH->SetBranchAddress("IndID",&hd.IndID);
// 		TH->SetBranchAddress("IndT",&hd.IndT);
// 		TH->SetBranchAddress("IndA",&hd.IndA);
// 		TH->SetBranchAddress("Indw",&hd.Indw);
// 		TH->SetBranchAddress("EventType",&hd.EventType);
// 		TT =  (TTree*) inroot->Get("Tracks");
// 		TT->SetBranchAddress("E",&td.E);
// 		TT->SetBranchAddress("StartEndColIndT",&td.StartEndColIndT);
// 		TT->SetBranchAddress("FitParams",&td.FitParams);
// 		TT->SetBranchAddress("FitNormChi2",&td.FitNormChi2);
// 		TT->SetBranchAddress("NHits",&td.NHits);
// 		TT->SetBranchAddress("Nexcl",&td.Nexcl);
// 		TT->SetBranchAddress("PMTIntegral",&td.PMTIntegral);
// 		TT->SetBranchAddress("ColTStartEnd",&td.ColTStartEnd);
// 		TT->SetBranchAddress("ColHitTStartEnd",&td.ColHitTStartEnd);
// 		
// 		sprintf(hname,"%s/Run_%d.root",NTupleFilePath,RunNo);
// // 		cout<<hname<<endl;
// 		inroot2=new TFile(hname);
// 		T =  (TTree*) inroot2->Get("T");
// 		T->SetBranchAddress("TPCWF",&ed.TPCWF);
// 		T->SetBranchAddress("PMTWF",&ed.PMTWF);
// 		
// 		for(int I=0;I<TH->GetEntries();I++)
// 		{
// 			TH->GetEntry(I);
// 			if(!(hd.EventType==0 || hd.EventType==1)) continue;
// 			T->GetEntry(hd.E);
// 			for(int i1=0;i1<3;i1++)
// 			{
// 				if(i1==1) continue;
// 				if(hd.PMTIntegral->at(i1)>=0) //unsaturated
// 				{
// 					qint=0.;
// 					for(int i2=0;i2<20;i2++){qintTOT[i2]=0;binmin[i2]=0;binmax[i2]=0;}
// 					for(int i2=1900;i2<2100;i2++)
// 					{
// 						sg=(-1.)*(((float)ed.PMTWF->at(i1)[i2])-PMTBaselines[i1][0]);
// 						qint+=sg;
// 						for(int i3=0;i3<20;i3++)
// 						{
// 							threshold=TOTThresholds[i3]+ThresholdOffset[i1];
// 							if(sg>threshold)
// 							{
// 								qintTOT[i3]+=sg;
// 								if(binmin[i3]==0) binmin[i3]=i2;
// 								binmax[i3]=i2;
// 							}
// 							else
// 							{
// 								if(binmin[i3]!=0) continue;
// 							}
// 						}
// 					}
// // 					Integral[i1]->Fill(qint);
// 					for(int i2=0;i2<20;i2++)
// 					{
// 						threshold=TOTThresholds[i2]+ThresholdOffset[i1];
// // 						IntegralforThreshold[i1][i2]->Fill(qintTOT[i2]);
// 						if((binmax[i2]-binmin[i2])>2)
// 						{
// 							if(((float)(binmax[i2]-binmin[i2]))>TOTRange[i1][i2][0] && ((float)(binmax[i2]-binmin[i2]))<TOTRange[i1][i2][1])
// 							{
// 	// 							DeltaTforThreshold[i1][i2]->Fill(binmax[i2]-binmin[i2]);
// 								TOT2[i1][i2]->Fill(binmax[i2]-binmin[i2],qintTOT[i2]);
// 							}
// 						}
// 					}
// 				}
// 			}
// 		}
// 		inroot->Close();
// 		inroot2->Close();
// 	}
	
// 	sprintf(hname,"TOTParameters_%d.txt",PRunNo);
// 	ofstream outfile(hname);
// 	TF1* tfpol4=new TF1("tfpol4","[0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4",0.,300.);
// 	for(int i1=0;i1<3;i1++)
// 	{
// 		if(i1==1) continue;
// 		for(int i2=0;i2<20;i2++)
// 		{
// 			TOT[i1][i2]->Fit(tfpol4,"q","q",TOTRange[i1][i2][0],TOTRange[i1][i2][1]);
// 			outfile<<i1<<" "<<(TOTThresholds[i2]+ThresholdOffset[i1])<<" "<<tfpol4->GetParameter(0)<<" "<<tfpol4->GetParameter(1)<<" "<<tfpol4->GetParameter(2)<<" "<<tfpol4->GetParameter(3)<<" "<<tfpol4->GetParameter(4)<<" "<<(tfpol4->GetChisquare()/tfpol4->GetNDF())<<endl;
// 		}
// 	}
// 	outfile.close();
	
	outroot->cd();
	for(int i1=0;i1<3;i1++)
	{
		if(i1==1) continue;
		Amplitude[i1]->Write();
		Ht0[i1]->Write();
		for(int i2=0;i2<RunNos.size();i2++)
		{
			Ht02[i1][i2]->Write();
		}
		MaxAmplitude[i1]->Write();
		Integral[i1]->Write();
		for(int i2=0;i2<20;i2++)
		{
			IntegralforThreshold[i1][i2]->Write();
			DeltaTforThreshold[i1][i2]->Write();
			DeltaTforThresholdSat[i1][i2]->Write();
			TOT[i1][i2]->Write();
			TOT2[i1][i2]->Write();
		}
		TOT3[i1]->Write();
		for(int i3=0;i3<RunNos.size();i3++)
		{
			AvgWFSat[i1][0][i3]->Divide(AvgWFSat[i1][1][i3]);
			AvgWFSat[i1][0][i3]->Write();
		}
	}
// 	tg1[0]->Write();
// 	tg1[2]->Write();
	
// 	TCanvas* cc=new TCanvas("cc","cc",600,600);
// 	TLegend* legend = new TLegend(0.6,0.7,0.9,0.9);
// 	for(int i1=0;i1<3;i1++)
// 	{
// 		if(i1==1) continue;
// 		legend->Clear();
// 		for(int i2=0;i2<20;i2++)
// 		{
// 			threshold=TOTThresholds[i2]+ThresholdOffset[i1];
// 			sprintf(hname,"Th. %d",((int)threshold));
// 			TOT[i1][i2]->SetName(hname);TOT[i1][i2]->SetTitle(hname);
// 			TOT[i1][i2]->SetLineColor(i2+1);TOT[i1][i2]->SetMarkerStyle(20);TOT[i1][i2]->SetMarkerColor(i2+1);
// 			TOT[i1][i2]->GetFunction("tfpol4")->SetLineColor(i2+1);
// 			if(i2==0)
// 			{
// 				TOT[i1][i2]->GetXaxis()->SetTitle("Time over Threshold");TOT[i1][i2]->GetXaxis()->CenterTitle();
// 				TOT[i1][i2]->GetYaxis()->SetTitle("Integral in Time over Threshold");TOT[i1][i2]->GetYaxis()->CenterTitle();
// 				TOT[i1][i2]->Draw();
// 				legend->AddEntry(TOT[i1][i2],TOT[i1][i2]->GetName(),"l");
// 			}
// 			else if(i2<19)
// 			{
// 				TOT[i1][i2]->Draw("same");
// 				legend->AddEntry(TOT[i1][i2],TOT[i1][i2]->GetName(),"l");
// 			}
// 			else
// 			{
// 				TOT[i1][i2]->Draw("same");
// 				legend->AddEntry(TOT[i1][i2],TOT[i1][i2]->GetName(),"l");
// 				legend->Draw("same");
// 				sprintf(hname,"IvsTOT_PMT_%d",i1);
// 				cc->SetName(hname);
// 				cc->Write();
// 			}
// 		}
// 	}
	
	
	outroot->Close();
// 	sprintf(hname,"cp PMTAnalysis_%d.root %s/Histos/PMTAnalysis_%d.root;wait;",PRunNo,AnalysisFilePath,RunNo);system(hname);
}

void TrackAndPMTAnalysis()//for light yield
{
	int PRunNo=RunNo;
	vector <int> RunNos;
	if(RunNo==1000)
	{
		//up to 7095 different PMT voltages
		RunNos.push_back(7074);//closed first
		RunNos.push_back(7076);//closed first
		RunNos.push_back(7078);//closed first
		RunNos.push_back(7082);//closed first
		RunNos.push_back(7084);//closed first
		RunNos.push_back(7092);//closed first
		RunNos.push_back(7094);//closed first
		RunNos.push_back(7096);//closed first
		RunNos.push_back(7104);//bump in slope//closed first
		RunNos.push_back(7106);//closed first
		
		RunNos.push_back(7107);
		RunNos.push_back(7110);
		RunNos.push_back(7121);
	}
	else if(RunNo==2000)
	{
		RunNos.push_back(7252);//closed first
		RunNos.push_back(7254);//closed first
		RunNos.push_back(7256);//closed first
		RunNos.push_back(7259);//closed first
		
		RunNos.push_back(7268);
		RunNos.push_back(7271);
		RunNos.push_back(7273);
	}
	else if(RunNo==3000)
	{
		RunNos.push_back(7308);
		RunNos.push_back(7312);
		RunNos.push_back(7314);
		RunNos.push_back(7316);
		RunNos.push_back(7318);
		RunNos.push_back(7320);
		RunNos.push_back(7322);
		RunNos.push_back(7325);//closed first
		RunNos.push_back(7327);//closed first
		RunNos.push_back(7385);//closed first
		RunNos.push_back(7388);
	}
	
	float SPE[3][5][2]={{{51.74,0.04},{38.53,0.06},{29.66,0.02},{32.04,0.01},{32.38,0.01}},{{0,0},{0,0},{0,0},{0,0},{0,0}},{{36.59,0.11},{28.46,0.07},{23.51,0.01},{22.65,0.01},{21.23,0.02}}};
	float tlims[5][2]={{1100,2040},{1100,2040},{1100,2040},{1150,2100},{1440,2350}};
	float DriftSpeeds[5]={0.149,0.149,0.149,0.147,0.150};//cm/us
	int runper=0;//0:7000-7087; 1:7090-7094; 2-7095-7144; 3:7200-7287; 4:7300-7403
	
	float PMTTOTFit[3][5]={{-3667.1,1036.12,-77.2925,2.88782,-0.0330686},{0,0,0,0,0},{4952.66,-1114.38,99.3245,-3.55144,0.0481073}};
	float PMTTOTTh[3]={110,0,40};
	float PMTintegral=0;
	float PMTintegrals[3]={0};
	float PMTTOTintegral=0;
	float PMTTOT=0;
	
	TH1F* hPMT[3];
	float PMTXY[3][2]={{95.,63.5},{63.5,32.},{63.5,95.}};
// 	string PMTNames[3]={"LAr w/TPB","LAr no TPB","LXe"};
	int PMTlc[3]={1,2,4};
	float TOTThresholds[20]={10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};
	float ThresholdOffset[3]={70,0,10};
	float threshold=0;
	
	TF1* tf1=new TF1("gaus","gaus",0.,3000.);
	TF1* tf2 = new TF1("fit",fitf3f,0.,1000,7);
	
	sprintf(hname,"TrackAndPMTAnalysis_%d.root",RunNo);
	TFile* outroot=new TFile(hname,"recreate");
	
	TH2F* NColHitsvsChi2Vertical=new TH2F("NColHitsvsChi2Vertical","NColHitsvsChi2Vertical",200,0,10,200,0,200);
	TH1F* AllHitT=new TH1F("AllHitT","AllHitT",500,0,5000);
	TH1F* VerticalChi2=new TH1F("VerticalChi2","VerticalChi2",100,0,10);
	TH1F* VerticalChi2Sel=new TH1F("VerticalChi2Sel","VerticalChi2Sel",100,0,10);
	TH1F* NColHits=new TH1F("NColHits","NColHits",1500,-0.5,1499.5);
	TH1F* NColHitsPerDeltaT=new TH1F("NColHitsPerDeltaT","NColHitsPerDeltaT",500,0,0.5);
	TH1F* NColHitsVertical=new TH1F("NColHitsVertical","NColHitsVertical",1500,-0.5,1499.5);
	
	TH1F* DistancetoPMT[2][2];
	DistancetoPMT[0][0]=new TH1F("DistancetoPMT_0_vertical","DistancetoPMT_0_vertical",150,0,60);
	DistancetoPMT[1][0]=new TH1F("DistancetoPMT_1_vertical","DistancetoPMT_1_vertical",150,0,60);
	DistancetoPMT[0][1]=new TH1F("DistancetoPMT_0_nonvertical","DistancetoPMT_0_nonvertical",150,0,60);
	DistancetoPMT[1][1]=new TH1F("DistancetoPMT_1_nonvertical","DistancetoPMT_1_nonvertical",150,0,60);
	
	TH1F* PMTIntegral[3][3][3];
	string s1[3]={"vertical","nonvertical","all"};
	string s2[3]={"nonsaturated","saturated","all"};
	for(int i1=0;i1<3;i1++)
	{
		if(i1==1) continue;
		for(int i2=0;i2<3;i2++)
		{
			for(int i3=0;i3<3;i3++)
			{
				sprintf(hname,"PMTIntegral_%d_%s_%s",i1,s1[i2].c_str(),s2[i3].c_str());
				PMTIntegral[i1][i2][i3]=new TH1F(hname,hname,300,0,60);
			}
		}
	}
	TProfile* PMTIntvsDist[3][3];
	for(int i1=0;i1<3;i1++)
	{
		if(i1==1) continue;
		for(int i2=0;i2<3;i2++)
		{
			sprintf(hname,"PMTIntvsDist_%d_%s",i1,s2[i2].c_str());
			PMTIntvsDist[i1][i2]=new TProfile(hname,hname,300,0,60,0,60);
		}
	}
	
	TH1F* TotalCollectionCharge=new TH1F("TotalCollectionCharge","TotalCollectionCharge",1200,-200000,1000000);
	TH1F* TotalCollectionChargeZS=new TH1F("TotalCollectionChargeZS","TotalCollectionChargeZS",1000,0,1000000);
	TProfile* DeltaTvsEvent=new TProfile("DeltaTvsEvent","DeltaTvsEvent",500,0,500,0,1000);
	TH1F* TrackColHitInt=new TH1F("TrackColHitInt","TrackColHitInt",2000,0,20000);
	TH1F* TrackColHitAmp=new TH1F("TrackColHitAmp","TrackColHitAmp",600,0,1500);
	TH1F* VerticalTrackDeltaR=new TH1F("VerticalTrackDeltaR","VerticalTrackDeltaR",100,0,200);
	TH1F* VerticalTrackLength=new TH1F("VerticalTrackLength","VerticalTrackLength",300,0,60);
	TH1F* NonVerticalTrackLength=new TH1F("NonVerticalTrackLength","NonVerticalTrackLength",300,0,60);
	TH1F* AllTrackLength=new TH1F("AllTrackLength","AllTrackLength",300,0,60);
	TH2F* VerticalTrackEntryXY=new TH2F("VerticalTrackEntryXY","VerticalTrackEntryXY",128,-0.5,127.5,128,-0.5,127.5);
	TH2F* VerticalTrackExitXY=new TH2F("VerticalTrackExitXY","VerticalTrackExitXY",128,-0.5,127.5,128,-0.5,127.5);
	TH1F* VerticalTrackColHitsDeltaT=new TH1F("VerticalTrackColHitsDeltaT","VerticalTrackColHitsDeltaT",1024,-0.5,4095.5);
	TH1F* VerticalTrackColHitsDeltaTSel=new TH1F("VerticalTrackColHitsDeltaTSel","VerticalTrackColHitsDeltaTSel",1024,-0.5,4095.5);
	TH1F* VerticalTrackColHitsDeltaTSel2=new TH1F("VerticalTrackColHitsDeltaTSel2","VerticalTrackColHitsDeltaTSel2",1024,-0.5,4095.5);
	TH1F* VerticalTrackTimeofFirstColHit=new TH1F("VerticalTrackTimeofFirstColHit","VerticalTrackTimeofFirstColHit",2000,0,4000);
	TH1F* VerticalTrackTimeofFirstColHits[15];
	TH1F* VerticalTrackColHitsDeltaTs[15];
	TProfile* DeltaTvsDeltaR=new TProfile("DeltaTvsDeltaR","DeltaTvsDeltaR",200,0,200,0,1000);
	for(int i1=0;i1<RunNos.size();i1++)
	{
		sprintf(hname,"VerticalTrackTimeofFirstColHits_R%d",RunNos[i1]);
		VerticalTrackTimeofFirstColHits[i1]=new TH1F(hname,hname,2000,0,4000);
		sprintf(hname,"VerticalTrackColHitsDeltaT_R%d",RunNos[i1]);
		VerticalTrackColHitsDeltaTs[i1]=new TH1F(hname,hname,2000,0,4000);
	}
// 	TH1F* InclinedTrackSlope;
	
	float x1(0),y1(0),t1(0),x2(0),y2(0),t2(0);
	float x1f(0),y1f(0),t1f(0),x2f(0),y2f(0),t2f(0);
	
	float tracklength=0;float pmtdist[3]={0};
	
	TFile* inroot;
	TFile* inroot2;
	TTree* TH;
	TTree* TT;
	TTree* T;
	bool PMTsat=false;
	float tmin=0;float tmax=0;
	int indtmin=0;int indtmax=0;
	float x=0;float y=0;
	float deltaR=0;
	float deltaT=0;
	bool vertical=true;
	int NEvtvertical=0;
	float sg=0;
	int TOTbinmin=0;int TOTbinmax=0;
	for(int J=0;J<RunNos.size();J++)
	{
		RunNo=RunNos[J];
		
		if(RunNo>=7000 && RunNo<=7087) runper=0;
		else if(RunNo>=7090 && RunNo<=7094) runper=1;
		else if(RunNo>=7095 && RunNo<=7144) runper=2;
		else if(RunNo>=7200 && RunNo<=7287) runper=3;
		else if(RunNo>=7300 && RunNo<=7403) runper=4;
		
		PMTBaselineRun=FindBaselineRun(RunNo);
		cout<<"Run: "<<RunNo<<" PMTBaselineRun: "<<PMTBaselineRun<<endl;
		ReadPMTBaselines();
		sprintf(hname,"%s/Histos/HitsAndTracks_%d.root",AnalysisFilePath,RunNo);
		inroot=new TFile(hname);
		TH =  (TTree*) inroot->Get("Hits");
		TH->SetBranchAddress("E",&hd.E);
		TH->SetBranchAddress("ColIndT",&hd.ColIndT);
		TH->SetBranchAddress("Int",&hd.Int);
		TH->SetBranchAddress("QColTot",&hd.QColTot);
		TH->SetBranchAddress("QColTotZS",&hd.QColTotZS);
		TH->SetBranchAddress("QHitTot",&hd.QHitTot);
		TH->SetBranchAddress("PMTIntegral",&hd.PMTIntegral);
		TH->SetBranchAddress("ColID",&hd.ColID);
		TH->SetBranchAddress("ColT",&hd.ColT);
		TH->SetBranchAddress("ColA",&hd.ColA);
		TH->SetBranchAddress("ColInt",&hd.ColInt);
		TH->SetBranchAddress("Colw",&hd.Colw);
		TH->SetBranchAddress("IndID",&hd.IndID);
		TH->SetBranchAddress("IndT",&hd.IndT);
		TH->SetBranchAddress("IndA",&hd.IndA);
		TH->SetBranchAddress("Indw",&hd.Indw);
		TH->SetBranchAddress("EventType",&hd.EventType);
		TT =  (TTree*) inroot->Get("Tracks");
		TT->SetBranchAddress("E",&td.E);
		TT->SetBranchAddress("StartEndColIndT",&td.StartEndColIndT);
		TT->SetBranchAddress("FitParams",&td.FitParams);
		TT->SetBranchAddress("FitNormChi2",&td.FitNormChi2);
		TT->SetBranchAddress("NHits",&td.NHits);
		TT->SetBranchAddress("Nexcl",&td.Nexcl);
		TT->SetBranchAddress("PMTIntegral",&td.PMTIntegral);
		TT->SetBranchAddress("ColTStartEnd",&td.ColTStartEnd);
		TT->SetBranchAddress("ColHitTStartEnd",&td.ColHitTStartEnd);
		
		sprintf(hname,"%s/Run_%d.root",NTupleFilePath,RunNo);
// 		cout<<hname<<endl;
		inroot2=new TFile(hname);
		T =  (TTree*) inroot2->Get("T");
		T->SetBranchAddress("TPCWF",&ed.TPCWF);
		T->SetBranchAddress("PMTWF",&ed.PMTWF);
		
		for(int I=0;I<TT->GetEntries();I++)
		{
			TT->GetEntry(I);
			
			for(int II=td.E;II<TH->GetEntries();II++)
			{
				TH->GetEntry(II);
				if(hd.E==td.E) break;
			}
			if(!(hd.EventType==0 || hd.EventType==1)) continue;
			T->GetEntry(hd.E);
			
			tmin=10000.;tmax=0.;vertical=true;
			for(int ik1=0;ik1<hd.ColID->size();ik1++)
			{
				x=td.FitParams->at(0)+td.FitParams->at(1)*(hd.ColT->at(ik1)-td.FitParams->at(4));
				if(fabs(x-((float)hd.ColID->at(ik1)))<=5)
				{
					TrackColHitInt->Fill(hd.ColInt->at(ik1));
					TrackColHitAmp->Fill(hd.ColA->at(ik1));
					if(hd.ColInt->at(ik1)>10)
					{
						if(hd.ColT->at(ik1)<tmin){tmin=hd.ColT->at(ik1);indtmin=ik1;}
						if(hd.ColT->at(ik1)>tmax){tmax=hd.ColT->at(ik1);indtmax=ik1;}
					}
				}
			}
			deltaT=tmax-tmin;
			NColHits->Fill(hd.ColID->size());
			if(tmin<tlims[runper][0] || tmax>tlims[runper][1]) continue;
			
			for(int ik1=0;ik1<hd.ColIndT->size();ik1++)
			{
				AllHitT->Fill(hd.ColIndT->at(ik1)[2]);
			}
			TotalCollectionCharge->Fill(hd.QColTot);
			TotalCollectionChargeZS->Fill(hd.QColTotZS);
			
			x1=((float)td.StartEndColIndT->at(0));
			y1=((float)td.StartEndColIndT->at(1));
			t1=((float)td.StartEndColIndT->at(2));
			x2=((float)td.StartEndColIndT->at(3));
			y2=((float)td.StartEndColIndT->at(4));
			t2=((float)td.StartEndColIndT->at(5));
			
			x1f=td.FitParams->at(0)+td.FitParams->at(1)*(tmin-td.FitParams->at(4));
			y1f=td.FitParams->at(2)+td.FitParams->at(3)*(tmin-td.FitParams->at(4));
			t1f=tmin;
			x2f=td.FitParams->at(0)+td.FitParams->at(1)*(tmax-td.FitParams->at(4));
			y2f=td.FitParams->at(2)+td.FitParams->at(3)*(tmax-td.FitParams->at(4));
			t2f=tmax;
			
			tracklength=sqrt(pow((x1f-x2f)*0.254,2)+pow((y1f-y2f)*0.254,2)+pow((tmax-tmin)*0.4*DriftSpeeds[runper],2));
// 			tracklength=sqrt(pow((x1-x2)*0.254,2)+pow((y1-y2)*0.254,2)+pow((tmax-tmin)*0.4*DriftSpeeds[runper],2));
			AllTrackLength->Fill(tracklength);
			
			for(int is1=0;is1<3;is1++)
			{
				if(is1==1) continue;
				pmtdist[is1]=sqrt(pow((x1f-PMTXY[is1][0])*0.254,2)+pow((y1f-PMTXY[is1][1])*0.254,2)+pow(((tmin-tlims[runper][0])*0.4*DriftSpeeds[runper])+4,2));
			}
			
			for(int i1=0;i1<3;i1++)
			{
				if(i1==1) continue;
				
				if(hd.PMTIntegral->at(i1)>0)//unsaturated
				{
					PMTsat=false;
				}
				else
				{
					PMTsat=true;
				}
				PMTintegral=0;
				if(!PMTsat)
				{
					for(int i2=1900;i2<10000;i2++)
					{
						PMTintegral+=(-1.)*(((float)ed.PMTWF->at(i1)[i2])-PMTBaselines[i1][0]);
					}
				}
				else//saturated
				{
					for(int i2=1900;i2<10000;i2++)
					{
						sg=(-1.)*(((float)ed.PMTWF->at(i1)[i2])-PMTBaselines[i1][0]);
						if(sg>PMTTOTTh[i1])
						{
							TOTbinmin=i2;
							break;
						}
						else
						{
							PMTintegral+=sg;
						}
					}
					for(int i2=TOTbinmin+1;i2<10000;i2++)
					{
						sg=(-1.)*(((float)ed.PMTWF->at(i1)[i2])-PMTBaselines[i1][0]);
						if(sg<PMTTOTTh[i1])
						{
							TOTbinmax=i2-1;
							break;
						}
					}
					PMTTOT=((float)TOTbinmax)-((float)TOTbinmin);
					PMTTOTintegral=PMTTOTFit[i1][0]+(PMTTOTFit[i1][1]*PMTTOT)+(PMTTOTFit[i1][1]*pow(PMTTOT,2))+(PMTTOTFit[i1][3]*pow(PMTTOT,3))+(PMTTOTFit[i1][4]*pow(PMTTOT,4));
					for(int i2=TOTbinmax+1;i2<10000;i2++)
					{
						sg=(-1.)*(((float)ed.PMTWF->at(i1)[i2])-PMTBaselines[i1][0]);
						PMTintegral+=sg;
					}
					PMTintegral+=PMTTOTintegral;
				}
				PMTintegrals[i1]=((PMTintegral/SPE[i1][runper][0])/tracklength)/2.1;//npe/MeV
			}
			
			for(int ip1=0;ip1<3;ip1++)
			{
				if(ip1==1) continue;
				PMTIntegral[ip1][2][2]->Fill(PMTintegrals[ip1]);
				PMTIntvsDist[ip1][2]->Fill(pmtdist[ip1],PMTintegrals[ip1]);
				if(hd.PMTIntegral->at(ip1)>0)//unsaturated
				{
					PMTIntegral[ip1][2][0]->Fill(PMTintegrals[ip1]);
					PMTIntvsDist[ip1][0]->Fill(pmtdist[ip1],PMTintegrals[ip1]);
				}
				else
				{
					PMTIntegral[ip1][2][1]->Fill(PMTintegrals[ip1]);
					PMTIntvsDist[ip1][1]->Fill(pmtdist[ip1],PMTintegrals[ip1]);
				}
			}
			
			deltaR=sqrt(pow(x1-x2,2)+pow(y1-y2,2));
			
// 			if(x1>2 && x1<125 && y1>2 && y1<125 && x2>2 && x2<125 && y2>2 && y2<125 && td.FitNormChi2<3)//vertical tracks
// 			if(x1>2 && x1<125 && y1>2 && y1<125 && x2>2 && x2<125 && y2>2 && y2<125)//vertical tracks
// 			if(x1>10 && x1<117 && y1>10 && y1<117 && x2>10 && x2<117 && y2>10 && y2<117)//vertical tracks//this was the original
			if(x1>10 && x1<117 && y1>10 && y1<117 && x2>10 && x2<117 && y2>10 && y2<117 && deltaR>15 && deltaR<35)//vertical tracks
			{
				NColHitsVertical->Fill(hd.ColID->size());
				VerticalChi2->Fill(td.FitNormChi2);
				NColHitsvsChi2Vertical->Fill(td.FitNormChi2,hd.ColID->size());
				NColHitsPerDeltaT->Fill(((float)hd.ColID->size())/(tmax-tmin));
				
				VerticalTrackEntryXY->Fill(x1,y1);
				VerticalTrackExitXY->Fill(x2,y2);
				VerticalTrackDeltaR->Fill(deltaR);
				
// 				tmin=td.ColHitTStartEnd->at(0);tmax=td.ColHitTStartEnd->at(1);
				
				VerticalTrackColHitsDeltaT->Fill(tmax-tmin);
// 				if(deltaR>10 && deltaR<30)
				if((((float)hd.ColID->size())/(tmax-tmin))>0. && (((float)hd.ColID->size())/(tmax-tmin))<0.1)
				{
					VerticalTrackColHitsDeltaTSel->Fill(tmax-tmin);
					if(runper==4 && RunNo>=7316) VerticalTrackColHitsDeltaTSel2->Fill(tmax-tmin);
				}
				VerticalTrackColHitsDeltaTs[J]->Fill(tmax-tmin);
				VerticalTrackTimeofFirstColHit->Fill(tmin);
				VerticalTrackTimeofFirstColHits[J]->Fill(tmin);
				if((tmax-tmin)<500 || (tmax-tmin)>1000)
				{
					VerticalChi2Sel->Fill(td.FitNormChi2);
// 					cout<<RunNos[J]<<" "<<td.E<<" "<<(tmax-tmin)<<" "<<tmin<<" "<<tmax<<endl;
				}
				
				DeltaTvsDeltaR->Fill(deltaR,tmax-tmin);
				DeltaTvsEvent->Fill(NEvtvertical/500,tmax-tmin);
				VerticalTrackLength->Fill(tracklength);
				
				DistancetoPMT[0][0]->Fill(pmtdist[0]);
				DistancetoPMT[1][0]->Fill(pmtdist[2]);
				
				for(int ip1=0;ip1<3;ip1++)
				{
					if(ip1==1) continue;
					PMTIntegral[ip1][0][2]->Fill(PMTintegrals[ip1]);
					if(hd.PMTIntegral->at(ip1)>0)//unsaturated
					{
						PMTIntegral[ip1][0][0]->Fill(PMTintegrals[ip1]);
					}
					else
					{
						PMTIntegral[ip1][0][1]->Fill(PMTintegrals[ip1]);
					}
				}
				
				NEvtvertical++;
			}
			else
			{
				NonVerticalTrackLength->Fill(tracklength);
				DistancetoPMT[0][1]->Fill(pmtdist[0]);
				DistancetoPMT[1][1]->Fill(pmtdist[2]);
				for(int ip1=0;ip1<3;ip1++)
				{
					if(ip1==1) continue;
					PMTIntegral[ip1][1][2]->Fill(PMTintegrals[ip1]);
					if(hd.PMTIntegral->at(ip1)>0)//unsaturated
					{
						PMTIntegral[ip1][1][0]->Fill(PMTintegrals[ip1]);
					}
					else
					{
						PMTIntegral[ip1][1][1]->Fill(PMTintegrals[ip1]);
					}
				}
			}
		}
		inroot->Close();
		inroot2->Close();
	}
	
	outroot->cd();
	
	VerticalTrackDeltaR->Scale(1./VerticalTrackDeltaR->Integral());
	VerticalTrackDeltaR->GetXaxis()->SetTitle("#Delta R (x 2.54 mm)");VerticalTrackDeltaR->GetXaxis()->CenterTitle();
// 	VerticalTrackDeltaR->GetYaxis()->SetTitle("Events / 1.27 mm");VerticalTrackDeltaR->GetYaxis()->CenterTitle();
	VerticalTrackDeltaR->Fit(tf1,"q","q",0.,100);
	VerticalTrackDeltaR->Write();
	
	
	tf2->SetParameter(0,VerticalTrackColHitsDeltaT->GetBinContent(VerticalTrackColHitsDeltaT->GetMaximumBin()));
	tf2->SetParameter(1,VerticalTrackColHitsDeltaT->GetBinCenter(VerticalTrackColHitsDeltaT->GetMaximumBin()));
	tf2->SetParameter(2,10);
	tf2->SetParameter(3,0.1);
	tf2->SetParameter(4,0.1);
	tf2->SetParameter(5,0.1);
// 	VerticalTrackColHitsDeltaT->Fit(tf2,"q","q",800,900);
// 	VerticalTrackColHitsDeltaT->Fit(tf2,"q","q",tf2->GetParameter(1)-2.5*tf2->GetParameter(2),tf2->GetParameter(1)+2*tf2->GetParameter(2));
	VerticalTrackColHitsDeltaT->Fit(tf1,"q","q",800,900);
	VerticalTrackColHitsDeltaT->Fit(tf2,"q","q",tf1->GetParameter(1)-2.*tf1->GetParameter(2),tf1->GetParameter(1)+2*tf1->GetParameter(2));
	
	
// 	VerticalTrackColHitsDeltaT->Fit(tf1,"q","q",800,900);
// 	VerticalTrackColHitsDeltaT->Fit(tf1,"q","q",tf1->GetParameter(1)-1.5*tf1->GetParameter(2),tf1->GetParameter(1)+1.5*tf1->GetParameter(2));
	VerticalTrackColHitsDeltaT->GetXaxis()->SetTitle("#Deltat (x 400 ns)");VerticalTrackColHitsDeltaT->GetXaxis()->CenterTitle();
	VerticalTrackColHitsDeltaT->Write();
	
// 	VerticalTrackColHitsDeltaTSel->Fit(tf2,"q","q",800,900);
// 	VerticalTrackColHitsDeltaTSel->Fit(tf2,"q","q",tf2->GetParameter(1)-2.5*tf2->GetParameter(2),tf2->GetParameter(1)+2*tf2->GetParameter(2));
	VerticalTrackColHitsDeltaTSel->Fit(tf1,"q","q",800,900);
	VerticalTrackColHitsDeltaTSel->Fit(tf2,"q","q",tf1->GetParameter(1)-2.*tf1->GetParameter(2),tf1->GetParameter(1)+2*tf1->GetParameter(2));
	
	
// 	VerticalTrackColHitsDeltaTSel->Fit(tf1,"q","q",800,900);
// 	VerticalTrackColHitsDeltaTSel->Fit(tf1,"q","q",tf1->GetParameter(1)-1.5*tf1->GetParameter(2),tf1->GetParameter(1)+1.5*tf1->GetParameter(2));
	VerticalTrackColHitsDeltaTSel->GetXaxis()->SetTitle("#Deltat (x 400 ns)");VerticalTrackColHitsDeltaTSel->GetXaxis()->CenterTitle();
	VerticalTrackColHitsDeltaTSel->Write();
	
	if(PRunNo==3000)
	{
// 		VerticalTrackColHitsDeltaTSel2->Fit(tf2,"q","q",800,900);
// 		VerticalTrackColHitsDeltaTSel2->Fit(tf2,"q","q",tf2->GetParameter(1)-2.5*tf2->GetParameter(2),tf2->GetParameter(1)+2*tf2->GetParameter(2));
		VerticalTrackColHitsDeltaTSel2->Fit(tf1,"q","q",800,900);
		VerticalTrackColHitsDeltaTSel2->Fit(tf2,"q","q",tf1->GetParameter(1)-2.*tf1->GetParameter(2),tf1->GetParameter(1)+2*tf1->GetParameter(2));
		
		
// 		VerticalTrackColHitsDeltaTSel2->Fit(tf1,"q","q",800,900);
// 		VerticalTrackColHitsDeltaTSel2->Fit(tf1,"q","q",tf1->GetParameter(1)-1.5*tf1->GetParameter(2),tf1->GetParameter(1)+1.5*tf1->GetParameter(2));
		VerticalTrackColHitsDeltaTSel2->GetXaxis()->SetTitle("#Deltat (x 400 ns)");VerticalTrackColHitsDeltaTSel2->GetXaxis()->CenterTitle();
		VerticalTrackColHitsDeltaTSel2->Write();
	}
	
	NColHitsvsChi2Vertical->Write();
	VerticalChi2->Write();
	VerticalChi2Sel->Write();
	NColHits->Write();
	NColHitsVertical->Write();
	NColHitsPerDeltaT->Write();
	TotalCollectionCharge->Write();
	TotalCollectionChargeZS->Write();
	DeltaTvsEvent->Write();
	AllHitT->Write();
	TrackColHitInt->Write();
	TrackColHitAmp->Write();
	
	for(int i1=0;i1<2;i1++)
	{
		for(int i2=0;i2<2;i2++)
		{
			DistancetoPMT[i1][i2]->Write();
		}
	}
	
	VerticalTrackLength->Write();
	NonVerticalTrackLength->Write();
	AllTrackLength->Write();
	VerticalTrackEntryXY->Write();
	VerticalTrackExitXY->Write();
	VerticalTrackTimeofFirstColHit->Write();
	for(int i1=0;i1<RunNos.size();i1++)
	{
		VerticalTrackTimeofFirstColHits[i1]->Write();
		VerticalTrackColHitsDeltaTs[i1]->Write();
	}
	DeltaTvsDeltaR->Write();
	
	for(int i1=0;i1<3;i1++)
	{
		if(i1==1) continue;
		for(int i2=0;i2<3;i2++)
		{
			for(int i3=0;i3<3;i3++)
			{
				PMTIntegral[i1][i2][i3]->Write();
			}
		}
	}
	for(int i1=0;i1<3;i1++)
	{
		if(i1==1) continue;
		for(int i2=0;i2<3;i2++)
		{
			PMTIntvsDist[i1][i2]->Write();
		}
	}
	
	outroot->Close();
// 	sprintf(hname,"cp PMTAnalysis_%d.root %s/Histos/PMTAnalysis_%d.root;wait;",PRunNo,AnalysisFilePath,RunNo);system(hname);
}

void GenerateMCInput()
{
	int PRunNo=RunNo;
	vector <int> RunNos;
	if(RunNo==1000)
	{
		//up to 7095 different PMT voltages
		RunNos.push_back(7074);//closed first
		RunNos.push_back(7076);//closed first
		RunNos.push_back(7078);//closed first
		RunNos.push_back(7082);//closed first
		RunNos.push_back(7084);//closed first
		RunNos.push_back(7092);//closed first
		RunNos.push_back(7094);//closed first
		RunNos.push_back(7096);//closed first
		RunNos.push_back(7104);//bump in slope//closed first
		RunNos.push_back(7106);//closed first
		
		RunNos.push_back(7107);
		RunNos.push_back(7110);
		RunNos.push_back(7121);
	}
	else if(RunNo==2000)
	{
		RunNos.push_back(7252);//closed first
		RunNos.push_back(7254);//closed first
		RunNos.push_back(7256);//closed first
		RunNos.push_back(7259);//closed first
		
		RunNos.push_back(7268);
		RunNos.push_back(7271);
		RunNos.push_back(7273);
	}
	else if(RunNo==3000)
	{
		RunNos.push_back(7308);
		RunNos.push_back(7312);
		RunNos.push_back(7314);
		RunNos.push_back(7316);
		RunNos.push_back(7318);
		RunNos.push_back(7320);
		RunNos.push_back(7322);
		RunNos.push_back(7325);//closed first
		RunNos.push_back(7327);//closed first
		RunNos.push_back(7385);//closed first
		RunNos.push_back(7388);
	}
	
	float SPE[3][5][2]={{{51.74,0.04},{38.53,0.06},{29.66,0.02},{32.04,0.01},{32.38,0.01}},{{0,0},{0,0},{0,0},{0,0},{0,0}},{{36.59,0.11},{28.46,0.07},{23.51,0.01},{22.65,0.01},{21.23,0.02}}};
	float tlims[5][2]={{1100,2040},{1100,2040},{1100,2040},{1150,2100},{1440,2350}};
	float DriftSpeeds[5]={0.149,0.149,0.149,0.147,0.150};//cm/us
	int runper=0;//0:7000-7087; 1:7090-7094; 2-7095-7144; 3:7200-7287; 4:7300-7403
	
	ifstream infile;
	sprintf(hname,"MCTrackParameters_%d_1.txt",PRunNo);
	ofstream outfile(hname);
	int aa=0;float bb[6]={0};
	int i=0;
	
	for(int J=0;J<RunNos.size();J++)
	{
		RunNo=RunNos[J];
		cout<<RunNo<<endl;
		
		if(RunNo>=7000 && RunNo<=7087) runper=0;
		else if(RunNo>=7090 && RunNo<=7094) runper=1;
		else if(RunNo>=7095 && RunNo<=7144) runper=2;
		else if(RunNo>=7200 && RunNo<=7287) runper=3;
		else if(RunNo>=7300 && RunNo<=7403) runper=4;
		
		sprintf(hname,"%s/Files/ValidTrackParameters_%d.txt",AnalysisFilePath,RunNo);
		cout<<hname<<endl;
		
		infile.open(hname);
		while(!infile.eof())
		{
// 			cout<<i<<endl;
			i++;
			infile>>aa>>bb[0]>>bb[1]>>bb[2]>>bb[3]>>bb[4]>>bb[5];
// 			cout<<aa<<" "<<bb[0]<<" "<<bb[1]<<" "<<bb[2]<<" "<<bb[3]<<" "<<bb[4]<<" "<<bb[5]<<endl;
			outfile<<RunNo<<" "<<aa<<" "<<((bb[0]-63.5)*2.54)<<" "<<((bb[1]-63.5)*2.54)<<" "<<(((bb[2]-tlims[runper][0])*DriftSpeeds[runper]*10.)+227.5)<<" "<<((bb[3]-63.5)*2.54)<<" "<<((bb[4]-63.5)*2.54)<<" "<<(((bb[5]-tlims[runper][0])*DriftSpeeds[runper]*10.)+227.5)<<endl;
		}
		infile.close();
	}
	outfile.close();
}

int main( int argc, const char* argv[] )
{
	if(argc==1)
	{
		cout<<"0 : GetTPCBaseline"<<endl;
		cout<<"1 : WriteTPCWF -- provide TPCBaselineRun"<<endl;
		cout<<"2 : FindHits -- provide TPCBaselineRun"<<endl;
		cout<<"3 : PlotHitsSummary"<<endl;
		cout<<"4 : PlotHits"<<endl;
		cout<<"5 : FitTracks"<<endl;
		cout<<"6 : GetPMTBaseline"<<endl;
		cout<<"7 : GetPMTCalibration -- provide PMTBaselineRun"<<endl;
		cout<<"8 : WritePMTWF -- provide PMTBaselineRun"<<endl;
		cout<<"9 : GetPMTCalibration2"<<endl;
		cout<<"10 : GetPMTIntegral -- provide PMTBaselineRun"<<endl;
		return -1;
	}
	opt=atoi(argv[1]);
	RunNo=atoi(argv[2]);
	
// 	FillICCH();
	
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
	
	TPCBaselineRun=FindBaselineRun(RunNo);
	PMTBaselineRun=TPCBaselineRun;
	
	if(opt==0)
	{
		GetBaseline();
	}
	if(opt==1)
	{
// 		TPCBaselineRun=atoi(argv[3]);
// 		PMTBaselineRun=TPCBaselineRun;
		WriteWF();
	}
	if(opt==2)
	{
// 		TPCBaselineRun=atoi(argv[3]);
// 		PMTBaselineRun=TPCBaselineRun;
// 		FindHits();
		FindHitsAndFitTracks();
	}
	if(opt==3)
	{
		FitTracks();
	}
	if(opt==4)
	{
// 		TPCBaselineRun=atoi(argv[3]);
		ColWireLines();
	}
	
	
	if(opt==6)
	{
// 		PMTBaselineRun=atoi(argv[3]);
		GetPMTCalibration3();
	}
	if(opt==7)
	{
// 		PMTBaselineRun=atoi(argv[3]);
		GetPMTCalibration();
	}
	if(opt==8)
	{
// 		PMTBaselineRun=atoi(argv[3]);
// 		WritePMTWF();
	}
	if(opt==9)
	{
// 		PMTBaselineRun=atoi(argv[3]);
// 		if(RunNo>7000){PMTBaselineRun=atoi(argv[3]);}
		GetPMTCalibration2();
	}
	if(opt==10)
	{
// 		PMTBaselineRun=atoi(argv[3]);
		GetPMTIntegral();
	}
	
	if(opt==100)
	{
		CompareTPCBaselines();
	}
	if(opt==101)
	{
		CombinedTPCPMTPlot();
	}
	if(opt==102)
	{
		ShowerAnalysis();
	}
	if(opt==103)
	{
// 		TPCBaselineRun=atoi(argv[3]);
		FindDeltaT();
	}
	if(opt==104)
	{
		TrackAnalysis();
	}
	if(opt==105)
	{
		MultiViewPlot();
	}
	if(opt==106)
	{
		PrintNevt();
	}
	if(opt==107)
	{
// 		TPCBaselineRun=atoi(argv[3]);
// 		PMTBaselineRun=TPCBaselineRun;
		HitAnalysis();
	}
	if(opt==108)
	{
// 		TPCBaselineRun=atoi(argv[3]);
		HitAnalysis2();
	}
	if(opt==109)
	{
// 		TPCBaselineRun=atoi(argv[3]);
		TPCCalibration();
	}
	if(opt==110)
	{
		TPCStats();
	}
	if(opt==111)
	{
		TrackAnalysis3();
	}
	if(opt==112)
	{
		PMTAnalysis();
	}
	if(opt==113)
	{
		PMTCalibrationsAll();
	}
	if(opt==114)
	{
		TrackAndPMTAnalysis();
	}
	if(opt==115)
	{
		PMTCalibrationsAll2();
	}
	if(opt==116)
	{
		GenerateMCInput();
	}
	if(opt==117)
	{
		GetPMTCalibration_CalibRuns();
	}
	if(opt==118)
	{
		GetPMTCalibration4();
	}
}

