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
//#include "Math/Vector3Dfwd.h"
#include "Fit/Fitter.h"
 
#include <netinet/in.h>

using namespace std;
using namespace ROOT::Math;

struct edata
{
	vector <vector<int> >* PMTWF;
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
int TimeRange[2]={750,2500};
float HitColIndDist=5.;
int NPMTs=8;//number of PMTs

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

// string PMTNames[8]={"LAr_wTPB","LAr_noTPB","LXe","LAr_wTPB_new","SC1","SC2"};
string PMTNames[8]={"LAr_noTPB","LAr_wTPB","LAr_wTPB_new","LXe","SC1","SiPM_Quartz","SiPM_TPB","S13371_Window"};
float baselines[8]={903.4,913.9,920.8,928.1,916.1,932.5,1,1};

int main( int argc, const char* argv[] )
{
	RunNo=atoi(argv[1]);
	
// 	//eos
 	sprintf(TPCinFilePath,"/eos/project/f/flic2019/Data/TPC/Runs");//where the Run00X.dat files reside
 	sprintf(PMTinFilePath,"/eos/project/f/flic2019/Data/PMT/WaveForms");//where the waveforms_chX.txt files reside
	sprintf(AnalysisFilePath,"/eos/project/f/flic2019/Analysis");//where Files and Histos will be placed
	sprintf(NTupleFilePath,"/eos/project/f/flic2019/Data/NTuple");//where the Run_XX.root files reside

	sprintf(hname,"%s/Run_%d.root",NTupleFilePath,RunNo);
	TFile* inroot=new TFile(hname);
	TTree* T =  (TTree*) inroot->Get("T");
	T->SetBranchAddress("PMTWF",&ed.PMTWF);
	
	sprintf(hname,"NitroXen_%d.root",RunNo);
	TFile* outroot=new TFile(hname,"recreate");
	
	TH1D* AllWF[8];
	TH1F* Baselines[8];
	TH1F* Integral_1000_10000[8];
	TH1F* Integral_1500_15000[8];
	TH1F* PE_1500_15000[8];
	TH1F* Integral_0_1800[8];
	TH1F* Integral_negative[8];
	TH1F* Frac_1800_2500_15000[8];
	TH1F* Integral_1500_8000[8];
	TH1F* PE_1500_8000[8];
	TH1F* Integral_1800_6000[8];
	TH1F* PE_1800_6000[8];
	TH1F* Integral_1500_15000_ZS[8];
	TH1F* Integral_10percent[8];
	TH2F* SC2_vs_SC1;
	TH2F* Int_neg_vs_pos[8];
	SC2_vs_SC1=new TH2F("SC2_vs_SC1","SC2_vs_SC1",4000,0,100000,4000,0,100000);
	TH1D* AvgWF[8][2];
	TH1D* AllUnsatAvgWF[8][2];
	TH1D* AvgWF_shifted[8][2];
	TH1D* AvgWF_NZintegral[8][2];
	
// 	float calibs[2][4]={{36.04,38.41,23.59,12.54},{36.04,88.21,23.59,21.37}};
// 	float calibs[2][4]={{35.01,36.16,23.92,16.92},{35.01,89.29,23.92,20.27}};
	float calibs[2][8]={{35.01,36.16,23.92,16.92},{35.01,89.29,23.92,20.27}};
	int calibID=0;
// 	if(RunNo>=299) calibID=1;
// 	else calibID=0;
	
	int mcolors[8]={3,1,4,2,1,2,3,4};
	
	TGraph* tg[8];TGraph* tgPE[8];
	for(int i1=0;i1<8;i1++)
	{
		sprintf(hname,"Signal_vs_Event_%s",PMTNames[i1].c_str());
		tg[i1]=new TGraph();
		tg[i1]->SetName(hname);tg[i1]->SetTitle(hname);
		tg[i1]->GetXaxis()->SetTitle("Event");tg[i1]->GetXaxis()->CenterTitle();
		tg[i1]->GetYaxis()->SetTitle("Charge (au)");tg[i1]->GetYaxis()->CenterTitle();
		tg[i1]->SetMarkerStyle(20);
		tg[i1]->SetMarkerColor(mcolors[i1]);
		
		sprintf(hname,"PE_vs_Event_%s",PMTNames[i1].c_str());
		tgPE[i1]=new TGraph();
		tgPE[i1]->SetName(hname);tg[i1]->SetTitle(hname);
		tgPE[i1]->GetXaxis()->SetTitle("Event");tg[i1]->GetXaxis()->CenterTitle();
		tgPE[i1]->GetYaxis()->SetTitle("Mean Number of Photoelectrons");tg[i1]->GetYaxis()->CenterTitle();
		tgPE[i1]->SetMarkerStyle(20);
		tgPE[i1]->SetMarkerColor(mcolors[i1]);
	}
	for(int i1=0;i1<8;i1++)
	{
		sprintf(hname,"AllWF_%s",PMTNames[i1].c_str());
		AllWF[i1]=new TH1D(hname,hname,15000,-0.5,14999.5);
		sprintf(hname,"Baselines_%s",PMTNames[i1].c_str());
		Baselines[i1]=new TH1F(hname,hname,2400,0,1200);
		sprintf(hname,"Integral_1000_10000_%s",PMTNames[i1].c_str());
		Integral_1000_10000[i1]=new TH1F(hname,hname,600,0,30000);
		sprintf(hname,"Integral_1500_15000_%s",PMTNames[i1].c_str());
		Integral_1500_15000[i1]=new TH1F(hname,hname,600,0,30000);
		sprintf(hname,"Integral_1500_8000_%s",PMTNames[i1].c_str());
		Integral_1500_8000[i1]=new TH1F(hname,hname,600,0,30000);
		sprintf(hname,"Integral_1800_6000_%s",PMTNames[i1].c_str());
		Integral_1800_6000[i1]=new TH1F(hname,hname,600,0,30000);
		sprintf(hname,"Integral_0_1800_%s",PMTNames[i1].c_str());
		Integral_0_1800[i1]=new TH1F(hname,hname,1000,-10000,10000);
		sprintf(hname,"Frac_1800_2500_15000_%s",PMTNames[i1].c_str());
		Frac_1800_2500_15000[i1]=new TH1F(hname,hname,2000,0,2);
		sprintf(hname,"Integral_1500_15000_ZS_%s",PMTNames[i1].c_str());
		Integral_1500_15000_ZS[i1]=new TH1F(hname,hname,600,0,30000);
		sprintf(hname,"Integral_negative_%s",PMTNames[i1].c_str());
		Integral_negative[i1]=new TH1F(hname,hname,1000,-20000,0);
		sprintf(hname,"Integral_10percent_%s",PMTNames[i1].c_str());
		Integral_10percent[i1]=new TH1F(hname,hname,600,0,30000);
		sprintf(hname,"AvgWF_%s",PMTNames[i1].c_str());
		AvgWF[i1][0]=new TH1D(hname,hname,15000,-0.5,14999.5);
		sprintf(hname,"AvgWF_%s_norm",PMTNames[i1].c_str());
		AvgWF[i1][1]=new TH1D(hname,hname,15000,-0.5,14999.5);
		sprintf(hname,"AllUnsatAvgWF_%s",PMTNames[i1].c_str());
		AllUnsatAvgWF[i1][0]=new TH1D(hname,hname,15000,-0.5,14999.5);
		sprintf(hname,"AllUnsatAvgWF_%s_norm",PMTNames[i1].c_str());
		AllUnsatAvgWF[i1][1]=new TH1D(hname,hname,15000,-0.5,14999.5);
		sprintf(hname,"PE_1500_15000_%s",PMTNames[i1].c_str());
		PE_1500_15000[i1]=new TH1F(hname,hname,400,0,1000);
		sprintf(hname,"PE_1500_8000_%s",PMTNames[i1].c_str());
		PE_1500_8000[i1]=new TH1F(hname,hname,400,0,1000);
		sprintf(hname,"PE_1800_6000_%s",PMTNames[i1].c_str());
		PE_1800_6000[i1]=new TH1F(hname,hname,400,0,1000);
		sprintf(hname,"AvgWFShifted_%s",PMTNames[i1].c_str());
		AvgWF_shifted[i1][0]=new TH1D(hname,hname,15000,-0.5,14999.5);
		sprintf(hname,"AvgWFShifted_%s_norm",PMTNames[i1].c_str());
		AvgWF_shifted[i1][1]=new TH1D(hname,hname,15000,-0.5,14999.5);
		sprintf(hname,"AvgWF_NZintegral_%s",PMTNames[i1].c_str());
		AvgWF_NZintegral[i1][0]=new TH1D(hname,hname,15000,-0.5,14999.5);
		sprintf(hname,"AvgWF_NZintegral_%s_norm",PMTNames[i1].c_str());
		AvgWF_NZintegral[i1][1]=new TH1D(hname,hname,15000,-0.5,14999.5);
		
		sprintf(hname,"Int_neg_vs_pos_%s",PMTNames[i1].c_str());
		Int_neg_vs_pos[i1]=new TH2F(hname,hname,1000,-20000,0,1000,0,20000);
	}
	
	TH1F* hh=new TH1F("PMTWF","PMTWF",15000,-0.5,14999.5);
	TH1F* hh_BS=new TH1F("PMTWF_BS","PMTWF_BS",15000,-0.5,14999.5);
	TH1F* hh1=new TH1F("PMTWF1","PMTWF1",15000,-0.5,14999.5);
	TH1F* hh_BS_shifted=new TH1F("PMTWF_BS_shifted","PMTWF_BS_shifted",15000,-0.5,14999.5);
	TH1F* hh_BS_shifted_norm=new TH1F("PMTWF_BS_shifted_norm","PMTWF_BS_shifted_norm",15000,-0.5,14999.5);
	
	for(int i1=1;i1<=hh1->GetNbinsX();i1++){hh1->SetBinContent(i1,1);}
	
	TF1* g1=new TF1("g1","gaus",-1000.,2000.);
// 	TF1* tf1=new TF1("tf1","[0]*(exp([0]+[1]*(x-[2]))+exp([3]+[4]*(x-[2]))+exp([5]+[8]*(x-[2]))+exp([7]+[8]*(x-[2])))",0.,15000.);
// 	TF1* tf1=new TF1("tf1","(exp([0]+[1]*(x-[2]))+exp([3]+[4]*(x-[2]))+exp([5]+[8]*(x-[2])))",0.,15000.);
// 	TF1* tf1=new TF1("tf1","gaus(0)+expo(3)+expo(5)",0.,15000.);
// 	TF1* tf1=new TF1("tf1","[0]*pow(x-[1],[2])*([3]*exp(-(x-[1])/[4])+[5]*exp(-(x-[1])/[8]))",0.,15000.);
// 	TF1* tf1=new TF1("tf1","[0]*pow(x-[1],[2])*(exp(-(x-[1])/[3])+exp(-(x-[1])/[4]))",0.,15000.);
// 	TF1* tf1=new TF1("tf1","[0]*exp(-0.5*(pow((log((x)/[2]))/[3],2)))",0.,15000.);
// 	TF1* tf1=new TF1("tf1","landau(0)+expo(3)",0.,15000.);
// 	TF1* tf1=new TF1("tf1","[0]*exp(-pow((x-[1])/[2],2))+[8]*exp(-x/[3])+[7]*exp(-x/[4])+[8]*exp(-x/[5])",0.,15000.);
// 	TF1 *tf1 = new TF1("fit",fitf2f,0.,15000,6);
// 	TF1 *tf1 = new TF1("fit",fitf3f,0.,15000,7);
//  	TF1* tf1=new TF1("tf1","([0]/[3])*exp(-(log(x-[1])-[2])^2/(2.*[3]^2))",0.,15000.);
	TF1* tf1=new TF1("tf1","expo",0.,15000.);
	TF1* tflin=new TF1("tflin","[0]",0.,15000.);
	TF1* tflin2=new TF1("tflin2","[0]+[1]*x",0.,15000.);
	
	vector <float> Amp;
	float hsum[20]={0.};
	float hped=0;float hsig=0;
	float hped2=0;float hsig2=0;
	float xlim=0;int fitstartbin=0;int xlimbin=0;
	float qints_1000_10000[8]={0};
	float qints_1500_15000[8]={0};
	float qints_1500_8000[8]={0};
	float qints_1800_6000[8]={0};
	float qints_1800_2500[8]={0};
	float qints_1800_15000[8]={0};
	float qints_0_1800[8]={0};
	float qints_negative[8]={0};
	float qints_positive[8]={0};
	float qints_1500_15000_ZS[8]={0};
	float qints_10percent[8]={0};
	float hpeds[8]={0};
	int sat[8]={0};int badBaseline[8]={0};
	int NSat[8]={0};
	
	float hmax=0;float hmax75=0;float hmax25=0;
	float xmax=0;float xmax75=0;float xmax25=0;
	int ibinmax=0;int ibinmax75=0;int ibinmax25=0;
	float xoffset=0;
	int binoffset=0;
	
	int sbin=0;int minbin=0;int maxbin=0;int peakbin=0;
	float maxamp=0;int nc=0;
	
	for(int I=0;I<T->GetEntries();I++)
	{
		if(RunNo==426 && I<1250) continue;//Xe doping to 10 ppm
		if(RunNo==435 && I<1050) continue;//Xe doping to 10 ppm
		if(RunNo==449 && I<1100) continue;//Xe doping to 10 ppm
		
		T->GetEntry(I);
		for(int i1=0;i1<8;i1++)
		{
			qints_1000_10000[i1]=0;
			qints_1500_15000[i1]=0;
			qints_1500_8000[i1]=0;
			qints_1800_6000[i1]=0;
			qints_1800_2500[i1]=0;
			qints_1800_15000[i1]=0;
			qints_0_1800[i1]=0;
			qints_1500_15000_ZS[i1]=0;
			qints_10percent[i1]=0;
			qints_negative[i1]=0;
			qints_positive[i1]=0;
			sat[i1]=0;
			badBaseline[i1]=0;
		}
		for(int i1=0;i1<8;i1++)
		{
			for(int i2=0;i2<ed.PMTWF->at(i1).size();i2++)
			{
				hh->SetBinContent(i2+1,ed.PMTWF->at(i1)[i2]);
				if(ed.PMTWF->at(i1)[i2]==0) sat[i1]=1;
			}
			if(hh->GetNbinsX()==0) continue;
			AllWF[i1]->Add(hh);
			hh->Fit(tflin,"q","q",0.,1500.);
			Baselines[i1]->Fill(tflin->GetParameter(0));
			hped=tflin->GetParameter(0);
			hpeds[i1]=hped;
			
			for(int i2=0;i2<ed.PMTWF->at(i1).size();i2++)
			{
				hsig=(-1.*(((float)ed.PMTWF->at(i1)[i2])-hpeds[i1]));
// 				hsig=(-1.*(((float)ed.PMTWF->at(i1)[i2])-baselines[i1]));
				hh_BS->SetBinContent(i2+1,hsig);
				if(i2>=1000 && i2<=10000) qints_1000_10000[i1]+=hsig;
				if(i2>=1500)
				{
					qints_1500_15000[i1]+=hsig;
					if(hsig>0)
					{
						qints_1500_15000_ZS[i1]+=hsig;
					}
					if(i2<=8000)
					{
						qints_1500_8000[i1]+=hsig;
					}
				}
				if(i2>=1800 && i2<=6000)
				{
					qints_1800_6000[i1]+=hsig;
				}
				if(i2>=1800 && i2<=2500)
				{
					qints_1800_2500[i1]+=hsig;
				}
				if(i2>=1800 && i2<=15000)
				{
					qints_1800_15000[i1]+=hsig;
				}
				if(i2<1800)
				{
					qints_0_1800[i1]+=hsig;
				}
				if(hsig<0){qints_negative[i1]+=hsig;}
				if(hsig>0){qints_positive[i1]+=hsig;}
			}
			Integral_0_1800[i1]->Fill(qints_0_1800[i1]);
			Integral_negative[i1]->Fill(qints_negative[i1]);
			Int_neg_vs_pos[i1]->Fill(qints_negative[i1],qints_positive[i1]);
			if(fabs(qints_0_1800[i1])>100) badBaseline[i1]=1;
			if(sat[i1]!=1 && badBaseline[i1]!=1)
			{
				Frac_1800_2500_15000[i1]->Fill(qints_1800_2500[i1]/qints_1800_15000[i1]);
				Integral_1000_10000[i1]->Fill(qints_1000_10000[i1]);
				Integral_1500_15000[i1]->Fill(qints_1500_15000[i1]);
				Integral_1500_8000[i1]->Fill(qints_1500_8000[i1]);
				Integral_1800_6000[i1]->Fill(qints_1800_6000[i1]);
// 				if(qints_1800_6000[i1]<500 && i1==0) cout<<I<<" "<<qints_1800_6000[i1]<<endl;
				Integral_1500_15000_ZS[i1]->Fill(qints_1500_15000_ZS[i1]);
				if(i1<4)
				{
					PE_1500_15000[i1]->Fill(qints_1500_15000[i1]/calibs[calibID][i1]);
					PE_1500_8000[i1]->Fill(qints_1500_8000[i1]/calibs[calibID][i1]);
					PE_1800_6000[i1]->Fill(qints_1800_6000[i1]/calibs[calibID][i1]);
					tgPE[i1]->SetPoint(tgPE[i1]->GetN(),I,qints_1500_15000[i1]/calibs[calibID][i1]);
// 					if(qints_1800_6000[i1]<500 && i1==0) cout<<I<<" "<<i1<<" "<<qints_1800_6000[i1]<<endl;
				}
			}
			
			peakbin=hh_BS->GetMaximumBin();
			maxamp=hh_BS->GetBinContent(peakbin);
			for(int i2=peakbin;i2>=1;i2--)
			{
				if((hh_BS->GetBinContent(i2)/maxamp)<0.1){minbin=i2;break;}
			}
			for(int i2=peakbin;i2<=hh_BS->GetNbinsX();i2++)
			{
				if((hh_BS->GetBinContent(i2)/maxamp)<0.1){maxbin=i2;break;}
			}
			for(int i2=minbin;i2<=maxbin;i2++)
			{
				qints_10percent[i1]+=hh_BS->GetBinContent(i2);
			}
			if(i1<4)
			{
				if(sat[i1]!=1 && badBaseline[i1]!=1)
				{
					Integral_10percent[i1]->Fill(qints_10percent[i1]);
				}
			}
			else
			{
				Integral_10percent[i1]->Fill(qints_10percent[i1]);
			}
			tg[i1]->SetPoint(tg[i1]->GetN(),I,qints_1500_15000[i1]);
// 			if(i1==0 && qints_1500_15000[i1]<-500) cout<<RunNo<<" "<<" "<<I<<" "<<qints_1500_15000[i1]<<endl;
			if(sat[i1]!=1 && badBaseline[i1]!=1)
			{
				AllUnsatAvgWF[i1][0]->Add(hh_BS);
				AllUnsatAvgWF[i1][1]->Add(hh1);
			}
			else
			{
				NSat[i1]++;
			}
			hh->Reset();
			hh_BS->Reset();
		}
// 		if(qints_10percent[4]>5000 && qints_10percent[4]<22000 && qints_10percent[5]>1500 && qints_10percent[5]<4000)
// 		if(qints_10percent[4]>2000)
		{
			for(int i1=0;i1<8;i1++)
			{
				if(sat[i1]==1 || badBaseline[i1]==1) continue;
				hh_BS->Reset();
				for(int i2=0;i2<ed.PMTWF->at(i1).size();i2++)
				{
					hh_BS->SetBinContent(i2+1,(-1.*(((float)ed.PMTWF->at(i1)[i2])-hpeds[i1])));
// 					hh_BS->SetBinContent(i2+1,(-1.*(((float)ed.PMTWF->at(i1)[i2])-baselines[i1])));
				}
				AvgWF[i1][0]->Add(hh_BS);
				AvgWF[i1][1]->Add(hh1);
				if(qints_1500_15000[i1]>0)
				{
					AvgWF_NZintegral[i1][0]->Add(hh_BS);
					AvgWF_NZintegral[i1][1]->Add(hh1);
				}
				
				hh_BS_shifted->Reset();
				hh_BS_shifted_norm->Reset();
				
				ibinmax=hh_BS->GetMaximumBin();
				
				binoffset=ibinmax;
				
// 				hmax=hh_BS->GetBinContent(ibinmax);
// 				xmax=hh_BS->GetBinCenter(ibinmax);
// 				for(int is1=ibinmax;is1>=1;is1--)
// 				{
// 					if((hh_BS->GetBinContent(is1)/hmax)<0.75){ibinmax75=is1;hmax75=hh_BS->GetBinContent(is1);xmax75=hh_BS->GetBinCenter(is1);break;}
// 				}
// 				for(int is1=ibinmax75;is1>=1;is1--)
// 				{
// 					if((hh_BS->GetBinContent(is1)/hmax)<0.25){ibinmax25=is1;hmax25=hh_BS->GetBinContent(is1);xmax25=hh_BS->GetBinCenter(is1);break;}
// 				}
// 				hh_BS->Fit(tflin2,"q","q",xmax25,xmax75);
// 				xoffset=-1.*tflin2->GetParameter(0)/tflin2->GetParameter(1);
// 				binoffset=hh_BS->FindBin(xoffset);
				
				for(int i2=0;i2<ed.PMTWF->at(i1).size();i2++)
				{
					if((binoffset-2000+i2)>=1 && (binoffset-2000+i2)<=hh_BS->GetNbinsX())
					{
						hh_BS_shifted->SetBinContent(binoffset-2000+i2,(-1.*(((float)ed.PMTWF->at(i1)[i2])-hpeds[i1])));
// 						hh_BS_shifted->SetBinContent(binoffset-2000+i2,(-1.*(((float)ed.PMTWF->at(i1)[i2])-baselines[i1])));
						hh_BS_shifted_norm->SetBinContent(binoffset-2000+i2,1);
					}
				}
				AvgWF_shifted[i1][0]->Add(hh_BS_shifted);
				AvgWF_shifted[i1][1]->Add(hh_BS_shifted_norm);
			}
		}
		SC2_vs_SC1->Fill(qints_10percent[5],qints_10percent[4]);
		if(I%1000==0) cout<<I<<" / "<<T->GetEntries()<<endl;
	}
	
	outroot->cd();
	for(int i1=0;i1<8;i1++)
	{
		AllWF[i1]->Write();
		Baselines[i1]->Write();
		Integral_1000_10000[i1]->Write();
		Integral_1500_15000[i1]->Write();
		Integral_1500_8000[i1]->Write();
		Integral_1800_6000[i1]->Write();
		Integral_0_1800[i1]->Write();
		Integral_1500_15000_ZS[i1]->Write();
		Integral_10percent[i1]->Write();
		Integral_negative[i1]->Write();
		Int_neg_vs_pos[i1]->Write();
		Frac_1800_2500_15000[i1]->Write();
		AvgWF[i1][0]->Divide(AvgWF[i1][1]);
		AvgWF[i1][0]->Write();
		AvgWF_NZintegral[i1][0]->Divide(AvgWF_NZintegral[i1][1]);
		AvgWF_NZintegral[i1][0]->Write();
		AvgWF_shifted[i1][0]->Divide(AvgWF_shifted[i1][1]);
		AvgWF_shifted[i1][0]->Write();
		AllUnsatAvgWF[i1][0]->Divide(AllUnsatAvgWF[i1][1]);
		AllUnsatAvgWF[i1][0]->Write();
// 		if(i1<4)
		{
			PE_1500_15000[i1]->Write();
			PE_1500_8000[i1]->Write();
			PE_1800_6000[i1]->Write();
			tg[i1]->Write();
			tgPE[i1]->Write();
		}
	}
	SC2_vs_SC1->Write();
	outroot->Close();
	
	cout<<"Run : "<<RunNo<<" Total events : "<<T->GetEntries()<<" Saturations :"<<endl;
	for(int i1=0;i1<8;i1++)
	{
		cout<<PMTNames[i1]<<" : "<<NSat[i1]<<endl;
	}
	
// 	sprintf(hname,"cp PMTCalibration3_%d.root %s/Histos/PMTCalibration3_%d.root",RunNo,AnalysisFilePath,RunNo);system(hname);
}


