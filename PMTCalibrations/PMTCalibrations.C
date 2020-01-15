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

int RunNo=0;
char hname[1000];
char hname2[1000];

// string PMTNames[6]={"LAr_wTPB","LAr_noTPB","LXe","LAr_wTPB_new","SC1","SC2"};
string PMTNames[8]={"LAr_noTPB","LAr_wTPB","LAr_wTPB_new","LXe","SC1","SiPM_Quartz","SiPM_TPB","S13371_Window"};

int main( int argc, const char* argv[] )
{
	sprintf(hname,"XArapucaXenPMTCalibrations.root");
	TFile* outroot=new TFile(hname,"recreate");
	TFile* inroot;
	
// 	int colors[8]={1,3,2,4,1,2,3,4};
	int colors[8]={3,1,4,2,1,2,3,4};
	
// 	int NRuns=25;
// 	int RunList[25]={408,409,410,411,412,419,420,424,425,428,429,433,434,440,441,445,447,451,452,456,458,466,467,472,473};
//	int NRuns=8;
//	int RunList[8]={425,428,429,441,451,452,456,458};
	int NRuns=13;
        int RunList[13]={408,409,411,412,419,420,424,433,434,466,467,472,473};	
	TH1F* hh[6];
	TGraphErrors* Calibrations[4];
	TGraphErrors* Calibrations_1;
	for(int i1=0;i1<4;i1++)
	{
		sprintf(hname,"Calibrations_%s",PMTNames[i1].c_str());
		Calibrations[i1]=new TGraphErrors();
		Calibrations[i1]->SetMarkerStyle(20+i1);Calibrations[i1]->SetMarkerColor(colors[i1]);
		Calibrations[i1]->GetXaxis()->SetTitle("Run ID");Calibrations[i1]->GetXaxis()->CenterTitle();
		Calibrations[i1]->GetYaxis()->SetTitle("Mean Charge/photoelectron (x20 fC)");Calibrations[i1]->GetYaxis()->CenterTitle();
	}
	Calibrations_1=new TGraphErrors();
	Calibrations_1->SetMarkerStyle(21);Calibrations_1->SetMarkerColor(colors[1]);
	Calibrations_1->GetXaxis()->SetTitle("Run ID");Calibrations_1->GetXaxis()->CenterTitle();
	Calibrations_1->GetYaxis()->SetTitle("Mean Charge/photoelectron (x20 fC)");Calibrations_1->GetYaxis()->CenterTitle();
	int np=0;
	
	TF1* tf1=new TF1("tf1","[0]*exp((x-[1])/[2])/([5]*exp((x-[1])/[3])+[6]*exp((x-[1])/[4]))",0.,15000.);
	TF1* tf2=new TF1("tf2","[0]*pow(x-[1],[2])*([3]*exp(-(x-[1])/[4])+[5]*exp(-(x-[1])/[6]))",0.,15000.);
	TF1* tf3=new TF1("tf3","[0]*exp(-(x-[1])/[2])+[3]*exp(-(x-[1])/[4])+[5]*exp(-(x-[1])/[6])",0.,15000.);
	TF1* tf4=new TF1("tf4","[0]*exp(-(x-[1])/[2])+[3]*exp(-(x-[1])/[4])",0.,15000.);
	TF1* tfl=new TF1("tfl","landau",0.,1000.);
	TF1* tflin=new TF1("tflin","[0]",0.,500.);
	TF1* tfh;
	
	TGraph* axisSC=new TGraph();
	axisSC->SetPoint(0,RunList[0]-1,0);
	axisSC->SetPoint(1,RunList[NRuns-1]+1,100);
	
	int binmax=0;int bin5percent=0;
	float valmax=0;float val5percent=0;
	float xmax=0;float x5percent=0;
	
	for(int i1=0;i1<NRuns;i1++)
	{
		RunNo=RunList[i1];
		cout<<RunNo<<endl;
		sprintf(hname,"/eos/project/f/flic2019/Analysis/Histos/PMTCalibration4_%d.root",RunNo);
		inroot=new TFile(hname);
		for(int i2=0;i2<4;i2++)
		{
			if(RunNo>210 && RunNo<297 && i2==3) continue;
			sprintf(hname,"Integral_%d",i2);
			inroot->GetObject(hname,hh[i2]);
			tfh=(TF1*)hh[i2]->GetFunction("g1");
			
			np=Calibrations[i2]->GetN();
			Calibrations[i2]->SetPoint(np,RunNo,tfh->GetParameter(1));
			Calibrations[i2]->SetPointError(np,0,tfh->GetParError(1));
			if(i2==1)
			{
				np=Calibrations_1->GetN();
				Calibrations_1->SetPoint(np,RunNo,tfh->GetParameter(1));
				Calibrations_1->SetPointError(np,0,tfh->GetParError(1));
			}
		}
		inroot->Close();
	}
	
	Calibrations[0]->Fit(tflin,"q","q",216.5,340);
	Calibrations[2]->Fit(tflin,"q","q",216.5,340);
	Calibrations[1]->Fit(tflin,"q","q",216.5,296);
// 	Calibrations[1]->Fit(tflin,"q+","q",296,340);
	Calibrations_1->Fit(tflin,"q","q",296,340);
	Calibrations[3]->Fit(tflin,"q","q",296,340);
	
	outroot->cd();
	
	TCanvas* cc=new TCanvas("cc","cc",600,600);
	axisSC->Draw("AP");
	axisSC->GetXaxis()->SetTitle("RunNo");axisSC->GetXaxis()->CenterTitle();
	axisSC->GetYaxis()->SetTitle("Mean Charge/photoelectron (x20 fC)");axisSC->GetYaxis()->CenterTitle();
	Calibrations[0]->Draw("same P");
	Calibrations[1]->Draw("same P");
	Calibrations_1->Draw("same P");
	Calibrations[2]->Draw("same P");
	Calibrations[3]->Draw("same P");
	cc->SetName("Calibrations");cc->SetTitle("Calibrations");
	cc->Write();
	
	for(int i1=0;i1<4;i1++)
	{
		Calibrations[i1]->Write();
	}
	Calibrations_1->Write();
	
	outroot->Close();
	
// 	sprintf(hname,"cp PMTCalibration3_%d.root %s/Histos/PMTCalibration3_%d.root",RunNo,AnalysisFilePath,RunNo);system(hname);
}


