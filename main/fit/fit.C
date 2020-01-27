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

string PMTNames[8]={"LAr_noTPB","LAr_wTPB","LAr_wTPB_new","LXe","SC1","SiPM_Quartz","SiPM_TPB","S13371_Window"};

int main( int argc, const char* argv[] )
{
	//RunNo=atoi(argv[1]);
	
 	sprintf(hname,"NitroXenFits_3-component.root");
	//sprintf(hname,"NitroXenFits_2-component.root");
	TFile* outroot=new TFile(hname,"recreate");
	TFile* inroot;
	
 	ofstream outfile("FitResults_3component.txt");
	//stream outfile("FitResults_2component.txt");
	
	int colors[8]={3,1,4,2,1,2,3,4};
	
	TH1D* hh[8];
	TH1F* pe;
	TH1F* integ;
	TGraphErrors* SlowComponents[8];
	TGraphErrors* Integrals[8];
	TGraphErrors* PEs[8];
	TGraphErrors* PEsFromIntegrals[8];
	for(int i1=0;i1<8;i1++)
	{
		sprintf(hname,"AllUnsatAvgWF_%s",PMTNames[i1].c_str());
		hh[i1]=new TH1D(hname,hname,15000,-0.5,14999.5);
		hh[i1]->GetXaxis()->SetTitle("Time (ns)");hh[i1]->GetXaxis()->CenterTitle();
		hh[i1]->GetYaxis()->SetTitle("Mean ADC Counts");hh[i1]->GetYaxis()->CenterTitle();
		
		sprintf(hname,"SlowComponent_%s",PMTNames[i1].c_str());
		SlowComponents[i1]=new TGraphErrors();
		SlowComponents[i1]->SetMarkerStyle(20+i1);SlowComponents[i1]->SetMarkerColor(colors[i1]);
		sprintf(hname,"Integral_%s",PMTNames[i1].c_str());
		Integrals[i1]=new TGraphErrors();
		Integrals[i1]->SetMarkerStyle(20+i1);Integrals[i1]->SetMarkerColor(colors[i1]);
		
		sprintf(hname,"PE_%s",PMTNames[i1].c_str());
		PEs[i1]=new TGraphErrors();
		PEs[i1]->SetMarkerStyle(20+i1);PEs[i1]->SetMarkerColor(colors[i1]);
		sprintf(hname,"PEFromIntegral_%s",PMTNames[i1].c_str());
		PEsFromIntegrals[i1]=new TGraphErrors();
		PEsFromIntegrals[i1]->SetMarkerStyle(20+i1);PEsFromIntegrals[i1]->SetMarkerColor(colors[i1]);
	}
	int np=0;
	
	TF1* tf1=new TF1("tf1","[0]*exp((x-[1])/[2])/([5]*exp((x-[1])/[3])+[8]*exp((x-[1])/[4]))",0.,15000.);
	TF1* tf2=new TF1("tf2","[0]*pow(x-[1],[2])*([3]*exp(-(x-[1])/[4])+[5]*exp(-(x-[1])/[8]))",0.,15000.);
	TF1* tf3=new TF1("tf3","[0]*exp(-(x-[1])/[2])+[3]*exp(-(x-[1])/[4])+[5]*exp(-(x-[1])/[8])",0.,15000.);
	TF1* tf4=new TF1("tf4","[0]*exp(-(x-[1])/[2])+[3]*exp(-(x-[1])/[4])",0.,15000.);
	TF1* tfl=new TF1("tfl","landau",0.,30000.);
	
	//int NRuns=34;
	//int RunList[34]={212,215,219,221,226,229,234,238,243,247,250,254,258,262,268,270,275,278,281,282,284,288,291,293,299,304,308,317,319,323,325,327,332,336};
	
	//int NRunsN2=15;
	//int RunListN2[15]={212,247,254,258,262,268,270,275,278,282,284,288,291,293,299};
	//int NRunsN2=13;
	//int RunListN2[13]={247,254,258,262,268,270,275,278,282,284,288,291,293};
	int NRunsN2=1;
	int RunListN2[1]={421};
	

	int NRunsXe=8;
	int RunListXe[8]={421,426,435,443,449,454,464,470};
	
	TGraph* axisSC=new TGraph();
	axisSC->SetPoint(0,RunListN2[0]-1,0);
	axisSC->SetPoint(1,RunListN2[NRunsN2-1]+1,2000);
	TGraph* axisI=new TGraph();
	axisI->SetPoint(0,RunListXe[0]-1,0);
	axisI->SetPoint(1,RunListXe[NRunsXe-1]+1,30000);
	TGraph* axisPE=new TGraph();
	axisPE->SetPoint(0,RunListXe[0]-1,0);
	axisPE->SetPoint(1,RunListXe[NRunsXe-1]+1,300);
	
	int binmax=0;int bin5percent=0;
	float valmax=0;float val5percent=0;
	float xmax=0;float x5percent=0;
	
	vector <float> nchi2;
	vector <float> ulim;
	float ul=0;
	int imin=0;
	float err=0;
	
	cout<<"slow components"<<endl;
	for(int i1=0;i1<NRunsN2;i1++)
	{
		RunNo=RunListN2[i1];
		cout<<RunNo<<endl;
		sprintf(hname,"../NitroXen_%d.root",RunNo);
		inroot=new TFile(hname);
		outfile<<RunNo;
		for(int i2=0;i2<8;i2++)
		{
			sprintf(hname,"AllUnsatAvgWF_%s",PMTNames[i2].c_str());
			//sprintf(hname,"AvgWF_%s",PMTNames[i2].c_str());//default after event selection
			//sprintf(hname,"AvgWFShifted_%s",PMTNames[i2].c_str());
			inroot->GetObject(hname,hh[i2]);
			
			// // //3-component fit
 			tf3->SetParameter(0,50);
 			tf3->SetParameter(1,2000);
 			tf3->SetParameter(2,6);
 			tf3->SetParameter(4,90);
 			tf3->SetParameter(6,1600);
 			
 			tf3->SetParLimits(0,0,100);
 			tf3->SetParLimits(1,1800,2100);
 			tf3->SetParLimits(2,2,20);
 			tf3->SetParLimits(3,0,100);
			//tf3->SetParLimits(4,10,400);
 			tf3->SetParLimits(4,10,800);
 			tf3->SetParLimits(5,0,100);
 			tf3->SetParLimits(6,500,2000);
			
			binmax=hh[i2]->GetMaximumBin();
			valmax=hh[i2]->GetBinContent(binmax);
			xmax=hh[i2]->GetBinCenter(binmax);
			
 			for(int is1=binmax;is1<=hh[i2]->GetNbinsX();is1++)
 			{
 				if((hh[i2]->GetBinContent(is1)/valmax)<0.05)
 				{
 					bin5percent=is1;
 					val5percent=hh[i2]->GetBinContent(is1);
 					x5percent=hh[i2]->GetBinCenter(is1);
 					break;
 				}
 			}
 			hh[i2]->Fit(tf3,"q","q",xmax,x5percent);
			
 			if(i2==0) hh[i2]->Fit(tf3,"q","q",hh[i2]->GetBinCenter(hh[i2]->GetMaximumBin()),4500);
 			else hh[i2]->Fit(tf3,"q","q",hh[i2]->GetBinCenter(hh[i2]->GetMaximumBin()),3500);
			
			nchi2.clear();ulim.clear();

			
			np=SlowComponents[i2]->GetN();
			SlowComponents[i2]->SetPoint(np,RunNo,tf3->GetParameter(6));
			SlowComponents[i2]->SetPointError(np,0,tf3->GetParError(6));
			
			outroot->cd();
			
			sprintf(hname2,"%s_%d",hh[i2]->GetName(),RunNo);
			hh[i2]->SetName(hname2);hh[i2]->SetTitle(hname2);
			hh[i2]->GetXaxis()->SetTitle("Time (ns)");hh[i2]->GetXaxis()->CenterTitle();
			hh[i2]->GetYaxis()->SetTitle("Mean ADC Counts");hh[i2]->GetYaxis()->CenterTitle();
			hh[i2]->SetLineColor(colors[i2]);
			hh[i2]->Write();
			
			if(i2==0 || i2==3)
			{
// 				if(i2>=0 && i2<=3) hh[i2]->Write();
				outfile<<" "<<i2<<" "<<tf3->GetParameter(0)<<" "<<tf3->GetParameter(1)<<" "<<tf3->GetParameter(2)<<" "<<tf3->GetParameter(3)<<" "<<tf3->GetParameter(4)<<" "<<tf3->GetParameter(5)<<" "<<tf3->GetParameter(6)<<" "<<tf3->GetChisquare()/tf3->GetNDF();//3-component
// 				outfile<<" "<<i2<<" "<<tf4->GetParameter(0)<<" "<<tf4->GetParameter(1)<<" "<<tf4->GetParameter(2)<<" "<<tf4->GetParameter(3)<<" "<<tf4->GetParameter(4)<<" "<<tf4->GetChisquare()/tf4->GetNDF();//2-component
			}
		}
		outfile<<endl;
		inroot->Close();
	}
	
// 	float calibs[2][4]={{36.04,38.41,23.59,12.54},{36.04,88.21,23.59,21.37}};
// 	float calibs[2][4]={{35.01,36.16,23.92,16.92},{35.01,89.29,23.92,20.27}};
// 	float caliberrs[2][4]={{0.02,0.06,0.02,0.06},{0.02,0.08,0.02,0.02}};
	//float calibs[2][8]={{88.9,31.6,18.7,22.7,1,1,1,1},{88.9,31.6,18.7,22.7,1,1,1,1}};
	float calibs[2][8]={{85.6,31.1,5.1,21.5,1,1,1,1},{85.6,31.1,5.1,21.5,1,1,1,1}};
	float caliberrs[2][8]={{0.06,0.02,0.02,0.02,0.01,0.01,0.01,0.01},{0.06,0.02,0.02,0.02,0.01,0.01,0.01,0.01}};
	int calibID=0;
// 	if(RunNo>=299) calibID=1;
// 	else calibID=0;
	int binoffset[8]={2,0,21,3,0,0,0,0};
	
	cout<<"integrals"<<endl;
	for(int i1=0;i1<NRunsXe;i1++)
	{
		RunNo=RunListXe[i1];
		cout<<RunNo<<endl;
		sprintf(hname,"../NitroXen_%d.root",RunNo);
		inroot=new TFile(hname);
		for(int i2=0;i2<8;i2++)
		{
			sprintf(hname,"AllUnsatAvgWF_%s",PMTNames[i2].c_str());
			//sprintf(hname,"AvgWF_%s",PMTNames[i2].c_str());//default after event selection
			//sprintf(hname,"AvgWFShifted_%s",PMTNames[i2].c_str());
			inroot->GetObject(hname,hh[i2]);
			//sprintf(hname,"PE_1500_15000_%s",PMTNames[i2].c_str());
			//sprintf(hname,"PE_1500_8000_%s",PMTNames[i2].c_str());
			sprintf(hname,"PE_1800_6000_%s",PMTNames[i2].c_str());
			inroot->GetObject(hname,pe);
			
			//sprintf(hname,"Integral_1500_15000_%s",PMTNames[i2].c_str());
			//sprintf(hname,"Integral_1500_8000_%s",PMTNames[i2].c_str());
			sprintf(hname,"Integral_1800_6000_%s",PMTNames[i2].c_str());
			inroot->GetObject(hname,integ);
			
			pe->Fit(tfl,"q","q",0.,300.);
			
			np=Integrals[i2]->GetN();
			Integrals[i2]->SetPoint(np,RunNo,hh[i2]->Integral());
			Integrals[i2]->SetPointError(np,0,0);
			np=PEs[i2]->GetN();
			PEs[i2]->SetPoint(np,RunNo,tfl->GetParameter(1));
			PEs[i2]->SetPointError(np,0,tfl->GetParError(1));
			
			integ->Fit(tfl,"q","q",0.,10000.);
			integ->Fit(tfl,"q","q",0.,tfl->GetParameter(1)+3*tfl->GetParameter(2));
			
			//special fits
			if(RunNo==270 && i2==3) integ->Fit(tfl,"q","q",0.,1000.);
			if(RunNo==304 && i2==2) integ->Fit(tfl,"q","q",0.,4000.);
			//if(RunNo==323 && i2==2) integ->Fit(tfl,"q","q",1500,5000.);
			if(RunNo==325 && i2==2) integ->Fit(tfl,"q","q",0,8000.);
			if(RunNo==327 && i2==2) integ->Fit(tfl,"q","q",0,8000.);
			if(RunNo==332 && i2==2) integ->Fit(tfl,"q","q",0,8000.);
			//if(RunNo==323 && i2==0) integ->Fit(tfl,"q","q",800,5000.);
			
			
				
			
			//np=PEsFromIntegrals[i2]->GetN();
			//PEsFromIntegrals[i2]->SetPoint(np,RunNo,hh[i2]->Integral()/calibs[calibID][i2]);
			//PEsFromIntegrals[i2]->SetPointError(np,0,0);
			
			np=PEsFromIntegrals[i2]->GetN();
			PEsFromIntegrals[i2]->SetPoint(np,RunNo,tfl->GetParameter(1)/calibs[calibID][i2]);
			PEsFromIntegrals[i2]->SetPointError(np,0,(tfl->GetParameter(1)/calibs[calibID][i2])*sqrt(pow(caliberrs[calibID][i2]/calibs[calibID][i2],2)+pow(tfl->GetParError(1)/tfl->GetParameter(1),2)));
			
			outroot->cd();
			//sprintf(hname2,"%s_%d",hh[i2]->GetName(),RunNo);
			//hh[i2]->SetName(hname2);hh[i2]->SetTitle(hname2);
			//hh[i2]->GetXaxis()->SetTitle("Time (ns)");hh[i2]->GetXaxis()->CenterTitle();
			//hh[i2]->GetYaxis()->SetTitle("Mean ADC Counts");hh[i2]->GetYaxis()->CenterTitle();
			//hh[i2]->SetLineColor(colors[i2]);
			//hh[i2]->Write();
			
			for(int is1=1;is1<=hh[i2]->GetNbinsX()-binoffset[i2];is1++)
			{
				hh[i2]->SetBinContent(is1,hh[i2]->GetBinContent(is1+binoffset[i2])/calibs[calibID][i2]);
			}
			//hh[i2]->Scale(hh[i2]->Integral()/calibs[calibID][i2]);
			sprintf(hname2,"%s_%d_PE",hh[i2]->GetName(),RunNo);
			hh[i2]->SetName(hname2);hh[i2]->SetTitle(hname2);
			hh[i2]->GetXaxis()->SetTitle("Time (ns)");hh[i2]->GetXaxis()->CenterTitle();
			hh[i2]->GetYaxis()->SetTitle("Average Signal");hh[i2]->GetYaxis()->CenterTitle();
			hh[i2]->SetLineColor(colors[i2]);
			hh[i2]->Write();
			
			sprintf(hname2,"%s_%d",pe->GetName(),RunNo);
			pe->SetName(hname2);pe->SetTitle(hname2);
			pe->GetXaxis()->SetTitle("Number of Photoelectrons");pe->GetXaxis()->CenterTitle();
			pe->GetYaxis()->SetTitle("Events / 2.5 photoelectrons");pe->GetYaxis()->CenterTitle();
			pe->SetLineColor(colors[i2]);
			pe->Write();
			
			sprintf(hname2,"%s_%d",integ->GetName(),RunNo);
			integ->SetName(hname2);integ->SetTitle(hname2);
			integ->GetXaxis()->SetTitle("Charge ");integ->GetXaxis()->CenterTitle();
			integ->GetYaxis()->SetTitle("Events / 200");integ->GetYaxis()->CenterTitle();
			integ->SetLineColor(colors[i2]);
			integ->Write();
			
		}
		inroot->Close();
	}
	
	TCanvas* cc=new TCanvas("cc","cc",600,600);
	axisSC->Draw("AP");
	axisSC->GetXaxis()->SetTitle("RunNo");axisSC->GetXaxis()->CenterTitle();
	axisSC->GetYaxis()->SetTitle("Slow Component (ns)");axisSC->GetYaxis()->CenterTitle();
	SlowComponents[0]->Draw("same P");
	SlowComponents[3]->Draw("same P");
	cc->SetName("SlowComponents");cc->SetTitle("SlowComponents");
	cc->Write();
	
	axisI->Draw("AP");
	axisI->GetXaxis()->SetTitle("RunNo");axisI->GetXaxis()->CenterTitle();
	axisI->GetYaxis()->SetTitle("Integral (au)");axisI->GetYaxis()->CenterTitle();
	Integrals[0]->Draw("same P");
	Integrals[1]->Draw("same P");
	Integrals[2]->Draw("same P");
	Integrals[3]->Draw("same P");
	cc->SetName("Integrals");cc->SetTitle("Integrals");
	cc->Write();
	
	axisPE->Draw("AP");
	axisPE->GetXaxis()->SetTitle("RunNo");axisPE->GetXaxis()->CenterTitle();
	axisPE->GetYaxis()->SetTitle("Number of Photoelectrons");axisPE->GetYaxis()->CenterTitle();
	PEs[0]->Draw("same P");
	PEs[1]->Draw("same P");
	PEs[2]->Draw("same P");
	PEs[3]->Draw("same P");
	cc->SetName("Photoelectrons");cc->SetTitle("Photoelectrons");
	cc->Write();
	
	axisPE->Draw("AP");
	axisPE->GetXaxis()->SetTitle("RunNo");axisPE->GetXaxis()->CenterTitle();
	axisPE->GetYaxis()->SetTitle("Number of Photoelectrons");axisPE->GetYaxis()->CenterTitle();
	PEsFromIntegrals[0]->Draw("same P");
	PEsFromIntegrals[1]->Draw("same P");
	PEsFromIntegrals[2]->Draw("same P");
	PEsFromIntegrals[3]->Draw("same P");
	cc->SetName("PhotoelectronsFromIntegrals");cc->SetTitle("PhotoelectronsFromIntegrals");
	cc->Write();
	
	outroot->Close();
	
		// 	sprintf(hname,"cp PMTCalibration3_%d.root %s/Histos/PMTCalibration3_%d.root",RunNo,AnalysisFilePath,RunNo);system(hname);
}


