Double_t fitf(Double_t *x,Double_t *par)
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

Double_t fitf2(Double_t *x,Double_t *par)
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
		sigma=par[2];
	}
	fitval=par[0]*exp(-0.5*pow(tt/sigma,2));
	return fitval;
}

void line(double t, const double *p, double &x, double &y, double &z)
{
	x = p[0] + p[1]*(t-p[4]);
	y = p[2] + p[3]*(t-p[4]);
	z = t;
}

struct SumDistance2
{
	TGraph2D *fGraph;
	SumDistance2(TGraph2D *g) : fGraph(g) {}
	double distance2(double x,double y,double z, const double *p)
	{
		XYZVector xp(x,y,z);
// 		XYZVector x0(p[0], p[2], 0. );
		XYZVector x0(p[0]-p[1]*p[4], p[2]-p[3]*p[4], 0. );
// 		XYZVector x1(p[0] + p[1], p[2] + p[3], 1. );
		XYZVector x1(p[0] + p[1]*(z-p[4]), p[2] + p[3]*(z-p[4]), z );
		XYZVector u = (x1-x0).Unit();
		double d2 = ((xp-x0).Cross(u)).Mag2();
		return d2;
	}
	double operator() (const double *par)
	{
		assert(fGraph != 0);
		double * x = fGraph->GetX();
		double * y = fGraph->GetY();
		double * z = fGraph->GetZ();
		int npoints = fGraph->GetN();
		double sum = 0;
		for (int i  = 0; i < npoints; ++i)
		{
			double d = distance2(x[i],y[i],z[i],par);
			sum += d;
		}
		return sum;
	}
};

float line3Dfit(TGraph2D * gr)
{
	ROOT::Fit::Fitter  fitter;
	SumDistance2 sdist(gr);
// 	ROOT::Math::Functor fcn(sdist,4);
	ROOT::Math::Functor fcn(sdist,5);
	fitter.SetFCN(fcn,pStart);
	for (int i = 0; i < 5; ++i) fitter.Config().ParSettings(i).SetStepSize(0.01);
	
// 	fitter.Config().ParSettings(0).SetLimits(-50,50);
// 	fitter.Config().ParSettings(1).SetLimits(-5,5);
// 	fitter.Config().ParSettings(2).SetLimits(-50,50);
// 	fitter.Config().ParSettings(3).SetLimits(-5,5);
// 	fitter.Config().ParSettings(4).SetLimits(0,4100);
	
// 	fitter.Config().ParSettings(0).SetStepSize(0.1);
// 	fitter.Config().ParSettings(1).SetStepSize(0.01);
// 	fitter.Config().ParSettings(2).SetStepSize(0.1);
// 	fitter.Config().ParSettings(3).SetStepSize(0.01);
// 	fitter.Config().ParSettings(4).SetStepSize(1);
	
// 	fitter.Config().SetMinimizer("Minuit","Migrad");
	
	bool ok = fitter.FitFCN();
// 	if (!ok)
// 	{
// // 		cout<<"Line3D Fit failed"<<endl;
// 		return -1;
// 	}
	const ROOT::Fit::FitResult & result = fitter.Result();
	
// 	parFit=result.GetParams();
	const double * parFit = result.GetParams();
	for(int i1=0;i1<5;i1++)
	{
		fitParams[i1]=parFit[i1];
	}
	return sqrt(result.MinFcnValue()/gr->GetN());
}

int FindBinAbove(TH1F* h,float th,int sb)
{
	for(int ip1=sb;ip1<=h->GetNbinsX();ip1++)
	{
		if(h->GetBinContent(ip1)>th){return ip1;}
	}
	return -1;
}

void GetBaseline()
{
	sprintf(hname,"%s/Run_%d.root",NTupleFilePath,RunNo);
	TFile* inroot=new TFile(hname);
	TTree* T =  (TTree*) inroot->Get("T");
	T->SetBranchAddress("TPCWF",&ed.TPCWF);
	T->SetBranchAddress("PMTWF",&ed.PMTWF);
	
	sprintf(hname,"Baselines_%d.root",RunNo);
	TFile* outroot=new TFile(hname,"recreate");
	
	TH1F* AllNormChi2=new TH1F("AllNormChi2","AllNormChi2",400,0.,2.);
	TH1D* AllAmp[3];
	AllAmp[0]=new TH1D("AllCollAmp","AllCollAmp",400,-100,100);
	AllAmp[1]=new TH1D("AllIndAmp","AllIndAmp",400,-100,100);
	AllAmp[2]=new TH1D("AllAmp","AllAmp",400,-100,100);
	
	TH1D* AllPMTAmp[3];
	for(int i1=0;i1<3;i1++)
	{
		sprintf(hname,"AllPMTAmp_%d",i1);
		AllPMTAmp[i1]=new TH1D(hname,hname,200,-100,100);
	}
	
	TGraphErrors* tg[2];
	for(int i1=0;i1<2;i1++)
	{
		tg[i1]=new TGraphErrors();
		sprintf(hname,"Baseline_%s",wt[i1].c_str());
		tg[i1]->SetName(hname);
		tg[i1]->SetTitle(hname);
		tg[i1]->SetMarkerStyle(24+i1);
		tg[i1]->SetMarkerColor(1+i1);
	}
	TGraphErrors* tgPMT;
	tgPMT=new TGraphErrors();
	sprintf(hname,"PMTBaseline");
	tgPMT->SetName(hname);
	tgPMT->SetTitle(hname);
	tgPMT->SetMarkerStyle(20);
	tgPMT->SetMarkerColor(1);
	
	TH1F* HB[2][128];
	for(int i1=0;i1<2;i1++)
	{
		for(int i2=0;i2<128;i2++)
		{
			sprintf(hname,"Amp_%s_%d",wt[i1].c_str(),i2);
			HB[i1][i2]=new TH1F(hname,hname,1000,-0.5,999.5);
		}
	}
	TH1F* PMTHB[3];
	for(int i1=0;i1<3;i1++)
	{
		sprintf(hname,"PMT_Amp_%d",i1);
		PMTHB[i1]=new TH1F(hname,hname,1500,-0.5,1499.5);
	}
	
	int ic=0;int gch=0;
	for(int I=0;I<T->GetEntries();I++)
	{
		T->GetEntry(I);
		for(int i1=0;i1<256;i1++)
		{
			ic=i1/128;gch=i1-((i1/128)*128);
			for(int i2=0;i2<4096;i2++)
			{
				HB[ic][gch]->Fill(ed.TPCWF->at(i1)[i2]);
			}
		}
		
		for(int i1=0;i1<3;i1++)
		{
			for(int i2=0;i2<ed.PMTWF->at(i1).size();i2++)
			{
				PMTHB[i1]->Fill(ed.PMTWF->at(i1)[i2]);
			}
		}
// 		if(I%10==0) cout<<"Event: "<<I<<endl;
	}
	
	float res[2][128][2]={{{0.}}};//mean sigma
	float resPMT[3][2]={{0.}};//mean sigma
	TF1* g=new TF1("gaus","gaus",-100.,3000.);
	for(int i1=0;i1<2;i1++)
	{
		for(int i2=0;i2<128;i2++)
		{
			g->SetParameter(0,HB[i1][i2]->GetMaximum());
			g->SetParameter(1,HB[i1][i2]->GetMean());
			g->SetParameter(2,HB[i1][i2]->GetRMS());
			HB[i1][i2]->Fit(g,"q","q",10.,1000.);
			res[i1][i2][0]=g->GetParameter(1);
			res[i1][i2][1]=g->GetParameter(2);
		}
	}
	for(int i1=0;i1<3;i1++)
	{
		g->SetParameter(0,PMTHB[i1]->GetMaximum());
		g->SetParameter(1,PMTHB[i1]->GetMean());
		g->SetParameter(2,PMTHB[i1]->GetRMS());
		PMTHB[i1]->Fit(g,"q","q",10.,2000.);
		resPMT[i1][0]=g->GetParameter(1);
		resPMT[i1][1]=g->GetParameter(2);
	}
	
	for(int I=0;I<T->GetEntries();I++)
	{
		T->GetEntry(I);
		for(int i1=0;i1<256;i1++)
		{
			ic=i1/128;gch=i1-((i1/128)*128);
			for(int i2=0;i2<4096;i2++)
			{
				AllAmp[ic]->Fill(ed.TPCWF->at(i1)[i2]-res[ic][gch][0]);
				AllAmp[2]->Fill(ed.TPCWF->at(i1)[i2]-res[ic][gch][0]);
			}
		}
		for(int i1=0;i1<3;i1++)
		{
			for(int i2=0;i2<ed.PMTWF->at(i1).size();i2++)
			{
				AllPMTAmp[i1]->Fill(ed.PMTWF->at(i1)[i2]-resPMT[i1][0]);
			}
		}
	}
	
	sprintf(hname,"TPCBaselines_%d.txt",RunNo);
	ofstream outfile(hname);
	
	for(int i1=0;i1<2;i1++)
	{
		for(int i2=0;i2<128;i2++)
		{
			outfile<<i1<<" "<<i2<<" "<<res[i1][i2][0]<<" "<<res[i1][i2][1]<<endl;//mean sigma
			tg[i1]->SetPoint(i2,i2,res[i1][i2][0]);
			tg[i1]->SetPointError(i2,0,res[i1][i2][1]);
		}
	}
	outfile.close();
	
	sprintf(hname,"PMTBaselines_%d.txt",RunNo);
	ofstream PMToutfile(hname);
	
	for(int i1=0;i1<3;i1++)
	{
		PMToutfile<<i1<<" "<<resPMT[i1][0]<<" "<<resPMT[i1][1]<<endl;
		tgPMT->SetPoint(i1,i1,resPMT[i1][0]);
		tgPMT->SetPointError(i1,0,resPMT[i1][1]);
	}
	
	outroot->cd();
	tgPMT->Write();
	tg[0]->Write();
	tg[1]->Write();
	for(int i1=0;i1<3;i1++)
	{
		g->SetParameter(0,AllAmp[i1]->GetMaximum());
		g->SetParameter(1,AllAmp[i1]->GetMean());
		g->SetParameter(2,AllAmp[i1]->GetRMS());
		AllAmp[i1]->Fit(g,"q","q",-50.,50.);
		AllAmp[i1]->Write();
		
		g->SetParameter(0,AllPMTAmp[i1]->GetMaximum());
		g->SetParameter(1,AllPMTAmp[i1]->GetMean());
		g->SetParameter(2,AllPMTAmp[i1]->GetRMS());
		AllPMTAmp[i1]->Fit(g,"q","q",-50.,50.);
		AllPMTAmp[i1]->Write();
	}
	for(int i1=0;i1<3;i1++)
	{
		PMTHB[i1]->Write();
	}
	
	for(int i1=0;i1<2;i1++)
	{
		for(int i2=0;i2<128;i2++)
		{
			HB[i1][i2]->Write();
		}
	}
	
	outroot->Close();
	inroot->Close();
	outfile.close();
	PMToutfile.close();
	
	sprintf(hname,"cp Baselines_%d.root %s/Histos/Baselines_%d.root;wait;",RunNo,AnalysisFilePath,RunNo);system(hname);
	sprintf(hname,"cp TPCBaselines_%d.txt %s/Files/TPCBaselines_%d.txt;wait;",RunNo,AnalysisFilePath,RunNo);system(hname);
	sprintf(hname,"cp PMTBaselines_%d.txt %s/Files/PMTBaselines_%d.txt",RunNo,AnalysisFilePath,RunNo);system(hname);
}

void ReadTPCBaselines()
{
	sprintf(hname,"cp %s/Files/TPCBaselines_%d.txt .;wait;",AnalysisFilePath,TPCBaselineRun);system(hname);
	sprintf(hname,"TPCBaselines_%d.txt",TPCBaselineRun);
	ifstream inTPCBaselines(hname);
	int a[2]={0};float b=0;float c=0;
	while(!inTPCBaselines.eof())
	{
		inTPCBaselines>>a[0]>>a[1]>>b>>c;
		TPCBaselines[a[0]][a[1]][0]=b;
		TPCBaselines[a[0]][a[1]][1]=c;
	}
	inTPCBaselines.close();
	sprintf(hname,"rm TPCBaselines_%d.txt",TPCBaselineRun);system(hname);
	
	sprintf(hname,"%s/Histos/Baselines_%d.root",AnalysisFilePath,TPCBaselineRun);
	TFile* inroot=new TFile(hname);
	TF1* g;TH1D* h1;
	inroot->GetObject("AllCollAmp",h1);g=(TF1*)h1->GetFunction("gaus");
	inroot->GetObject("AllIndAmp",h1);g=(TF1*)h1->GetFunction("gaus");
// 	ThColInd[0]=g->GetParameter(2)*4.;//3.5 sigma
// 	ThColInd[1]=g->GetParameter(2)*4.;//3.5 sigma
	
	ThColInd[0]=g->GetParameter(2)*ThSigma[0];
	ThColInd[1]=g->GetParameter(2)*ThSigma[1];
	
	inroot->Close();
	cout<<"Collection Threshold = "<<ThColInd[0]<<" Induction Threshold = "<<ThColInd[1]<<endl;
}

void ReadPMTBaselines()
{
	sprintf(hname,"cp %s/Files/PMTBaselines_%d.txt .;wait;",AnalysisFilePath,PMTBaselineRun);system(hname);
	sprintf(hname,"PMTBaselines_%d.txt",PMTBaselineRun);
	ifstream inPMTBaselines(hname);
	int a[2]={0};float c(0);float b(0);
	while(!inPMTBaselines.eof())
	{
		inPMTBaselines>>a[0]>>b>>c;
		PMTBaselines[a[0]][0]=b;
		PMTBaselines[a[0]][1]=c;
	}
	inPMTBaselines.close();
	sprintf(hname,"rm PMTBaselines_%d.txt",PMTBaselineRun);system(hname);
	
	sprintf(hname,"%s/Histos/Baselines_%d.root",AnalysisFilePath,PMTBaselineRun);
	TFile* inroot=new TFile(hname);
	TF1* g;TH1D* h1;
	for(int i1=0;i1<3;i1++)
	{
		sprintf(hname,"AllPMTAmp_%d",i1);
		inroot->GetObject(hname,h1);g=(TF1*)h1->GetFunction("gaus");ThPMT[i1]=g->GetParameter(2)*3.;//3 sigma
		cout<<"PMT "<<i1<<" Threshold = "<<ThPMT[i1]<<endl;
	}
	inroot->Close();
}

int FindBaselineRun(int r)
{
	ifstream inRunList("RunList.txt");
	int a[2]={0};
	while(!inRunList.eof())
	{
		inRunList>>a[0]>>a[1];
		if(a[0]==r) return a[1];
	}
	inRunList.close();
	return -1;
}

void GetPMTCalibration()
{
// 	ReadPMTBaselines();
// 	sprintf(hname,"%s/Histos/PMTCalibration_%d.root",AnalysisFilePath,RunNo);
// 	TFile* outroot=new TFile(hname,"recreate");
// 	
// 	TH1F* hh[4][3];
// 	
// 	for(int i1=0;i1<3;i1++)
// 	{
// 		sprintf(hname,"BL_amp_%d",i1);
// 		hh[0][i1]=new TH1F(hname,hname,1500,-1000,2000);
// 		sprintf(hname,"BL_amp2_%d",i1);
// 		hh[1][i1]=new TH1F(hname,hname,1500,-1000,2000);
// 		sprintf(hname,"BL_amp3_%d",i1);
// 		hh[2][i1]=new TH1F(hname,hname,1500,-1000,2000);
// 		sprintf(hname,"BL_amp4_%d",i1);
// 		hh[3][i1]=new TH1F(hname,hname,1500,-1000,2000);
// 	}
// 	
// 	TF1* lin=new TF1("lin","[0]",0.,4095.);
// 	TH1F* hb=new TH1F("hb","hb",15000,-0.5,14999.5);
// // 	TF1* tf=new TF1("GL","gaus(0)+landau(3)",-100.,200.);
// // 	TF1* tf=new TF1("GLG","gaus(0)+landau(3)+gaus(6)",-100.,200.);
// 	TF1* tf=new TF1("GGG","gaus(0)+gaus(3)+gaus(6)",-100.,200.);
// 	
// 	char cNum[10];
// 	int ii=0;int Ev=0;
// 	ifstream pmtinfile;
// 	for(int i1=0;i1<3;i1++)
// 	{
// 		Ev=0;
// 		sprintf(hname,"%s/waveforms_ch%d_run%d.txt",PMTinFilePath,i1,RunNo);
// 		cout<<hname<<endl;
// 		pmtinfile.open(hname);
// 		string strarr;
// 		while(getline(pmtinfile,strarr,'\n'))
// 		{
// 			istringstream iss(strarr);
// 			ii=0;
// // 			float qint=0;
// 			float qint[4]={0};
// 			for(std::string s; iss >> s; )
// 			{
// 				if(ii>0)
// 				{
// 					hb->SetBinContent(ii,atof(s.c_str())-PMTBaselines[i1][0]);
// // 					qint+=((-1.)*(atof(s.c_str())-PMTBaselines[i1][0]));
// 					for(int is1=0;is1<4;is1++)
// 					{
// 						qint[is1]+=((-1.)*(atof(s.c_str())-PMTBaselines[i1][0]));
// 					}
// 				}
// 				ii++;
// 				if(ii%20==0){hh[0][i1]->Fill(qint[0]);qint[0]=0;}
// 				if(ii%10==0){hh[1][i1]->Fill(qint[1]);qint[1]=0;}
// 				if(ii%100==0){hh[2][i1]->Fill(qint[2]);qint[2]=0;}
// 				hh[3][i1]->Fill(qint[3]);qint[3]=0;
// // 				hh[i1]->Fill(qint);qint=0;
// 			}
// // 			hh[i1]->Fill(qint);
// 			Ev++;
// 		}
// 		pmtinfile.close();
// 	}
// 	
// 	for(int i1=0;i1<3;i1++)
// 	{
// 		for(int i2=0;i2<4;i2++)
// 		{
// 			tf->SetParameter(0,hh[i2][i1]->GetBinContent(hh[i2][i1]->GetMaximumBin()));
// 			tf->SetParameter(1,0);
// 			tf->SetParameter(2,4);
// 			tf->SetParameter(3,hh[i2][i1]->GetBinContent(hh[i2][i1]->GetMaximumBin())/1000);
// 			tf->SetParameter(4,30);
// 			tf->SetParameter(5,10);
// 			tf->SetParameter(6,hh[i2][i1]->GetBinContent(hh[i2][i1]->GetMaximumBin())/1e6);
// 			tf->SetParameter(7,100);
// 			tf->SetParameter(8,50);
// 			
// 			hh[i2][i1]->Fit(tf,"q","q",-100.,200.);
// 		}
// 		
// // 		cout<<"PMT "<<i1<<" "<<(tf->GetParameter(4)-tf->GetParameter(1))<<" ("<<sqrt(pow(tf->GetParError(4),2)+pow(tf->GetParError(1),2))<<")"<<endl;
// 	}
// 	
// 	
// // 	sprintf(hname,"Files/Calibrations_%d.txt",RunNo);
// // 	ofstream outfile(hname);
// 	
// 	outroot->cd();
// 	for(int i1=0;i1<3;i1++)
// 	{
// 		for(int i2=0;i2<4;i2++)
// 		{
// // 		outfile<<i1<<" "<<hh[i1][1]->GetBinCenter(hh[i1][1]->GetMaximumBin())<<endl;
// 			hh[i2][i1]->Write();
// 		}
// 	}
// // 	outfile.close();
// 	outroot->Close();
}

void GetPMTCalibration2()
{
	sprintf(hname,"PMTCalibration_%d.root",RunNo);
	TFile* outroot=new TFile(hname,"recreate");
	
	sprintf(hname,"%s/Run_%d.root",NTupleFilePath,RunNo);
	TFile* inroot=new TFile(hname);
	TTree* T =  (TTree*) inroot->Get("T");
	T->SetBranchAddress("TPCWF",&ed.TPCWF);
	T->SetBranchAddress("PMTWF",&ed.PMTWF);
	
	TH1D* hh[11][3];TH1F* htr[3];TH1F* hs[3];TH1F* hs2[3];
	TH1F* hfitmean[3];TH1F* hfitrms[3];
	TH1D* Ampvst[3];
	TH1F* BLs[3];
	
	if(RunNo>7000)
	{
		ReadPMTBaselines();
	}
	
	for(int i1=0;i1<3;i1++)
	{
		for(int i2=0;i2<10;i2++)
		{
			sprintf(hname,"Int_%d_%d",(i2+1)*5,i1);
			hh[i2][i1]=new TH1D(hname,hname,3000,-200.,2800.);
// 			hh[i2][i1]=new TH1D(hname,hname,6000,-200.,2800.);
		}
		
		sprintf(hname,"Amp_vs_t_%d",i1);
		Ampvst[i1]=new TH1D(hname,hname,15000,-0.5,14999.5);
		sprintf(hname,"BLs_%d",i1);
		BLs[i1]=new TH1F(hname,hname,2000,0.,1000.);
		
		sprintf(hname,"TransientsWindow_%d",i1);
// 		htr[i1]=new TH1F(hname,hname,2000,-100,900);
		htr[i1]=new TH1F(hname,hname,3000,-200.,2800.);
		sprintf(hname,"SignalWindow_%d",i1);
// 		hs[i1]=new TH1F(hname,hname,2000,-100,900);
		hs[i1]=new TH1F(hname,hname,3000,-200.,2800.);
		sprintf(hname,"SignalWindow2_%d",i1);
		hs2[i1]=new TH1F(hname,hname,3000,-200.,2800.);
		sprintf(hname,"BLfitmean_%d",i1);
		hfitmean[i1]=new TH1F(hname,hname,2000,-100,900);
		sprintf(hname,"BLfitrms_%d",i1);
		hfitrms[i1]=new TH1F(hname,hname,2000,-100,900);
	}
// 	TH2F* Evvst[3];
// 	for(int i1=0;i1<3;i1++)
// 	{
// 		sprintf(hname,"Ev_vs_t_%d",i1);
// 		Evvst[i1]=new TH2F(hname,hname,15000,-0.5,14999.5,1000,-0.5,999.5);
// 	}
	
	TGraph* BLvsE[3];
	for(int i1=0;i1<3;i1++)
	{
		sprintf(hname,"BL_vs_E_%d",i1);
		BLvsE[i1]=new TGraph();
		BLvsE[i1]->SetName(hname);BLvsE[i1]->SetTitle(hname);
		BLvsE[i1]->SetMarkerStyle(24+i1);
		BLvsE[i1]->SetMarkerColor(1+i1);
		BLvsE[i1]->GetXaxis()->SetTitle("Event ID");BLvsE[i1]->GetXaxis()->CenterTitle();
		BLvsE[i1]->GetYaxis()->SetTitle("Mean Amplitude of First 500 Time Bins");BLvsE[i1]->GetYaxis()->CenterTitle();
	}
	
	TF1* g1=new TF1("g1","gaus",-1000.,2000.);
	
	TF1* lin=new TF1("lin","[0]",0.,4095.);
	TH1F* hb=new TH1F("hb","hb",15000,-0.5,14999.5);
// 	TF1* tf=new TF1("GL","gaus(0)+landau(3)",-100.,200.);
// 	TF1* tf=new TF1("GLG","gaus(0)+landau(3)+gaus(6)",-100.,200.);
	TF1* tf=new TF1("GGG","gaus(0)+gaus(3)+gaus(6)",-50.,100.);
	
// 	TF1* tf=new TF1("GG","gaus(0)+gaus(3)",-50.,100.);
// 	TF1* tf=new TF1("GL","gaus(0)+landau(3)",-50.,100.);
// 	TF1* tf=new TF1("GLG","gaus(0)+landau(3)+gaus(6)",-50.,100.);
	TF1* tfg=new TF1("G","gaus",-50.,100.);
	TF1* tfgg=new TF1("GG","gaus(0)+gaus(3)",-50.,100.);
	TF1* tfggg=new TF1("GGG","gaus(0)+gaus(3)+gaus(6)",-50.,100.);
	
	vector <float> Amp;
	float hsum[20]={0.};
	float hped=0;
	float hped2=0;float hsig2=0;
	float s=0;float str=0;float ss=0.;
	
	for(int I=0;I<T->GetEntries();I++)
	{
		T->GetEntry(I);
		for(int i1=0;i1<3;i1++)
		{
			hped=0;
			for(int i2=0;i2<500;i2++){hped+=ed.PMTWF->at(i1)[i2];}
			hped/=500.;
			BLvsE[i1]->SetPoint(BLvsE[i1]->GetN(),I,hped);
			if(I<200) continue;
			
			for(int is1=0;is1<20;is1++){hsum[is1]=0;}
			for(int i2=0;i2<ed.PMTWF->at(i1).size();i2++)
// 			for(int i2=2500;i2<10000;i2++)
			{
				Ampvst[i1]->Fill(i2,-1.*(ed.PMTWF->at(i1)[i2]-hped));
// 				if(I<1000) Evvst[i1]->Fill(i2,I,ed.PMTWF->at(i1)[i2]);
				
				for(int i3=0;i3<10;i3++)
				{
					if(i2>0 && i2%((i3+1)*5)==0){hh[i3][i1]->Fill(hsum[i3]);hsum[i3]=0;}
				}
				if(RunNo>7000)
				{
					for(int is1=0;is1<10;is1++){hsum[is1]+=(((-1.*(ed.PMTWF->at(i1)[i2]-PMTBaselines[i1][0]))>ThPMT[i1])?(-1.*(ed.PMTWF->at(i1)[i2]-PMTBaselines[i1][0])):0);}
				}
				else
				{
					for(int is1=0;is1<10;is1++){hsum[is1]+=(-1.*(ed.PMTWF->at(i1)[i2]-hped));}
				}
			}
			for(int i2=0;i2<ed.PMTWF->at(i1).size();i2+=60)
			{
				str=0;ss=0;
				for(int i3=i2;i3<(i2+15);i3++)
				{
					s=-1.*(ed.PMTWF->at(i1)[i3]-hped);
					str+=s;
				}
				for(int i3=i2+45;i3<i2+60;i3++)
				{
					s=-1.*(ed.PMTWF->at(i1)[i3]-hped);
					str+=s;
				}
				for(int i3=i2+15;i3<i2+45;i3++)
				{
					s=-1.*(ed.PMTWF->at(i1)[i3]-hped);
					ss+=s;
				}
				htr[i1]->Fill(str);
			}
		}
		if(I%1000==0) cout<<I<<" / "<<T->GetEntries()<<endl;
	}
// 	sprintf(hname,"Files/Calibrations_%d.txt",RunNo);
// 	ofstream outfile(hname);
	
	outroot->cd();
	
	float transientcuts[3]={0};
	for(int i1=0;i1<3;i1++)
	{
		htr[i1]->Fit(tfg,"q","q",-50,100);
		htr[i1]->Write();
		transientcuts[i1]=tfg->GetParameter(1)+3.5*tfg->GetParameter(2);
	}
	for(int I=0;I<T->GetEntries();I++)
	{
		T->GetEntry(I);
		for(int i1=0;i1<3;i1++)
		{
			hped=0;
			for(int i2=0;i2<500;i2++){hped+=ed.PMTWF->at(i1)[i2];}
			hped/=500.;
			if(I<200) continue;
			
			for(int i2=0;i2<ed.PMTWF->at(i1).size();i2+=60)
			{
				str=0;ss=0;
				for(int i3=i2;i3<(i2+15);i3++)
				{
					s=-1.*(ed.PMTWF->at(i1)[i3]-hped);
					str+=s;
				}
				for(int i3=i2+45;i3<i2+60;i3++)
				{
					s=-1.*(ed.PMTWF->at(i1)[i3]-hped);
					str+=s;
				}
				for(int i3=i2+15;i3<i2+45;i3++)
				{
					s=-1.*(ed.PMTWF->at(i1)[i3]-hped);
					ss+=s;
				}
				if(str<transientcuts[i1])
				{
					hs[i1]->Fill(ss);
				}
			}
		}
		if(I%1000==0) cout<<I<<" / "<<T->GetEntries()<<endl;
	}
	for(int i1=0;i1<3;i1++)
	{
		hs[i1]->Fit(tfg,"q","q",-50.,100.);
		tfgg->SetParameter(0,tfg->GetParameter(0));
		tfgg->SetParameter(1,tfg->GetParameter(1));
		tfgg->SetParameter(2,tfg->GetParameter(2));
		tfgg->SetParameter(3,tfg->GetParameter(0)/100.);
		tfgg->SetParameter(4,tfg->GetParameter(1)+30.);
		tfgg->SetParameter(5,tfg->GetParameter(2)*2);
		hs[i1]->Fit(tfgg,"q","q",-50.,60.);
		float chi2b=1000;float xchi2b=0.;float xchi2bsel=0.;
		xchi2b=40.;
		while(xchi2b<60.)
		{
			hs[i1]->Fit(tfgg,"q","q",-50.,xchi2b);
			if((tfgg->GetChisquare()/tfgg->GetNDF())<chi2b)
			{
				chi2b=(tfgg->GetChisquare()/tfgg->GetNDF());
				xchi2bsel=xchi2b;
			}
			xchi2b+=0.1;
		}
		hs[i1]->Fit(tfgg,"q","q",-50.,xchi2bsel);
		hs[i1]->Write();
		
		hs2[i1]=(TH1F*) hs[i1]->Clone();
		sprintf(hname,"SignalWindow2_%d",i1);
		hs2[i1]->SetName(hname);hs2[i1]->SetTitle(hname);
		hs2[i1]->Fit(tfgg,"q","q",-50.,xchi2bsel);
		tfggg->SetParameter(0,tfgg->GetParameter(0));
		tfggg->SetParameter(1,tfgg->GetParameter(1));
		tfggg->SetParameter(2,tfgg->GetParameter(2));
		tfggg->SetParameter(3,tfgg->GetParameter(3));
		tfggg->SetParameter(4,tfgg->GetParameter(4));
		tfggg->SetParameter(5,tfgg->GetParameter(5));
		tfggg->SetParameter(6,tfgg->GetParameter(3)/10);
		tfggg->SetParameter(7,tfgg->GetParameter(4)*2);
		tfggg->SetParameter(8,tfgg->GetParameter(5)*4);
		hs2[i1]->Fit(tfggg,"q","q",-50.,tfgg->GetParameter(4)+1.5*tfgg->GetParameter(5));
		hs2[i1]->Write();
	}
	
	for(int i1=0;i1<3;i1++)
	{
// 		htr[i1]->Fit(tfg,"q","q",-50,100);
// 		htr[i1]->Write();
// 		hs[i1]->Write();
		
		for(int i2=0;i2<10;i2++)
		{
			if(hh[i2][i1]->GetEntries()>0)
			{
				hh[i2][i1]->Fit(tfg,"q","q",-50.,100.);
				tfgg->SetParameter(0,tfg->GetParameter(0));
				tfgg->SetParameter(1,tfg->GetParameter(1));
				tfgg->SetParameter(2,tfg->GetParameter(2));
				tfgg->SetParameter(3,tfg->GetParameter(0)/100.);
				tfgg->SetParameter(4,tfg->GetParameter(1)+30.);
				tfgg->SetParameter(5,tfg->GetParameter(2)*2);
				
				tfgg->SetParLimits(3,0,1e6);
				tfgg->SetParLimits(4,20,60);
				tfgg->SetParLimits(5,0,100);
				
				hh[i2][i1]->Fit(tfgg,"q","q",-50.,60.);
				
				tfggg->SetParameter(0,tfgg->GetParameter(0));
				tfggg->SetParameter(1,tfgg->GetParameter(1));
				tfggg->SetParameter(2,tfgg->GetParameter(2));
				tfggg->SetParameter(3,tfgg->GetParameter(3));
				tfggg->SetParameter(4,tfgg->GetParameter(4));
				tfggg->SetParameter(5,tfgg->GetParameter(5));
				tfggg->SetParameter(6,tfgg->GetParameter(3)/10);
				tfggg->SetParameter(7,tfgg->GetParameter(4)*2);
				tfggg->SetParameter(8,tfgg->GetParameter(5)*4);
				
				tfggg->SetParLimits(6,0,1e6);
				tfggg->SetParLimits(7,30,100);
				tfggg->SetParLimits(8,0,200);
				
				hh[i2][i1]->Fit(tfggg,"q","q",-50.,100.);
			}
		}
	}
// 	
	for(int i1=0;i1<3;i1++)
	{
		for(int i2=0;i2<10;i2++)
		{
			if(hh[i2][i1]->GetEntries()>0) hh[i2][i1]->Write();
		}
		Ampvst[i1]->Write();
// 		hh_5[i1]->Write();
// 		hh_10[i1]->Write();
// 		hh_20[i1]->Write();
// 		Evvst[i1]->Write();
		BLvsE[i1]->Write();
		BLs[i1]->Write();
	}
	
// 	outfile.close();
	outroot->Close();
	
	sprintf(hname,"cp PMTCalibration_%d.root %s/Histos/PMTCalibration_%d.root",RunNo,AnalysisFilePath,RunNo);system(hname);
// 	for(int i1=0;i1<3;i1++)
// 	{
// 		sprintf(hname,"rm waveforms_ch%d_run%d.txt",i1,RunNo);system(hname);
// 	}
}

void GetPMTCalibration3()
{
	sprintf(hname,"%s/Run_%d.root",NTupleFilePath,RunNo);
	TFile* inroot=new TFile(hname);
	TTree* T =  (TTree*) inroot->Get("T");
	T->SetBranchAddress("TPCWF",&ed.TPCWF);
	T->SetBranchAddress("PMTWF",&ed.PMTWF);
	
	ReadPMTBaselines();
	
	sprintf(hname,"PMTCalibration3_%d.root",RunNo);
	TFile* outroot=new TFile(hname,"recreate");
	
// 	string PMTNames[3]={"LAr w/TPB","LAr no TPB","LXe"};
	
	TH1F* Integral_30[3];TH2F* Integral_vs_width[3];TH1F* Integral_w[3][30];
	TH1F* MaxAmplitude[3];TH1F* Integral_MAgt2[3];TH1F* Integral_MAgt2_wgt10[3];
	for(int i1=0;i1<3;i1++)
	{
// 		sprintf(hname,"Integral_30_%s",PMTNames[i1].c_str());
		sprintf(hname,"Integral_30_%d",i1);
		Integral_30[i1]=new TH1F(hname,hname,200,0,200);
		sprintf(hname,"MaxAmplitude_%d",i1);
		MaxAmplitude[i1]=new TH1F(hname,hname,200,0,200);
		sprintf(hname,"Integral_MAgt2_%d",i1);
		Integral_MAgt2[i1]=new TH1F(hname,hname,200,0,200);
		sprintf(hname,"Integral_MAgt2_wgt10_%d",i1);
		Integral_MAgt2_wgt10[i1]=new TH1F(hname,hname,200,0,200);
// 		sprintf(hname,"Integral_vs_width_%s",PMTNames[i1].c_str());
		sprintf(hname,"Integral_vs_width_%d",i1);
		Integral_vs_width[i1]=new TH2F(hname,hname,50,0,50,200,0,200);
		for(int i2=0;i2<30;i2++)
		{
// 			sprintf(hname,"Integral_w_%s_%d",PMTNames[i1].c_str(),i2);
			sprintf(hname,"Integral_w_%d_%d",i1,i2);
			Integral_w[i1][i2]=new TH1F(hname,hname,200,0,200);
		}
	}
	
	TH1F* hh=new TH1F("PMTWF","PMTWF",15000,-0.5,14999.5);
	
	TF1* g1=new TF1("g1","gaus",-1000.,2000.);
// 	TF1* tf1=new TF1("tf1","[0]*(exp([0]+[1]*(x-[2]))+exp([3]+[4]*(x-[2]))+exp([5]+[6]*(x-[2]))+exp([7]+[8]*(x-[2])))",0.,15000.);
// 	TF1* tf1=new TF1("tf1","(exp([0]+[1]*(x-[2]))+exp([3]+[4]*(x-[2]))+exp([5]+[6]*(x-[2])))",0.,15000.);
// 	TF1* tf1=new TF1("tf1","gaus(0)+expo(3)+expo(5)",0.,15000.);
// 	TF1* tf1=new TF1("tf1","[0]*pow(x-[1],[2])*([3]*exp(-(x-[1])/[4])+[5]*exp(-(x-[1])/[6]))",0.,15000.);
// 	TF1* tf1=new TF1("tf1","[0]*pow(x-[1],[2])*(exp(-(x-[1])/[3])+exp(-(x-[1])/[4]))",0.,15000.);
// 	TF1* tf1=new TF1("tf1","[0]*exp(-0.5*(pow((log((x)/[2]))/[3],2)))",0.,15000.);
// 	TF1* tf1=new TF1("tf1","landau(0)+expo(3)",0.,15000.);
// 	TF1* tf1=new TF1("tf1","[0]*exp(-pow((x-[1])/[2],2))+[6]*exp(-x/[3])+[7]*exp(-x/[4])+[8]*exp(-x/[5])",0.,15000.);
// 	TF1 *tf1 = new TF1("fit",fitf2f,0.,15000,6);
// 	TF1 *tf1 = new TF1("fit",fitf3f,0.,15000,7);
//  	TF1* tf1=new TF1("tf1","([0]/[3])*exp(-(log(x-[1])-[2])^2/(2.*[3]^2))",0.,15000.);
	TF1* tf1=new TF1("tf1","expo",0.,15000.);
	
	vector <float> Amp;
	float hsum[20]={0.};
	float hped=0;
	float hped2=0;float hsig2=0;
	float xlim=0;int fitstartbin=0;int xlimbin=0;
	float qint=0;int sbin=0;int minbin=0;int maxbin=0;int peakbin=0;
	float maxamp=0;int nc=0;
	
	for(int I=0;I<T->GetEntries();I++)
	{
		T->GetEntry(I);
		for(int i1=0;i1<3;i1++)
		{
			for(int i2=0;i2<ed.PMTWF->at(i1).size();i2++)
			{
				hh->SetBinContent(i2+1,-1.*(ed.PMTWF->at(i1)[i2]-PMTBaselines[i1][0]));
			}
			if(hh->GetNbinsX()==0) continue;
// 			if(I<10)
// 			{
				sprintf(hname,"PMTWF_%d_%d",i1,I);
				hh->SetName(hname);hh->SetTitle(hname);
				
				for(int ik1=hh->GetMaximumBin();ik1<hh->GetNbinsX();ik1++)
				{
					if(hh->GetBinContent(ik1)>0.1*hh->GetBinContent(hh->GetMaximumBin()))
// 					if(hh->GetBinContent(ik1)>10)
					{
						fitstartbin=ik1;
					}
					else break;
				}
				
				hh->Fit(tf1,"q","q",hh->GetBinCenter(fitstartbin),10000.);
				xlim=tf1->GetX(1,0,15000);
				xlimbin=hh->FindBin(xlim);
				
				sbin=xlimbin;
// 				cout<<I<<" "<<i1<<" "<<sbin<<endl;
				while(sbin<10000)
				{
					nc=0;
					for(int ik1=sbin;ik1<hh->GetNbinsX();ik1++)
					{
						if(hh->GetBinContent(ik1+1)>0 && hh->GetBinContent(ik1)<0)
						{
							nc++;
							minbin=ik1+1;
							for(int ik2=ik1+1;ik2<hh->GetNbinsX();ik2++)
							{
								if(hh->GetBinContent(ik2+1)<0 && hh->GetBinContent(ik2)>0)
								{
									maxbin=ik2;
									break;
								}
							}
							qint=0;maxamp=0;
							for(int ik2=minbin;ik2<=maxbin;ik2++)
							{
								qint+=hh->GetBinContent(ik2);
								if(hh->GetBinContent(ik2)>maxamp) maxamp=hh->GetBinContent(ik2);
							}
							MaxAmplitude[i1]->Fill(maxamp);
							if(maxamp>2) Integral_MAgt2[i1]->Fill(qint);
							if(maxamp>2 && (maxbin-minbin+1)>10) Integral_MAgt2_wgt10[i1]->Fill(qint);
							if((maxbin-minbin+1)>=5)
							{
								Integral_vs_width[i1]->Fill(maxbin-minbin+1,qint);
								for(int ik3=0;ik3<30;ik3++)
								{
									if((maxbin-minbin+1)>=(ik3+5))
									{
										Integral_w[i1][ik3]->Fill(qint);
									}
								}
							}
	// 						if(i1==2) cout<<I<<" "<<i1<<" "<<(maxbin-minbin+1)<<endl;
// 							nc=0;
							sbin=maxbin+1;
							break;
						}
// 						if(sbin>10000) break;
					}
					if(nc==0) break;
				}
// 				for(int ik1=sbin;ik1<hh->GetNbinsX();ik1++)
// 				{
// 					if(hh->GetBinContent(ik1+1)>0 && hh->GetBinContent(ik1)<0)
// 					{
// 						minbin=ik1;
// 						for(int ik2=ik1;ik2<hh->GetNbinsX();ik2++)
// 						{
// 							if(hh->GetBinContent(ik2+1)<0 && hh->GetBinContent(ik2)>0)
// 							{
// 								maxbin=ik2;
// 								break;
// 							}
// 						}
// 						qint=0;maxamp=0;
// 						for(int ik2=minbin;ik2<=maxbin;ik2++)
// 						{
// 							qint+=hh->GetBinContent(ik2);
// 							if(hh->GetBinContent(ik2)>maxamp) maxamp=hh->GetBinContent(ik2);
// 						}
// 						if(maxamp>2) Integral_MAgt2[i1]->Fill(qint);
// 						if(maxamp>2 && (maxbin-minbin+1)>10) Integral_MAgt2_wgt10[i1]->Fill(qint);
// 						if((maxbin-minbin+1)>=5)
// 						{
// 							Integral_vs_width[i1]->Fill(maxbin-minbin+1,qint);
// 							for(int ik3=0;ik3<30;ik3++)
// 							{
// 								if((maxbin-minbin+1)>=(ik3+5))
// 								{
// 									Integral_w[i1][ik3]->Fill(qint);
// 								}
// 							}
// 						}
// // 						if(i1==2) cout<<I<<" "<<i1<<" "<<(maxbin-minbin+1)<<endl;
// 					}
// 					sbin=maxbin+1;
// 					if(sbin>10000) break;
// 				}
// 				cout<<I<<" "<<i1<<" "<<xlim<<endl;
				for(int ik1=xlimbin;ik1<hh->GetNbinsX();ik1+=30)
				{
					qint=0;
					for(int ik2=ik1;ik2<(ik1+30);ik2++)
					{
						qint+=hh->GetBinContent(ik2);
					}
					Integral_30[i1]->Fill(qint);
				}
				
				
// 				hh->Write();
// 			}
			hh->Reset();
		}
		if(I%1000==0) cout<<I<<" / "<<T->GetEntries()<<endl;
// 		if(I>1000) break;
	}
	
	outroot->cd();
	for(int i1=0;i1<3;i1++)
	{
		Integral_30[i1]->Write();
		Integral_vs_width[i1]->Write();
		MaxAmplitude[i1]->Write();
		Integral_MAgt2[i1]->Write();
		Integral_MAgt2_wgt10[i1]->Write();
		for(int i2=0;i2<30;i2++)
		{
			Integral_w[i1][i2]->Write();
		}
	}
	outroot->Close();
	
	sprintf(hname,"cp PMTCalibration3_%d.root %s/Histos/PMTCalibration3_%d.root",RunNo,AnalysisFilePath,RunNo);system(hname);
}

void GetPMTCalibration4()//for PMT calibration runs with VME
{
	sprintf(hname,"%s/Run_%d.root",NTupleFilePath,RunNo);
	TFile* inroot=new TFile(hname);
	TTree* T =  (TTree*) inroot->Get("T");
	T->SetBranchAddress("TPCWF",&ed.TPCWF);
	T->SetBranchAddress("PMTWF",&ed.PMTWF);
	
	sprintf(hname,"PMTCalibration4_%d.root",RunNo);
	TFile* outroot=new TFile(hname,"recreate");
	
	if(RunNo>=208) NPMTs=4;
	
	TH1F* Integral[NPMTs];TH1F* IntegralSel[NPMTs];
	TH1F* MaxAmplitude[NPMTs];TH1F* PulseWidth[NPMTs];TH1F* PulseWidthSel[NPMTs];
	TH1F* Baseline[NPMTs];TH1F* SelAmplitude[NPMTs];TH1F* SelAmplitudeOverRMS[NPMTs];
	for(int i1=0;i1<NPMTs;i1++)
	{
		sprintf(hname,"Integral_%d",i1);
// 		Integral[i1]=new TH1F(hname,hname,1200,0,200);
		Integral[i1]=new TH1F(hname,hname,2400,0,200);
		sprintf(hname,"IntegralSel_%d",i1);
		IntegralSel[i1]=new TH1F(hname,hname,1200,0,200);
		sprintf(hname,"MaxAmplitude_%d",i1);
		MaxAmplitude[i1]=new TH1F(hname,hname,1000,0,200);
		sprintf(hname,"SelAmplitude_%d",i1);
		SelAmplitude[i1]=new TH1F(hname,hname,1000,0,200);
		sprintf(hname,"SelAmplitudeOverRMS_%d",i1);
// 		SelAmplitudeOverRMS[i1]=new TH1F(hname,hname,1200,0,200);
		SelAmplitudeOverRMS[i1]=new TH1F(hname,hname,1600,0,200);
		sprintf(hname,"PulseWidth_%d",i1);
		PulseWidth[i1]=new TH1F(hname,hname,200,-0.5,199.5);
		sprintf(hname,"PulseWidthSel_%d",i1);
		PulseWidthSel[i1]=new TH1F(hname,hname,200,-0.5,199.5);
		sprintf(hname,"Baseline_%d",i1);
		Baseline[i1]=new TH1F(hname,hname,4500,0,1500);
	}
	
	TH1F* hh=new TH1F("PMTWF","PMTWF",15000,-0.5,14999.5);
	TH1F* hbl=new TH1F("BaselineHist","BaselineHist",1500,-0.5,1499.5);
	
	TF1* g1=new TF1("g1","gaus",-1000.,2000.);
// 	TF1* tf1=new TF1("tf1","[0]*(exp([0]+[1]*(x-[2]))+exp([3]+[4]*(x-[2]))+exp([5]+[6]*(x-[2]))+exp([7]+[8]*(x-[2])))",0.,15000.);
// 	TF1* tf1=new TF1("tf1","(exp([0]+[1]*(x-[2]))+exp([3]+[4]*(x-[2]))+exp([5]+[6]*(x-[2])))",0.,15000.);
// 	TF1* tf1=new TF1("tf1","gaus(0)+expo(3)+expo(5)",0.,15000.);
// 	TF1* tf1=new TF1("tf1","[0]*pow(x-[1],[2])*([3]*exp(-(x-[1])/[4])+[5]*exp(-(x-[1])/[6]))",0.,15000.);
// 	TF1* tf1=new TF1("tf1","[0]*pow(x-[1],[2])*(exp(-(x-[1])/[3])+exp(-(x-[1])/[4]))",0.,15000.);
// 	TF1* tf1=new TF1("tf1","[0]*exp(-0.5*(pow((log((x)/[2]))/[3],2)))",0.,15000.);
// 	TF1* tf1=new TF1("tf1","landau(0)+expo(3)",0.,15000.);
// 	TF1* tf1=new TF1("tf1","[0]*exp(-pow((x-[1])/[2],2))+[6]*exp(-x/[3])+[7]*exp(-x/[4])+[8]*exp(-x/[5])",0.,15000.);
// 	TF1 *tf1 = new TF1("fit",fitf2f,0.,15000,6);
// 	TF1 *tf1 = new TF1("fit",fitf3f,0.,15000,7);
//  	TF1* tf1=new TF1("tf1","([0]/[3])*exp(-(log(x-[1])-[2])^2/(2.*[3]^2))",0.,15000.);
	TF1* tf1=new TF1("tf1","expo",0.,15000.);
	TF1* tflin=new TF1("tflin","[0]",0.,2000.);
	
	vector <float> Amp;
	float hsum[20]={0.};
	float hped=0;
	float hped2=0;float hsig2=0;
	float xlim=0;int fitstartbin=0;int xlimbin=0;
	float maxamp=0;int nc=0;
	
	float baseline[NPMTs];baseline[NPMTs]={0.};
	float baselineRMS[NPMTs];baselineRMS[NPMTs]={0.};
	float lmax=0;float pulseamp=0;
	float qint=0;int sbin=0;int minbin=0;int maxbin=0;int peakbin=0;
	
	for(int I=0;I<T->GetEntries();I++)
	{
		T->GetEntry(I);
		if(I<100) continue;
		for(int i1=0;i1<NPMTs;i1++)
		{
			hbl->Reset();
			for(int i2=0;i2<1000;i2++)
			{
				hbl->Fill(ed.PMTWF->at(i1)[i2]);
			}
			hbl->Fit(g1,"q","q",hbl->GetMean()-3*hbl->GetRMS(),hbl->GetMean()+3*hbl->GetRMS());
			baseline[i1]=g1->GetParameter(1);
			baselineRMS[i1]=g1->GetParameter(2);
			Baseline[i1]->Fill(baseline[i1]);
			hh->Reset();
			lmax=0;
			for(int i2=0;i2<ed.PMTWF->at(i1).size();i2++)
			{
				pulseamp=-1.*(ed.PMTWF->at(i1)[i2]-baseline[i1]);
				hh->SetBinContent(i2+1,pulseamp);
				if(pulseamp>lmax){lmax=pulseamp;}
				if(i2>0 && i2%1000==0)
				{
					MaxAmplitude[i1]->Fill(lmax);
					lmax=0;
				}
			}
			sbin=0;
			while(sbin<14000)
			{
				minbin=-1;
// 				minbin=FindBinAbove(hh,3*baselineRMS[i1],sbin);
				minbin=FindBinAbove(hh,3.5*baselineRMS[i1],sbin);
				if(minbin==-1) break;
				maxbin=-1;
				for(int ik1=minbin;ik1<=hh->GetNbinsX();ik1++)
				{
					if(hh->GetBinContent(ik1)<0) {maxbin=ik1-1;break;}
				}
				if(maxbin==-1) break;
				for(int ik1=minbin;ik1>=1;ik1--)
				{
					if(hh->GetBinContent(ik1)<0) {minbin=ik1+1;break;}
				}
				if(minbin==2) break;
				qint=0.;
				lmax=0.;
				for(int ik1=minbin;ik1<=maxbin;ik1++)
				{
					qint+=hh->GetBinContent(ik1);
					if(hh->GetBinContent(ik1)>lmax) lmax=hh->GetBinContent(ik1);
				}
				Integral[i1]->Fill(qint);
				PulseWidth[i1]->Fill(maxbin-minbin+1);
				SelAmplitude[i1]->Fill(lmax);
				SelAmplitudeOverRMS[i1]->Fill(lmax/baselineRMS[i1]);
// 				if((lmax/baselineRMS[i1])>4.)
// 				{
// 					IntegralSel[i1]->Fill(qint);
// 					PulseWidthSel[i1]->Fill(maxbin-minbin+1);
// 				}
				sbin=maxbin+1;
			}
		}
		if(I%1000==0) cout<<I<<" / "<<T->GetEntries()<<endl;
	}
	
// 	float maxampth[NPMTs];
// 	for(int i1=0;i1<NPMTs;i1++)
// 	{
// 		maxbin=SelAmplitudeOverRMS[i1]->GetMaximumBin();
// 		for(int is1=maxbin;is1<SelAmplitudeOverRMS[i1]->GetNbinsX();is1++)
// 		{
// 			if(SelAmplitudeOverRMS[i1]->GetBinContent(is1+1)>SelAmplitudeOverRMS[i1]->GetBinContent(is1))
// 			{
// 				sbin=is1;
// 				break;
// 			}
// 		}
// 		
// 		
// // 		maxamp=SelAmplitudeOverRMS[i1]->GetBinContent(maxbin);
// // 		for(int is1=maxbin;is1<=SelAmplitudeOverRMS[i1]->GetNbinsX();is1++)
// // 		{
// // 			if((SelAmplitudeOverRMS[i1]->GetBinContent(is1)/maxamp)<0.01){minbin=is1;break;}
// // 		}
// // 		sbin=minbin;
// // 		for(int is1=minbin;is1<SelAmplitudeOverRMS[i1]->GetNbinsX();is1++)
// // 		{
// // 			if(SelAmplitudeOverRMS[i1]->GetBinContent(is1+1)>SelAmplitudeOverRMS[i1]->GetBinContent(is1))
// // 			{
// // 				sbin=is1;
// // 				break;
// // 			}
// // 		}
// 		maxampth[i1]=SelAmplitudeOverRMS[i1]->GetBinCenter(sbin);
// 		
// // 		g1->SetParameter(0,SelAmplitudeOverRMS[i1]->GetBinContent(SelAmplitudeOverRMS[i1]->GetMaximumBin()));
// // 		g1->SetParameter(1,SelAmplitudeOverRMS[i1]->GetBinCenter(SelAmplitudeOverRMS[i1]->GetMaximumBin()));
// // 		g1->SetParameter(2,0.25);
// // 		SelAmplitudeOverRMS[i1]->Fit(g1,"q","q",SelAmplitudeOverRMS[i1]->GetBinCenter(SelAmplitudeOverRMS[i1]->GetMaximumBin())-1,SelAmplitudeOverRMS[i1]->GetBinCenter(SelAmplitudeOverRMS[i1]->GetMaximumBin())+1);
// 	}
// 	
// 	for(int I=0;I<T->GetEntries();I++)
// 	{
// 		T->GetEntry(I);
// 		if(I<100) continue;
// 		for(int i1=0;i1<NPMTs;i1++)
// 		{
// 			hbl->Reset();
// 			for(int i2=0;i2<1000;i2++)
// 			{
// 				hbl->Fill(ed.PMTWF->at(i1)[i2]);
// 			}
// 			hbl->Fit(g1,"q","q",hbl->GetMean()-3*hbl->GetRMS(),hbl->GetMean()+3*hbl->GetRMS());
// 			baseline[i1]=g1->GetParameter(1);
// 			baselineRMS[i1]=g1->GetParameter(2);
// // 			Baseline[i1]->Fill(baseline[i1]);
// 			hh->Reset();
// 			lmax=0;
// 			for(int i2=0;i2<ed.PMTWF->at(i1).size();i2++)
// 			{
// 				pulseamp=-1.*(ed.PMTWF->at(i1)[i2]-baseline[i1]);
// 				hh->SetBinContent(i2+1,pulseamp);
// // 				if(pulseamp>lmax){lmax=pulseamp;}
// // 				if(i2>0 && i2%1000==0)
// // 				{
// // 					MaxAmplitude[i1]->Fill(lmax);
// // 					lmax=0;
// // 				}
// 			}
// 			sbin=0;
// 			while(sbin<14000)
// 			{
// 				minbin=-1;
// 				minbin=FindBinAbove(hh,maxampth[i1],sbin);
// // 				minbin=FindBinAbove(hh,3*baselineRMS[i1],sbin);
// // 				minbin=FindBinAbove(hh,2*baselineRMS[i1],sbin);
// 				if(minbin==-1) break;
// 				maxbin=-1;
// 				for(int ik1=minbin;ik1<=hh->GetNbinsX();ik1++)
// 				{
// 					if(hh->GetBinContent(ik1)<0) {maxbin=ik1-1;break;}
// 				}
// 				if(maxbin==-1) break;
// 				for(int ik1=minbin;ik1>=1;ik1--)
// 				{
// 					if(hh->GetBinContent(ik1)<0) {minbin=ik1+1;break;}
// 				}
// 				if(minbin==2) break;
// 				qint=0.;
// 				lmax=0.;
// 				for(int ik1=minbin;ik1<=maxbin;ik1++)
// 				{
// 					qint+=hh->GetBinContent(ik1);
// 					if(hh->GetBinContent(ik1)>lmax) lmax=hh->GetBinContent(ik1);
// 				}
// // 				Integral[i1]->Fill(qint);
// // 				PulseWidth[i1]->Fill(maxbin-minbin+1);
// // 				SelAmplitude[i1]->Fill(lmax);
// // 				SelAmplitudeOverRMS[i1]->Fill(lmax/baselineRMS[i1]);
// // 				if((lmax/baselineRMS[i1])>4.)
// // 				{
// 					IntegralSel[i1]->Fill(qint);
// 					PulseWidthSel[i1]->Fill(maxbin-minbin+1);
// // 				}
// 				sbin=maxbin+1;
// 			}
// 		}
// 		if(I%1000==0) cout<<I<<" / "<<T->GetEntries()<<endl;
// 	}
	
	float PMTMeans[4]={36,37,24,12};
	float PMTSigmaF[4]={1.5,1.5,1,1.5};
	if(RunNo>=297 && RunNo<=338){PMTMeans[1]=88;PMTMeans[3]=20;}
	if(RunNo==310){PMTMeans[1]=37;PMTMeans[3]=12;}
	
	outroot->cd();
	for(int i1=0;i1<NPMTs;i1++)
	{
		g1->SetParameter(0,500);
		g1->SetParameter(1,PMTMeans[i1]);
		g1->SetParameter(2,10);
		//The problematic runs first
		if(RunNo==208 && i1==3) Integral[i1]->Fit(g1,"q","q",9,20);
		else if(RunNo==210 && i1==3) Integral[i1]->Fit(g1,"q","q",8,20);
		else if(RunNo==217 && i1==1) Integral[i1]->Fit(g1,"q","q",9,18);
		else if(RunNo==245 && i1==2) Integral[i1]->Fit(g1,"q","q",18,35);
		else if(RunNo==260 && i1==2) Integral[i1]->Fit(g1,"q","q",18,35);
		else if(RunNo==273 && i1==1) Integral[i1]->Fit(g1,"q","q",25,50);
		else if(RunNo==302 && i1==1) Integral[i1]->Fit(g1,"q","q",40,120);
		else if(RunNo==329 && i1==3) Integral[i1]->Fit(g1,"q","q",14,28);
		else
		{
			Integral[i1]->Fit(g1,"q","q",PMTMeans[i1]-10,PMTMeans[i1]+10);
			Integral[i1]->Fit(g1,"q","q",g1->GetParameter(1)-PMTSigmaF[i1]*g1->GetParameter(2),g1->GetParameter(1)+PMTSigmaF[i1]*g1->GetParameter(2));
			Integral[i1]->Fit(g1,"q","q",g1->GetParameter(1)-PMTSigmaF[i1]*g1->GetParameter(2),g1->GetParameter(1)+PMTSigmaF[i1]*g1->GetParameter(2));
		}
		Integral[i1]->Write();
		
// 		g1->SetParameter(0,IntegralSel[i1]->GetBinContent(IntegralSel[i1]->GetMaximumBin()));
// 		g1->SetParameter(1,IntegralSel[i1]->GetBinCenter(IntegralSel[i1]->GetMaximumBin()));
// 		g1->SetParameter(2,10);
// 		IntegralSel[i1]->Fit(g1,"q","q",IntegralSel[i1]->GetBinContent(IntegralSel[i1]->GetMaximumBin())-10,IntegralSel[i1]->GetBinContent(IntegralSel[i1]->GetMaximumBin())+10);
// 		IntegralSel[i1]->Fit(g1,"q","q",g1->GetParameter(1)-1.5*g1->GetParameter(2),g1->GetParameter(1)+1.5*g1->GetParameter(2));
// 		IntegralSel[i1]->Fit(g1,"q","q",g1->GetParameter(1)-1.5*g1->GetParameter(2),g1->GetParameter(1)+1.5*g1->GetParameter(2));
// 		IntegralSel[i1]->Write();
		
		MaxAmplitude[i1]->Write();
		SelAmplitude[i1]->Write();
		SelAmplitudeOverRMS[i1]->Write();
		PulseWidth[i1]->Write();
// 		PulseWidthSel[i1]->Write();
		Baseline[i1]->Write();
	}
	outroot->Close();
	
	sprintf(hname,"cp PMTCalibration4_%d.root %s/Histos/PMTCalibration4_%d.root",RunNo,AnalysisFilePath,RunNo);system(hname);
}

void GetPMTCalibration_CalibRuns()
{
	sprintf(hname,"PMTCalibration_CalibRun_%d.root",RunNo);
	TFile* outroot=new TFile(hname,"recreate");
	
	sprintf(hname,"%s/Run_%d.root",NTupleFilePath,RunNo);
	TFile* inroot=new TFile(hname);
	TTree* T =  (TTree*) inroot->Get("T");
	T->SetBranchAddress("TPCWF",&ed.TPCWF);
	T->SetBranchAddress("PMTWF",&ed.PMTWF);
	
	TH1F* hs2[NPMTs];
	TH1F* MeanPedestalsAll[NPMTs];
	TH1F* MeanPedestals[NPMTs];TH1F* LinFitPedestals[NPMTs];TH1F* GausFitPedestals[NPMTs];TH1F* GausFitPedestalRMS[NPMTs];
	TH1F* BLFitOverMean[NPMTs];
	TH1D* Ampvst[NPMTs];
	TH1F* BLs[NPMTs];
	TH1F* MaxAmpover3sigma[NPMTs];
	TH2F* MaxAmpover3sigmavsWidth[NPMTs];
	TH1F* Integral[NPMTs];
	TH1F* Amplitudes[NPMTs];
	TH1F* MaxAmplitudes[NPMTs];
	TH1F* PulseWidth[NPMTs];
	TH1F* IntegralSel[NPMTs][6];
	TH1F* IntegralSel2[NPMTs];
	TH1F* IntegralSel3[NPMTs];
	
	for(int i1=0;i1<NPMTs;i1++)
	{
		sprintf(hname,"Amp_vs_t_%d",i1);
		Ampvst[i1]=new TH1D(hname,hname,15000,-0.5,14999.5);
		sprintf(hname,"BLs_%d",i1);
		BLs[i1]=new TH1F(hname,hname,2000,0.,1000.);
		
		sprintf(hname,"MeanPedestalsAll_%d",i1);
		MeanPedestalsAll[i1]=new TH1F(hname,hname,4000,700.,1100.);
		sprintf(hname,"MeanPedestals_%d",i1);
		MeanPedestals[i1]=new TH1F(hname,hname,4000,700.,1100.);
		sprintf(hname,"LinFitPedestals_%d",i1);
		LinFitPedestals[i1]=new TH1F(hname,hname,4000,700.,1100.);
		sprintf(hname,"GausFitPedestals_%d",i1);
		GausFitPedestals[i1]=new TH1F(hname,hname,4000,700.,1100.);
		sprintf(hname,"GausFitPedestalRMS_%d",i1);
		GausFitPedestalRMS[i1]=new TH1F(hname,hname,400,0.,2.);
		sprintf(hname,"MaxAmpover3sigma_%d",i1);
		MaxAmpover3sigma[i1]=new TH1F(hname,hname,10000,0.,10.);
		sprintf(hname,"MaxAmpover3sigmavsWidth_%d",i1);
		MaxAmpover3sigmavsWidth[i1]=new TH2F(hname,hname,20,-0.5,19.5,10000,0.,10.);
		
		sprintf(hname,"BLFitOverMean_%d",i1);
		BLFitOverMean[i1]=new TH1F(hname,hname,5000,0.9,1.1);
		
		sprintf(hname,"Integral_%d",i1);
		Integral[i1]=new TH1F(hname,hname,1000,0,200);
		
		sprintf(hname,"PulseWidth_%d",i1);
		PulseWidth[i1]=new TH1F(hname,hname,30,-0.5,29.5);
		sprintf(hname,"Amplitudes_%d",i1);
		Amplitudes[i1]=new TH1F(hname,hname,200,0,20);
		sprintf(hname,"MaxAmplitudes_%d",i1);
		MaxAmplitudes[i1]=new TH1F(hname,hname,200,0,20);
		
		for(int i2=0;i2<6;i2++)
		{
			sprintf(hname,"IntegralSel_%d_%d",i1,i2);
			IntegralSel[i1][i2]=new TH1F(hname,hname,1000,0,200);
		}
		
		sprintf(hname,"IntegralSel2_%d",i1);
		IntegralSel2[i1]=new TH1F(hname,hname,1000,0,200);
		sprintf(hname,"IntegralSel3_%d",i1);
		IntegralSel3[i1]=new TH1F(hname,hname,1000,0,200);
	}
	TGraph* BLvsE[NPMTs];
	for(int i1=0;i1<NPMTs;i1++)
	{
		sprintf(hname,"BL_vs_E_%d",i1);
		BLvsE[i1]=new TGraph();
		BLvsE[i1]->SetName(hname);BLvsE[i1]->SetTitle(hname);
		BLvsE[i1]->SetMarkerStyle(24+i1);
		BLvsE[i1]->SetMarkerColor(1+i1);
		BLvsE[i1]->GetXaxis()->SetTitle("Event ID");BLvsE[i1]->GetXaxis()->CenterTitle();
		BLvsE[i1]->GetYaxis()->SetTitle("Mean Amplitude of First 1000 Time Bins");BLvsE[i1]->GetYaxis()->CenterTitle();
	}
	
	TF1* g1=new TF1("g1","gaus",-1000.,2000.);
	TF1* lin=new TF1("lin","[0]",0.,15000);
	
	TH1F* hb=new TH1F("hb","hb",15000,-0.5,14999.5);
	TH1F* hp=new TH1F("hp","hp",2000,-0.5,1999.5);
	
// 	TF1* tf=new TF1("GL","gaus(0)+landau(3)",-100.,200.);
// 	TF1* tf=new TF1("GLG","gaus(0)+landau(3)+gaus(6)",-100.,200.);
	TF1* tf=new TF1("GGG","gaus(0)+gaus(3)+gaus(6)",-50.,100.);
	
// 	TF1* tf=new TF1("GG","gaus(0)+gaus(3)",-50.,100.);
// 	TF1* tf=new TF1("GL","gaus(0)+landau(3)",-50.,100.);
// 	TF1* tf=new TF1("GLG","gaus(0)+landau(3)+gaus(6)",-50.,100.);
	TF1* tfg=new TF1("G","gaus",-50.,100.);
	TF1* tfgg=new TF1("GG","gaus(0)+gaus(3)",-50.,100.);
	TF1* tfggg=new TF1("GGG","gaus(0)+gaus(3)+gaus(6)",-50.,100.);
	
	vector <float> Amp;
	float hsum[20]={0.};
	float hped=0;
	float hped2=0;
	float hped3=0;float hsig3=0;
	float s=0;float str=0;float ss=0.;
	int minbin=0;int maxbin=0;
	float maxamp=0;
	float maxampratio=0;
	int isel=0;
	
	vector <float> peds[NPMTs][4];//mean, lin fit, gaus fit mean, gaus fit sigma
	
// 	cout<<T->GetEntries()<<endl;s
	for(int I=0;I<T->GetEntries();I++)
	{
		T->GetEntry(I);
		for(int i1=0;i1<NPMTs;i1++)
		{
			hped=0;
			hb->Reset();hp->Reset();
			for(int i2=0;i2<1000;i2++)
			{
				hped+=ed.PMTWF->at(i1)[i2];
				hb->SetBinContent(i2+1,ed.PMTWF->at(i1)[i2]);
				hp->Fill(ed.PMTWF->at(i1)[i2]);
			}
			hped/=1000.;
			hb->Fit(lin,"q","q",0.,1000.);
			hped2=lin->GetParameter(0);
			hp->Fit(g1,"q","q",0.,2000.);
			hped3=g1->GetParameter(1);
			hsig3=g1->GetParameter(2);
			MeanPedestalsAll[i1]->Fill(hped);
			BLvsE[i1]->SetPoint(BLvsE[i1]->GetN(),I,hped);
			
			peds[i1][0].push_back(hped);
			peds[i1][1].push_back(hped2);
			peds[i1][2].push_back(hped3);
			peds[i1][3].push_back(hsig3);
			
			if(I<200) continue;
			if(hped>0) BLFitOverMean[i1]->Fill(hped2/hped);
			MeanPedestals[i1]->Fill(hped);
			LinFitPedestals[i1]->Fill(hped2);
			GausFitPedestals[i1]->Fill(hped3);
			GausFitPedestalRMS[i1]->Fill(hsig3);
		}
	}
	for(int I=0;I<T->GetEntries();I++)
	{
// 		cout<<I<<endl;
		if(I<200) continue;
		T->GetEntry(I);
		for(int i1=0;i1<NPMTs;i1++)
		{
			hped3=peds[i1][2][I];
			hsig3=peds[i1][3][I];
			for(int i2=50;i2<ed.PMTWF->at(i1).size()-50;i2++)
			{
				s=-1*(ed.PMTWF->at(i1)[i2]-hped3);
				if(s>3.*hsig3)
				{
					maxbin=i2-1;
					for(int i3=i2;i3<ed.PMTWF->at(i1).size();i3++)
					{
						s=-1*(ed.PMTWF->at(i1)[i3]-hped3);
						if(s<=0)
						{
							maxbin=i3-1;
							break;
						}
					}
					minbin=i2+1;
					for(int i3=i2-1;i3>=1;i3--)
					{
						s=-1*(ed.PMTWF->at(i1)[i3]-hped3);
						if(s<=0)
						{
							minbin=i3+1;
							break;
						}
					}
					PulseWidth[i1]->Fill((maxbin-minbin+1));
					if((maxbin-minbin+1)>=8)
					{
						ss=0;maxamp=0;
						for(int i3=minbin;i3<=maxbin;i3++)
						{
							s=-1*(ed.PMTWF->at(i1)[i3]-hped3);
							ss+=s;
							if(s>maxamp) maxamp=s;
							Amplitudes[i1]->Fill(s);
						}
						MaxAmplitudes[i1]->Fill(maxamp);
						maxampratio=maxamp/(3.*hsig3);
						isel=0;
						if(i1==0)
						{
							if(maxampratio<1.15) isel=0;
							else if(maxampratio<1.61) isel=1;
							else if(maxampratio<2.08) isel=2;
							else if(maxampratio<2.6) isel=3;
							else if(maxampratio<3.05) isel=4;
							else isel=5;
						}
						else if(i1==1)
						{
							if(maxampratio<1.2) isel=0;
							else if(maxampratio<1.7) isel=1;
							else isel=2;
						}
						else if(i1==2)
						{
							if(maxampratio<1.3) isel=0;
							else if(maxampratio<1.8) isel=1;
							else if(maxampratio<2.35) isel=2;
							else if(maxampratio<2.95) isel=3;
							else if(maxampratio<3.5) isel=4;
							else isel=5;
						}
						
						
// 						if((i1==0 && (maxamp/(3.*hsig3))>1.15 && (maxamp/(3.*hsig3))<1.55)||(i1==1 && (maxamp/(3.*hsig3))>1.2 && (maxamp/(3.*hsig3))<1.6)||(i1==2 && (maxamp/(3.*hsig3))>1.3 && (maxamp/(3.*hsig3))<1.8))
// 						{
							Integral[i1]->Fill(ss);
							IntegralSel[i1][isel]->Fill(ss);
							if(maxampratio>2.) IntegralSel2[i1]->Fill(ss);
							if(maxampratio>2. && maxampratio<3.) IntegralSel3[i1]->Fill(ss);
// 						}
						MaxAmpover3sigma[i1]->Fill(maxamp/(3.*hsig3));
						MaxAmpover3sigmavsWidth[i1]->Fill(maxbin-minbin+1,maxamp/(3.*hsig3));
					}
					i2=maxbin+1;
					if(i2>(ed.PMTWF->at(i1).size()-50)) break;
				}
			}
		}
	}
	outroot->cd();
	
	for(int i1=0;i1<NPMTs;i1++)
	{
		Ampvst[i1]->Write();
		BLvsE[i1]->Write();
		MeanPedestalsAll[i1]->Write();
		MeanPedestals[i1]->Fit(g1,"q","q",700.,1100.);
		MeanPedestals[i1]->Write();
		LinFitPedestals[i1]->Fit(g1,"q","q",700.,1100.);
		LinFitPedestals[i1]->Write();
		GausFitPedestals[i1]->Fit(g1,"q","q",700.,1100.);
		GausFitPedestals[i1]->Write();
		GausFitPedestalRMS[i1]->Fit(g1,"q","q",0.5,1.5);
		GausFitPedestalRMS[i1]->Write();
		BLFitOverMean[i1]->Write();
		MaxAmpover3sigma[i1]->Write();
		MaxAmpover3sigmavsWidth[i1]->Write();
		PulseWidth[i1]->Write();
		Amplitudes[i1]->Write();
		MaxAmplitudes[i1]->Write();
		Integral[i1]->Write();
		IntegralSel2[i1]->Write();
		IntegralSel3[i1]->Write();
		for(int i2=0;i2<6;i2++)
		{
			IntegralSel[i1][i2]->SetLineColor(i2+1);
			IntegralSel[i1][i2]->Write();
		}
	}
	
	outroot->Close();
	
// 	sprintf(hname,"cp PMTCalibration_%d.root %s/Histos/PMTCalibration_%d.root",RunNo,AnalysisFilePath,RunNo);system(hname);
// 	for(int i1=0;i1<3;i1++)
// 	{
// 		sprintf(hname,"rm waveforms_ch%d_run%d.txt",i1,RunNo);system(hname);
// 	}
}

void GetPMTIntegral()
{
// 	ReadPMTBaselines();
// 	sprintf(hname,"PMTIntegrals_%d.root",RunNo);
// 	TFile* outroot=new TFile(hname,"recreate");
// 	
// 	vector <int> rr[3];
// 	
// 	TH1F* hh[3];
// 	
// 	for(int i1=0;i1<3;i1++)
// 	{
// 		sprintf(hname,"Integrals_%d",i1);
// 		hh[i1]=new TH1F(hname,hname,20000,-1000000,1000000);
// 	}
// 	
// // 	TF1* lin=new TF1("lin","[0]",0.,4095.);
// // 	TH1F* hb=new TH1F("hb","hb",15000,-0.5,14999.5);
// // // 	TF1* tf=new TF1("GL","gaus(0)+landau(3)",-100.,200.);
// // // 	TF1* tf=new TF1("GLG","gaus(0)+landau(3)+gaus(6)",-100.,200.);
// // 	TF1* tf=new TF1("GGG","gaus(0)+gaus(3)+gaus(6)",-100.,200.);
// 	
// 	int NSat[3]={0};
// 	
// 	for(int i1=0;i1<3;i1++)
// 	{
// 		sprintf(hname,"cp %s/waveforms_ch%d_run%d.txt .;wait;",PMTinFilePath,i1,RunNo);
// 		system(hname);
// 	}
// 	
// 	char cNum[10];
// 	int ii=0;int Ev=0;
// 	ifstream pmtinfile;
// 	bool isSaturated=false;
// 	float sg=0;
// 	float qint=0;
// 	for(int i1=0;i1<3;i1++)
// 	{
// 		Ev=0;
// 		sprintf(hname,"waveforms_ch%d_run%d.txt",i1,RunNo);
// 		cout<<hname<<endl;
// 		pmtinfile.open(hname);
// 		string strarr;
// 		while(getline(pmtinfile,strarr,'\n'))
// 		{
// 			istringstream iss(strarr);
// 			ii=0;
// 			qint=0;
// 			isSaturated=false;
// 			for(std::string s; iss >> s; )
// 			{
// 				if(ii>1500 && ii<7500)
// 				{
// 					sg=((-1.)*(atof(s.c_str())-PMTBaselines[i1][0]));
// 					if(sg>ThPMT)
// 					{
// 						qint+=sg;
// 					}
// 					if(atoi(s.c_str())==0){isSaturated=true;}
// 				}
// 				ii++;
// 			}
// 			if(isSaturated){qint*=-1.;NSat[i1]++;}
// 			rr[i1].push_back(qint);
// 			hh[i1]->Fill(qint);
// 			Ev++;
// 		}
// 		pmtinfile.close();
// 	}
// // 	sprintf(hname,"%s/Files/PMTIntegrals_%d.txt",AnalysisFilePath,RunNo);
// 	sprintf(hname,"PMTIntegrals_%d.txt",RunNo);
// 	ofstream outfile(hname);
// 	for(int i1=0;i1<rr[0].size();i1++)
// 	{
// 		outfile<<i1<<" "<<rr[0][i1]<<" "<<rr[1][i1]<<" "<<rr[2][i1]<<endl;
// 	}
// 	outfile.close();
// 	
// 	outroot->cd();
// 	for(int i1=0;i1<3;i1++)
// 	{
// 		hh[i1]->Write();
// 	}
// 	outroot->Close();
// 	
// 	sprintf(hname,"mv PMTIntegrals_%d.root %s/Histos/PMTIntegrals_%d.root",RunNo,AnalysisFilePath,RunNo);system(hname);
// 	sprintf(hname,"mv PMTIntegrals_%d.txt %s/Files/PMTIntegrals_%d.txt",RunNo,AnalysisFilePath,RunNo);system(hname);
// 	for(int i1=0;i1<3;i1++)
// 	{
// 		sprintf(hname,"rm waveforms_ch%d_run%d.txt",i1,RunNo);system(hname);
// 	}
// 	
// 	sprintf(hname,"echo \"Run %d saturations %d %d %d out of %d events\" >> /afs/cern.ch/work/f/flic2019/flicanalyzer/logs/PMTIntLog.txt",RunNo,NSat[0],NSat[1],NSat[2],Ev);
// 	system(hname);
}

bool compareByHitTime(const hits &a, const hits &b)
{
    return a.t < b.t;
}

bool compareByTLTime(const tlist &a, const tlist &b)
{
    return a.t < b.t;
}

void FitTracks()
{
	sprintf(hname,"ValidTrackParameters_%d.txt",RunNo);
	ofstream outfile(hname);
	vector <int> eventlist;
	
	sprintf(hname,"Tracks_%d.root",RunNo);
	TFile* outroot=new TFile(hname,"recreate");
	TTree* T1 =  new TTree("T","T");
	T1->Branch("E",&td.E);
	T1->Branch("StartEndColIndT",&td.StartEndColIndT);
	T1->Branch("FitParams",&td.FitParams);
	T1->Branch("FitNormChi2",&td.FitNormChi2);
	T1->Branch("QColTot",&td.QColTot);
	T1->Branch("NHits",&td.NHits);
	T1->Branch("Nexcl",&td.Nexcl);
	T1->Branch("PMTIntegral",&td.PMTIntegral);
	
	TH1F* NClperEvent=new TH1F("NClperEvent","NClperEvent",60,-0.5,59.5);
	TH1F* ClSize=new TH1F("ClSize","ClSize",400,-0.5,399.5);
	TH1F* MaxClSize=new TH1F("MaxClSize","MaxClSize",400,-0.5,399.5);
	TH1F* FitNormChi2=new TH1F("FitNormChi2","FitNormChi2",200,0.,100.);
	TH1F* FitClSize=new TH1F("FitClSize","FitClSize",200,-0.5,199.5);
	TH1F* NExclHits=new TH1F("NExclHits","NExclHits",200,-0.5,199.5);
	TH1F* TrackDeltaT=new TH1F("TrackDeltaT","TrackDeltaT",500,0,5000);
	TH1F* ColTrackDeltaT=new TH1F("ColTrackDeltaT","ColTrackDeltaT",500,0,5000);
	TH1F* TrackDeltaTSel=new TH1F("TrackDeltaTSel","TrackDeltaTSel",500,0,5000);
	TH1F* AllDeltaT=new TH1F("AllDeltaT","AllDeltaT",500,0,5000);
	TH2F* AllXY=new TH2F("AllXY","AllXY",128,-0.5,127.5,128,-0.5,127.5);
	TH2F* TrackXY=new TH2F("TrackXY","TrackXY",128,-0.5,127.5,128,-0.5,127.5);
	TH2F* TrackEntryXY=new TH2F("TrackEntryXY","TrackEntryXY",128,-0.5,127.5,128,-0.5,127.5);
	TH2F* TrackExitXY=new TH2F("TrackExitXY","TrackExitXY",128,-0.5,127.5,128,-0.5,127.5);
	TH2F* DeltaTvsTrackAngle=new TH2F("DeltaTvsTrackAngle","DeltaTvsTrackAngle",314,0,1.57,500,0,1000);
	TProfile* DeltaTvsTrackAnglePr=new TProfile("DeltaTvsTrackAnglePr","DeltaTvsTrackAnglePr",250,0,1.57,0,1000);
	TH1F* Tfirst=new TH1F("Tfirst","Tfirst",500,0,5000);
	TH1F* Tlast=new TH1F("Tlast","Tlast",500,0,5000);
	TH1F* Trackt0=new TH1F("Trackt0","Trackt0",500,0,5000);
	TH2F* DeltaTvsTrackt0=new TH2F("DeltaTvsTrackt0","DeltaTvsTrackt0",500,0,5000,500,0,1000);
	TH1F* DeltacolperDeltat=new TH1F("DeltacolperDeltat","DeltacolperDeltat",500,0,50);
	TH1F* DeltaindperDeltat=new TH1F("DeltaindperDeltat","DeltaindperDeltat",500,0,50);
	TH1F* DeltaRperDeltat=new TH1F("DeltaRperDeltat","DeltaRperDeltat",2000,0,100);
	TH1F* ColDeltaRR=new TH1F("ColDeltaR","ColDeltaR",2000,0,2000);
	TH1F* HitDeltaRR=new TH1F("HitDeltaRR","HitDeltaRR",2000,0,2000);
	TH1F* HitDeltaT=new TH1F("HitDeltaT","HitDeltaT",2000,0,1000);
	TH1F* HitDeltaR=new TH1F("HitDeltaR","HitDeltaR",200,0,100);
	TH2F* HitDeltaR_vs_DeltaT=new TH2F("HitDeltaR_vs_DeltaT","HitDeltaR_vs_DeltaT",200,-0.5,199.5,200,0,20);
	
	TH2F* DeltaTLengthvsAngle=new TH2F("DeltaTLengthvsAngle","DeltaTLengthvsAngle",314,0,1.57,100,0,2000);
	TH2F* DeltaTLengthvsAngle_norm=new TH2F("DeltaTLengthvsAngle_norm","DeltaTLengthvsAngle_norm",314,0,1.57,100,0,2000);
	
	TH1F* FitNormChi2_Col_vs_t=new TH1F("FitNormChi2_Col_vs_t","FitNormChi2_Col_vs_t",1000,0.,4000.);
	TH1F* FitNormChi2_Ind_vs_t=new TH1F("FitNormChi2_Ind_vs_t","FitNormChi2_Ind_vs_t",1000,0.,4000.);
	
	TH2F* NColHits_vs_fitnormchi2=new TH2F("NColHits_vs_fitnormchi2","NColHits_vs_fitnormchi2",200,-0.5,199.5,500,0,2000);
	TH2F* NIndHits_vs_fitnormchi2=new TH2F("NIndHits_vs_fitnormchi2","NIndHits_vs_fitnormchi2",200,-0.5,199.5,500,0,2000);
	
	TGraph2D* tg2d=new TGraph2D();
	sprintf(hname,"Ev3d");
	tg2d->SetName(hname);tg2d->SetTitle(hname);
	tg2d->SetMarkerStyle(20);tg2d->SetMarkerColor(2);
	
	TGraph2D* tg2dA=new TGraph2D();
	sprintf(hname,"Ev3dout");
	tg2dA->SetName(hname);tg2dA->SetTitle(hname);
	tg2dA->SetMarkerStyle(20);tg2dA->SetMarkerColor(4);
	
	TGraph* tg1dcol=new TGraph();
	sprintf(hname,"Col_vs_t");
	tg1dcol->SetName(hname);tg1dcol->SetTitle(hname);
	
	TGraph* tg1dind=new TGraph();
	sprintf(hname,"Ind_vs_t");
	tg1dind->SetName(hname);tg1dind->SetTitle(hname);
	
	TF1* tflin=new TF1("lin","[0]+[1]*x",750,2500);
	TF1* tflin2=new TF1("lin2","[0]+[1]*x",750,2500);
	
	TCanvas* cc=new TCanvas("cc","cc",600,600);
	TGraph2D* tg2daxis=new TGraph2D();
	tg2daxis->SetName("Axis");tg2daxis->SetTitle("Axis");
	tg2daxis->SetPoint(0,0,0,0);tg2daxis->SetPoint(0,127,127,4095);
	
	sprintf(hname,"%s/Histos/Hits_%d.root",AnalysisFilePath,RunNo);
	TFile* inroot=new TFile(hname);
	TTree* T =  (TTree*) inroot->Get("T");
	T->SetBranchAddress("E",&hd.E);
	T->SetBranchAddress("ColIndT",&hd.ColIndT);
	T->SetBranchAddress("Int",&hd.Int);
	T->SetBranchAddress("QColTot",&hd.QColTot);
	T->SetBranchAddress("QColTotZS",&hd.QColTotZS);
	T->SetBranchAddress("PMTIntegral",&hd.PMTIntegral);
	T->SetBranchAddress("ColID",&hd.ColID);
	T->SetBranchAddress("ColT",&hd.ColT);
	T->SetBranchAddress("ColA",&hd.ColA);
	T->SetBranchAddress("ColInt",&hd.ColInt);
	T->SetBranchAddress("Colw",&hd.Colw);
	T->SetBranchAddress("IndID",&hd.IndID);
	T->SetBranchAddress("IndT",&hd.IndT);
	T->SetBranchAddress("IndA",&hd.IndA);
	T->SetBranchAddress("Indw",&hd.Indw);
	
// 	sprintf(hname,"%s/Run_%d.root",NTupleFilePath,RunNo);
// 	TFile* inroot=new TFile(hname);
// 	TTree* TN =  (TTree*) inroot->Get("T");
// 	TN->SetBranchAddress("TPCWF",&ed.TPCWF);
// 	TN->SetBranchAddress("PMTWF",&ed.PMTWF);
	
// 	ReadTPCBaselines();
	
	bool found=false;
	bool skipEvent=false;
	vector <hits> H;
	vector <int> inCL;
	int b[7]={0};
	int Ev=0;
	bool foundcl=false;
	int nid=0;
	int indCLHP=0;int maxCLsize=0;
	float fitnormchi2=0;float fitnormchi2col=0;float fitnormchi2ind=0;
	TPolyLine3D *l;
	int tmin=5000;int tmax=0;
	int mincol=0;int maxcol=0;
	float tmincol=5000;float tmaxcol=0;
	float tmintrack=5000;float tmaxtrack=0;
	float tminall=5000;float tmaxall=0;
	int indoftmin=0;int indoftmax=0;
	TH3F *frame3d=new TH3F("frame3d","frame3d",1,-0.5,127.5,1,-0.5,127.5,1,-0.5,4095.5);
	TPolyMarker3D *pm3d1;
	TPolyMarker3D *pm3d2;
	double x1(0),y1(0),z1(0),t1(0);
	double x2(0),y2(0),z2(0),t2(0);
	double tflinpars[2][2]={{0}};//col, ind - const, slope
	int Nvalid=0;
	float QColtot=0;
	
// 	for(int I=0;I<T->GetEntries();I++)
// 	{
// 		T->GetEntry(I);
// 		
// 		if(hd.ColIndT->size()>300) continue;
// 		for(int i1=0;i1<hd.ColIndT->size();i1++)
// 		{
// 			for(int i2=i1+1;i2<hd.ColIndT->size();i2++)
// 			{
// 				DeltaRperDeltat->Fill(fabs(sqrt(pow(hd.ColIndT->at(i1)[0]-hd.ColIndT->at(i2)[0],2)+pow(hd.ColIndT->at(i1)[1]-hd.ColIndT->at(i2)[1],2))/((float)(hd.ColIndT->at(i1)[2]-hd.ColIndT->at(i2)[2]))));
// 				HitDeltaT->Fill(((float)(hd.ColIndT->at(i1)[2]-hd.ColIndT->at(i2)[2])));
// 				HitDeltaR_vs_DeltaT->Fill(((float)(hd.ColIndT->at(i1)[2]-hd.ColIndT->at(i2)[2])),fabs(sqrt(pow(hd.ColIndT->at(i1)[0]-hd.ColIndT->at(i2)[0],2)+pow(hd.ColIndT->at(i1)[1]-hd.ColIndT->at(i2)[1],2))));
// 				HitDeltaR->Fill(fabs(sqrt(pow(hd.ColIndT->at(i1)[0]-hd.ColIndT->at(i2)[0],2)+pow(hd.ColIndT->at(i1)[1]-hd.ColIndT->at(i2)[1],2))));
// 				HitDeltaRR->Fill(fabs(sqrt(pow(hd.ColIndT->at(i1)[0]-hd.ColIndT->at(i2)[0],2)+pow(hd.ColIndT->at(i1)[1]-hd.ColIndT->at(i2)[1],2)+pow(hd.ColIndT->at(i1)[2]-hd.ColIndT->at(i2)[2],2))));
// 			}
// 		}
// 		for(int i1=0;i1<hd.ColID->size();i1++)
// 		{
// 			for(int i2=i1+1;i2<hd.ColID->size();i2++)
// 			{
// 				ColDeltaRR->Fill(fabs(sqrt(pow(hd.ColID->at(i1)-hd.ColID->at(i2),2)+pow(hd.ColT->at(i1)-hd.ColT->at(i2),2))));
// 			}
// 		}
// 	}
	
	
// 	for(int I=0;I<T->GetEntries();I++)
// 	{
// 		T->GetEntry(I);
// 		
// 		if(hd.ColIndT->size()<5) continue;
// 		if(hd.ColIndT->size()>300) continue;
// 		
// 		vector <clusters> CL;
// 		cl.hi.clear();
// 		inCL.clear();
// 		for(int i1=0;i1<hd.ColIndT->size();i1++)
// 		{
// 			inCL.push_back(0);
// 		}
// 		tmin=5000;tmax=0;
// 		for(int i1=0;i1<hd.ColIndT->size();i1++)
// 		{
// 			if(hd.ColIndT->at(i1)[2]>tmax) tmax=hd.ColIndT->at(i1)[2];
// 			if(hd.ColIndT->at(i1)[2]<tmin) tmin=hd.ColIndT->at(i1)[2];
// 			if(inCL[i1]!=0) continue;
// 			cl.hi.clear();
// 			cl.hi.push_back(i1);
// 			inCL[i1]=1;
// 			nid=1;
// 			while(nid>0)
// 			{
// 				nid=0;
// 				for(int i2=0;i2<hd.ColIndT->size();i2++)
// 				{
// 					if(inCL[i2]!=0) continue;
// 					for(int i3=0;i3<cl.hi.size();i3++)
// 					{
// // 						if(sqrt(pow(H[cl.hi[i3]].col-H[i2].col,2)+pow(H[cl.hi[i3]].ind-H[i2].ind,2))<10 && (H[cl.hi[i3]].t-H[i2].t)<20)
// // 						if(sqrt(pow(hd.ColIndT->at(i1)[0]-hd.ColIndT->at(i2)[0],2)+pow(hd.ColIndT->at(i1)[1]-hd.ColIndT->at(i2)[1],2))<10 && (hd.ColIndT->at(i1)[2]-hd.ColIndT->at(i2)[2])<20)
// // 						if(fabs(sqrt(pow(hd.ColIndT->at(i1)[0]-hd.ColIndT->at(i2)[0],2)+pow(hd.ColIndT->at(i1)[1]-hd.ColIndT->at(i2)[1],2)))<5 && abs(hd.ColIndT->at(i1)[2]-hd.ColIndT->at(i2)[2])<10)
// 						if(fabs(sqrt(pow(hd.ColIndT->at(i1)[0]-hd.ColIndT->at(i2)[0],2)+pow(hd.ColIndT->at(i1)[1]-hd.ColIndT->at(i2)[1],2)))<3)
// // 						if((sqrt(pow(hd.ColIndT->at(i1)[0]-hd.ColIndT->at(i2)[0],2)+pow(hd.ColIndT->at(i1)[1]-hd.ColIndT->at(i2)[1],2))/(hd.ColIndT->at(i1)[2]-hd.ColIndT->at(i2)[2]))<=3)
// // 						if(fabs(sqrt(pow(hd.ColIndT->at(i1)[0]-hd.ColIndT->at(i2)[0],2)+pow(hd.ColIndT->at(i1)[1]-hd.ColIndT->at(i2)[1],2))/((float)(hd.ColIndT->at(i1)[2]-hd.ColIndT->at(i2)[2])))<1.5)
// 						{
// 							cl.hi.push_back(i2);
// 							nid++;
// 							inCL[i2]=1;
// 							break;
// 						}
// 					}
// 				}
// 			}
// 			CL.push_back(cl);
// 			cl.hi.clear();
// 		}
// 		AllDeltaT->Fill(tmax-tmin);
// 		
// 		int indbestcl=-1;int nhitsintr=0;
// 		for(int i2=0;i2<CL.size();i2++)
// 		{
// 			if(CL[i2].hi.size()<5) continue;
// 			tg2d->Set(0);
// 			tmin=5000;tmax=0;
// 			for(int i1=0;i1<CL[i2].hi.size();i1++)
// 			{
// 				tg2d->SetPoint(tg2d->GetN(),hd.ColIndT->at(CL[i2].hi[i1])[0],hd.ColIndT->at(CL[i2].hi[i1])[1],hd.ColIndT->at(CL[i2].hi[i1])[2]);
// 				if(hd.ColIndT->at(CL[i2].hi[i1])[2]<tmin){tmin=hd.ColIndT->at(CL[i2].hi[i1])[2];indoftmin=i1;}
// 				if(hd.ColIndT->at(CL[i2].hi[i1])[2]>tmax){tmax=hd.ColIndT->at(CL[i2].hi[i1])[2];indoftmax=i1;}
// 			}
// 		
// 			pStart[0]=hd.ColIndT->at(CL[i2].hi[indoftmin])[0];
// 			pStart[1]=((float)(hd.ColIndT->at(CL[i2].hi[indoftmax])[0]-hd.ColIndT->at(CL[i2].hi[indoftmin])[0]))/((float)(hd.ColIndT->at(CL[i2].hi[indoftmax])[2]-hd.ColIndT->at(CL[i2].hi[indoftmin])[2]));
// 			pStart[2]=hd.ColIndT->at(CL[i2].hi[indoftmin])[1];
// 			pStart[3]=((float)(hd.ColIndT->at(CL[i2].hi[indoftmax])[1]-hd.ColIndT->at(CL[i2].hi[indoftmin])[1]))/((float)(hd.ColIndT->at(CL[i2].hi[indoftmax])[2]-hd.ColIndT->at(CL[i2].hi[indoftmin])[2]));
// 			pStart[4]=1500;
// 		
// 			fitnormchi2=line3Dfit(tg2d);
// 		
// 			for(int i1=0;i1<CL.size();i1++)
// 			{
// 				if(i1==i2) continue;
// 				for(int i3=0;i3<CL[i1].hi.size();i3++)
// 				{
// 					x1=fitParams[0]+fitParams[1]*(hd.ColIndT->at(CL[i1].hi[i3])[2]-fitParams[4]);
// 					y1=fitParams[2]+fitParams[3]*(hd.ColIndT->at(CL[i1].hi[i3])[2]-fitParams[4]);
// 					if(fabs(sqrt(pow(x1-hd.ColIndT->at(CL[i1].hi[i3])[0],2.)+pow(y1-hd.ColIndT->at(CL[i1].hi[i3])[1],2.)))<10.)
// 					{
// 						tg2d->SetPoint(tg2d->GetN(),hd.ColIndT->at(CL[i1].hi[i3])[0],hd.ColIndT->at(CL[i1].hi[i3])[1],hd.ColIndT->at(CL[i1].hi[i3])[2]);
// 					}
// 				}
// 			}
// 			if(tg2d->GetN()>nhitsintr)
// 			{
// 				nhitsintr=tg2d->GetN();
// 				indbestcl=i2;
// 			}
// 		}
// 		if(indbestcl<0) continue;
// 		tg2d->Set(0);
// 		tmin=5000;tmax=0;
// 		for(int i1=0;i1<CL[indbestcl].hi.size();i1++)
// 		{
// 			tg2d->SetPoint(tg2d->GetN(),hd.ColIndT->at(CL[indbestcl].hi[i1])[0],hd.ColIndT->at(CL[indbestcl].hi[i1])[1],hd.ColIndT->at(CL[indbestcl].hi[i1])[2]);
// 			if(hd.ColIndT->at(CL[indbestcl].hi[i1])[2]<tmin){tmin=hd.ColIndT->at(CL[indbestcl].hi[i1])[2];indoftmin=i1;}
// 			if(hd.ColIndT->at(CL[indbestcl].hi[i1])[2]>tmax){tmax=hd.ColIndT->at(CL[indbestcl].hi[i1])[2];indoftmax=i1;}
// 		}
// 	
// 		pStart[0]=hd.ColIndT->at(CL[indbestcl].hi[indoftmin])[0];
// 		pStart[1]=((float)(hd.ColIndT->at(CL[indbestcl].hi[indoftmax])[0]-hd.ColIndT->at(CL[indbestcl].hi[indoftmin])[0]))/((float)(hd.ColIndT->at(CL[indbestcl].hi[indoftmax])[2]-hd.ColIndT->at(CL[indbestcl].hi[indoftmin])[2]));
// 		pStart[2]=hd.ColIndT->at(CL[indbestcl].hi[indoftmin])[1];
// 		pStart[3]=((float)(hd.ColIndT->at(CL[indbestcl].hi[indoftmax])[1]-hd.ColIndT->at(CL[indbestcl].hi[indoftmin])[1]))/((float)(hd.ColIndT->at(CL[indbestcl].hi[indoftmax])[2]-hd.ColIndT->at(CL[indbestcl].hi[indoftmin])[2]));
// 		pStart[4]=1500;
// 	
// 		fitnormchi2=line3Dfit(tg2d);
// 	
// 		for(int i1=0;i1<CL.size();i1++)
// 		{
// 			if(i1==indbestcl) continue;
// 			for(int i2=0;i2<CL[i1].hi.size();i2++)
// 			{
// 				x1=fitParams[0]+fitParams[1]*(hd.ColIndT->at(CL[i1].hi[i2])[2]-fitParams[4]);
// 				y1=fitParams[2]+fitParams[3]*(hd.ColIndT->at(CL[i1].hi[i2])[2]-fitParams[4]);
// 				if(fabs(sqrt(pow(x1-hd.ColIndT->at(CL[i1].hi[i2])[0],2.)+pow(y1-hd.ColIndT->at(CL[i1].hi[i2])[1],2.)))<10.)
// 				{
// 					CL[indbestcl].hi.push_back(CL[i1].hi[i2]);
// 					
// 					tg2d->SetPoint(tg2d->GetN(),hd.ColIndT->at(CL[i1].hi[i2])[0],hd.ColIndT->at(CL[i1].hi[i2])[1],hd.ColIndT->at(CL[i1].hi[i2])[2]);
// 					if(hd.ColIndT->at(CL[i1].hi[i2])[2]<tmin){tmin=hd.ColIndT->at(CL[i1].hi[i2])[2];indoftmin=CL[indbestcl].hi.size()-1;}
// 					if(hd.ColIndT->at(CL[i1].hi[i2])[2]>tmax){tmax=hd.ColIndT->at(CL[i1].hi[i2])[2];indoftmax=CL[indbestcl].hi.size()-1;}
// 					TrackXY->Fill(hd.ColIndT->at(CL[i1].hi[i2])[0],hd.ColIndT->at(CL[i1].hi[i2])[1]);
// 					QColtot+=hd.QColTot;
// 					
// 					CL[i1].hi[i2]=CL[i1].hi[CL[i1].hi.size()-1];
// 					CL[i1].hi.pop_back();
// 					i2--;
// 				}
// 			}
// 		}
// 		
// // 		tmin=5000;tmax=0;
// // 		QColtot=0;
// // 		tg2d->Set(0);
// // 		
// // 		for(int i1=0;i1<CL[indCLHP].hi.size();i1++)
// // 		{
// // 			tg2d->SetPoint(tg2d->GetN(),hd.ColIndT->at(CL[indCLHP].hi[i1])[0],hd.ColIndT->at(CL[indCLHP].hi[i1])[1],hd.ColIndT->at(CL[indCLHP].hi[i1])[2]);
// // 			if(hd.ColIndT->at(CL[indCLHP].hi[i1])[2]<tmin){tmin=hd.ColIndT->at(CL[indCLHP].hi[i1])[2];indoftmin=i1;}
// // 			if(hd.ColIndT->at(CL[indCLHP].hi[i1])[2]>tmax){tmax=hd.ColIndT->at(CL[indCLHP].hi[i1])[2];indoftmax=i1;}
// // 			TrackXY->Fill(hd.ColIndT->at(CL[indCLHP].hi[i1])[0],hd.ColIndT->at(CL[indCLHP].hi[i1])[1]);
// // 			QColtot+=hd.QColTot;
// // 		}
// // 		
// // 		pStart[0]=hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[0];
// // 		pStart[1]=((float)(hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[0]-hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[0]))/((float)(hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[2]-hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[2]));
// // 		pStart[2]=hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[1];
// // 		pStart[3]=((float)(hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[1]-hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[1]))/((float)(hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[2]-hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[2]));
// // 		pStart[4]=1500;
// // 		
// // 		fitnormchi2=line3Dfit(tg2d);
// 		
// // 		cout<<I<<endl;
// // 		cout<<"before"<<endl;
// // 		for(int i1=0;i1<CL.size();i1++)
// // 		{
// // 			cout<<i1<<" "<<CL[i1].hi.size()<<endl;
// // 		}
// 		
// // 		if(I<10)
// // 		{
// // 			sprintf(hname,"TrackHits_%d",I);
// // 			tg2d->SetName(hname);tg2d->SetTitle(hname);
// // 			outroot->cd();
// // 			tg2d->Write();
// // 		}
// 		
// // 		for(int i1=0;i1<CL.size();i1++)
// // 		{
// // 			if(i1==indCLHP) continue;
// // 			for(int i2=0;i2<CL[i1].hi.size();i2++)
// // 			{
// // 				x1=fitParams[0]+fitParams[1]*(hd.ColIndT->at(CL[i1].hi[i2])[2]-fitParams[4]);
// // 				y1=fitParams[2]+fitParams[3]*(hd.ColIndT->at(CL[i1].hi[i2])[2]-fitParams[4]);
// // 				if(fabs(sqrt(pow(x1-hd.ColIndT->at(CL[i1].hi[i2])[0],2.)+pow(y1-hd.ColIndT->at(CL[i1].hi[i2])[1],2.)))<10.)
// // 				{
// // 					CL[indCLHP].hi.push_back(CL[i1].hi[i2]);
// // 					
// // 					tg2d->SetPoint(tg2d->GetN(),hd.ColIndT->at(CL[i1].hi[i2])[0],hd.ColIndT->at(CL[i1].hi[i2])[1],hd.ColIndT->at(CL[i1].hi[i2])[2]);
// // 					if(hd.ColIndT->at(CL[i1].hi[i2])[2]<tmin){tmin=hd.ColIndT->at(CL[i1].hi[i2])[2];indoftmin=CL[indCLHP].hi.size()-1;}
// // 					if(hd.ColIndT->at(CL[i1].hi[i2])[2]>tmax){tmax=hd.ColIndT->at(CL[i1].hi[i2])[2];indoftmax=CL[indCLHP].hi.size()-1;}
// // 					TrackXY->Fill(hd.ColIndT->at(CL[i1].hi[i2])[0],hd.ColIndT->at(CL[i1].hi[i2])[1]);
// // 					QColtot+=hd.QColTot;
// // 					
// // 					CL[i1].hi[i2]=CL[i1].hi[CL[i1].hi.size()-1];
// // 					CL[i1].hi.pop_back();
// // 					i2--;
// // 				}
// // 			}
// // 		}
// 		
// 		
// // 		cout<<"after"<<endl;
// // 		for(int i1=0;i1<CL.size();i1++)
// // 		{
// // 			cout<<i1<<" "<<CL[i1].hi.size()<<endl;
// // 		}
// 		for(int i1=0;i1<CL.size();i1++)
// 		{
// 			if(CL[i1].hi.size()==0)
// 			{
// 				CL[i1].hi=CL[CL.size()-1].hi;
// 				CL.pop_back();
// 				i1--;
// 			}
// 		}
// 		
// 		int HitsinCLs1=0;int HitsinCLs2=0;
// 		indCLHP=-1;maxCLsize=0;
// 		for(int i1=0;i1<CL.size();i1++)
// 		{
// 			if(CL[i1].hi.size()>maxCLsize){maxCLsize=CL[i1].hi.size();indCLHP=i1;}
// 			HitsinCLs1+=CL[i1].hi.size();
// 		}
// 		
// 		indCLHP=-1;maxCLsize=0;
// 		for(int i1=0;i1<CL.size();i1++)
// 		{
// 			ClSize->Fill(CL[i1].hi.size());
// 			if(CL[i1].hi.size()>maxCLsize){maxCLsize=CL[i1].hi.size();indCLHP=i1;}
// 			HitsinCLs2+=CL[i1].hi.size();
// 		}
// // 		cout<<I<<" "<<hd.ColIndT->size()<<" "<<HitsinCLs1<<" "<<HitsinCLs2<<endl;
// 		
// // 		cout<<"after 2"<<endl;
// // 		for(int i1=0;i1<CL.size();i1++)
// // 		{
// // 			cout<<i1<<" "<<CL[i1].hi.size()<<endl;
// // 		}
// // 		cout<<endl;
// 		
// 		maxCLsize=CL[indCLHP].hi.size();
// 		MaxClSize->Fill(maxCLsize);
// 		NClperEvent->Fill(CL.size());
// 		
// 		pStart[0]=hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[0];
// 		pStart[1]=((float)(hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[0]-hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[0]))/((float)(hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[2]-hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[2]));
// 		pStart[2]=hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[1];
// 		pStart[3]=((float)(hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[1]-hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[1]))/((float)(hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[2]-hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[2]));
// 		pStart[4]=1500;
// 		
// 		fitnormchi2=line3Dfit(tg2d);
// 		
// 		TrackDeltaT->Fill(tmax-tmin);
// 		FitNormChi2->Fill(fitnormchi2);
// 		
// // 		td.E=I;
// 		td.E=hd.E;
// 		td.StartEndColIndT->clear();
// 		td.StartEndColIndT->push_back(hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[0]);
// 		td.StartEndColIndT->push_back(hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[1]);
// 		td.StartEndColIndT->push_back(hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[2]);
// 		td.StartEndColIndT->push_back(hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[0]);
// 		td.StartEndColIndT->push_back(hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[1]);
// 		td.StartEndColIndT->push_back(hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[2]);
// 		td.FitParams->clear();
// 		for(int is2=0;is2<5;is2++)
// 		{
// 			td.FitParams->push_back(fitParams[is2]);
// 		}
// 		td.FitNormChi2=fitnormchi2;
// 		td.QColTot=QColtot;
// 		td.NHits=CL[indCLHP].hi.size();
// 		td.Nexcl=hd.ColIndT->size()-CL[indCLHP].hi.size();
// 		td.PMTIntegral->clear();
// 		for(int is2=0;is2<3;is2++)
// 		{
// 			td.PMTIntegral->push_back(hd.PMTIntegral->at(is2));
// 		}
// 		
// 		T1->Fill();
// 		
// 		if(fitnormchi2<5.)
// 		{
// // // 			cc->cd(); 
// // // 			frame3d->Draw();
// // // 			
// // // 			pm3d1 = new TPolyMarker3D();
// // // 			pm3d2 = new TPolyMarker3D();
// // // 			
// // // // 			pm3d1 = new TPolyMarker3D(CL[indCLHP].hi.size());
// // // // // 			pm3d2 = new TPolyMarker3D(H.size()-CL[indCLHP].hi.size());
// // // // 			pm3d2 = new TPolyMarker3D(hd.ColIndT->size()-CL[indCLHP].hi.size());
// // // 			
// // // 			for(int i1=0;i1<CL[indCLHP].hi.size();i1++)
// // // 			{
// // // // 				pm3d1->SetPoint(i1,H[CL[indCLHP].hi[i1]].col,H[CL[indCLHP].hi[i1]].ind,H[CL[indCLHP].hi[i1]].t);
// // // 				pm3d1->SetPoint(i1,hd.ColIndT->at(CL[indCLHP].hi[i1])[0],hd.ColIndT->at(CL[indCLHP].hi[i1])[1],hd.ColIndT->at(CL[indCLHP].hi[i1])[2]);
// // // 			}
// // // 			
// // // 			int na=0;
// // // 			for(int i1=0;i1<CL.size();i1++)
// // // 			{
// // // 				if(i1==indCLHP) continue;
// // // 				for(int i2=0;i2<CL[i1].hi.size();i2++)
// // // 				{
// // // // 					pm3d2->SetPoint(na,H[CL[i1].hi[i2]].col,H[CL[i1].hi[i2]].ind,H[CL[i1].hi[i2]].t);
// // // 					pm3d2->SetPoint(na,hd.ColIndT->at(CL[i1].hi[i2])[0],hd.ColIndT->at(CL[i1].hi[i2])[1],hd.ColIndT->at(CL[i1].hi[i2])[2]);
// // // 					na++;
// // // 				}
// // // 			}
// // // // 					pm3d1->SetMarkerSize(2);
// // // 			pm3d1->SetMarkerColor(kRed);
// // // 			pm3d1->SetMarkerStyle(24);   
// // // 			pm3d1->Draw();
// // // // 					pm3d2->SetMarkerSize(2);
// // // 			pm3d2->SetMarkerColor(kBlue);
// // // 			pm3d2->SetMarkerStyle(24);   
// // // 			pm3d2->Draw();
// // // 			
// // // 			l = new TPolyLine3D(2);
// 			
// 			x1=fitParams[0]+fitParams[1]*(tmin-fitParams[4]);
// 			y1=fitParams[2]+fitParams[3]*(tmin-fitParams[4]);
// 			z1=tmin;
// // // 			l->SetPoint(0,x1,y1,z1);
// 			
// 			x2=fitParams[0]+fitParams[1]*(tmax-fitParams[4]);
// 			y2=fitParams[2]+fitParams[3]*(tmax-fitParams[4]);
// 			z2=tmax;
// // // 			l->SetPoint(1,x2,y2,z2);
// // // 			l->SetLineColor(kBlack);
// // // 			l->SetLineWidth(2);
// // // 			l->Draw("same");
// 			
// 			DeltaTvsTrackAngle->Fill(atan(sqrt(pow(x2-x1,2)+pow(y2-y1,2))/(tmax-tmin)),(tmax-tmin));
// 			DeltaTvsTrackAnglePr->Fill(atan(sqrt(pow(x2-x1,2)+pow(y2-y1,2))/(tmax-tmin)),(tmax-tmin));
// 			Tfirst->Fill(tmin);
// 			Tlast->Fill(tmax);
// //					if(tmin<=1160) Tlast->Fill(tmax);
// 			Trackt0->Fill(fitParams[4]);
// 			TrackEntryXY->Fill(x1,y1);
// 			TrackExitXY->Fill(x2,y2);
// 			if(tmin<=1160 && tmax>=1920) TrackDeltaTSel->Fill(tmax-tmin);
// 			DeltaTvsTrackt0->Fill(fitParams[4],tmax-tmin);
// 			DeltaTLengthvsAngle->Fill(atan(sqrt(pow(x2-x1,2)+pow(y2-y1,2))/(tmax-tmin)),sqrt(pow(x2-x1,2)+pow(y2-y1,2)+pow(tmax-tmin,2)),tmax-tmin);
// 			DeltaTLengthvsAngle_norm->Fill(atan(sqrt(pow(x2-x1,2)+pow(y2-y1,2))/(tmax-tmin)),sqrt(pow(x2-x1,2)+pow(y2-y1,2)+pow(tmax-tmin,2)));
// 			
// //					cout<<Ev<<" "<<tmin<<" "<<tmax<<" "<<(tmax-tmin)<<endl;
// 			
// // 			sprintf(hname,"Ev_%d_%5.3f",td.E,fitnormchi2);
// // 			cc->SetName(hname);cc->SetTitle(hname);
// // 			outroot->cd();
// // 			cc->Write();
// // 			tg2dA->Set(0);
// // 			delete l;
// 			
// 			FitClSize->Fill(CL[indCLHP].hi.size());
// 			NExclHits->Fill(hd.ColIndT->size()-CL[indCLHP].hi.size());
// // 			if(CL[indCLHP].hi.size()>15 && (H.size()-CL[indCLHP].hi.size())<10)
// // 			{
// // 				Nvalid++;
// // 				outfile<<Ev<<" "<<fitParams[0]<<" "<<fitParams[1]<<" "<<fitParams[2]<<" "<<fitParams[3]<<" "<<fitParams[4]<<" "<<fitnormchi2<<" "<<CL[indCLHP].hi.size()<<endl;
// // 			}
// // 			if(td.NHits>15 && td.Nexcl<10)
// 			if(td.NHits>10 && td.Nexcl<10)
// 			{
// 				outfile<<td.E<<" "<<hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[0]<<" "<<hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[1]<<" "<<hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[2]<<" "<<hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[0]<<" "<<hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[1]<<" "<<hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[2]<<endl;
// 			}
// 		}
// 		tg2d->Set(0);
// 		inCL.clear();
// 		cl.hi.clear();
// 		CL.clear();
// // 		if(I==10) break;
// 		if(I%100==0) cout<<I<<" / "<<T->GetEntries()<<endl;
// 	}
	
// 	for(int I=0;I<T->GetEntries();I++)
// 	{
// 		T->GetEntry(I);
// 		
// 		if(hd.ColIndT->size()>300) continue;
// 		for(int i1=0;i1<hd.ColIndT->size();i1++)
// 		{
// 			for(int i2=i1+1;i2<hd.ColIndT->size();i2++)
// 			{
// 				DeltaRperDeltat->Fill(fabs(sqrt(pow(hd.ColIndT->at(i1)[0]-hd.ColIndT->at(i2)[0],2)+pow(hd.ColIndT->at(i1)[1]-hd.ColIndT->at(i2)[1],2))/((float)(hd.ColIndT->at(i1)[2]-hd.ColIndT->at(i2)[2]))));
// 			}
// 		}
// 	}
	
// 	for(int I=0;I<T->GetEntries();I++)
// 	{
// 		T->GetEntry(I);
// 		if(hd.ColIndT->size()<5) continue;
// 		
// 		for(int i1=0;i1<hd.ColIndT->size();i1++)
// 		{
// 			for(int i2=i1+1;i2<hd.ColIndT->size();i2++)
// 			{
// 				DeltaRperDeltat->Fill(fabs((sqrt(pow(hd.ColIndT->at(i1)[0]-hd.ColIndT->at(i2)[0],2.)+pow(hd.ColIndT->at(i1)[1]-hd.ColIndT->at(i2)[1],2.))/((float)(hd.ColIndT->at(i1)[2]-hd.ColIndT->at(i2)[2])))));
// 			}
// 		}
// 		if(I%100==0) cout<<I<<" / "<<T->GetEntries()<<endl;
// 	}
	
	
	for(int I=0;I<T->GetEntries();I++)
	{
		T->GetEntry(I);
		if(hd.ColIndT->size()<5) continue;
		if(hd.ColIndT->size()>300) continue;
		
		vector <clusters> CL;
		cl.hi.clear();
		inCL.clear();
		for(int i1=0;i1<hd.ColID->size();i1++)
		{
			inCL.push_back(0);
		}
		tminall=5000;tmaxall=0;
		for(int i1=0;i1<hd.ColID->size();i1++)
		{
			if(hd.ColT->at(i1)>tmaxall) tmaxall=hd.ColT->at(i1);
			if(hd.ColT->at(i1)<tminall) tminall=hd.ColT->at(i1);
			if(inCL[i1]!=0) continue;
			cl.hi.clear();
			cl.hi.push_back(i1);
			inCL[i1]=1;
			nid=1;
			while(nid>0)
			{
				nid=0;
				for(int i2=0;i2<hd.ColID->size();i2++)
				{
					if(inCL[i2]!=0) continue;
					for(int i3=0;i3<cl.hi.size();i3++)
					{
// 						if(abs(hd.ColID->at(i2)-hd.ColID->at(cl.hi[i3]))<3 && abs(hd.ColT->at(i2)-hd.ColT->at(cl.hi[i3]))<10)
// 						if(abs(hd.ColID->at(i2)-hd.ColID->at(cl.hi[i3]))<=1 && abs(hd.ColT->at(i2)-hd.ColT->at(cl.hi[i3]))<=100)
						if(abs(hd.ColID->at(i2)-hd.ColID->at(cl.hi[i3]))<=1)
// 						if(fabs(sqrt(pow(hd.ColID->at(i2)-hd.ColID->at(cl.hi[i3]),2)+pow(hd.ColT->at(i2)-hd.ColT->at(cl.hi[i3]),2)))<140)
// 						if(fabs(((float)hd.ColIndT->at(i1)[0]-hd.ColIndT->at(i2)[0])/((float)hd.ColIndT->at(i1)[2]-hd.ColIndT->at(i2)[2]))<1.5)
						{
							cl.hi.push_back(i2);
							nid++;
							inCL[i2]=1;
							break;
						}
					}
				}
			}
			CL.push_back(cl);
			cl.hi.clear();
		}
		AllDeltaT->Fill(tmaxall-tminall);
		int HitsinCLs1=0;int HitsinCLs2=0;
		indCLHP=-1;maxCLsize=0;
		for(int i1=0;i1<CL.size();i1++)
		{
			if(CL[i1].hi.size()>maxCLsize){maxCLsize=CL[i1].hi.size();indCLHP=i1;}
			HitsinCLs1+=CL[i1].hi.size();
		}
		tg1dcol->Set(0);tg1dind->Set(0);
		
		tmincol=5000;tmaxcol=0;
		mincol=200;maxcol=0;
		for(int i1=0;i1<CL[indCLHP].hi.size();i1++)
		{
			tg1dcol->SetPoint(tg1dcol->GetN(),hd.ColT->at(CL[indCLHP].hi[i1]),hd.ColID->at(CL[indCLHP].hi[i1]));
// 			tg1dind->SetPoint(tg1dind->GetN(),hd.ColIndT->at(CL[indCLHP].hi[i1])[2],hd.ColIndT->at(CL[indCLHP].hi[i1])[1]);
			if(hd.ColT->at(CL[indCLHP].hi[i1])>tmaxcol) {tmaxcol=hd.ColT->at(CL[indCLHP].hi[i1]);indoftmax=i1;}
			if(hd.ColT->at(CL[indCLHP].hi[i1])<tmincol) {tmincol=hd.ColT->at(CL[indCLHP].hi[i1]);indoftmin=i1;}
			if(hd.ColID->at(CL[indCLHP].hi[i1])>maxcol) {maxcol=hd.ColID->at(CL[indCLHP].hi[i1]);}
			if(hd.ColID->at(CL[indCLHP].hi[i1])<mincol) {mincol=hd.ColID->at(CL[indCLHP].hi[i1]);}
		}
		ColTrackDeltaT->Fill(tmaxcol-tmincol);
		
		tflin->SetParameter(0,mincol);
		tflin->SetParameter(1,((float)(maxcol-mincol))/(tmaxcol-tmincol));
		tg1dcol->Fit(tflin,"q","q",tmincol-1,tmaxcol+1);
		fitnormchi2col=tflin->GetChisquare()/tflin->GetNDF();
		tflinpars[0][0]=tflin->GetParameter(0);tflinpars[0][1]=tflin->GetParameter(1);
// 		tg1dind->Fit(tflin,"q","q",tmin-1,tmax+1);
// 		fitnormchi2ind=tflin->GetChisquare()/tflin->GetNDF();
// 		tflinpars[1][0]=tflin->GetParameter(0);tflinpars[1][1]=tflin->GetParameter(1);
		
		if(I<10)
		{
			sprintf(hname,"Col_vs_t_%d",I);
			tg1dcol->SetName(hname);tg1dcol->SetTitle(hname);
			tg1dcol->SetMarkerStyle(20);
			outroot->cd();
			tg1dcol->Write();
		}
		
// 		for(int i1=0;i1<tg1dcol->GetN();i1++)
// 		{
// 			tg1dcol->GetPoint(i1,x1,y1);
// 			if(fabs(tflinpars[0][0]+tflinpars[0][1]*x1-y1)>10)
// 			{
// 				tg1dcol->RemovePoint(i1);
// 				i1--;
// 			}
// 		}
// 		for(int i1=0;i1<tg1dind->GetN();i1++)
// 		{
// 			tg1dind->GetPoint(i1,x1,y1);
// 			if(fabs(tflinpars[1][0]+tflinpars[1][1]*x1-y1)>10)
// 			{
// 				tg1dind->RemovePoint(i1);
// 				i1--;
// 			}
// 		}
// 		
// 		tg1dcol->Fit(tflin,"q","q",tmin-1,tmax+1);
// 		fitnormchi2col=tflin->GetChisquare()/tflin->GetNDF();
// 		tflinpars[0][0]=tflin->GetParameter(0);tflinpars[0][1]=tflin->GetParameter(1);
// 		tg1dind->Fit(tflin,"q","q",tmin-1,tmax+1);
// 		fitnormchi2ind=tflin->GetChisquare()/tflin->GetNDF();
// 		tflinpars[1][0]=tflin->GetParameter(0);tflinpars[1][1]=tflin->GetParameter(1);
		
		
// 		bool nd=true;
// 		while(nd)
// 		{
// 			nd=false;
// 			tg1dcol->Fit(tflin,"q","q",tmin-1,tmax+1);
// 			tflinpars[0][0]=tflin->GetParameter(0);tflinpars[0][1]=tflin->GetParameter(1);
// 			for(int i1=0;i1<tg1dcol->GetN();i1++)
// 			{
// 				tg1dcol->GetPoint(i1,x1,y1);
// 				if(fabs(tflinpars[0][0]+tflinpars[0][1]*x1-y1)>10){tg1dcol->RemovePoint(i1);nd=true;break;}
// 			}
// 		}
// 		fitnormchi2col=tflin->GetChisquare()/tflin->GetNDF();
		
		FitNormChi2_Col_vs_t->Fill(fitnormchi2col);
		NColHits_vs_fitnormchi2->Fill(fitnormchi2col,tg1dcol->GetN());
		
		FitNormChi2_Ind_vs_t->Fill(fitnormchi2ind);
		NIndHits_vs_fitnormchi2->Fill(fitnormchi2ind,tg1dind->GetN());
		
// 		cout<<I<<" "<<hd.ColIndT->at(nmaxind)[0]<<" "<<hd.ColIndT->at(nmaxind)[1]<<" "<<hd.ColIndT->at(nmaxind)[2]<<" "<<nmax<<" "<<fitnormchi2col<<" "<<fitnormchi2ind<<endl;
		
// 		if(fitnormchi2col>100 || fitnormchi2ind>100) continue;
		
		
		
// 		sprintf(hname,"[0]+%f*x",ameancol);
// 		tflin2=new TF1("lin2",hname,750,2500);
// 		tg1dcol->Fit(tflin2,"q","q",tmin-1,tmax+1);
		
		
		
// 		tmin=5000;tmax=0;
// 		indoftmin=0;indoftmax=0;
// 		for(int i1=0;i1<hd.IndID->size();i1++)
// 		{
// 			tg1dind->SetPoint(tg1dind->GetN(),hd.IndT->at(i1),hd.IndID->at(i1));
// 			if(hd.IndT->at(i1)<tmin) {tmin=hd.IndT->at(i1);indoftmin=i1;}
// 			if(hd.IndT->at(i1)>tmax) {tmax=hd.IndT->at(i1);indoftmax=i1;}
// 		}
// 		tg1dind->Fit(tflin,"q","q",tmin-1,tmax+1);
// // 		tg1dind->Fit(tflin,"q","q",tmin-1,tmin+500);
// 		fitnormchi2ind=tflin->GetChisquare()/tflin->GetNDF();
// 		tflinpars[1][0]=tflin->GetParameter(0);tflinpars[1][1]=tflin->GetParameter(1);
// 		
// 		nd=true;
// 		while(nd)
// 		{
// 			nd=false;
// 			tg1dind->Fit(tflin,"q","q",tmin-1,tmax+1);
// 			tflinpars[1][0]=tflin->GetParameter(0);tflinpars[1][1]=tflin->GetParameter(1);
// 			for(int i1=0;i1<tg1dind->GetN();i1++)
// 			{
// 				tg1dind->GetPoint(i1,x1,y1);
// 				if(fabs(tflinpars[1][0]+tflinpars[1][1]*x1-y1)>10){tg1dind->RemovePoint(i1);nd=true;break;}
// 			}
// 		}
// 		fitnormchi2col=tflin->GetChisquare()/tflin->GetNDF();
		
		
		
// 		if(fitnormchi2col>100 || fitnormchi2ind>100) continue;
		
		
		//for the track
		tmintrack=5000;tmaxtrack=0;
		indoftmin=0;indoftmax=0;Nvalid=0;
		
		QColtot=0;
		tg2d->Set(0);
		for(int i1=0;i1<hd.ColIndT->size();i1++)
		{
			if(float(hd.ColIndT->at(i1)[2])<tmincol || float(hd.ColIndT->at(i1)[2])>tmaxcol) continue;
			x1=tflinpars[0][0]+tflinpars[0][1]*((float)hd.ColIndT->at(i1)[2]);
// 			y1=tflinpars[1][0]+tflinpars[1][1]*hd.ColIndT->at(i1)[2];
// 			if(sqrt(pow(x1-hd.ColIndT->at(i1)[0],2)+pow(y1-hd.ColIndT->at(i1)[1],2))<=10.)
			if(fabs(x1-hd.ColIndT->at(i1)[0])<=10.)
			{
				tg2d->SetPoint(tg2d->GetN(),hd.ColIndT->at(i1)[0],hd.ColIndT->at(i1)[1],hd.ColIndT->at(i1)[2]);
				Nvalid++;
				if(((float)hd.ColIndT->at(i1)[2])>tmaxtrack) {tmaxtrack=hd.ColIndT->at(i1)[2];indoftmax=i1;}
				if(((float)hd.ColIndT->at(i1)[2])<tmintrack) {tmintrack=hd.ColIndT->at(i1)[2];indoftmin=i1;}
			}
		}
		pStart[0]=hd.ColIndT->at(indoftmin)[0];
		pStart[1]=((float)(hd.ColIndT->at(indoftmax)[0]-hd.ColIndT->at(indoftmin)[0]))/((float)(hd.ColIndT->at(indoftmax)[2]-hd.ColIndT->at(indoftmin)[2]));
		pStart[2]=hd.ColIndT->at(indoftmin)[1];
		pStart[3]=((float)(hd.ColIndT->at(indoftmax)[1]-hd.ColIndT->at(indoftmin)[1]))/((float)(hd.ColIndT->at(indoftmax)[2]-hd.ColIndT->at(indoftmin)[2]));
		pStart[4]=750;
		
		fitnormchi2=line3Dfit(tg2d);
		
		TrackDeltaT->Fill(tmaxtrack-tmintrack);
		FitNormChi2->Fill(fitnormchi2);
		
// 		if(fitnormchi2<5.)
		{
			td.E=hd.E;
			td.StartEndColIndT->clear();
			td.StartEndColIndT->push_back(hd.ColIndT->at(indoftmin)[0]);
			td.StartEndColIndT->push_back(hd.ColIndT->at(indoftmin)[1]);
			td.StartEndColIndT->push_back(hd.ColIndT->at(indoftmin)[2]);
			td.StartEndColIndT->push_back(hd.ColIndT->at(indoftmax)[0]);
			td.StartEndColIndT->push_back(hd.ColIndT->at(indoftmax)[1]);
			td.StartEndColIndT->push_back(hd.ColIndT->at(indoftmax)[2]);
			td.FitParams->clear();
			for(int is2=0;is2<5;is2++)
			{
				td.FitParams->push_back(fitParams[is2]);
			}
			td.FitNormChi2=fitnormchi2;
			td.QColTot=QColtot;
			td.NHits=Nvalid;
			td.Nexcl=hd.ColIndT->size()-Nvalid;
			td.PMTIntegral->clear();
			for(int is2=0;is2<3;is2++)
			{
				td.PMTIntegral->push_back(hd.PMTIntegral->at(is2));
			}
			
			T1->Fill();
			
			x1=fitParams[0]+fitParams[1]*(tmintrack-fitParams[4]);
			y1=fitParams[2]+fitParams[3]*(tmintrack-fitParams[4]);
			z1=tmintrack;
			
			x2=fitParams[0]+fitParams[1]*(tmaxtrack-fitParams[4]);
			y2=fitParams[2]+fitParams[3]*(tmaxtrack-fitParams[4]);
			z2=tmaxtrack;
			
// 			DeltaTvsTrackAngle->Fill(atan(sqrt(pow(x2-x1,2)+pow(y2-y1,2))/(tmax-tmin)),(tmax-tmin));
// 			DeltaTvsTrackAnglePr->Fill(atan(sqrt(pow(x2-x1,2)+pow(y2-y1,2))/(tmax-tmin)),(tmax-tmin));
			Tfirst->Fill(tmincol);
			Tlast->Fill(tmaxcol);
			Trackt0->Fill(fitParams[4]);
			TrackEntryXY->Fill(x1,y1);
			TrackExitXY->Fill(x2,y2);
// 			if(tmin<=1160 && tmax>=1920) TrackDeltaTSel->Fill(tmax-tmin);
			DeltaTvsTrackt0->Fill(fitParams[4],tmaxcol-tmincol);
// 			DeltaTLengthvsAngle->Fill(atan(sqrt(pow(x2-x1,2)+pow(y2-y1,2))/(tmax-tmin)),sqrt(pow(x2-x1,2)+pow(y2-y1,2)+pow(tmax-tmin,2)),tmax-tmin);
// 			DeltaTLengthvsAngle_norm->Fill(atan(sqrt(pow(x2-x1,2)+pow(y2-y1,2))/(tmax-tmin)),sqrt(pow(x2-x1,2)+pow(y2-y1,2)+pow(tmax-tmin,2)));
			NExclHits->Fill(td.Nexcl);
// // 			if(td.NHits>15 && td.Nexcl<10)
// 			if(td.NHits>10 && td.Nexcl<10)
// 			{
// 				outfile<<td.E<<" "<<hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[0]<<" "<<hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[1]<<" "<<hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[2]<<" "<<hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[0]<<" "<<hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[1]<<" "<<hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[2]<<endl;
// 			}
		}
		if(I%100==0) cout<<I<<" / "<<T->GetEntries()<<endl;
		if(I>50) break;
	}
	
// 	for(int I=0;I<T->GetEntries();I++)
// 	{
// 		T->GetEntry(I);
// 		
// 		if(hd.ColIndT->size()<5) continue;
// 		
// // 		if(hd.ColIndT->size()<20) continue;
// 		vector <clusters> CL;
// 		cl.hi.clear();
// 		inCL.clear();
// 		for(int i1=0;i1<hd.ColIndT->size();i1++)
// 		{
// 			inCL.push_back(0);
// 		}
// 		tmin=5000;tmax=0;
// 		for(int i1=0;i1<hd.ColIndT->size();i1++)
// 		{
// 			if(hd.ColIndT->at(i1)[2]>tmax) tmax=hd.ColIndT->at(i1)[2];
// 			if(hd.ColIndT->at(i1)[2]<tmin) tmin=hd.ColIndT->at(i1)[2];
// 			if(inCL[i1]!=0) continue;
// 			cl.hi.clear();
// 			cl.hi.push_back(i1);
// 			inCL[i1]=1;
// 			nid=1;
// 			while(nid>0)
// 			{
// 				nid=0;
// 				for(int i2=0;i2<hd.ColIndT->size();i2++)
// 				{
// 					if(inCL[i2]!=0) continue;
// 					for(int i3=0;i3<cl.hi.size();i3++)
// 					{
// // 						if(sqrt(pow(H[cl.hi[i3]].col-H[i2].col,2)+pow(H[cl.hi[i3]].ind-H[i2].ind,2))<10 && (H[cl.hi[i3]].t-H[i2].t)<20)
// // 						if(sqrt(pow(hd.ColIndT->at(i1)[0]-hd.ColIndT->at(i2)[0],2)+pow(hd.ColIndT->at(i1)[1]-hd.ColIndT->at(i2)[1],2))<10 && (hd.ColIndT->at(i1)[2]-hd.ColIndT->at(i2)[2])<20)
// 						if(sqrt(pow(hd.ColIndT->at(i1)[0]-hd.ColIndT->at(i2)[0],2)+pow(hd.ColIndT->at(i1)[1]-hd.ColIndT->at(i2)[1],2))<5 && (hd.ColIndT->at(i1)[2]-hd.ColIndT->at(i2)[2])<30)
// 						{
// 							cl.hi.push_back(i2);
// 							nid++;
// 							inCL[i2]=1;
// 							break;
// 						}
// 					}
// 				}
// 			}
// 			CL.push_back(cl);
// 			cl.hi.clear();
// 		}
// 		AllDeltaT->Fill(tmax-tmin);
// 		
// 		int HitsinCLs1=0;int HitsinCLs2=0;
// 		indCLHP=-1;maxCLsize=0;
// 		for(int i1=0;i1<CL.size();i1++)
// 		{
// 			if(CL[i1].hi.size()>maxCLsize){maxCLsize=CL[i1].hi.size();indCLHP=i1;}
// 			HitsinCLs1+=CL[i1].hi.size();
// 		}
// 		
// 		tmin=5000;tmax=0;
// 		QColtot=0;
// 		tg2d->Set(0);
// 		
// 		for(int i1=0;i1<CL[indCLHP].hi.size();i1++)
// 		{
// 			tg2d->SetPoint(tg2d->GetN(),hd.ColIndT->at(CL[indCLHP].hi[i1])[0],hd.ColIndT->at(CL[indCLHP].hi[i1])[1],hd.ColIndT->at(CL[indCLHP].hi[i1])[2]);
// 			if(hd.ColIndT->at(CL[indCLHP].hi[i1])[2]<tmin){tmin=hd.ColIndT->at(CL[indCLHP].hi[i1])[2];indoftmin=i1;}
// 			if(hd.ColIndT->at(CL[indCLHP].hi[i1])[2]>tmax){tmax=hd.ColIndT->at(CL[indCLHP].hi[i1])[2];indoftmax=i1;}
// 			TrackXY->Fill(hd.ColIndT->at(CL[indCLHP].hi[i1])[0],hd.ColIndT->at(CL[indCLHP].hi[i1])[1]);
// 			QColtot+=hd.QColTot;
// 		}
// 		
// 		pStart[0]=hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[0];
// 		pStart[1]=((float)(hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[0]-hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[0]))/((float)(hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[2]-hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[2]));
// 		pStart[2]=hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[1];
// 		pStart[3]=((float)(hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[1]-hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[1]))/((float)(hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[2]-hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[2]));
// 		pStart[4]=1500;
// 		
// 		fitnormchi2=line3Dfit(tg2d);
// 		
// // 		cout<<I<<endl;
// // 		cout<<"before"<<endl;
// // 		for(int i1=0;i1<CL.size();i1++)
// // 		{
// // 			cout<<i1<<" "<<CL[i1].hi.size()<<endl;
// // 		}
// 		
// 		for(int i1=0;i1<CL.size();i1++)
// 		{
// 			if(i1==indCLHP) continue;
// 			for(int i2=0;i2<CL[i1].hi.size();i2++)
// 			{
// 				x1=fitParams[0]+fitParams[1]*(hd.ColIndT->at(CL[i1].hi[i2])[2]-fitParams[4]);
// 				y1=fitParams[2]+fitParams[3]*(hd.ColIndT->at(CL[i1].hi[i2])[2]-fitParams[4]);
// 				if(sqrt(pow(x1-hd.ColIndT->at(CL[i1].hi[i2])[0],2.)+pow(y1-hd.ColIndT->at(CL[i1].hi[i2])[1],2.))<10.)
// 				{
// 					CL[indCLHP].hi.push_back(CL[i1].hi[i2]);
// 					
// 					tg2d->SetPoint(tg2d->GetN(),hd.ColIndT->at(CL[i1].hi[i2])[0],hd.ColIndT->at(CL[i1].hi[i2])[1],hd.ColIndT->at(CL[i1].hi[i2])[2]);
// 					if(hd.ColIndT->at(CL[i1].hi[i2])[2]<tmin){tmin=hd.ColIndT->at(CL[i1].hi[i2])[2];indoftmin=CL[indCLHP].hi.size()-1;}
// 					if(hd.ColIndT->at(CL[i1].hi[i2])[2]>tmax){tmax=hd.ColIndT->at(CL[i1].hi[i2])[2];indoftmax=CL[indCLHP].hi.size()-1;}
// 					TrackXY->Fill(hd.ColIndT->at(CL[i1].hi[i2])[0],hd.ColIndT->at(CL[i1].hi[i2])[1]);
// 					QColtot+=hd.QColTot;
// 					
// 					CL[i1].hi[i2]=CL[i1].hi[CL[i1].hi.size()-1];
// 					CL[i1].hi.pop_back();
// 					i2--;
// 				}
// 			}
// 		}
// // 		cout<<"after"<<endl;
// // 		for(int i1=0;i1<CL.size();i1++)
// // 		{
// // 			cout<<i1<<" "<<CL[i1].hi.size()<<endl;
// // 		}
// 		for(int i1=0;i1<CL.size();i1++)
// 		{
// 			if(CL[i1].hi.size()==0)
// 			{
// 				CL[i1].hi=CL[CL.size()-1].hi;
// 				CL.pop_back();
// 				i1--;
// 			}
// 		}
// 		indCLHP=-1;maxCLsize=0;
// 		for(int i1=0;i1<CL.size();i1++)
// 		{
// 			ClSize->Fill(CL[i1].hi.size());
// 			if(CL[i1].hi.size()>maxCLsize){maxCLsize=CL[i1].hi.size();indCLHP=i1;}
// 			HitsinCLs2+=CL[i1].hi.size();
// 		}
// // 		cout<<I<<" "<<hd.ColIndT->size()<<" "<<HitsinCLs1<<" "<<HitsinCLs2<<endl;
// 		
// // 		cout<<"after 2"<<endl;
// // 		for(int i1=0;i1<CL.size();i1++)
// // 		{
// // 			cout<<i1<<" "<<CL[i1].hi.size()<<endl;
// // 		}
// // 		cout<<endl;
// 		
// 		maxCLsize=CL[indCLHP].hi.size();
// 		MaxClSize->Fill(maxCLsize);
// 		NClperEvent->Fill(CL.size());
// 		
// 		pStart[0]=hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[0];
// 		pStart[1]=((float)(hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[0]-hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[0]))/((float)(hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[2]-hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[2]));
// 		pStart[2]=hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[1];
// 		pStart[3]=((float)(hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[1]-hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[1]))/((float)(hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[2]-hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[2]));
// 		pStart[4]=1500;
// 		
// 		fitnormchi2=line3Dfit(tg2d);
// 		
// 		TrackDeltaT->Fill(tmax-tmin);
// 		FitNormChi2->Fill(fitnormchi2);
// 		
// // 		td.E=I;
// 		td.E=hd.E;
// 		td.StartEndColIndT->clear();
// 		td.StartEndColIndT->push_back(hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[0]);
// 		td.StartEndColIndT->push_back(hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[1]);
// 		td.StartEndColIndT->push_back(hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[2]);
// 		td.StartEndColIndT->push_back(hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[0]);
// 		td.StartEndColIndT->push_back(hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[1]);
// 		td.StartEndColIndT->push_back(hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[2]);
// 		td.FitParams->clear();
// 		for(int is2=0;is2<5;is2++)
// 		{
// 			td.FitParams->push_back(fitParams[is2]);
// 		}
// 		td.FitNormChi2=fitnormchi2;
// 		td.QColTot=QColtot;
// 		td.NHits=CL[indCLHP].hi.size();
// 		td.Nexcl=hd.ColIndT->size()-CL[indCLHP].hi.size();
// 		td.PMTIntegral->clear();
// 		for(int is2=0;is2<3;is2++)
// 		{
// 			td.PMTIntegral->push_back(hd.PMTIntegral->at(is2));
// 		}
// 		
// 		T1->Fill();
// 		
// 		if(fitnormchi2<5.)
// 		{
// // // 			cc->cd(); 
// // // 			frame3d->Draw();
// // // 			
// // // 			pm3d1 = new TPolyMarker3D();
// // // 			pm3d2 = new TPolyMarker3D();
// // // 			
// // // // 			pm3d1 = new TPolyMarker3D(CL[indCLHP].hi.size());
// // // // // 			pm3d2 = new TPolyMarker3D(H.size()-CL[indCLHP].hi.size());
// // // // 			pm3d2 = new TPolyMarker3D(hd.ColIndT->size()-CL[indCLHP].hi.size());
// // // 			
// // // 			for(int i1=0;i1<CL[indCLHP].hi.size();i1++)
// // // 			{
// // // // 				pm3d1->SetPoint(i1,H[CL[indCLHP].hi[i1]].col,H[CL[indCLHP].hi[i1]].ind,H[CL[indCLHP].hi[i1]].t);
// // // 				pm3d1->SetPoint(i1,hd.ColIndT->at(CL[indCLHP].hi[i1])[0],hd.ColIndT->at(CL[indCLHP].hi[i1])[1],hd.ColIndT->at(CL[indCLHP].hi[i1])[2]);
// // // 			}
// // // 			
// // // 			int na=0;
// // // 			for(int i1=0;i1<CL.size();i1++)
// // // 			{
// // // 				if(i1==indCLHP) continue;
// // // 				for(int i2=0;i2<CL[i1].hi.size();i2++)
// // // 				{
// // // // 					pm3d2->SetPoint(na,H[CL[i1].hi[i2]].col,H[CL[i1].hi[i2]].ind,H[CL[i1].hi[i2]].t);
// // // 					pm3d2->SetPoint(na,hd.ColIndT->at(CL[i1].hi[i2])[0],hd.ColIndT->at(CL[i1].hi[i2])[1],hd.ColIndT->at(CL[i1].hi[i2])[2]);
// // // 					na++;
// // // 				}
// // // 			}
// // // // 					pm3d1->SetMarkerSize(2);
// // // 			pm3d1->SetMarkerColor(kRed);
// // // 			pm3d1->SetMarkerStyle(24);   
// // // 			pm3d1->Draw();
// // // // 					pm3d2->SetMarkerSize(2);
// // // 			pm3d2->SetMarkerColor(kBlue);
// // // 			pm3d2->SetMarkerStyle(24);   
// // // 			pm3d2->Draw();
// // // 			
// // // 			l = new TPolyLine3D(2);
// 			
// 			x1=fitParams[0]+fitParams[1]*(tmin-fitParams[4]);
// 			y1=fitParams[2]+fitParams[3]*(tmin-fitParams[4]);
// 			z1=tmin;
// // // 			l->SetPoint(0,x1,y1,z1);
// 			
// 			x2=fitParams[0]+fitParams[1]*(tmax-fitParams[4]);
// 			y2=fitParams[2]+fitParams[3]*(tmax-fitParams[4]);
// 			z2=tmax;
// // // 			l->SetPoint(1,x2,y2,z2);
// // // 			l->SetLineColor(kBlack);
// // // 			l->SetLineWidth(2);
// // // 			l->Draw("same");
// 			
// 			DeltaTvsTrackAngle->Fill(atan(sqrt(pow(x2-x1,2)+pow(y2-y1,2))/(tmax-tmin)),(tmax-tmin));
// 			DeltaTvsTrackAnglePr->Fill(atan(sqrt(pow(x2-x1,2)+pow(y2-y1,2))/(tmax-tmin)),(tmax-tmin));
// 			Tfirst->Fill(tmin);
// 			Tlast->Fill(tmax);
// //					if(tmin<=1160) Tlast->Fill(tmax);
// 			Trackt0->Fill(fitParams[4]);
// 			TrackEntryXY->Fill(x1,y1);
// 			TrackExitXY->Fill(x2,y2);
// 			if(tmin<=1160 && tmax>=1920) TrackDeltaTSel->Fill(tmax-tmin);
// 			DeltaTvsTrackt0->Fill(fitParams[4],tmax-tmin);
// 			DeltaTLengthvsAngle->Fill(atan(sqrt(pow(x2-x1,2)+pow(y2-y1,2))/(tmax-tmin)),sqrt(pow(x2-x1,2)+pow(y2-y1,2)+pow(tmax-tmin,2)),tmax-tmin);
// 			DeltaTLengthvsAngle_norm->Fill(atan(sqrt(pow(x2-x1,2)+pow(y2-y1,2))/(tmax-tmin)),sqrt(pow(x2-x1,2)+pow(y2-y1,2)+pow(tmax-tmin,2)));
// 			
// //					cout<<Ev<<" "<<tmin<<" "<<tmax<<" "<<(tmax-tmin)<<endl;
// 			
// // 			sprintf(hname,"Ev_%d_%5.3f",td.E,fitnormchi2);
// // 			cc->SetName(hname);cc->SetTitle(hname);
// // 			outroot->cd();
// // 			cc->Write();
// // 			tg2dA->Set(0);
// // 			delete l;
// 			
// 			FitClSize->Fill(CL[indCLHP].hi.size());
// 			NExclHits->Fill(hd.ColIndT->size()-CL[indCLHP].hi.size());
// // 			if(CL[indCLHP].hi.size()>15 && (H.size()-CL[indCLHP].hi.size())<10)
// // 			{
// // 				Nvalid++;
// // 				outfile<<Ev<<" "<<fitParams[0]<<" "<<fitParams[1]<<" "<<fitParams[2]<<" "<<fitParams[3]<<" "<<fitParams[4]<<" "<<fitnormchi2<<" "<<CL[indCLHP].hi.size()<<endl;
// // 			}
// // 			if(td.NHits>15 && td.Nexcl<10)
// 			if(td.NHits>10 && td.Nexcl<10)
// 			{
// 				outfile<<td.E<<" "<<hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[0]<<" "<<hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[1]<<" "<<hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[2]<<" "<<hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[0]<<" "<<hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[1]<<" "<<hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[2]<<endl;
// 			}
// 		}
// 		tg2d->Set(0);
// 		inCL.clear();
// 		cl.hi.clear();
// 		CL.clear();
// // 		if(I==10) break;
// 		if(I%100==0) cout<<I<<" / "<<T->GetEntries()<<endl;
// 	}
// 	cout<<Nvalid<<" / "<<Ev<<" "<<(((float)Nvalid)/((float)Ev))<<endl;
	DeltaTLengthvsAngle->Divide(DeltaTLengthvsAngle_norm);
	outroot->Write();
	outroot->Close();
	outfile.close();
	inroot->Close();
	
// 	sprintf(hname,"cp TrackParameters_%d.txt %s/Files/TrackParameters_%d.txt;wait;",RunNo,AnalysisFilePath,RunNo);
// 	system(hname);
	sprintf(hname,"cp Tracks_%d.root %s/Histos/Tracks_%d.root",RunNo,AnalysisFilePath,RunNo);system(hname);
	sprintf(hname,"cp ValidTrackParameters_%d.txt %s/Files/ValidTrackParameters_%d.txt;wait;",RunNo,AnalysisFilePath,RunNo);system(hname);
}

void FindDeltaT()
{
// 	sprintf(hname,"%s/Histos/DeltaTFinding_%d.root",AnalysisFilePath,RunNo);
// 	TFile* outroot=new TFile(hname,"recreate");
// 	TH1F* Tfirst=new TH1F("Tfirst","Tfirst",4096,-0.5,4095.5);
// 	TH1F* Tlast=new TH1F("Tlast","Tlast",4096,-0.5,4095.5);
// 	TH1F* DeltaT=new TH1F("DeltaT","DeltaT",4096,-0.5,4095.5);
// 	
// 	sprintf(hname,"%s/Run00%d.dat",TPCinFilePath,RunNo);
// 	file.open(hname, ios::in|ios::binary|ios::ate);
// 	
// 	ReadTPCBaselines();
// 	
// 	int Ev=0;
// 	int integral=0;int integral2=0;int integral3=0;
// 	int peakbin=0;
// 	int tmin=-1;int tmax=0;
// 	float ss=0;
// 	file.seekg(0, ios::beg);
// 	while(!file.eof())
// 	{
// 		ReadEvent();
// 		
// 		tmin=-1;
// 		for(int i2=0;i2<4096;i2++)
// 		{
// 			for(int i1=0;i1<256;i1++)
// 			{
// 				if(IC[i1]!=0) continue;
// 				ss=((float)e.signal[i1][i2])-TPCBaselines[IC[i1]][CH[i1]][0];
// 				if(ss>ThCol)
// 				{
// 					if(tmin==-1) tmin=i2;
// 					tmax=i2;
// 				}
// 			}
// 		}
// 		Tfirst->Fill(tmin);
// 		Tlast->Fill(tmax);
// 		DeltaT->Fill(tmax-tmin);
// 			
// 		Ev++;
// 		for(int i1=0;i1<256;i1++)
// 		{
// 			e.signal[i1].clear();
// 		}
// 		if(Ev%10==0) cout<<"Event: "<<Ev<<endl;
// 	}
// 	file.close();
// 	outroot->cd();
// 	outroot->Write();
// 	outroot->Close();
}

void FindHitsAndFitTracks()
{
	if(RunNo==7329) HitColIndDist=5.43;
	if(RunNo==7331 || RunNo==7283) HitColIndDist=4.54;
	if(RunNo==7333 || RunNo==7334) HitColIndDist=4.27;
	if(RunNo==7336) HitColIndDist=6.27;
	if(RunNo==7338 || RunNo==7339 || RunNo==7340 || RunNo==7341 || RunNo==7342) HitColIndDist=4.05;
	if(RunNo==7344 || RunNo==7345 || RunNo==7346) HitColIndDist=4.16;
	if(RunNo==7348) HitColIndDist=4.26;
	if(RunNo==7349) HitColIndDist=4.27;
	if(RunNo==7351) HitColIndDist=4.54;
	if(RunNo>=7353 && RunNo<=7371) HitColIndDist=4.42;
	
// 	sprintf(hname,"Hits_%d.root",RunNo);
	sprintf(hname,"HitsAndTracks_%d.root",RunNo);
	TFile* outroot=new TFile(hname,"recreate");
// 	TTree* T1 =  new TTree("T","T");
	TTree* T1 =  new TTree("Hits","Hits");
	T1->Branch("E",&hd.E);
	T1->Branch("ColIndT",&hd.ColIndT);
	T1->Branch("Int",&hd.Int);
	T1->Branch("QColTot",&hd.QColTot);
	T1->Branch("QColTotZS",&hd.QColTotZS);
	T1->Branch("QHitTot",&hd.QHitTot);
	T1->Branch("PMTIntegral",&hd.PMTIntegral);
	
	T1->Branch("ColID",&hd.ColID);
	T1->Branch("ColT",&hd.ColT);
	T1->Branch("ColA",&hd.ColA);
	T1->Branch("ColInt",&hd.ColInt);
	T1->Branch("Colw",&hd.Colw);
	T1->Branch("IndID",&hd.IndID);
	T1->Branch("IndT",&hd.IndT);
	T1->Branch("IndA",&hd.IndA);
	T1->Branch("Indw",&hd.Indw);
	T1->Branch("EventType",&hd.EventType);
	
	outroot->mkdir("Histos");
	
	sprintf(hname,"ValidTrackParameters_%d.txt",RunNo);
	ofstream outfile(hname);
	vector <int> eventlist;
	
// 	sprintf(hname,"Tracks_%d.root",RunNo);
// 	TFile* outroottr=new TFile(hname,"recreate");
// 	TTree* Ttr =  new TTree("T","T");
	TTree* Ttr =  new TTree("Tracks","Tracks");
	Ttr->Branch("E",&td.E);
	Ttr->Branch("StartEndColIndT",&td.StartEndColIndT);
	Ttr->Branch("FitParams",&td.FitParams);
	Ttr->Branch("FitNormChi2",&td.FitNormChi2);
	Ttr->Branch("QColTot",&td.QColTot);
	Ttr->Branch("NHits",&td.NHits);
	Ttr->Branch("Nexcl",&td.Nexcl);
	Ttr->Branch("PMTIntegral",&td.PMTIntegral);
	Ttr->Branch("ColTStartEnd",&td.ColTStartEnd);
	Ttr->Branch("ColHitTStartEnd",&td.ColHitTStartEnd);
	
	sprintf(hname,"%s/Run_%d.root",NTupleFilePath,RunNo);
	TFile* inroot=new TFile(hname);
	TTree* T =  (TTree*) inroot->Get("T");
	T->SetBranchAddress("TPCWF",&ed.TPCWF);
	T->SetBranchAddress("PMTWF",&ed.PMTWF);
	
	TH1F* HitCol=new TH1F("HitCol","HitCol",128,-0.5,127.5);
	TH1F* HitInd=new TH1F("HitInd","HitInd",128,-0.5,127.5);
	TH1F* HitTime=new TH1F("HitTime","HitTime",4096,-0.5,4095.5);
	TH1F* HitIntegral1=new TH1F("HitIntegral1","HitIntegral1",500,0,3000);
	TH1F* HitIntegral2=new TH1F("HitIntegral2","HitIntegral2",500,0,3000);
	TH1F* HitIntegral3=new TH1F("HitIntegral3","HitIntegral3",500,0,3000);
	TH1F* NHits=new TH1F("NHits","NHits",2000,-0.5,1999.5);
	TH1F* NHitsPerColPeak=new TH1F("NHitsPerColPeak","NHitsPerColPeak",200,-0.5,199.5);
	TH1F* HQTotCol=new TH1F("QTotCol","QTotCol",4000,0,1000000);
	TH1F* NPoints=new TH1F("NPoints","NPoints",2240,0,224000);//224001,-0.5,224000.5);
	TH1F* ClusterNPointsperArea=new TH1F("ClusterNPointsperArea","ClusterNPointsperArea",1020,-0.01,1.01);
	TH1F* NPointsperArea=new TH1F("NPointsperArea","NPointsperArea",1020,-0.01,1.01);
	TH1F* TotalArea=new TH1F("TotalArea","TotalArea",2240,0,224000);//224001,-0.5,224000.5);
	TH1F* ClmaxPointFraction=new TH1F("ClmaxPointFraction","ClmaxPointFraction",510,-0.01,1.01);
	TH1F* DeltaTTrackColT=new TH1F("DeltaTTrackColT","DeltaTTrackColT",500,0,5000);
	TH1F* DeltaTTrackColHitT=new TH1F("DeltaTTrackColHitT","DeltaTTrackColHitT",500,0,5000);
	TH1F* NClperEvent=new TH1F("NClperEvent","NClperEvent",60,-0.5,59.5);
	TH1F* ClSize=new TH1F("ClSize","ClSize",400,-0.5,399.5);
	TH1F* MaxClSize=new TH1F("MaxClSize","MaxClSize",400,-0.5,399.5);
	
	TH1F* FitNormChi2=new TH1F("FitNormChi2","FitNormChi2",200,0.,100.);
	TH1F* FitClSize=new TH1F("FitClSize","FitClSize",200,-0.5,199.5);
	TH1F* NExclHits=new TH1F("NExclHits","NExclHits",200,-0.5,199.5);
	TH1F* TrackDeltaT=new TH1F("TrackDeltaT","TrackDeltaT",500,0,5000);
	TH1F* ColTrackDeltaT=new TH1F("ColTrackDeltaT","ColTrackDeltaT",500,0,5000);
	
	TH1F* TrackDeltaTSel=new TH1F("TrackDeltaTSel","TrackDeltaTSel",500,0,5000);
	TH1F* AllDeltaT=new TH1F("AllDeltaT","AllDeltaT",500,0,5000);
	TH2F* AllXY=new TH2F("AllXY","AllXY",128,-0.5,127.5,128,-0.5,127.5);
	TH2F* TrackXY=new TH2F("TrackXY","TrackXY",128,-0.5,127.5,128,-0.5,127.5);
	TH2F* TrackEntryXY=new TH2F("TrackEntryXY","TrackEntryXY",128,-0.5,127.5,128,-0.5,127.5);
	TH2F* TrackExitXY=new TH2F("TrackExitXY","TrackExitXY",128,-0.5,127.5,128,-0.5,127.5);
	TH2F* DeltaTvsTrackAngle=new TH2F("DeltaTvsTrackAngle","DeltaTvsTrackAngle",314,0,1.57,500,0,1000);
	TProfile* DeltaTvsTrackAnglePr=new TProfile("DeltaTvsTrackAnglePr","DeltaTvsTrackAnglePr",250,0,1.57,0,1000);
	TH1F* Tfirst=new TH1F("Tfirst","Tfirst",500,0,5000);
	TH1F* Tlast=new TH1F("Tlast","Tlast",500,0,5000);
	TH1F* TColfirst=new TH1F("TColfirst","TColfirst",500,0,5000);
	TH1F* TCollast=new TH1F("TCollast","TCollast",500,0,5000);
	TH1F* Trackt0=new TH1F("Trackt0","Trackt0",500,0,5000);
	TH2F* DeltaTvsTrackt0=new TH2F("DeltaTvsTrackt0","DeltaTvsTrackt0",500,0,5000,500,0,1000);
// 	TH1F* DeltacolperDeltat=new TH1F("DeltacolperDeltat","DeltacolperDeltat",500,0,50);
// 	TH1F* DeltaindperDeltat=new TH1F("DeltaindperDeltat","DeltaindperDeltat",500,0,50);
// 	TH1F* DeltaRperDeltat=new TH1F("DeltaRperDeltat","DeltaRperDeltat",2000,0,100);
// 	TH1F* ColDeltaRR=new TH1F("ColDeltaR","ColDeltaR",2000,0,2000);
// 	TH1F* HitDeltaRR=new TH1F("HitDeltaRR","HitDeltaRR",2000,0,2000);
// 	TH1F* HitDeltaT=new TH1F("HitDeltaT","HitDeltaT",2000,0,1000);
// 	TH1F* HitDeltaR=new TH1F("HitDeltaR","HitDeltaR",200,0,100);
// 	TH2F* HitDeltaR_vs_DeltaT=new TH2F("HitDeltaR_vs_DeltaT","HitDeltaR_vs_DeltaT",200,-0.5,199.5,200,0,20);
	
	TH2F* DeltaTLengthvsAngle=new TH2F("DeltaTLengthvsAngle","DeltaTLengthvsAngle",314,0,1.57,100,0,2000);
	TH2F* DeltaTLengthvsAngle_norm=new TH2F("DeltaTLengthvsAngle_norm","DeltaTLengthvsAngle_norm",314,0,1.57,100,0,2000);
	
	TH1F* FitNormChi2_Col_vs_t=new TH1F("FitNormChi2_Col_vs_t","FitNormChi2_Col_vs_t",1000,0.,4000.);
	TH1F* FitNormChi2_Ind_vs_t=new TH1F("FitNormChi2_Ind_vs_t","FitNormChi2_Ind_vs_t",1000,0.,4000.);
	
	TH2F* NColHits_vs_fitnormchi2=new TH2F("NColHits_vs_fitnormchi2","NColHits_vs_fitnormchi2",200,-0.5,199.5,500,0,2000);
	TH2F* NIndHits_vs_fitnormchi2=new TH2F("NIndHits_vs_fitnormchi2","NIndHits_vs_fitnormchi2",200,-0.5,199.5,500,0,2000);
	
	
	
	TH1F* StartTime=new TH1F("StartTime","StartTime",512,-0.5,4095.5);
	TH1F* EndTime=new TH1F("EndTime","EndTime",512,-0.5,4095.5);
	TH1F* StartCollection=new TH1F("StartCollection","StartCollection",128,-0.5,127.5);
	TH1F* EndCollection=new TH1F("EndCollection","EndCollection",128,-0.5,127.5);
	TH1F* StartInduction=new TH1F("StartInduction","StartInduction",128,-0.5,127.5);
	TH1F* EndInduction=new TH1F("EndInduction","EndInduction",128,-0.5,127.5);
	
	TH1F* ColTWidth=new TH1F("ColTWidth","ColTWidth",128,-0.5,127.5);
	TH1F* IndTWidth=new TH1F("IndTWidth","IndTWidth",128,-0.5,127.5);
	
	TH2F* NHits_T_vs_E=new TH2F("NHits_T_vs_E","NHits_T_vs_E",T->GetEntries(),-0.5,T->GetEntries()-0.5,4096,-0.5,4095.5);
	TH2F* NHits_Col_vs_Ind=new TH2F("NHits_Col_vs_Ind","NHits_Col_vs_Ind",20,-0.5,19.5,20,-0.5,19.5);
	TH2F* QCol_vs_QInd=new TH2F("QCol_vs_QInd","QCol_vs_QInd",50,0,200,50,0,200);
	
	TH1F* ColIndDeltaT=new TH1F("ColIndDeltaT","ColIndDeltaT",2000,-5.,5.);
	TH1F* ColIndMinDeltaT5=new TH1F("ColIndMinDeltaT5","ColIndMinDeltaT5",4000,-10.,10.);
	TH1F* ColWidth=new TH1F("ColWidth","ColWidth",128,-0.5,127.5);
	TH1F* IndWidth=new TH1F("IndWidth","IndWidth",128,-0.5,127.5);
	TH1F* ColAmp=new TH1F("ColAmp","ColAmp",300,0.,300.);
	TH1F* IndAmp=new TH1F("IndAmp","IndAmp",300,0.,300.);
	TH2F* Amp_Col_vs_Ind=new TH2F("Amp_Col_vs_Ind","Amp_Col_vs_Ind",300,0,300,300,0,300);
	TH2F* Width_Col_vs_Ind=new TH2F("Width_Col_vs_Ind","Width_Col_vs_Ind",128,-0.5,127.5,128,-0.5,127.5);
	TH1F* HitIntegralLost=new TH1F("HitIntegralLost","HitIntegralLost",4000,0,1000000);
	TH1F* TrueTotalIntegral=new TH1F("TrueTotalIntegral","TrueTotalIntegral",4000,0,1000000);
	TH2F* Col_Int_vs_Amp=new TH2F("Col_Int_vs_Amp","Col_Int_vs_Amp",300,0,300,500,0,3000);
	
	
	TH1D* QperT[3];
	QperT[0]=new TH1D("QperT","QperT",3500,750,2500);
	QperT[1]=new TH1D("QperT_norm","QperT_norm",3500,750,2500);
	QperT[2]=new TH1D("QperT_div","QperT_div",3500,750,2500);
	
	ReadTPCBaselines();
	ReadPMTBaselines();
	
	TSpectrum *s;
	Double_t *noxpeaks;
	vector <int> xpeaks;
	TH1F* hh[2][128];TH1F* hh2[2][128];TH1F* hh3[2][128];
	TH1F* ht[4096];TH1F* ht2[4096];TH1F* ht3[4096];
	vector <int> inCL;
	double x1(0),y1(0),z1(0),t1(0);
	double x2(0),y2(0),z2(0),t2(0);
	int nid=0;
	
	for(int i1=0;i1<2;i1++)
	{
		for(int i2=0;i2<128;i2++)
		{
			sprintf(hname,"WF_%s_%d",wt[i1].c_str(),i2);
			hh[i1][i2]=new TH1F(hname,hname,4096,-0.5,4095.5);
			sprintf(hname,"WFZS_%s_%d",wt[i1].c_str(),i2);
			hh2[i1][i2]=new TH1F(hname,hname,4096,-0.5,4095.5);
			sprintf(hname,"WFAbs_%s_%d",wt[i1].c_str(),i2);
			hh3[i1][i2]=new TH1F(hname,hname,4096,-0.5,4095.5);
		}
	}
	for(int i1=0;i1<4096;i1++)
	{	
		sprintf(hname,"WFt_%d",i1);
		ht[i1]=new TH1F(hname,hname,256,-0.5,255.5);
		sprintf(hname,"WFtZS_%d",i1);
		ht2[i1]=new TH1F(hname,hname,256,-0.5,255.5);
		sprintf(hname,"WFtIndInvZS_%d",i1);
		ht3[i1]=new TH1F(hname,hname,256,-0.5,255.5);
	}
	
	TGraph* W_vs_t[2];
	W_vs_t[0]=new TGraph();
	W_vs_t[0]->SetName("Col_vs_t");
	W_vs_t[0]->SetTitle("Col_vs_t");
	W_vs_t[0]->SetMarkerStyle(20);
	W_vs_t[1]=new TGraph();
	W_vs_t[1]->SetName("Ind_vs_t");
	W_vs_t[1]->SetTitle("Ind_vs_t");
	W_vs_t[1]->SetMarkerStyle(20);
	TGraph* Col_vs_t_fit=new TGraph();
	Col_vs_t_fit->SetName("Col_vs_t_fit");
	Col_vs_t_fit->SetTitle("Col_vs_t_fit");
	Col_vs_t_fit->SetMarkerStyle(20);
	
	TGraph2D* tg2d=new TGraph2D();
	sprintf(hname,"Ev3d");
	tg2d->SetName(hname);tg2d->SetTitle(hname);
	tg2d->SetMarkerStyle(20);tg2d->SetMarkerColor(4);
	int nh=0;
	vector <hits> H;
	TH3F *frame3d=new TH3F("frame3d","frame3d",1,-0.5,127.5,1,-0.5,127.5,1,-0.5,4095.5);
	TPolyMarker3D *pm3d1;
	TCanvas* cc=new TCanvas("cc","cc",600,600);
	cc->cd();
	gStyle->SetOptStat(0);
	
	bool isSaturated=false;
	float sg=0.;
	float qint=0.;
	float QtotTrue=0;
	float QtotHits=0;
	float QtotCol=0;
	
	vector <int> bb;
	vector <float> bc;
	int Ev=0;
	int integral=0;int integral2=0;int integral3=0;
	int peakbin=0;
	int peakbin2=0;
	int colminbin=0;int colmaxbin=0;
	int indminbin=0;int indmaxbin=0;
	int xmin=0;int xmax=0;
	int xminpre=0;int xmaxpre=0;
	int ic=0;int gch=0;
	float QColTot=0;
	float QColTotZS=0;
	float colindmindeltaT=0;
	float tmincl=0;float tmaxcl=0;
	float colmincl=0;float colmaxcl=0;
	float tminallpt=0;float tmaxallpt=0;
	float colminallpt=0;float colmaxallpt=0;
	float colhittmin=0;float colhittmax=0;
	
	float tmintrack=5000;float tmaxtrack=0;
	int indoftmin=0;int indoftmax=0;
	int Nvalid=0;
	float QColtot=0;
	float fitnormchi2=0;
	
	TF1* gf=new TF1("gaus","gaus",0.,5000.);
	TF1* lf=new TF1("lin","[0]-[1]*x",0.,5000.);
	TF1* tflin=new TF1("lin","[0]+[1]*x",750,2500);
	
	float minQ=0;
	int coltsel=0;int coltseldiff=0;
	int indtsel=0;int indtseldiff=0;int indisel=0;
	float undershootoffset=0;
	int nhitspercolpeak=0;
	int nhitsprecol=0;
	int nhitspostcol=0;
	int maxclsize=0;
	int indCLHP=0;
	double tflinpars[2][2]={{0}};//col, ind - const, slope
	float fitnormchi2col=0;
	vector <float> trackXY[2];
	
	for(int I=0;I<T->GetEntries();I++)
// 	for(int I=8;I<=8;I++)
	{
		T->GetEntry(I);
		int npeaks = 200;
		int nfound=0;
		float ss=0;
		bool found=true;
		bool HitOK=true;
		
		QColTot=0;QColTotZS=0;
		for(int i1=0;i1<256;i1++)
		{
			ic=i1/128;gch=i1-((i1/128)*128);
// 			for(int i2=0;i2<4096;i2++)
			for(int i2=TimeRange[0];i2<TimeRange[1];i2++)
			{
				ss=((float)ed.TPCWF->at(i1)[i2])-TPCBaselines[ic][gch][0];
				if(i1<128)
				{
					QColTot+=ss;
					if(ss>0) QColTotZS+=ss;
					hh[ic][gch]->SetBinContent(hh[ic][gch]->FindBin(i2),ss);
// 					if(ss>ThColInd[ic])
					if(ss>(TPCBaselines[ic][gch][1]*ThSigma[ic]))
					{
// 						ht3[i2]->SetBinContent(ht3[i2]->FindBin(ic*128+gch),ss);
						hh2[ic][gch]->SetBinContent(hh2[ic][gch]->FindBin(i2),ss);
					}
				}
				else
				{
					hh[ic][gch]->SetBinContent(hh[ic][gch]->FindBin(i2),ss);
					hh3[ic][gch]->SetBinContent(hh3[ic][gch]->FindBin(i2),abs(ss));
					ss*=-1;
// 					if(ss>ThColInd[ic])
					if(ss>(TPCBaselines[ic][gch][1]*ThSigma[ic]))
					{
// // 						ht3[i2]->SetBinContent(ht3[i2]->FindBin(ic*128+gch),ss);
						hh2[ic][gch]->SetBinContent(hh2[ic][gch]->FindBin(i2),ss);
					}
				}
			}
		}
		QtotCol=0.;
		tminallpt=5000;tmaxallpt=0;
		colminallpt=200;colmaxallpt=0;
		for(int i1=0;i1<2;i1++)
		{
			W_vs_t[i1]->Set(0);
			for(int i2=0;i2<128;i2++)
			{
				for(int is1=1;is1<=hh2[i1][i2]->GetNbinsX();is1++)
				{
					if(hh2[i1][i2]->GetBinContent(is1)>0)
					{
						if(i1==0)
						{
							QtotCol+=hh2[i1][i2]->GetBinContent(is1);
							if(i2>colmaxallpt) colmaxallpt=i2;
							if(i2<colminallpt) colminallpt=i2;
							if(hh2[i1][i2]->GetBinCenter(is1)>tmaxallpt) tmaxallpt=hh2[i1][i2]->GetBinCenter(is1);
							if(hh2[i1][i2]->GetBinCenter(is1)<tminallpt) tminallpt=hh2[i1][i2]->GetBinCenter(is1);
						}
						W_vs_t[i1]->SetPoint(W_vs_t[i1]->GetN(),hh2[i1][i2]->GetBinCenter(is1),i2);
					}
				}
			}
		}
		HQTotCol->Fill(QtotCol);
		NPoints->Fill(W_vs_t[0]->GetN());
		NPointsperArea->Fill(((float)W_vs_t[0]->GetN())/((tmaxallpt-tminallpt)*(colmaxallpt-colminallpt)));
		TotalArea->Fill((tmaxallpt-tminallpt)*(colmaxallpt-colminallpt));
		
		
		
		
// // 		Fill the Hit Candidates collection then induction
		vector <hitcandidateshort> HCS[2][128];
		vector <hitcandidateshort> HCSCand;
		int startbin=0;int nextstartbin=0;int endbin=0;int minbin=0;int maxbin=0;
		int prepeakbin=0;int midpeakbin=0;
		float fitminx=0;float fitmaxx=0;
		QtotTrue=0;
		for(int i2=0;i2<128;i2++)
		{
			QtotTrue+=hh2[0][i2]->Integral();
// 			cout<<"0 "<<i2<<endl;
			startbin=FindBinAbove(hh2[0][i2],0,1);
			while(startbin>=0)
			{
				peakbin=startbin;
				prepeakbin=startbin;
// 				for(int ik1=startbin;ik1<hh2[0][i2]->GetNbinsX();ik1++)
				for(int ik1=startbin;ik1<hh2[0][i2]->FindBin(TimeRange[1]);ik1++)
				{
					if(hh2[0][i2]->GetBinContent(ik1+1)>=hh2[0][i2]->GetBinContent(ik1))
					{
						peakbin=ik1+1;
					}
					else
					{
						for(int ik2=peakbin;ik2>=startbin;ik2--)
						{
							if(hh2[0][i2]->GetBinContent(ik2)==hh2[0][i2]->GetBinContent(peakbin)){prepeakbin=ik2;}
							else{break;}
						}
						maxbin=peakbin;
// 						midpeakbin=(peakbin+prepeakbin)/2;
// 						for(int ik2=peakbin;ik2<hh2[0][i2]->GetNbinsX();ik2++)
						for(int ik2=peakbin;ik2<hh2[0][i2]->FindBin(TimeRange[1]);ik2++)
						{
// 							if(hh2[0][i2]->GetBinContent(ik2+1)>hh2[0][i2]->GetBinContent(ik2))
							if(hh[0][i2]->GetBinContent(ik2+1)>hh[0][i2]->GetBinContent(ik2))
							{
// 								endbin=ik2;
								maxbin=ik2;
// 								nextstartbin=ik2+1;
								break;
							}
// 							else if(hh2[0][i2]->GetBinContent(ik2)==0)
// 							{
// 								endbin=ik2-1;
// 								nextstartbin=FindBinAbove(hh2[0][i2],0,ik2);
// 								break;
// 							}
						}
						minbin=peakbin;
// 						for(int ik2=prepeakbin-1;ik2>1;ik2--)
						for(int ik2=prepeakbin-1;ik2>hh2[0][i2]->FindBin(TimeRange[0]);ik2--)
						{
// 							if(hh[0][i2]->GetBinContent(ik2-1)>=hh[0][i2]->GetBinContent(ik2))
							if(hh[0][i2]->GetBinContent(ik2-1)>hh[0][i2]->GetBinContent(ik2))
							{
								minbin=ik2;
								break;
							}
						}
						for(int ik2=minbin;ik2<peakbin;ik2++)
						{
							if(hh[0][i2]->GetBinContent(ik2+1)>hh[0][i2]->GetBinContent(ik2))
							{
								minbin=ik2+1;
								break;
							}
						}
// 						if((endbin-startbin)>3)
						if((maxbin-minbin)>3)
						{
// 							fitminx=(((prepeakbin-2)>=1)?(hh2[0][i2]->GetBinCenter(prepeakbin-2)-0.1):0.);
							fitminx=(((prepeakbin-2)>=minbin)?(hh2[0][i2]->GetBinCenter(prepeakbin-2)-0.1):(hh2[0][i2]->GetBinCenter(minbin)-0.1));
// 							fitmaxx=(((peakbin+2)<=endbin)?(hh2[0][i2]->GetBinCenter(peakbin+2)+0.1):(hh2[0][i2]->GetBinCenter(endbin)+0.1));
							fitmaxx=(((peakbin+2)<=maxbin)?(hh2[0][i2]->GetBinCenter(peakbin+2)+0.1):(hh2[0][i2]->GetBinCenter(maxbin)+0.1));
							hh2[0][i2]->Fit(gf,"q","q",fitminx,fitmaxx);
							hcs.t=gf->GetParameter(1);
							hcs.Int=0;
							hcs.ci=0;
							hcs.wid=i2;
// 							hcs.w=endbin-startbin+1;
							hcs.w=maxbin-minbin+1;
							
							undershootoffset=0.;
// 							if(hh[0][i2]->GetBinContent(minbin)<0 || hh[0][i2]->GetBinContent(maxbin)<0)
// 							{
// 								undershootoffset=hh[0][i2]->GetBinContent(minbin);
// 								if(hh[0][i2]->GetBinContent(maxbin)<undershootoffset)
// 								{
// 									undershootoffset=hh[0][i2]->GetBinContent(maxbin);
// 								}
// 							}
							
// 							hcs.a=gf->GetParameter(0)+fabs(undershootoffset);
							hcs.a=gf->GetParameter(0);
// 							for(int ik2=startbin;ik2<=endbin;ik2++)
							for(int ik2=minbin;ik2<=maxbin;ik2++)
							{
// 								hcs.Int+=hh2[0][i2]->GetBinContent(ik2);
// 								hcs.Int+=(hh[0][i2]->GetBinContent(ik2)+fabs(undershootoffset));
								if(hh[0][i2]->GetBinContent(ik2)>0) hcs.Int+=(hh[0][i2]->GetBinContent(ik2));
							}
							HCS[0][i2].push_back(hcs);
						}
						nextstartbin=maxbin+1;
						break;
					}
				}
// 				cout<<I<<" 0 "<<i2<<" "<<startbin<<" "<<peakbin<<" "<<minbin<<" "<<maxbin<<" "<<nextstartbin<<" "<<FindBinAbove(hh2[0][i2],0,nextstartbin)<<endl;
// 				startbin=nextstartbin;
				startbin=FindBinAbove(hh2[0][i2],0,nextstartbin);
// 				if((hh2[0][i2]->FindBin(TimeRange[1])-startbin)<10) break;
// 				if(startbin==hh2[0][i2]->FindBin(2500)) startbin=-1;
			}
		}
		for(int i2=0;i2<128;i2++)
		{
// 			cout<<"1 "<<i2<<endl;
			startbin=FindBinAbove(hh3[1][i2],ThColInd[1],1);
			while(startbin>=0)
			{
// 				cout<<I<<" 1 "<<i2<<" "<<startbin<<endl;
				
// 				minbin=((startbin-1)>=1)?(startbin-1):1;
// 				maxbin=((startbin+1)<=hh[1][i2]->GetNbinsX())?(startbin+1):hh[1][i2]->GetNbinsX();
				while(!((hh[1][i2]->GetBinContent(startbin)<hh[1][i2]->GetBinContent(startbin-1))&&(hh[1][i2]->GetBinContent(startbin)>hh[1][i2]->GetBinContent(startbin+1))))
				{
					startbin=FindBinAbove(hh3[1][i2],ThColInd[1],startbin+1);
					if(startbin==-1) break;
				}
				if(startbin==-1) break;
				for(int ik1=startbin;ik1<hh[1][i2]->GetNbinsX();ik1++)
				{
					if(hh[1][i2]->GetBinContent(ik1+1)>=hh[1][i2]->GetBinContent(ik1))
					{
						maxbin=ik1;
						break;
					}
				}
				for(int ik1=startbin;ik1>1;ik1--)
				{
					if(hh[1][i2]->GetBinContent(ik1-1)<=hh[1][i2]->GetBinContent(ik1))
					{
						minbin=ik1;
						break;
					}
				}
				if((maxbin-minbin)>3 && ((hh[1][i2]->GetBinContent(minbin)-hh[1][i2]->GetBinContent(maxbin))>(ThColInd[1]*2.)))
				{
					midpeakbin=(minbin+maxbin)/2;
					fitminx=(((midpeakbin-2)>=minbin)?(hh[1][i2]->GetBinCenter(midpeakbin-2)-0.1):(hh[1][i2]->GetBinCenter(minbin)-0.1));
					fitmaxx=(((midpeakbin+2)<=maxbin)?(hh[1][i2]->GetBinCenter(midpeakbin+2)+0.1):(hh[1][i2]->GetBinCenter(maxbin)+0.1));
					hh[1][i2]->Fit(lf,"q","q",fitminx,fitmaxx);
					hcs.t=(lf->GetParameter(0)/lf->GetParameter(1));
					hcs.Int=0;
					hcs.ci=0;
					hcs.wid=i2;
					hcs.w=maxbin-minbin+1;
					hcs.a=hh[1][i2]->GetBinContent(minbin)-hh[1][i2]->GetBinContent(maxbin);
					HCS[1][i2].push_back(hcs);
				}
				startbin=FindBinAbove(hh3[1][i2],ThColInd[1],maxbin);
			}
		}
		vector <hits> Hits;
// 		for(int i1=0;i1<128;i1++)
// 		{
// 			for(int ik1=0;ik1<HCS[0][i1].size();ik1++)
// 			{
// 				colindmindeltaT=10000.;
// 				nhitspercolpeak=0;
// 				for(int i2=0;i2<128;i2++)
// 				{
// 					for(int ik2=0;ik2<HCS[1][i2].size();ik2++)
// 					{
// 						if(fabs(HCS[0][i1][ik1].t-(HCS[1][i2][ik2].t+5.))<fabs(colindmindeltaT))
// 						{
// 							colindmindeltaT=HCS[0][i1][ik1].t-(HCS[1][i2][ik2].t+5.);
// 						}
// 						if(fabs(HCS[0][i1][ik1].t-(HCS[1][i2][ik2].t+5.))<=2.)
// 						{
// 							hit.col=i1;
// 							hit.ind=i2;
// 							hit.t=((int)HCS[0][i1][ik1].t);
// 							hit.Int=HCS[0][i1][ik1].Int+HCS[1][i2][ik2].Int;
// 							hit.Int2=HCS[1][i2][ik2].Int;
// 							hit.Int3=HCS[0][i1][ik1].Int;
// 							HitCol->Fill(hit.col);
// 							HitInd->Fill(hit.ind);
// 							HitTime->Fill(hit.t);
// 							HitIntegral1->Fill(hit.Int);
// 							HitIntegral2->Fill(hit.Int2);
// 							HitIntegral3->Fill(hit.Int3);
// 							HitOK=true;
// // 							if(hit.col==0 || hit.col==127 || hit.ind==0 || hit.ind==127) HitOK=false;
// 	// 						if(RunNo==7074 && hit.col==93) HitOK=false;
// 	// 						if(RunNo==7076 && hit.col==93) HitOK=false;
// 	// 						if(RunNo==7084 && hit.col==88) HitOK=false;
// 							
// 							hit.colT=HCS[0][i1][ik1].t;
// 							hit.indT=HCS[1][i2][ik2].t;
// 							hit.colWidth=HCS[0][i1][ik1].w;
// 							hit.indWidth=HCS[1][i2][ik2].w;
// 							hit.colAmp=HCS[0][i1][ik1].a;
// 							hit.indAmp=HCS[1][i2][ik2].a;
// 							
// 							if(HitOK)
// 							{
// 								Hits.push_back(hit);
// 							}
// 							nhitspercolpeak++;
// 						}
// 					}
// 				}
// 				ColIndMinDeltaT5->Fill(colindmindeltaT);
// 				NHitsPerColPeak->Fill(nhitspercolpeak);
// 			}
// 		}
		for(int i1=0;i1<128;i1++)
		{
			for(int ik1=0;ik1<HCS[0][i1].size();ik1++)
			{
				colindmindeltaT=10000.;
				nhitspercolpeak=0;
				
				if(HCS[0][i1][ik1].a<(ThColInd[0]*1.5)) continue;
				nhitsprecol=((int)(HCS[0][i1][ik1].Int/400.))+1;
				HCSCand.clear();
				
				for(int i2=0;i2<128;i2++)
				{
					for(int ik2=0;ik2<HCS[1][i2].size();ik2++)
					{
						if(fabs(HCS[0][i1][ik1].t-(HCS[1][i2][ik2].t+HitColIndDist))<fabs(colindmindeltaT))
						{
							colindmindeltaT=HCS[0][i1][ik1].t-(HCS[1][i2][ik2].t+HitColIndDist);
						}
						if(fabs(HCS[0][i1][ik1].t-(HCS[1][i2][ik2].t+HitColIndDist))>2.) continue;
						if(HCSCand.size()==0)
						{
							HCSCand.push_back(HCS[1][i2][ik2]);
							continue;
						}
						for(int ip1=0;ip1<HCSCand.size();ip1++)
						{
							if(fabs(HCS[0][i1][ik1].t-(HCS[1][i2][ik2].t+HitColIndDist))<fabs(HCS[0][i1][ik1].t-(HCSCand[ip1].t+HitColIndDist)))
							{
								HCSCand.push_back(HCS[1][i2][ik2]);
								for(int ip2=ip1+1;ip2<HCSCand.size();ip2++)
								{
									HCSCand[ip2]=HCSCand[ip2-1];
								}
								HCSCand[ip1]=HCS[1][i2][ik2];
								break;
							}
						}
					}
				}
// 				if(nhitsprecol<HCSCand.size()) nhitspostcol=nhitsprecol;
// 				else nhitspostcol=HCSCand.size();
				nhitspostcol=HCSCand.size();
				NHitsPerColPeak->Fill(nhitspostcol);
				if(nhitspostcol==0) continue;
				for(int ip1=0;ip1<HCSCand.size();ip1++)
				{
					hit.col=i1;
					hit.ind=HCSCand[ip1].wid;
					hit.t=((int)HCS[0][i1][ik1].t);
// 					hit.Int=(HCS[0][i1][ik1].Int/((float)nhitspostcol))+HCSCand[ip1].Int;
// 					hit.Int2=HCSCand[ip1].Int;
// 					hit.Int3=(HCS[0][i1][ik1].Int/((float)nhitspostcol));
					hit.Int=(HCS[0][i1][ik1].Int/((float)nhitspostcol));
					hit.Int2=HCS[0][i1][ik1].Int;
					hit.Int3=(HCS[0][i1][ik1].Int/((float)nhitspostcol));
					HitCol->Fill(hit.col);
					HitInd->Fill(hit.ind);
					HitTime->Fill(hit.t);
					HitIntegral1->Fill(hit.Int);
					HitIntegral2->Fill(hit.Int2);
					HitIntegral3->Fill(hit.Int3);
					HitOK=true;
// 					if(hit.col==0 || hit.col==127 || hit.ind==0 || hit.ind==127) HitOK=false;
// 					if(RunNo==7074 && hit.col==93) HitOK=false;
// 					if(RunNo==7076 && hit.col==93) HitOK=false;
// 					if(RunNo==7084 && hit.col==88) HitOK=false;
					
					hit.colT=HCS[0][i1][ik1].t;
					hit.indT=HCSCand[ip1].t;
					hit.colWidth=HCS[0][i1][ik1].w;
					hit.indWidth=HCSCand[ip1].w;
					hit.colAmp=HCS[0][i1][ik1].a;
					hit.indAmp=HCSCand[ip1].a;
					
					if(HitOK)
					{
						Hits.push_back(hit);
					}
					nhitspercolpeak++;
					if(nhitspercolpeak==nhitspostcol) break;
				}
				ColIndMinDeltaT5->Fill(colindmindeltaT);
			}
		}
		
		NHits->Fill(Hits.size());
		
		hd.E=I;
		hd.QColTot=QColTot;
		hd.QColTotZS=QColTotZS;
		hd.ColIndT->clear();
		hd.Int->clear();
		int mins[3]={5000,5000,5000};int maxs[3]={0,0,0};//Col, ind, t
		
// 		cc->cd();
// 		frame3d->Draw();
// // 		pm3d1 = new TPolyMarker3D(H.size());
// 		pm3d1 = new TPolyMarker3D();
		
		QtotHits=0;
		
		for(int i1=0;i1<Hits.size();i1++)
		{
			bb.clear();bc.clear();
			bb.push_back(Hits[i1].col);bb.push_back(Hits[i1].ind);bb.push_back(Hits[i1].t);
			bc.push_back(Hits[i1].Int);bc.push_back(Hits[i1].Int2);bc.push_back(Hits[i1].Int3);
			hd.ColIndT->push_back(bb);
			hd.Int->push_back(bc);
			if(Hits[i1].col<mins[0]) mins[0]=Hits[i1].col;if(Hits[i1].col>maxs[0]) maxs[0]=Hits[i1].col;
			if(Hits[i1].ind<mins[1]) mins[1]=Hits[i1].ind;if(Hits[i1].ind>maxs[1]) maxs[1]=Hits[i1].ind;
			if(Hits[i1].t<mins[2]) mins[2]=Hits[i1].t;if(Hits[i1].t>maxs[2]) maxs[2]=Hits[i1].t;
			
			ColIndDeltaT->Fill(Hits[i1].colT-(Hits[i1].indT+5.));
			ColWidth->Fill(Hits[i1].colWidth);
			IndWidth->Fill(Hits[i1].indWidth);
			ColAmp->Fill(Hits[i1].colAmp);
			IndAmp->Fill(Hits[i1].indAmp);
			Amp_Col_vs_Ind->Fill(Hits[i1].indAmp,Hits[i1].colAmp);
			Width_Col_vs_Ind->Fill(Hits[i1].indWidth,Hits[i1].colWidth);
			QtotHits+=Hits[i1].Int3;
			Col_Int_vs_Amp->Fill(Hits[i1].colAmp,Hits[i1].Int3);
			
			AllXY->Fill(Hits[i1].col,Hits[i1].ind);
			
// 			if(Hits[i1].Int3<30 || Hits[i1].Int3>2000)
// 			{
// 				cout<<I<<" "<<Hits[i1].col<<" "<<Hits[i1].ind<<" "<<Hits[i1].t<<" "<<Hits[i1].Int3<<endl;
// 			}
			
// 			cout<<i1<<" "<<Hits[i1].col<<" "<<Hits[i1].ind<<" "<<Hits[i1].t<<endl;
// 			pm3d1->SetPoint(i1,Hits[i1].col,Hits[i1].ind,Hits[i1].t);
		}
		
		hd.QHitTot=QtotHits;
		
		HitIntegralLost->Fill(QtotTrue-QtotHits);
		TrueTotalIntegral->Fill(QtotTrue);
		AllDeltaT->Fill(maxs[2]-mins[2]);
		
// 		pm3d1->SetMarkerColor(kRed);
// 		pm3d1->SetMarkerStyle(24);   
// 		pm3d1->Draw();
// 		sprintf(hname,"Ev_%d",I);
// 		cc->SetName(hname);cc->SetTitle(hname);
// 		outroot->cd();
// 		cc->Write();
		
		StartTime->Fill(mins[2]);EndTime->Fill(maxs[2]);
		StartCollection->Fill(mins[0]);EndCollection->Fill(maxs[0]);
		StartInduction->Fill(mins[1]);EndInduction->Fill(maxs[1]);
		
		hd.PMTIntegral->clear();
		for(int i1=0;i1<3;i1++)
		{
			isSaturated=false;
			qint=0.;
// 			for(int i2=0;i2<ed.PMTWF->at(i1).size();i2++)
			for(int i2=1500;i2<7500;i2++)
			{
				sg=(-1.)*(((float)ed.PMTWF->at(i1)[i2])-PMTBaselines[i1][0]);
				if(sg>ThPMT[i1]) qint+=sg;
				if(ed.PMTWF->at(i1)[i2]==0){isSaturated=true;}
			}
			if(isSaturated) qint*=-1.;
			hd.PMTIntegral->push_back(qint);
		}
		
		hd.ColID->clear();
		hd.ColT->clear();
		hd.ColA->clear();
		hd.ColInt->clear();
		hd.Colw->clear();
		hd.IndID->clear();
		hd.IndT->clear();
		hd.IndA->clear();
		hd.Indw->clear();
		for(int i1=0;i1<128;i1++)
		{
			for(int ik1=0;ik1<HCS[0][i1].size();ik1++)
			{
				hd.ColID->push_back(HCS[0][i1][ik1].wid);
				hd.ColT->push_back(HCS[0][i1][ik1].t);
				hd.ColA->push_back(HCS[0][i1][ik1].a);
				hd.ColInt->push_back(HCS[0][i1][ik1].Int);
				hd.Colw->push_back(HCS[0][i1][ik1].w);
			}
			for(int ik1=0;ik1<HCS[1][i1].size();ik1++)
			{
				hd.IndID->push_back(HCS[1][i1][ik1].wid);
				hd.IndT->push_back(HCS[1][i1][ik1].t);
				hd.IndA->push_back(HCS[1][i1][ik1].a);
				hd.Indw->push_back(HCS[1][i1][ik1].w);
			}
		}
		
		if(W_vs_t[0]->GetN()<1000) hd.EventType=3;
		else if(W_vs_t[0]->GetN()>4000) hd.EventType=2;
		else	//going into track fitting
		{
			hd.EventType=1;//assume partial traversing, then identify
			vector <clusters> CL;
			cl.hi.clear();
			inCL.clear();
			for(int i1=0;i1<W_vs_t[0]->GetN();i1++)
			{
				inCL.push_back(0);
			}
			for(int i1=0;i1<W_vs_t[0]->GetN();i1++)
			{
				if(inCL[i1]!=0) continue;
				cl.hi.clear();
				cl.hi.push_back(i1);
				inCL[i1]=1;
				nid=1;
				while(nid>0)
				{
					nid=0;
					for(int i2=0;i2<W_vs_t[0]->GetN();i2++)
					{
						if(inCL[i2]!=0) continue;
						W_vs_t[0]->GetPoint(i2,x2,y2);
						for(int i3=0;i3<cl.hi.size();i3++)
						{
							W_vs_t[0]->GetPoint(cl.hi[i3],x1,y1);
							if(fabs(x1-x2)<=1 && fabs(y1-y2)<=1)
							{
								cl.hi.push_back(i2);
								nid++;
								inCL[i2]=1;
								break;
							}
						}
					}
				}
				CL.push_back(cl);
				cl.hi.clear();
			}
			maxclsize=0;
			for(int i1=0;i1<CL.size();i1++)
			{
				tmincl=5000;tmaxcl=0;
				colmincl=200;colmaxcl=0;
				for(int i2=0;i2<CL[i1].hi.size();i2++)
				{
					W_vs_t[0]->GetPoint(CL[i1].hi[i2],x1,y1);
					if(x1<tmincl) tmincl=x1;if(x1>tmaxcl) tmaxcl=x1;
					if(y1<colmincl) colmincl=y1;if(y1>colmaxcl) colmaxcl=y1;
	// 				if(x1<tminallpt) tminallpt=x1;if(x1>tmaxallpt) tmaxallpt=x1;
	// 				if(y1<colminallpt) colminallpt=y1;if(y1>colmaxallpt) colmaxallpt=y1;
				}
				ClusterNPointsperArea->Fill(((float)CL[i1].hi.size())/((tmaxcl-tmincl)*(colmaxcl-colmincl)));
				ClSize->Fill(CL[i1].hi.size());
				if(CL[i1].hi.size()>maxclsize) {maxclsize=CL[i1].hi.size();indCLHP=i1;}
			}
			MaxClSize->Fill(maxclsize);
			NClperEvent->Fill(CL.size());
			Col_vs_t_fit->Set(0);
			tmincl=5000;tmaxcl=0;
			colmincl=200;colmaxcl=0;
			for(int i1=0;i1<CL[indCLHP].hi.size();i1++)
			{
				W_vs_t[0]->GetPoint(CL[indCLHP].hi[i1],x1,y1);
				Col_vs_t_fit->SetPoint(Col_vs_t_fit->GetN(),x1,y1);
				if(x1>tmaxcl) tmaxcl=x1;if(x1<tmincl) tmincl=x1;
				if(y1>colmaxcl) colmaxcl=y1;if(y1<colmincl) colmincl=y1;
			}
			tflin->SetParameter(0,colmincl);
			tflin->SetParameter(1,((float)(colmaxcl-colmincl))/(tmaxcl-tmincl));
			Col_vs_t_fit->Fit(tflin,"q","q",tmincl-1,tmaxcl+1);
			fitnormchi2col=((float)tflin->GetChisquare())/((float)tflin->GetNDF());
			tflinpars[0][0]=tflin->GetParameter(0);tflinpars[0][1]=tflin->GetParameter(1);
			
			FitNormChi2_Col_vs_t->Fill(fitnormchi2col);
			ColTrackDeltaT->Fill(tmaxcl-tmincl);
			ClmaxPointFraction->Fill(((float)CL[indCLHP].hi.size())/((float)W_vs_t[0]->GetN()));
			
// 			cout<<I<<" "<<(((float)CL[indCLHP].hi.size())/((float)W_vs_t[0]->GetN()))<<endl;
			if((((float)CL[indCLHP].hi.size())/((float)W_vs_t[0]->GetN()))>0.8)
// 			if(fitnormchi2col<20)
			{
// 				cout<<I<<" passed point fraction"<<endl;
				//for the track
				tmintrack=5000;tmaxtrack=0;
				indoftmin=0;indoftmax=0;Nvalid=0;
				
				tg2d->Set(0);
				trackXY[0].clear();trackXY[1].clear();
				for(int i1=0;i1<hd.ColIndT->size();i1++)
				{
					if(float(hd.ColIndT->at(i1)[2])<tmincl || float(hd.ColIndT->at(i1)[2])>tmaxcl) continue;
					x1=tflinpars[0][0]+tflinpars[0][1]*((float)hd.ColIndT->at(i1)[2]);
					if(fabs(x1-((float)hd.ColIndT->at(i1)[0]))<=5.)
					{
						tg2d->SetPoint(tg2d->GetN(),hd.ColIndT->at(i1)[0],hd.ColIndT->at(i1)[1],hd.ColIndT->at(i1)[2]);
						Nvalid++;
						if(((float)hd.ColIndT->at(i1)[2])>tmaxtrack) {tmaxtrack=hd.ColIndT->at(i1)[2];indoftmax=i1;}
						if(((float)hd.ColIndT->at(i1)[2])<tmintrack) {tmintrack=hd.ColIndT->at(i1)[2];indoftmin=i1;}
						trackXY[0].push_back(hd.ColIndT->at(i1)[0]);
						trackXY[1].push_back(hd.ColIndT->at(i1)[1]);
					}
				}
				FitClSize->Fill(Nvalid);
				if(Nvalid>=5)
				{
// 					cout<<I<<" passed nvalid"<<endl;
					pStart[0]=hd.ColIndT->at(indoftmin)[0];
					pStart[1]=((float)(hd.ColIndT->at(indoftmax)[0]-hd.ColIndT->at(indoftmin)[0]))/((float)(hd.ColIndT->at(indoftmax)[2]-hd.ColIndT->at(indoftmin)[2]));
					pStart[2]=hd.ColIndT->at(indoftmin)[1];
					pStart[3]=((float)(hd.ColIndT->at(indoftmax)[1]-hd.ColIndT->at(indoftmin)[1]))/((float)(hd.ColIndT->at(indoftmax)[2]-hd.ColIndT->at(indoftmin)[2]));
					pStart[4]=750;
					
					fitnormchi2=line3Dfit(tg2d);
					
					if(fitnormchi2<10)
					{
// 						cout<<I<<" passed normchi2"<<endl;
						TrackDeltaT->Fill(tmaxtrack-tmintrack);
						FitNormChi2->Fill(fitnormchi2);
						
						
						td.E=hd.E;
						td.StartEndColIndT->clear();
						td.StartEndColIndT->push_back(hd.ColIndT->at(indoftmin)[0]);
						td.StartEndColIndT->push_back(hd.ColIndT->at(indoftmin)[1]);
						td.StartEndColIndT->push_back(hd.ColIndT->at(indoftmin)[2]);
						td.StartEndColIndT->push_back(hd.ColIndT->at(indoftmax)[0]);
						td.StartEndColIndT->push_back(hd.ColIndT->at(indoftmax)[1]);
						td.StartEndColIndT->push_back(hd.ColIndT->at(indoftmax)[2]);
						td.FitParams->clear();
						for(int is2=0;is2<5;is2++)
						{
							td.FitParams->push_back(fitParams[is2]);
						}
						td.FitNormChi2=fitnormchi2;
						td.QColTot=QtotCol;
						td.NHits=Nvalid;
						td.Nexcl=hd.ColIndT->size()-Nvalid;
						td.PMTIntegral->clear();
						for(int is2=0;is2<3;is2++)
						{
							td.PMTIntegral->push_back(hd.PMTIntegral->at(is2));
						}
						
						td.ColTStartEnd->clear();
						td.ColHitTStartEnd->clear();
						
						td.ColTStartEnd->push_back(tmincl);
						td.ColTStartEnd->push_back(tmaxcl);
						
						colhittmin=5000;colhittmax=0;
						for(int i1=0;i1<128;i1++)
						{
							for(int ik1=0;ik1<HCS[0][i1].size();ik1++)
							{
								x1=tflinpars[0][0]+tflinpars[0][1]*HCS[0][i1][ik1].t;
								if(fabs(x1-((float)i1))<=5)
								{
									if(HCS[0][i1][ik1].t<colhittmin) colhittmin=HCS[0][i1][ik1].t;
									if(HCS[0][i1][ik1].t>colhittmax) colhittmax=HCS[0][i1][ik1].t;
								}
							}
						}
						td.ColHitTStartEnd->push_back(colhittmin);
						td.ColHitTStartEnd->push_back(colhittmax);
						
						Ttr->Fill();
						
						for(int ip1=0;ip1<trackXY[0].size();ip1++)
						{
							TrackXY->Fill(trackXY[0][ip1],trackXY[1][ip1]);
						}
						
// 						//track timing
// 						x1=fitParams[0]+fitParams[1]*(tmintrack-fitParams[4]);
// 						y1=fitParams[2]+fitParams[3]*(tmintrack-fitParams[4]);
// 						z1=tmintrack;
// 						
// 						x2=fitParams[0]+fitParams[1]*(tmaxtrack-fitParams[4]);
// 						y2=fitParams[2]+fitParams[3]*(tmaxtrack-fitParams[4]);
// 						z2=tmaxtrack;
						
						//collection hit timing
						x1=fitParams[0]+fitParams[1]*(colhittmin-fitParams[4]);
						y1=fitParams[2]+fitParams[3]*(colhittmin-fitParams[4]);
						z1=colhittmin;
						
						x2=fitParams[0]+fitParams[1]*(colhittmax-fitParams[4]);
						y2=fitParams[2]+fitParams[3]*(colhittmax-fitParams[4]);
						z2=colhittmax;
						
						
			// 			DeltaTvsTrackAngle->Fill(atan(sqrt(pow(x2-x1,2)+pow(y2-y1,2))/(tmax-tmin)),(tmax-tmin));
			// 			DeltaTvsTrackAnglePr->Fill(atan(sqrt(pow(x2-x1,2)+pow(y2-y1,2))/(tmax-tmin)),(tmax-tmin));
						Tfirst->Fill(tmintrack);
						Tlast->Fill(tmaxtrack);
						TColfirst->Fill(colhittmin);
						TCollast->Fill(colhittmax);
						Trackt0->Fill(fitParams[4]);
						TrackEntryXY->Fill(x1,y1);
						TrackExitXY->Fill(x2,y2);
			// 			if(tmin<=1160 && tmax>=1920) TrackDeltaTSel->Fill(tmax-tmin);
						DeltaTvsTrackt0->Fill(fitParams[4],tmaxtrack-tmintrack);
			// 			DeltaTLengthvsAngle->Fill(atan(sqrt(pow(x2-x1,2)+pow(y2-y1,2))/(tmax-tmin)),sqrt(pow(x2-x1,2)+pow(y2-y1,2)+pow(tmax-tmin,2)),tmax-tmin);
			// 			DeltaTLengthvsAngle_norm->Fill(atan(sqrt(pow(x2-x1,2)+pow(y2-y1,2))/(tmax-tmin)),sqrt(pow(x2-x1,2)+pow(y2-y1,2)+pow(tmax-tmin,2)));
						NExclHits->Fill(td.Nexcl);
						
						if((x1>1 && x1<127) && (x2>1 && x2<127) && (y1>1 && y1<127) && (y2>1 && y2<127))
						{
							hd.EventType=0;
							DeltaTTrackColT->Fill(tmaxcl-tmincl);
							DeltaTTrackColHitT->Fill(colhittmax-colhittmin);
							
							for(int i1=0;i1<128;i1++)
							{
								for(int ik1=0;ik1<HCS[0][i1].size();ik1++)
								{
									x1=tflinpars[0][0]+tflinpars[0][1]*HCS[0][i1][ik1].t;
									if(fabs(x1-((float)i1))<=5)
									{
										QperT[0]->Fill(HCS[0][i1][ik1].t,HCS[0][i1][ik1].Int);
										QperT[2]->Fill(HCS[0][i1][ik1].t,HCS[0][i1][ik1].Int);
										QperT[1]->Fill(HCS[0][i1][ik1].t);
									}
								}
							}
						}
// 						if(td.NHits>10 && td.Nexcl<10)
						if(td.Nexcl<10)
						{
// 							outfile<<td.E<<" "<<hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[0]<<" "<<hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[1]<<" "<<hd.ColIndT->at(CL[indCLHP].hi[indoftmin])[2]<<" "<<hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[0]<<" "<<hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[1]<<" "<<hd.ColIndT->at(CL[indCLHP].hi[indoftmax])[2]<<endl;
							
							outfile<<td.E<<" "<<x1<<" "<<y1<<" "<<z1<<" "<<x2<<" "<<y2<<" "<<z2<<endl;
							
						}
					}
					else hd.EventType=3;
				}
				else hd.EventType=3;
			}
			else hd.EventType=3;
			
// 			if(I<50)
// 			{
// // 				sprintf(hname,"Col_vs_t_%d",I);W_vs_t[0]->SetName(hname);W_vs_t[0]->SetTitle(hname);
// // 				sprintf(hname,"Ind_vs_t_%d",I);W_vs_t[1]->SetName(hname);W_vs_t[1]->SetTitle(hname);
// // 				outroot->cd();
// // 				W_vs_t[0]->Write();
// // 				W_vs_t[1]->Write();
// 				
// 				sprintf(hname,"Col_vs_t_%d",I);Col_vs_t_fit->SetName(hname);Col_vs_t_fit->SetTitle(hname);
// 				outroot->cd();
// 				Col_vs_t_fit->Write();
// 			}
		}
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		T1->Fill();
		
// 		if((RunNo==7268 && I==77) || (RunNo==7285 && I==2))
// 		{
// 			outroot->cd("Histos");
// 			for(int i1=0;i1<2;i1++)
// 			{
// 				for(int i2=0;i2<128;i2++)
// 				{
// 					hh2[i1][i2]->Write();
// 				}
// 			}
// 			for(int i1=750;i1<2500;i1++)
// 			{
// 				ht3[i1]->Write();
// 			}
// 		}
		
// 		for(int i1=0;i1<4096;i1++)
		for(int i1=750;i1<2500;i1++)
		{
			ht[i1]->Reset();
			ht2[i1]->Reset();
			ht3[i1]->Reset();
		}
		for(int i1=0;i1<2;i1++)
		{
			for(int i2=0;i2<128;i2++)
			{
				hh[i1][i2]->Reset();
				hh2[i1][i2]->Reset();
				hh3[i1][i2]->Reset();
			}
		}
		Hits.clear();
// 		for(int i1=0;i1<256;i1++)
// 		{
// 			hcw[i1].clear();
// 		}
		if(I%100==0) cout<<"Event: "<<I<<endl;
// 		if(I==50) break;
	}
	outroot->cd();
	T1->Write();
	Ttr->Write();
	gStyle->SetOptStat(0);
	outroot->cd("Histos");
	HitCol->Write();
	HitInd->Write();
	HitTime->Write();
	HitIntegral1->Write();
	HitIntegral2->Write();
	HitIntegral3->Write();
	NHits->Write();
	StartTime->Write();
	EndTime->Write();
	StartCollection->Write();
	EndCollection->Write();
	StartInduction->Write();
	EndInduction->Write();
	NHits_T_vs_E->Write();
	NHits_Col_vs_Ind->Write();
	ColTWidth->Write();
	IndTWidth->Write();
	FitNormChi2_Col_vs_t->Write();
	ClmaxPointFraction->Write();
	TrackDeltaT->Write();
	FitNormChi2->Write();
	FitClSize->Write();
	NExclHits->Write();
	
	Tfirst->Write();
	Tlast->Write();
	TColfirst->Write();
	TCollast->Write();
	Trackt0->Write();
	TrackEntryXY->Write();
	TrackExitXY->Write();
	TrackXY->Write();
	DeltaTvsTrackt0->Write();
	DeltaTTrackColT->Write();
	DeltaTTrackColHitT->Write();
	
	AllDeltaT->Write();
	AllXY->Write();
	
	
	HQTotCol->Write();
	NPoints->Write();
	ClusterNPointsperArea->Write();
	NPointsperArea->Write();
	TotalArea->Write();
	
	NClperEvent->Write();
	ClSize->Write();
	MaxClSize->Write();
	
// 	TrackDeltaT->Write();
	ColTrackDeltaT->Write();
	QperT[0]->Write();
	QperT[1]->Write();
	QperT[2]->Divide(QperT[1]);
	QperT[2]->Write();
	
// 	TH1F* TrackDeltaTSel=new TH1F("TrackDeltaTSel","TrackDeltaTSel",500,0,5000);
// 	TH1F* AllDeltaT=new TH1F("AllDeltaT","AllDeltaT",500,0,5000);
// 	TH2F* AllXY=new TH2F("AllXY","AllXY",128,-0.5,127.5,128,-0.5,127.5);
// 	TH2F* TrackXY=new TH2F("TrackXY","TrackXY",128,-0.5,127.5,128,-0.5,127.5);
// 	TH2F* DeltaTvsTrackt0=new TH2F("DeltaTvsTrackt0","DeltaTvsTrackt0",500,0,5000,500,0,1000);
// 	TH1F* DeltacolperDeltat=new TH1F("DeltacolperDeltat","DeltacolperDeltat",500,0,50);
// 	TH1F* DeltaindperDeltat=new TH1F("DeltaindperDeltat","DeltaindperDeltat",500,0,50);
// 	TH1F* DeltaRperDeltat=new TH1F("DeltaRperDeltat","DeltaRperDeltat",2000,0,100);
// 	TH1F* ColDeltaRR=new TH1F("ColDeltaR","ColDeltaR",2000,0,2000);
// 	TH1F* HitDeltaRR=new TH1F("HitDeltaRR","HitDeltaRR",2000,0,2000);
// 	TH1F* HitDeltaT=new TH1F("HitDeltaT","HitDeltaT",2000,0,1000);
// 	TH1F* HitDeltaR=new TH1F("HitDeltaR","HitDeltaR",200,0,100);
// 	TH2F* HitDeltaR_vs_DeltaT=new TH2F("HitDeltaR_vs_DeltaT","HitDeltaR_vs_DeltaT",200,-0.5,199.5,200,0,20);
	
// 	TH2F* DeltaTLengthvsAngle=new TH2F("DeltaTLengthvsAngle","DeltaTLengthvsAngle",314,0,1.57,100,0,2000);
// 	TH2F* DeltaTLengthvsAngle_norm=new TH2F("DeltaTLengthvsAngle_norm","DeltaTLengthvsAngle_norm",314,0,1.57,100,0,2000);
	
// 	TH1F* FitNormChi2_Ind_vs_t=new TH1F("FitNormChi2_Ind_vs_t","FitNormChi2_Ind_vs_t",1000,0.,4000.);
	
// 	TH2F* NColHits_vs_fitnormchi2=new TH2F("NColHits_vs_fitnormchi2","NColHits_vs_fitnormchi2",200,-0.5,199.5,500,0,2000);
// 	TH2F* NIndHits_vs_fitnormchi2=new TH2F("NIndHits_vs_fitnormchi2","NIndHits_vs_fitnormchi2",200,-0.5,199.5,500,0,2000);
	
	
// 	TH2F* NHits_T_vs_E=new TH2F("NHits_T_vs_E","NHits_T_vs_E",T->GetEntries(),-0.5,T->GetEntries()-0.5,4096,-0.5,4095.5);
// 	TH2F* NHits_Col_vs_Ind=new TH2F("NHits_Col_vs_Ind","NHits_Col_vs_Ind",20,-0.5,19.5,20,-0.5,19.5);
// 	TH2F* QCol_vs_QInd=new TH2F("QCol_vs_QInd","QCol_vs_QInd",50,0,200,50,0,200);
	
	
	ColIndDeltaT->Write();
	ColIndMinDeltaT5->Write();
	ColWidth->Write();
	IndWidth->Write();
	ColAmp->Write();
	IndAmp->Write();
	Amp_Col_vs_Ind->Write();
	Width_Col_vs_Ind->Write();
	HitIntegralLost->Write();
	TrueTotalIntegral->Write();
	Col_Int_vs_Amp->Write();
	NHitsPerColPeak->Write();
	DeltaTTrackColT->Write();
	DeltaTTrackColHitT->Write();
	
	outroot->Close();
// 	outroottr->cd();
// 	Ttr->Write();
// 	outroottr->Close();
	inroot->Close();
	outfile.close();
	
// 	sprintf(hname,"cp Hits_%d.root %s/Histos/Hits_%d.root;wait;",RunNo,AnalysisFilePath,RunNo);system(hname);
// 	sprintf(hname,"cp Tracks_%d.root %s/Histos/Tracks_%d.root",RunNo,AnalysisFilePath,RunNo);system(hname);
	sprintf(hname,"cp HitsAndTracks_%d.root %s/Histos/HitsAndTracks_%d.root;wait;",RunNo,AnalysisFilePath,RunNo);system(hname);
	sprintf(hname,"cp ValidTrackParameters_%d.txt %s/Files/ValidTrackParameters_%d.txt;wait;",RunNo,AnalysisFilePath,RunNo);system(hname);
}

void HitAnalysis()
{
	sprintf(hname,"%s/Run_%d.root",NTupleFilePath,RunNo);
	TFile* inroot=new TFile(hname);
	TTree* T =  (TTree*) inroot->Get("T");
	T->SetBranchAddress("TPCWF",&ed.TPCWF);
	T->SetBranchAddress("PMTWF",&ed.PMTWF);
	
	ReadTPCBaselines();
	ReadPMTBaselines();
	
	sprintf(hname,"HitAnalysis_%d.root",RunNo);
	TFile* outroot=new TFile(hname,"recreate");
	TH1F* FirstColIndTimeDiff=new TH1F("FirstColIndTimeDiff","FirstColIndTimeDiff",39,-19.5,19.5);
	TH1F* InterIndTimeDiff=new TH1F("InterIndTimeDiff","InterIndTimeDiff",99,-49.5,49.5);
	
	TH1F* NColperWT=new TH1F("NColperWt","NColperWt",150,-0.5,149.5);
	TH1F* NColperTW=new TH1F("NColpertW","NColpertW",150,-0.5,149.5);
	TH1F* ColWidthZSWT=new TH1F("ColWidthZSWt","ColWidthZSWt",250,-0.5,249.5);
	TH1F* ColWidthZSTW=new TH1F("ColWidthZStW","ColWidthZStW",250,-0.5,249.5);
	TH1F* ColIntZS=new TH1F("ColIntZS","ColIntZS",1000,0,100000);
	TH1F* ColIntZSWT=new TH1F("ColIntZSWT","ColIntZSWT",1000,0,100000);
	TH1F* ColIntZSTW=new TH1F("ColIntZSTW","ColIntZSTW",1000,0,100000);
	TH1F* QTotColZSWT=new TH1F("QTotColZSWt","QTotColZSWt",1000,0.,10000.);
	TH1F* QTotColZSTW=new TH1F("QTotColZStW","QTotColZStW",1000,0.,10000.);
	
	TGraph* tgint1=new TGraph();
	tgint1->SetName("ColIntegralWt");tgint1->SetTitle("ColIntegralWt");
	tgint1->SetMarkerStyle(25);tgint1->SetMarkerColor(1);
	tgint1->GetXaxis()->SetTitle("Event ID");tgint1->GetXaxis()->CenterTitle();
	tgint1->GetYaxis()->SetTitle("Total Charge (arbitrary units)");tgint1->GetYaxis()->CenterTitle();
	
	TGraph* tgint2=new TGraph();
	tgint2->SetName("ColIntegraltW");tgint2->SetTitle("ColIntegraltW");
	tgint2->SetMarkerStyle(26);tgint2->SetMarkerColor(2);
	tgint2->GetXaxis()->SetTitle("Event ID");tgint2->GetXaxis()->CenterTitle();
	tgint2->GetYaxis()->SetTitle("Total Charge (arbitrary units)");tgint2->GetYaxis()->CenterTitle();
	
	TGraph* tgint3=new TGraph();
	tgint3->SetName("ColIntegralTrue");tgint3->SetTitle("ColIntegralTrue");
	tgint3->SetMarkerStyle(24);tgint3->SetMarkerColor(4);
	tgint3->GetXaxis()->SetTitle("Event ID");tgint3->GetXaxis()->CenterTitle();
	tgint3->GetYaxis()->SetTitle("Total Charge (arbitrary units)");tgint3->GetYaxis()->CenterTitle();
	
	TCanvas* cc1=new TCanvas("ColWires","ColWires",600,600);
	TCanvas* cc2=new TCanvas("IndWiresInv","IndWires",600,600);
	TCanvas* cc3=new TCanvas("IndWires","IndWires",600,600);
	
	TSpectrum *s;
	Double_t *xpeaks;Double_t *xpeaks2;
	TH1F* hh[2][128];TH1F* hh2[2][128];TH1F* hh3[2][128];
	TH1F* ht[4096];TH1F* ht2[4096];TH1F* ht3[4096];
	int wcolors[5]={1,2,3,4,6};
	
	for(int i1=0;i1<2;i1++)
	{
		for(int i2=0;i2<128;i2++)
		{
			sprintf(hname,"WF_%s_%d",wt[i1].c_str(),i2);
			hh[i1][i2]=new TH1F(hname,hname,4096,-0.5,4095.5);
			hh[i1][i2]->SetLineColor(wcolors[i2%5]);
			sprintf(hname,"WFZS_%s_%d",wt[i1].c_str(),i2);
			hh2[i1][i2]=new TH1F(hname,hname,4096,-0.5,4095.5);
			hh2[i1][i2]->SetLineColor(wcolors[i2%5]);
			sprintf(hname,"WFInvZS_%s_%d",wt[i1].c_str(),i2);
			hh3[i1][i2]=new TH1F(hname,hname,4096,-0.5,4095.5);
			hh3[i1][i2]->SetLineColor(wcolors[i2%5]);
		}
	}
	for(int i1=0;i1<4096;i1++)
	{	
		sprintf(hname,"WFt_%d",i1);
		ht[i1]=new TH1F(hname,hname,256,-0.5,255.5);
		sprintf(hname,"WFtZS_%d",i1);
		ht2[i1]=new TH1F(hname,hname,256,-0.5,255.5);
		sprintf(hname,"WFtZSIndInv_%d",i1);
		ht3[i1]=new TH1F(hname,hname,256,-0.5,255.5);
	}
	
	TGraph2D* tg2d=new TGraph2D();
	sprintf(hname,"Ev3d");
	tg2d->SetName(hname);tg2d->SetTitle(hname);
	tg2d->SetMarkerStyle(20);tg2d->SetMarkerColor(4);
	int nh=0;
	vector <hits> H;
	TH3F *frame3d=new TH3F("frame3d","frame3d",1,-0.5,127.5,1,-0.5,127.5,1,-0.5,4095.5);
	TPolyMarker3D *pm3d1;
	TCanvas* cc=new TCanvas("cc","cc",600,600);
	cc->cd();
	gStyle->SetOptStat(0);
	
	bool isSaturated=false;
	float sg=0.;
	float qint=0.;
	
	TGraph* tgacol=new TGraph();tgacol->SetPoint(tgacol->GetN(),-0.5,-20);tgacol->SetPoint(tgacol->GetN(),4095.5,100);
	TGraph* tgaind=new TGraph();tgaind->SetPoint(tgaind->GetN(),-0.5,-30);tgaind->SetPoint(tgaind->GetN(),4095.5,30);
	
	vector <int> bb;
	vector <float> bc;
	int Ev=0;
	int integral=0;int integral2=0;int integral3=0;
	int peakbin=0;int binmin=0;int binmax=0;
	int ic=0;int gch=0;
	float QColTot=0;
	float QColTotZSWt=0;
	float QColTotZStW=0;
	int minTcol=0;int minTind=0;int minTindDiff=0;
	int ncolpeaks=0;
	for(int I=0;I<T->GetEntries();I++)
	{
		T->GetEntry(I);
		int npeaks = 200;
		int nfound=0;
		int nfound2=0;
		float ss=0;
		bool found=true;
		bool HitOK=true;
		
		for(int i1=0;i1<256;i1++)
		{
			ic=i1/128;gch=i1-((i1/128)*128);
			for(int i2=0;i2<4096;i2++)
			{
				ss=((float)ed.TPCWF->at(i1)[i2])-TPCBaselines[ic][gch][0];
				
				hh[ic][gch]->SetBinContent(hh[ic][gch]->FindBin(i2),ss);
				if(ss>ThColInd[ic]) {hh2[ic][gch]->SetBinContent(hh2[ic][gch]->FindBin(i2),ss);}
				
				ht[i2]->SetBinContent(ht[i2]->FindBin(ic*128+gch),ss);
				if(ss>ThColInd[ic]) {ht2[i2]->SetBinContent(ht2[i2]->FindBin(ic*128+gch),ss);}
				
				if(i1<128)
				{
					if(ss>ThColInd[ic])
					{
						ht3[i2]->SetBinContent(ht3[i2]->FindBin(ic*128+gch),ss);
					}
				}
				else
				{
					ss*=-1;
					if(ss>ThColInd[ic])
					{
						hh3[ic][gch]->SetBinContent(hh3[ic][gch]->FindBin(i2),ss);
						ht3[i2]->SetBinContent(ht3[i2]->FindBin(ic*128+gch),ss);
					}
				}
			}
		}
		QColTot=0;
		for(int i1=0;i1<128;i1++)
		{
			QColTot+=hh2[0][i1]->Integral();
		}
		tgint3->SetPoint(tgint3->GetN(),I,QColTot);
		ColIntZS->Fill(QColTot);
		
		minTcol=5000;minTind=5000;
		QColTotZSWt=0;
		for(int i1=0;i1<128;i1++)
		{
			s = new TSpectrum(2*npeaks);
			nfound = s->Search(hh2[0][i1],3.,"",0.1);
			xpeaks = s->GetPositionX();
			NColperWT->Fill(nfound);
			for(int is1=0;is1<nfound;is1++)
			{
				if(xpeaks[is1]<minTcol){minTcol=xpeaks[is1];}
				peakbin=hh2[0][i1]->FindBin(xpeaks[is1]);
				integral3=0;
				for(int ip1=peakbin;ip1<=hh2[0][i1]->GetNbinsX();ip1++)
				{
					if(hh2[0][i1]->GetBinContent(ip1)>ThColInd[0])
					{
						integral3+=hh2[0][i1]->GetBinContent(ip1);
						binmax=ip1;
					}
					else break;
				}
				for(int ip1=peakbin-1;ip1>=1;ip1--)
				{
					if(hh2[0][i1]->GetBinContent(ip1)>ThColInd[0])
					{
						integral3+=hh2[0][i1]->GetBinContent(ip1);
						binmin=ip1;
					}
					else break;
				}
				ColWidthZSWT->Fill(binmax-binmin+1);
				QTotColZSWT->Fill(integral3);
				QColTotZSWt+=integral3;
			}
			delete s;
		}
		ColIntZSWT->Fill(QColTotZSWt);
		tgint1->SetPoint(tgint1->GetN(),I,QColTotZSWt);
		for(int i2=0;i2<128;i2++)
		{
			s = new TSpectrum(2*npeaks);
			nfound = s->Search(hh3[1][i2],3.,"",0.1);
			xpeaks = s->GetPositionX();
			nfound2 = s->Search(hh2[1][i2],3.,"",0.1);
			xpeaks2 = s->GetPositionX();
			minTindDiff=1000;
			for(int is1=0;is1<nfound;is1++)
			{
				if(xpeaks[is1]<minTind){minTind=xpeaks[is1];}
				for(int is2=0;is2<nfound2;is2++)
				{
					if(abs(xpeaks[is1]-xpeaks2[is2])<minTindDiff)
					{
						minTindDiff=abs(xpeaks[is1]-xpeaks2[is2]);
					}
				}
				if(minTindDiff<1000)
				{
					InterIndTimeDiff->Fill(minTindDiff);
				}
			}
			delete s;
		}
		FirstColIndTimeDiff->Fill(minTcol-minTind);
		
		QColTotZStW=0;
		for(int i1=0;i1<4096;i1++)
		{
			s = new TSpectrum(2*npeaks);
			nfound = s->Search(ht3[i1],3.,"",0.1);
			xpeaks = s->GetPositionX();
			ncolpeaks=0;
			for(int is1=0;is1<nfound;is1++)
			{
				if(xpeaks[is1]>=128) continue;
				ncolpeaks++;
				peakbin=ht3[i1]->FindBin(xpeaks[is1]);
				integral3=0;
				for(int ip1=peakbin;ip1<=128;ip1++)
				{
					if(ht3[i1]->GetBinContent(ip1)>ThColInd[0])
					{
						integral3+=ht3[i1]->GetBinContent(ip1);
						binmax=ip1;
					}
					else break;
				}
				for(int ip1=peakbin-1;ip1>=1;ip1--)
				{
					if(ht3[i1]->GetBinContent(ip1)>ThColInd[0])
					{
						integral3+=ht3[i1]->GetBinContent(ip1);
						binmin=ip1;
					}
					else break;
				}
				ColWidthZSTW->Fill(binmax-binmin+1);
				QTotColZSTW->Fill(integral3);
				QColTotZStW+=integral3;
			}
			NColperTW->Fill(ncolpeaks);
			delete s;
		}
		ColIntZSTW->Fill(QColTotZStW);
		tgint2->SetPoint(tgint2->GetN(),I,QColTotZStW);
		
// 		outroot->cd();
// 		cc1->cd();
// 		tgacol->Draw("AP");
// 		for(int i2=0;i2<128;i2++)
// 		{
// 			hh2[0][i2]->Draw("same");
// 		}
// 		sprintf(hname,"Ev_%d_Col_%d",I,minTcol);
// 		cc1->SetName(hname);
// 		cc1->Write();
// 		
// 		cc2->cd();
// 		tgaind->Draw("AP");
// 		for(int i2=0;i2<128;i2++)
// 		{
// 			hh2[1][i2]->Draw("same");
// 		}
// 		sprintf(hname,"Ev_%d_Ind_%d",I,minTind);
// 		cc2->SetName(hname);
// 		cc2->Write();
// 		
// 		cc3->cd();
// 		tgaind->Draw("AP");
// 		for(int i2=0;i2<128;i2++)
// 		{
// 			hh[1][i2]->Draw("same");
// 		}
// 		sprintf(hname,"Ev_%d_IndFull_%d",I,minTind);
// 		cc3->SetName(hname);
// 		cc3->Write();
// 		
// 		FirstColIndTimeDiff->Fill(minTcol-minTind);
		
		for(int i1=0;i1<2;i1++)
		{
			for(int i2=0;i2<128;i2++)
			{
				hh[i1][i2]->Reset();
				hh2[i1][i2]->Reset();
				hh3[i1][i2]->Reset();
			}
		}
		for(int i1=0;i1<4096;i1++)
		{
			ht[i1]->Reset();
			ht2[i1]->Reset();
			ht3[i1]->Reset();
		}
		if(I%10==0) cout<<"Event: "<<I<<endl;
// 		if(I==200) break;
	}
	outroot->cd();
	FirstColIndTimeDiff->Write();
	InterIndTimeDiff->Write();
	NColperWT->Write();
	NColperTW->Write();
	ColWidthZSWT->Write();
	ColWidthZSTW->Write();
	ColIntZS->Write();
	QTotColZSWT->Write();
	QTotColZSTW->Write();
	tgint1->Write();
	tgint2->Write();
	tgint3->Write();
	outroot->Close();
	inroot->Close();
	
	sprintf(hname,"cp HitAnalysis_%d.root %s/Histos/HitAnalysis_%d.root;wait;",RunNo,AnalysisFilePath,RunNo);system(hname);
}



void HitAnalysis2()
{
	sprintf(hname,"%s/Run_%d.root",NTupleFilePath,RunNo);
	TFile* inroot=new TFile(hname);
	TTree* T =  (TTree*) inroot->Get("T");
	T->SetBranchAddress("TPCWF",&ed.TPCWF);
	T->SetBranchAddress("PMTWF",&ed.PMTWF);
	
	ReadTPCBaselines();
	
	sprintf(hname,"HitAnalysis2_%d.root",RunNo);
	TFile* outroot=new TFile(hname,"recreate");
	
	TGraph* Col_vs_t=new TGraph();
	Col_vs_t->SetName("Col_vs_t");Col_vs_t->SetTitle("Col_vs_t");
	Col_vs_t->SetMarkerStyle(25);Col_vs_t->SetMarkerColor(1);
	
	TGraph* Ind_vs_t=new TGraph();
	Ind_vs_t->SetName("Ind_vs_t");Ind_vs_t->SetTitle("Ind_vs_t");
	Ind_vs_t->SetMarkerStyle(26);Ind_vs_t->SetMarkerColor(2);
	
	TH2F* Col_vs_Ind=new TH2F("Col_vs_Ind","Col_vs_Ind",128,-0.5,127.5,128,-0.5,127.5);
	
	int ic=0;int gch=0;
	float ss=0;
	outroot->cd();
	for(int I=0;I<T->GetEntries();I++)
	{
		T->GetEntry(I);
		Col_vs_t->Set(0);Ind_vs_t->Set(0);
		for(int i1=0;i1<256;i1++)
		{
			ic=i1/128;gch=i1-((i1/128)*128);
			for(int i2=750;i2<2500;i2++)
			{
				ss=((float)ed.TPCWF->at(i1)[i2])-TPCBaselines[ic][gch][0];
				if(i1<128)
				{
					if(ss>(TPCBaselines[ic][gch][1]*ThSigma[ic]))
					{
						Col_vs_t->SetPoint(Col_vs_t->GetN(),i2,gch);
					}
				}
				else
				{
					ss*=-1;
					if(ss>(TPCBaselines[ic][gch][1]*ThSigma[ic]))
					{
						Ind_vs_t->SetPoint(Ind_vs_t->GetN(),i2,gch);
					}
				}
			}
		}
// 		for(int i2=750;i2<2500;i2++)
// 		{
// 			for(int i1=0;i1<256;i1++)
// 			{
// 				ic=i1/128;gch=i1-((i1/128)*128);
// 				ss=((float)ed.TPCWF->at(i1)[i2])-TPCBaselines[ic][gch][0];
				
				
				
				
		
		
		sprintf(hname,"Col_vs_t_%d",I);
		Col_vs_t->SetName(hname);Col_vs_t->SetTitle(hname);
		Col_vs_t->Write();
		sprintf(hname,"Ind_vs_t_%d",I);
		Ind_vs_t->SetName(hname);Ind_vs_t->SetTitle(hname);
		Ind_vs_t->Write();
		if(I==10) break;
	}
	outroot->Close();
}












void WriteWF()
{
	sprintf(hname,"%s/Run_%d.root",NTupleFilePath,RunNo);
	TFile* inroot=new TFile(hname);
	TTree* T =  (TTree*) inroot->Get("T");
	T->SetBranchAddress("TPCWF",&ed.TPCWF);
	T->SetBranchAddress("PMTWF",&ed.PMTWF);
	
	TH1F* Amp[2];
	Amp[0]=new TH1F("Amp_Col","Amp_Col",1000,-500,500);
	Amp[1]=new TH1F("Amp_Ind","Amp_Ind",1000,-500,500);
	
	ReadTPCBaselines();
	ReadPMTBaselines();

	TFile* outroot;
	
	sprintf(hname,"Waveforms_%d.root",RunNo);
	outroot=new TFile(hname,"recreate");
	outroot->mkdir("Ch0");
	outroot->mkdir("Ch1");
	outroot->mkdir("Ch2");
	
	TH1F* hb;
	TH1F* hh[2][128];
	for(int i1=0;i1<2;i1++)
	{
		for(int i2=0;i2<128;i2++)
		{
			sprintf(hname,"TPCWF_%s_%d",wt[i1].c_str(),i2);
			hh[i1][i2]=new TH1F(hname,hname,4096,-0.5,4095.5);
		}
	}
	sprintf(hname,"Ev2D");
	TH2F* h2=new TH2F(hname,hname,4096,-0.5,4095.5,256,-0.5,255.5);
	
	TH1F* hhPMT[3];
	for(int i1=0;i1<3;i1++)
	{
		sprintf(hname,"Int_%d",i1);
		hhPMT[i1]=new TH1F(hname,hname,1200,-20000,4000);
	}
	
	int ic=0;int gch=0;
	float ss=0.;
	int NSat[3]={0};
	bool isSaturated=false;
	float qint=0.;
	float sg=0.;
	
	for(int I=0;I<T->GetEntries();I++)
	{
		T->GetEntry(I);
		
		for(int i1=0;i1<256;i1++)
		{
			ic=i1/128;gch=i1-((i1/128)*128);
			for(int i2=0;i2<4096;i2++)
			{
// 				ss=((float)ed.TPCWF->at(i1)[i2])-TPCBaselines[ic][gch][1];
				ss=((float)ed.TPCWF->at(i1)[i2])-TPCBaselines[ic][gch][0];
				hh[ic][gch]->SetBinContent(hh[ic][gch]->FindBin(i2),ss);
				h2->Fill(i2,i1,ss);
				Amp[ic]->Fill(ss);
			}
		}
		outroot->cd();
		
		sprintf(hname2,"Ev2D_%d",I);
		h2->SetName(hname2);h2->SetTitle(hname2);
		h2->SetMinimum(-50);h2->SetMaximum(50);
		h2->Write();
		h2->SetName("Ev2D");h2->SetTitle("Ev2D");
		h2->Reset();
		
// 		for(int i1=0;i1<256;i1++)
// 		{
// 			for(int i2=0;i2<4096;i2++)
// 			{
// 				sprintf(hname2,"%s",hh[i1][i2]->GetName());
// 				sprintf(hname,"Ev_%d_%s",I,hh[i1][i2]->GetName());
// 				hh[i1][i2]->SetName(hname);hh[i1][i2]->SetTitle(hname);
// 				hh[i1][i2]->Write();
// 				hh[i1][i2]->SetName(hname2);hh[i1][i2]->SetTitle(hname2);
// 			}
// 		}
		
		for(int i1=0;i1<2;i1++)
		{
			for(int i2=0;i2<128;i2++)
			{
				hh[i1][i2]->Reset();
			}
		}
		
		for(int i1=0;i1<3;i1++)
		{
			sprintf(hname,"WF_PMT_%d_Ev_%d",i1,I);
			hb=new TH1F(hname,hname,15000,-0.5,14999.5);
			isSaturated=false;
			qint=0.;
			for(int i2=0;i2<ed.PMTWF->at(i1).size();i2++)
			{
				sg=(-1.)*(((float)ed.PMTWF->at(i1)[i2])-PMTBaselines[i1][0]);
				hb->SetBinContent(i2+1,sg);
				if(sg>ThPMT[i1]) qint+=sg;
				if(ed.PMTWF->at(i1)[i2]==0){isSaturated=true;}
			}
			if(isSaturated){NSat[i1]++;}
			outroot->cd();
			
			sprintf(hname2,"Int=%5.3f",qint);
			hb->SetTitle(hname2);
			if(isSaturated)
			{
				sprintf(hname,"%s Saturated",hname2);
				hb->SetTitle(hname);
			}
			sprintf(hname,"Ch%d",i1);
			outroot->cd(hname);
			
			hb->Write();
			hb->Reset();
			hhPMT[i1]->Fill(qint);
		}
		if(I%100==0) cout<<"Ev : "<<I<<" / "<<T->GetEntries()<<endl;
// 		if(I==100) break;
	}
	
	outroot->cd();
	Amp[0]->Write();
	Amp[1]->Write();
	
	for(int i1=0;i1<3;i1++)
	{
		hhPMT[i1]->Write();
	}
	outroot->Close();
	file.close();
	
	cout<<"Run "<<RunNo<<" saturations: "<<NSat[0]<<" "<<NSat[1]<<" "<<NSat[2]<<" out of "<<T->GetEntries()<<endl; 
	
	sprintf(hname,"cp Waveforms_%d.root %s/Histos/Waveforms_%d.root",RunNo,AnalysisFilePath,RunNo);system(hname);
}

void TPCCalibration()
{
	sprintf(hname,"TPCCalibration.root");
	TFile* outroot=new TFile(hname,"recreate");
	
	ReadTPCBaselines();
	
	int RunNos[2]={7384,7383};//coll ind
	
	TH1F* NegPulseArea[2][128];
	TH1F* PosPulseArea[2][128];
	for(int i1=0;i1<2;i1++)
	{
		for(int i2=0;i2<128;i2++)
		{
			sprintf(hname,"NegPulseArea_%s_%d",wt[i1].c_str(),i2);
			NegPulseArea[i1][i2]=new TH1F(hname,hname,200,0,2000);
			sprintf(hname,"PosPulseArea_%s_%d",wt[i1].c_str(),i2);
			PosPulseArea[i1][i2]=new TH1F(hname,hname,200,0,2000);
		}
	}
	
	TH1F* PeakAmp[2];
	TH1D* AllAmp[2];
	
	PeakAmp[0]=new TH1F("ColPeakAmp","ColPeakAmp",400,0,200);
	PeakAmp[1]=new TH1F("IndPeakAmp","IndPeakAmp",400,0,200);
	
	AllAmp[0]=new TH1D("ColAllAmp","ColAllAmp",800,-200,200);
	AllAmp[1]=new TH1D("IndAllAmp","IndAllAmp",800,-200,200);
	
	TH1F* hh[2][128];TH1F* hh2[2][128];TH1F* hh3[2][128];TH1F* hh4[2][128];
	float QtotCol=0;
	int startbin=0;
	float ss=0;
	int ic=0;int gch=0;
	
	for(int i1=0;i1<2;i1++)
	{
		for(int i2=0;i2<128;i2++)
		{
			sprintf(hname,"WFPos_%s_%d",wt[i1].c_str(),i2);
			hh[i1][i2]=new TH1F(hname,hname,4096,-0.5,4095.5);
			sprintf(hname,"WFPosZS_%s_%d",wt[i1].c_str(),i2);
			hh2[i1][i2]=new TH1F(hname,hname,4096,-0.5,4095.5);
			sprintf(hname,"WFNegZS_%s_%d",wt[i1].c_str(),i2);
			hh3[i1][i2]=new TH1F(hname,hname,4096,-0.5,4095.5);
			sprintf(hname,"WFNeg_%s_%d",wt[i1].c_str(),i2);
			hh4[i1][i2]=new TH1F(hname,hname,4096,-0.5,4095.5);
		}
	}
	
	int cbin=0;
	int peakbin=0;
	TFile* inroot;
	TTree* T;
	for(int J=0;J<2;J++)
	{
		RunNo=RunNos[J];
		sprintf(hname,"%s/Run_%d.root",NTupleFilePath,RunNo);
		inroot=new TFile(hname);
		T =  (TTree*) inroot->Get("T");
		T->SetBranchAddress("TPCWF",&ed.TPCWF);
		T->SetBranchAddress("PMTWF",&ed.PMTWF);
		
		for(int I=0;I<T->GetEntries();I++)
		{
			T->GetEntry(I);
			if(I%100==0) cout<<I<<endl;
			for(int i1=0;i1<256;i1++)
			{
				ic=i1/128;gch=i1-((i1/128)*128);
				for(int i2=0;i2<4096;i2++)
// 				for(int i2=TimeRange[0];i2<TimeRange[1];i2++)
				{
					ss=((float)ed.TPCWF->at(i1)[i2])-TPCBaselines[ic][gch][0];
					
					AllAmp[J]->Fill(ss);
					
					hh[ic][gch]->SetBinContent(hh[ic][gch]->FindBin(i2),ss);
// 					if(ss>(TPCBaselines[ic][gch][1]*ThSigma[ic]))
					if(ss>20)
					{
						hh2[ic][gch]->SetBinContent(hh2[ic][gch]->FindBin(i2),ss);
					}
					ss*=-1;
					hh4[ic][gch]->SetBinContent(hh4[ic][gch]->FindBin(i2),ss);
					if(ss>20)
					{
						hh3[ic][gch]->SetBinContent(hh3[ic][gch]->FindBin(i2),ss);
					}
				}
			}
// 			for(int i1=0;i1<2;i1++)
			{
				int i1=J;
				for(int i2=0;i2<128;i2++)
				{
					startbin=FindBinAbove(hh2[i1][i2],0,1);
					while(startbin>=0)
					{
// 						cout<<i1<<" "<<i2<<" pos "<<startbin<<endl;
						for(int ik1=startbin;ik1<hh2[i1][i2]->GetNbinsX();ik1++)
						{
							if(hh2[i1][i2]->GetBinContent(ik1+1)<hh2[i1][i2]->GetBinContent(ik1))
							{
								PeakAmp[i1]->Fill(hh2[i1][i2]->GetBinContent(ik1));
								cbin=ik1;
								break;
							}
						}
						peakbin=cbin;
						for(int ik1=cbin;ik1<=hh2[i1][i2]->GetNbinsX();ik1++)
						{
							if(hh2[i1][i2]->GetBinContent(ik1)==0)
							{
								cbin=ik1;
								break;
							}
						}
						if(cbin>100 && cbin<4000)
						{
							QtotCol=0.;
							for(int ik1=cbin;ik1>1;ik1--)
							{
								if(hh[i1][i2]->GetBinContent(ik1)>0)
								{
									QtotCol+=hh[i1][i2]->GetBinContent(ik1);
								}
								else break;
							}
							for(int ik1=cbin+1;ik1<=hh[i1][i2]->GetNbinsX();ik1++)
							{
								if(hh[i1][i2]->GetBinContent(ik1)>0)
								{
									QtotCol+=hh[i1][i2]->GetBinContent(ik1);
								}
								else break;
							}
							PosPulseArea[i1][i2]->Fill(QtotCol);
						}
						else break;
						startbin=FindBinAbove(hh2[i1][i2],0,cbin+1);
						if(startbin>4000) break;
					}
					
					startbin=FindBinAbove(hh3[i1][i2],0,1);
					while(startbin>=0)
					{
// 						cout<<i1<<" "<<i2<<" pos "<<startbin<<endl;
						for(int ik1=startbin;ik1<hh3[i1][i2]->GetNbinsX();ik1++)
						{
							if(hh3[i1][i2]->GetBinContent(ik1+1)<hh3[i1][i2]->GetBinContent(ik1))
							{
								PeakAmp[i1]->Fill(hh3[i1][i2]->GetBinContent(ik1));
								cbin=ik1;
								break;
							}
						}
						peakbin=cbin;
						for(int ik1=cbin;ik1<=hh3[i1][i2]->GetNbinsX();ik1++)
						{
							if(hh3[i1][i2]->GetBinContent(ik1)==0)
							{
								cbin=ik1;
								break;
							}
						}
						if(cbin>100 && cbin<4000)
						{
							QtotCol=0.;
							for(int ik1=cbin;ik1>1;ik1--)
							{
								if(hh4[i1][i2]->GetBinContent(ik1)>0)
								{
									QtotCol+=hh4[i1][i2]->GetBinContent(ik1);
								}
								else break;
							}
							for(int ik1=cbin+1;ik1<=hh4[i1][i2]->GetNbinsX();ik1++)
							{
								if(hh4[i1][i2]->GetBinContent(ik1)>0)
								{
									QtotCol+=hh4[i1][i2]->GetBinContent(ik1);
								}
								else break;
							}
							NegPulseArea[i1][i2]->Fill(QtotCol);
						}
						else break;
						startbin=FindBinAbove(hh3[i1][i2],0,cbin+1);
						if(startbin>4000) break;
					}
				}
			}
			for(int i1=0;i1<2;i1++)
			{
				for(int i2=0;i2<128;i2++)
				{
					hh[i1][i2]->Reset();
					hh2[i1][i2]->Reset();
					hh3[i1][i2]->Reset();
					hh4[i1][i2]->Reset();
				}
			}
		}
		inroot->Close();
	}
	TF1* tf1=new TF1("landau","landau",0.,2000.);
	TF1* tf2=new TF1("g+l","gaus(0)+landau(3)",0.,2000.);
	TF1* tf3=new TF1("g+g","gaus(0)+gaus(3)",0.,2000.);
	TGraphErrors* tg[2][2];
	tg[0][0]=new TGraphErrors();
	tg[0][0]->SetName("NegInt_Col");tg[0][0]->SetTitle("NegInt_Col");
	tg[0][0]->SetMarkerStyle(20);tg[0][0]->SetMarkerColor(1);
	tg[0][1]=new TGraphErrors();
	tg[0][1]->SetName("PosInt_Col");tg[0][1]->SetTitle("PosInt_Col");
	tg[0][1]->SetMarkerStyle(21);tg[0][1]->SetMarkerColor(2);
	tg[1][0]=new TGraphErrors();
	tg[1][0]->SetName("NegInt_Ind");tg[1][0]->SetTitle("NegInt_Ind");
	tg[1][0]->SetMarkerStyle(22);tg[1][0]->SetMarkerColor(3);
	tg[1][1]=new TGraphErrors();
	tg[1][1]->SetName("PosInt_Ind");tg[1][1]->SetTitle("PosInt_Ind");
	tg[1][1]->SetMarkerStyle(23);tg[1][1]->SetMarkerColor(4);
	
	float maxfit=0;
	outroot->cd();
	for(int i1=0;i1<2;i1++){AllAmp[i1]->Write();}
	for(int i1=0;i1<2;i1++){PeakAmp[i1]->Write();}
	for(int i1=0;i1<2;i1++)
	{
		for(int i2=0;i2<128;i2++)
		{
			for(int ik1=1;ik1<NegPulseArea[i1][i2]->GetNbinsX();ik1++)
			{
				if((NegPulseArea[i1][i2]->Integral(1,ik1)/NegPulseArea[i1][i2]->Integral())>0.95){maxfit=NegPulseArea[i1][i2]->GetBinCenter(ik1);break;}
			}
// 			NegPulseArea[i1][i2]->Fit(tf1,"q","q",0.,maxfit);
// 			tg[i1][0]->SetPoint(tg[i1][0]->GetN(),i2,tf1->GetParameter(1));
// 			tg[i1][0]->SetPointError(tg[i1][0]->GetN()-1,0,tf1->GetParError(1));
			
// 			tf2->SetParameter(0,NegPulseArea[i1][i2]->GetBinContent(NegPulseArea[i1][i2]->GetMaximumBin()));
// 			tf2->SetParameter(1,NegPulseArea[i1][i2]->GetBinCenter(NegPulseArea[i1][i2]->GetMaximumBin()));
// 			tf2->SetParameter(3,NegPulseArea[i1][i2]->GetBinContent(NegPulseArea[i1][i2]->GetMaximumBin())/10);
// 			tf2->SetParameter(4,NegPulseArea[i1][i2]->GetBinCenter(NegPulseArea[i1][i2]->GetMaximumBin()));
// 			NegPulseArea[i1][i2]->Fit(tf2,"q","q",0.,maxfit);
// 			tg[i1][0]->SetPoint(tg[i1][0]->GetN(),i2,tf2->GetParameter(1));
// 			tg[i1][0]->SetPointError(tg[i1][0]->GetN()-1,0,tf2->GetParError(1));
			
			NegPulseArea[i1][i2]->Fit(tf1,"q","q",0.,maxfit);
			tf3->SetParameter(0,tf1->GetParameter(0));
			tf3->SetParameter(1,tf1->GetParameter(1));
			tf3->SetParameter(3,tf1->GetParameter(0)/10);
			tf3->SetParameter(4,tf1->GetParameter(1)+5);
			
// 			tf3->SetParameter(0,NegPulseArea[i1][i2]->GetBinContent(NegPulseArea[i1][i2]->GetMaximumBin()));
// 			tf3->SetParameter(1,NegPulseArea[i1][i2]->GetBinCenter(NegPulseArea[i1][i2]->GetMaximumBin()));
// 			tf3->SetParameter(3,NegPulseArea[i1][i2]->GetBinContent(NegPulseArea[i1][i2]->GetMaximumBin())/10);
// 			tf3->SetParameter(4,NegPulseArea[i1][i2]->GetBinCenter(NegPulseArea[i1][i2]->GetMaximumBin()));
			NegPulseArea[i1][i2]->Fit(tf3,"q","q",0.,maxfit);
			tg[i1][0]->SetPoint(tg[i1][0]->GetN(),i2,tf3->GetParameter(1));
			tg[i1][0]->SetPointError(tg[i1][0]->GetN()-1,0,tf3->GetParError(1));
			
			NegPulseArea[i1][i2]->Write();
			
			for(int ik1=1;ik1<PosPulseArea[i1][i2]->GetNbinsX();ik1++)
			{
				if((PosPulseArea[i1][i2]->Integral(1,ik1)/PosPulseArea[i1][i2]->Integral())>0.95){maxfit=PosPulseArea[i1][i2]->GetBinCenter(ik1);break;}
			}
			PosPulseArea[i1][i2]->Fit(tf1,"q","q",0.,maxfit);
			tg[i1][1]->SetPoint(tg[i1][1]->GetN(),i2,tf1->GetParameter(1));
			tg[i1][1]->SetPointError(tg[i1][1]->GetN()-1,0,tf1->GetParError(1));
			PosPulseArea[i1][i2]->Write();
		}
	}
	tg[0][0]->Write();
	tg[0][1]->Write();
	tg[1][0]->Write();
	tg[1][1]->Write();
	
	outroot->Close();
}

void TPCStats()
{
	sprintf(hname,"TPCStats.root");
	TFile* outroot=new TFile(hname,"recreate");
	
// 	int RunNos[32]={7074,7076,7078,7082,7084,7092,7094,7096,7104,7106,7107,7110,7121,7252,7254,7256,7259,7268,7271,7273,7275,7308,7312,7314,7316,7318,7320,7322,7325,7327,7385,7388};
	int RunNos[31]={7074,7076,7078,7082,7084,7092,7094,7096,7104,7106,7107,7110,7121,7252,7254,7256,7259,7268,7271,7273,7308,7312,7314,7316,7318,7320,7322,7325,7327,7385,7388};
	
	int NEvt=0;
	int NPMTSat[3]={0};
	int NSat=0;
	int NEvtType[4]={0};
	
	TFile* inroot;
	TTree* TH;
	TTree* TT;
	bool PMTsat=false;
	for(int J=0;J<31;J++)
	{
		RunNo=RunNos[J];
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
		
		NEvt=TH->GetEntries();
		NSat=0;
		NPMTSat[0]=0;NPMTSat[1]=0;NPMTSat[2]=0;
		NEvtType[0]=0;NEvtType[1]=0;NEvtType[2]=0;NEvtType[3]=0;
		for(int I=0;I<TH->GetEntries();I++)
		{
			TH->GetEntry(I);
			
			PMTsat=false;
			for(int i1=0;i1<3;i1++)
			{
				if(hd.PMTIntegral->at(i1)<0){NPMTSat[i1]++;PMTsat=true;}
			}
			if(!PMTsat) NEvtType[hd.EventType]++;
			else{NSat++;}
		}
		cout<<"Run: "<<RunNo<<" "<<NEvt<<" Events "<<NSat<<" with PMT saturation. Event Types: "<<NEvtType[0]<<" "<<NEvtType[1]<<" "<<NEvtType[2]<<" "<<NEvtType[3]<<endl;
		inroot->Close();
	}
}

void PMTCalibrationsAll()
{
// 	int Run1List[85]={50,7073,51,7074,7075,52,7076,57,7077,7078,7079,59,7080,7081,60,7082,61,7083,7084,7085,7086,7087,7088,68,7089,7090,70,7091,71,74,7092,7093,75,76,7094,7095,77,7096,7098,78,7099,7100,7101,79,7102,7103,7104,7105,80,7106,7107,7108,81,7109,7110,7120,82,7121,7122,83,7123,7124,84,7125,7126,7127,85,7128,7129,7130,86,7131,7132,7133,7134,7135,7136,7137,7138,7139,7140,7141,7142,7143,7144};
// 	int Run2List[73]={88,91,92,93,94,7240,7241,95,96,7242,7243,7244,7245,7246,97,98,7247,7248,7249,7250,7251,99,7252,7253,100,7254,7255,101,7256,7257,7258,102,7259,7260,103,7261,7262,7263,7264,7265,104,7266,7267,105,7268,7269,7270,106,7271,7272,107,7273,7274,108,7275,7276,109,7277,7278,110,7279,7280,111,7281,7282,112,7283,7284,113,7285,114,7286,7287};
// 	int Run3List[123]={115,7307,7308,116,7309,7310,117,7311,7312,118,7313,7314,119,7315,7316,120,7317,7318,121,7319,7320,122,7321,7322,7323,7324,123,7325,7326,124,7327,7328,125,7329,7330,126,7331,7332,127,7333,7334,7335,128,7336,7337,129,7338,7339,7340,7341,7342,7343,7344,7345,7346,7347,7348,7349,7350,7351,7352,130,7353,7354,7355,7356,7357,7358,7359,7360,7361,7362,7363,7364,7365,7366,7367,7368,7369,7370,131,132,7371,7372,7373,7374,7375,7376,7377,7378,133,7379,7380,7381,7382,7383,7384,134,7385,7386,135,7387,7388,7389,7390,7391,136,7392,7393,137,7394,7395,7396,138,7397,7398,7399,7400,7401,139,7402,7403,140};
	
	int NR1L=29;
	int NR2L=20;
	int NR3L=24;
	
	int Run1List[29]={50,51,7074,7076,57,59,60,61,7084,70,71,74,75,76,78,80,7106,7110,82,83,84,7126,85,7131,7132,7136,7141,7142,7143};//7087,7081,7082,7096,7099,7100,7123,7128,7129,7104,7107,7121,7140,7135,7133,7134,7144,7137,7138,7139,7092,77,86,79,81,68,7078,7089,7094,7125  
	int Run2List[20]={7241,7245,99,100,7254,101,102,7259,103,7264,104,105,7268,106,7271,107,7273,7281,7283,7285};//7252,7279,7287,92,94,96,97,93,95,98,7257,7263,108,7275,109,110,112,113,114,7266,7277 weird ,7256
	int Run3List[24]={115,116,117,118,119,120,7318,121,7320,122,7322,7323,123,7325,124,7327,125,126,127,128,129,7369,7394,7397};//,130,132,133,134,135,136,137,138,139 rate ,7371,7374,7375 something wrong with LXe ,7312,7314,7316 /// ,7331,7334,7336,7385,7389,7392
	
	TGraphErrors* tge[3][3][3];//pmt runperiod method
	for(int i1=0;i1<3;i1++)
	{
		if(i1==1) continue;
		for(int i2=0;i2<3;i2++)
		{
			for(int i3=0;i3<2;i3++)
			{
				sprintf(hname,"PMT_%s_Run%d_M%d",PMTNames[i1].c_str(),i2+1,i3+1);
				tge[i1][i2][i3]=new TGraphErrors();
				tge[i1][i2][i3]->SetName(hname);tge[i1][i2][i3]->SetTitle(hname);
				tge[i1][i2][i3]->SetMarkerStyle(20+i3);
				tge[i1][i2][i3]->SetMarkerColor(1+i3*3+i1);
				tge[i1][i2][i3]->GetXaxis()->SetTitle("Run ID");tge[i1][i2][i3]->GetXaxis()->CenterTitle();
				tge[i1][i2][i3]->GetYaxis()->SetTitle("Single Photoelectron Charge (arbitrary units)");tge[i1][i2][i3]->GetYaxis()->CenterTitle();
			}
			sprintf(hname,"PMT_%s_Run%d_All",PMTNames[i1].c_str(),i2+1);
			tge[i1][i2][2]=new TGraphErrors();
			tge[i1][i2][2]->SetName(hname);tge[i1][i2][2]->SetTitle(hname);
			tge[i1][i2][2]->SetMarkerStyle(20);
			tge[i1][i2][2]->SetMarkerColor(1);
			tge[i1][i2][2]->GetXaxis()->SetTitle("Run ID");tge[i1][i2][2]->GetXaxis()->CenterTitle();
			tge[i1][i2][2]->GetYaxis()->SetTitle("Single Photoelectron Charge (arbitrary units)");tge[i1][i2][2]->GetYaxis()->CenterTitle();
		}
	}
	
	TFile* outroot=new TFile("PMTCalibrationsAll.root","recreate");
	
	TH1F* hh;TF1* tf1;
	TF1* tf2=new TF1("GG","gaus(0)+gaus(3)",0.,200.);
	TF1* tf3=new TF1("gaus","gaus",0.,100.);
	TFile* inroot;
	float fitulim[3]={45,0,35};
	int maxbin1=0;int maxbin2=0;int minbin1=0;
	int minbin=0;int maxbin=0;
	float xmin=0;float xmax=0;
	float fitchi21=0;float fitchi22=0;
	for(int i1=0;i1<NR1L;i1++)
	{
// 		if(Run1List[i1]==52) continue;
		cout<<i1<<" "<<Run1List[i1];
		if(Run1List[i1]<7000)
		{
			sprintf(hname,"%s/Histos/PMTCalibration_%d.root",AnalysisFilePath,Run1List[i1]);
			inroot=new TFile(hname);
			for(int i2=0;i2<3;i2++)
			{
				if(i2==1) continue;
				sprintf(hname,"SignalWindow2_%d",i2);
				inroot->GetObject(hname,hh);
				
// 				tf1=(TF1*)hh->GetFunction("GGG");
// 				tge[i2][0][0]->SetPoint(tge[i2][0][0]->GetN(),i1,tf1->GetParameter(4));
// 				tge[i2][0][0]->SetPointError(tge[i2][0][0]->GetN()-1,0,tf1->GetParError(4));
// 				cout<<" "<<i2<<" "<<tf1->GetParameter(4)<<" "<<tf1->GetParError(4)<<endl;
				
				maxbin1=hh->GetMaximumBin();
				hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(maxbin1),30);
				minbin1=hh->GetMinimumBin();
				hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(minbin1),60);
				maxbin2=hh->GetMaximumBin();
				hh->GetXaxis()->SetRangeUser(-50,500);
				
				for(int ik1=maxbin1;ik1>=1;ik1--)
				{
					if((hh->GetBinContent(maxbin1)/hh->GetBinContent(ik1))>100){minbin=ik1;break;}
				}
				for(int ik1=maxbin2;ik1<=hh->GetNbinsX();ik1++)
				{
					if((hh->GetBinContent(maxbin2)/hh->GetBinContent(ik1))>5){maxbin=ik1;break;}
				}
				tf2->SetParameter(0,hh->GetBinContent(maxbin1));
				tf2->SetParameter(1,hh->GetBinCenter(maxbin1));
				tf2->SetParameter(2,5);
				tf2->SetParameter(3,hh->GetBinContent(maxbin2));
				tf2->SetParameter(4,hh->GetBinCenter(maxbin2));
				tf2->SetParameter(5,20);
				tf2->SetParLimits(5,0,50);
				hh->Fit(tf2,"q","q",hh->GetBinCenter(minbin),hh->GetBinCenter(maxbin));
				tge[i2][0][0]->SetPoint(tge[i2][0][0]->GetN(),i1,tf2->GetParameter(4)-tf2->GetParameter(1));
				tge[i2][0][0]->SetPointError(tge[i2][0][0]->GetN()-1,0,sqrt(pow(tf2->GetParError(4),2)+pow(tf2->GetParError(1),2)));
				tge[i2][0][2]->SetPoint(tge[i2][0][2]->GetN(),i1,tf2->GetParameter(4)-tf2->GetParameter(1));
				tge[i2][0][2]->SetPointError(tge[i2][0][2]->GetN()-1,0,sqrt(pow(tf2->GetParError(4),2)+pow(tf2->GetParError(1),2)));
				cout<<" "<<i2<<" "<<(tf2->GetParameter(4)-tf2->GetParameter(1))<<" "<<sqrt(pow(tf2->GetParError(4),2)+pow(tf2->GetParError(1),2))<<endl;
				
				
				
				sprintf(hname,"%d_%d",Run1List[i1],i2);
				hh->SetName(hname);hh->SetTitle(hname);
				outroot->cd();hh->Write();
			}
			inroot->Close();
		}
		else
		{
			sprintf(hname,"%s/Histos/PMTCalibration3_%d.root",AnalysisFilePath,Run1List[i1]);
			inroot=new TFile(hname);
			for(int i2=0;i2<3;i2++)
			{
				if(i2==1) continue;
				if(Run1List[i1]<7095 && i2==2) sprintf(hname,"Integral_w_%d_20",i2);
				else sprintf(hname,"Integral_w_%d_6",i2);
// 				sprintf(hname,"Integral_w_%d_6",i2);
				inroot->GetObject(hname,hh);
				
// 				tf2->SetParameter(0,500);
// 				tf2->SetParameter(1,8);
// 				tf2->SetParameter(2,2);
// 				tf2->SetParameter(3,1000);
// 				tf2->SetParameter(4,30);
// 				tf2->SetParameter(5,10);
// 				hh->Fit(tf2,"q","q",0.,fitulim[i2]);
				
				if(i2==0)
				{
					hh->GetXaxis()->SetRangeUser(0,20);
					maxbin1=hh->GetMaximumBin();
					hh->GetXaxis()->SetRangeUser(20,80);
					maxbin2=hh->GetMaximumBin();
					hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(maxbin1),hh->GetBinCenter(maxbin2));
					minbin1=hh->GetMinimumBin();
					hh->GetXaxis()->SetRangeUser(0,200);
					if(minbin1==maxbin1)
					{
						minbin1-=5;
					}
					tf3->SetParameter(0,hh->GetBinContent(maxbin2));
					tf3->SetParameter(1,hh->GetBinCenter(maxbin2));
					tf3->SetParameter(2,20);
					xmin=hh->GetBinCenter(minbin1);
					xmax=2*hh->GetBinCenter(maxbin2)-hh->GetBinCenter(minbin1);
					hh->Fit(tf3,"q","q",xmin,xmax);
					xmin=tf3->GetParameter(1)-1.5*tf3->GetParameter(2);
					xmax=tf3->GetParameter(1)+1.5*tf3->GetParameter(2);
					if(xmin<hh->GetBinCenter(minbin1)) xmin=hh->GetBinCenter(minbin1);
					hh->Fit(tf3,"q","q",xmin,xmax);
				}
				else
				{
// 					hh->GetXaxis()->SetRangeUser(0,20);
// 					maxbin1=hh->GetMaximumBin();
// 					hh->GetXaxis()->SetRangeUser(20,80);
// 					maxbin2=hh->GetMaximumBin();
// 					hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(maxbin1),hh->GetBinCenter(maxbin2));
// 					minbin1=hh->GetMinimumBin();
// 					hh->GetXaxis()->SetRangeUser(0,200);
// 					if(minbin1==maxbin1) minbin1-=5;
// 					if(minbin1==maxbin2) maxbin2+=5;
// 					tf3->SetParameter(0,hh->GetBinContent(maxbin2));
// 					tf3->SetParameter(1,hh->GetBinCenter(maxbin2));
// 					tf3->SetParameter(2,20);
// 					xmin=hh->GetBinCenter(minbin1);
// 					xmax=2*hh->GetBinCenter(maxbin2)-hh->GetBinCenter(minbin1);
// 					hh->Fit(tf3,"q","q",xmin,xmax);
// 					xmin=tf3->GetParameter(1)-1.5*tf3->GetParameter(2);
// 					xmax=tf3->GetParameter(1)+1.5*tf3->GetParameter(2);
// 					if(xmin<hh->GetBinCenter(minbin1)) xmin=hh->GetBinCenter(minbin1);
// 					hh->Fit(tf3,"q","q",xmin,xmax);
					
					
					hh->Fit(tf3,"q","q",0,100);
					xmin=tf3->GetParameter(1)-1.5*tf3->GetParameter(2);
					xmax=tf3->GetParameter(1)+1.5*tf3->GetParameter(2);
					hh->Fit(tf3,"q","q",xmin,xmax);
					
					if(Run1List[i1]<7095)
					{
						hh->GetXaxis()->SetRangeUser(0,20);
						maxbin1=hh->GetMaximumBin();
						hh->GetXaxis()->SetRangeUser(30,80);
						maxbin2=hh->GetMaximumBin();
						hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(maxbin1),hh->GetBinCenter(maxbin2));
						minbin1=hh->GetMinimumBin();
						hh->GetXaxis()->SetRangeUser(0,200);
						tf3->SetParameter(0,hh->GetBinContent(maxbin2));
						tf3->SetParameter(1,hh->GetBinCenter(maxbin2));
						tf3->SetParameter(2,20);
						xmin=hh->GetBinCenter(minbin1);
						xmax=2*hh->GetBinCenter(maxbin2)-hh->GetBinCenter(minbin1);
						hh->Fit(tf3,"q","q",xmin,xmax);
						xmin=tf3->GetParameter(1)-1.5*tf3->GetParameter(2);
						xmax=tf3->GetParameter(1)+1.5*tf3->GetParameter(2);
						if(xmin<hh->GetBinCenter(minbin1)) xmin=hh->GetBinCenter(minbin1);
						hh->Fit(tf3,"q","q",xmin,xmax);
					}
					else if(Run1List[i1]>7131)
					{
						hh->GetXaxis()->SetRangeUser(15,80);
						maxbin2=hh->GetMaximumBin();
						hh->GetXaxis()->SetRangeUser(0,200);
						tf3->SetParameter(0,hh->GetBinContent(maxbin2));
						tf3->SetParameter(1,hh->GetBinCenter(maxbin2));
						tf3->SetParameter(2,20);
						hh->Fit(tf3,"q","q",hh->GetBinCenter(maxbin2)-5,hh->GetBinCenter(maxbin2)+5);
						xmin=tf3->GetParameter(1)-1.5*tf3->GetParameter(2);
						xmax=tf3->GetParameter(1)+1.5*tf3->GetParameter(2);
						hh->Fit(tf3,"q","q",xmin,xmax);
					}
				}
// 				hh->GetXaxis()->SetRangeUser(0,20);
// 				maxbin1=hh->GetMaximumBin();
// 				hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(maxbin1),30);
// 				minbin1=hh->GetMinimumBin();
// 				hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(minbin1),80);
// 				maxbin2=hh->GetMaximumBin();
// 				hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(minbin1),hh->GetBinCenter(maxbin2));
// 				minbin1=hh->GetMinimumBin();
// 				hh->GetXaxis()->SetRangeUser(0,200);
// 				tf3->SetParameter(0,hh->GetBinContent(maxbin2));
// 				tf3->SetParameter(1,hh->GetBinCenter(maxbin2));
// 				tf3->SetParameter(2,20);
// 				
// 				for(int ik1=maxbin2;ik1>=1;ik1--)
// 				{
// 					if((hh->GetBinContent(ik1)/hh->GetBinContent(maxbin2))<0.7){minbin=ik1;break;}
// 				}
// 				
// 				for(int ik1=maxbin2;ik1<=hh->GetNbinsX();ik1++)
// 				{
// 					if((hh->GetBinContent(ik1)/hh->GetBinContent(maxbin2))<0.7){maxbin=ik1;break;}
// 				}
// 				
// 				if(minbin<minbin1) minbin=minbin1;
// 				
// 				hh->Fit(tf3,"q","q",hh->GetBinCenter(minbin),hh->GetBinCenter(maxbin));
// // 				fitchi21=tf3->GetChisquare()/tf3->GetNDF();
// 				xmin=tf3->GetParameter(1)-1.5*tf3->GetParameter(2);
// 				xmax=tf3->GetParameter(1)+1.5*tf3->GetParameter(2);
// 				if(xmin<hh->GetBinCenter(minbin1)) xmin=hh->GetBinCenter(minbin1);
// 				hh->Fit(tf3,"q","q",xmin,xmax);
// // 				fitchi22=tf3->GetChisquare()/tf3->GetNDF();
// // 				if(fitchi21<fitchi22) hh->Fit(tf3,"q","q",hh->GetBinCenter(minbin),hh->GetBinCenter(maxbin));
// // 				hh->Fit(tf3,"q","q",tf3->GetParameter(1)-1.5*tf3->GetParameter(2),tf3->GetParameter(1)+1.5*tf3->GetParameter(2));
// 				
				
// 				tf3->SetParameter(1,50);
// // 				tf3->SetParLimits(1,20,60);
// 				hh->Fit(tf3,"q","q",20,80);
// 				if(tf3->GetParameter(1)<0)
// 				{
// 					tf3->SetParameter(1,50);
// 					hh->Fit(tf3,"q","q",15,45);
// 				}
// 				hh->Fit(tf3,"q","q",tf3->GetParameter(1)-1.5*tf3->GetParameter(2),tf3->GetParameter(1)+1.5*tf3->GetParameter(2));
				tge[i2][0][1]->SetPoint(tge[i2][0][1]->GetN(),i1,tf3->GetParameter(1));
				tge[i2][0][1]->SetPointError(tge[i2][0][1]->GetN()-1,0,tf3->GetParError(1));
				tge[i2][0][2]->SetPoint(tge[i2][0][2]->GetN(),i1,tf3->GetParameter(1));
				tge[i2][0][2]->SetPointError(tge[i2][0][2]->GetN()-1,0,tf3->GetParError(1));
				cout<<" "<<i2<<" "<<tf3->GetParameter(1)<<" "<<tf3->GetParError(1)<<endl;
				
// 				tge[i2][0][1]->SetPoint(tge[i2][0][1]->GetN(),i1,tf2->GetParameter(4));
// 				tge[i2][0][1]->SetPointError(tge[i2][0][1]->GetN()-1,0,tf2->GetParError(4));
// 				cout<<" "<<i2<<" "<<tf2->GetParameter(4)<<" "<<tf2->GetParError(4)<<endl;
				sprintf(hname,"%d_%d",Run1List[i1],i2);
				hh->SetName(hname);hh->SetTitle(hname);
				outroot->cd();hh->Write();
			}
			inroot->Close();
		}
	}
	for(int i1=0;i1<NR2L;i1++)
	{
// 		if(Run2List[i1]==88 || Run2List[i1]==91) continue;
		cout<<i1<<" "<<Run2List[i1];
		if(Run2List[i1]<7000)
		{
			sprintf(hname,"%s/Histos/PMTCalibration_%d.root",AnalysisFilePath,Run2List[i1]);
			inroot=new TFile(hname);
			for(int i2=0;i2<3;i2++)
			{
				if(i2==1) continue;
				sprintf(hname,"SignalWindow2_%d",i2);
				inroot->GetObject(hname,hh);
				
// 				tf1=(TF1*)hh->GetFunction("GGG");
// 				tge[i2][0][0]->SetPoint(tge[i2][0][0]->GetN(),i1,tf1->GetParameter(4));
// 				tge[i2][0][0]->SetPointError(tge[i2][0][0]->GetN()-1,0,tf1->GetParError(4));
// 				cout<<" "<<i2<<" "<<tf1->GetParameter(4)<<" "<<tf1->GetParError(4)<<endl;
				
				maxbin1=hh->GetMaximumBin();
				hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(maxbin1),30);
				minbin1=hh->GetMinimumBin();
				hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(minbin1),60);
				maxbin2=hh->GetMaximumBin();
				hh->GetXaxis()->SetRangeUser(-50,500);
				
				for(int ik1=maxbin1;ik1>=1;ik1--)
				{
					if((hh->GetBinContent(maxbin1)/hh->GetBinContent(ik1))>100){minbin=ik1;break;}
				}
				for(int ik1=maxbin2;ik1<=hh->GetNbinsX();ik1++)
				{
					if((hh->GetBinContent(maxbin2)/hh->GetBinContent(ik1))>5){maxbin=ik1;break;}
				}
				tf2->SetParameter(0,hh->GetBinContent(maxbin1));
				tf2->SetParameter(1,hh->GetBinCenter(maxbin1));
				tf2->SetParameter(2,5);
				tf2->SetParameter(3,hh->GetBinContent(maxbin2));
				tf2->SetParameter(4,hh->GetBinCenter(maxbin2));
				tf2->SetParameter(5,20);
				tf2->SetParLimits(5,0,50);
				hh->Fit(tf2,"q","q",hh->GetBinCenter(minbin),hh->GetBinCenter(maxbin));
				tge[i2][1][0]->SetPoint(tge[i2][1][0]->GetN(),i1,tf2->GetParameter(4)-tf2->GetParameter(1));
				tge[i2][1][0]->SetPointError(tge[i2][1][0]->GetN()-1,0,sqrt(pow(tf2->GetParError(4),2)+pow(tf2->GetParError(1),2)));
				tge[i2][1][2]->SetPoint(tge[i2][1][2]->GetN(),i1,tf2->GetParameter(4)-tf2->GetParameter(1));
				tge[i2][1][2]->SetPointError(tge[i2][1][2]->GetN()-1,0,sqrt(pow(tf2->GetParError(4),2)+pow(tf2->GetParError(1),2)));
				cout<<" "<<i2<<" "<<(tf2->GetParameter(4)-tf2->GetParameter(1))<<" "<<sqrt(pow(tf2->GetParError(4),2)+pow(tf2->GetParError(1),2))<<endl;
				
				
				
				sprintf(hname,"%d_%d",Run2List[i1],i2);
				hh->SetName(hname);hh->SetTitle(hname);
				outroot->cd();hh->Write();
			}
			inroot->Close();
// 			for(int i2=0;i2<3;i2++)
// 			{
// 				if(i2==1) continue;
// 				sprintf(hname,"SignalWindow2_%d",i2);
// 				inroot->GetObject(hname,hh);
// 				tf1=(TF1*)hh->GetFunction("GGG");
// 				tge[i2][1][0]->SetPoint(tge[i2][1][0]->GetN(),i1,tf1->GetParameter(4));
// 				tge[i2][1][0]->SetPointError(tge[i2][1][0]->GetN()-1,0,tf1->GetParError(4));
// 				cout<<" "<<i2<<" "<<tf1->GetParameter(4)<<" "<<tf1->GetParError(4)<<endl;
// 				sprintf(hname,"%d_%d",Run2List[i1],i2);
// 				hh->SetName(hname);hh->SetTitle(hname);
// 				outroot->cd();hh->Write();
// 			}
// 			inroot->Close();
		}
		else
		{
			sprintf(hname,"%s/Histos/PMTCalibration3_%d.root",AnalysisFilePath,Run2List[i1]);
			inroot=new TFile(hname);
			for(int i2=0;i2<3;i2++)
			{
				if(i2==1) continue;
				sprintf(hname,"Integral_w_%d_6",i2);
				inroot->GetObject(hname,hh);
				
// 				tf2->SetParameter(0,500);
// 				tf2->SetParameter(1,8);
// 				tf2->SetParameter(2,2);
// 				tf2->SetParameter(3,1000);
// 				tf2->SetParameter(4,30);
// 				tf2->SetParameter(5,10);
// 				hh->Fit(tf2,"q","q",0.,fitulim[i2]);
				
				if(i2==0)
				{
					hh->GetXaxis()->SetRangeUser(0,20);
					maxbin1=hh->GetMaximumBin();
					hh->GetXaxis()->SetRangeUser(25,80);
					maxbin2=hh->GetMaximumBin();
					hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(maxbin1),hh->GetBinCenter(maxbin2));
					minbin1=hh->GetMinimumBin();
					hh->GetXaxis()->SetRangeUser(0,200);
					if(minbin1==maxbin1)
					{
						minbin1-=5;
					}
					tf3->SetParameter(0,hh->GetBinContent(maxbin2));
					tf3->SetParameter(1,hh->GetBinCenter(maxbin2));
					tf3->SetParameter(2,20);
					xmin=hh->GetBinCenter(minbin1);
					xmax=2*hh->GetBinCenter(maxbin2)-hh->GetBinCenter(minbin1);
					hh->Fit(tf3,"q","q",xmin,xmax);
					xmin=tf3->GetParameter(1)-1.5*tf3->GetParameter(2);
					xmax=tf3->GetParameter(1)+1.5*tf3->GetParameter(2);
					if(xmin<hh->GetBinCenter(minbin1)) xmin=hh->GetBinCenter(minbin1);
					hh->Fit(tf3,"q","q",xmin,xmax);
				}
				else
				{
					hh->Fit(tf3,"q","q",0,50);
					xmin=tf3->GetParameter(1)-1.5*tf3->GetParameter(2);
					xmax=tf3->GetParameter(1)+1.5*tf3->GetParameter(2);
					hh->Fit(tf3,"q","q",xmin,xmax);
// 					
					if(Run2List[i1]==7241 || Run2List[i1]==7271 || Run2List[i1]==7273)
					{
						hh->GetXaxis()->SetRangeUser(0,15);
						maxbin1=hh->GetMaximumBin();
						hh->GetXaxis()->SetRangeUser(15,50);
						maxbin2=hh->GetMaximumBin();
						hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(maxbin1),hh->GetBinCenter(maxbin2));
						minbin1=hh->GetMinimumBin();
						hh->GetXaxis()->SetRangeUser(0,200);
						tf3->SetParameter(0,hh->GetBinContent(maxbin2));
						tf3->SetParameter(1,hh->GetBinCenter(maxbin2));
						tf3->SetParameter(2,20);
						xmin=hh->GetBinCenter(minbin1);
						xmax=2*hh->GetBinCenter(maxbin2)-hh->GetBinCenter(minbin1);
						hh->Fit(tf3,"q","q",xmin,xmax);
						xmin=tf3->GetParameter(1)-1.5*tf3->GetParameter(2);
						xmax=tf3->GetParameter(1)+1.5*tf3->GetParameter(2);
						if(xmin<hh->GetBinCenter(minbin1)) xmin=hh->GetBinCenter(minbin1);
						hh->Fit(tf3,"q","q",xmin,xmax);
					}
				}
				tge[i2][1][1]->SetPoint(tge[i2][1][1]->GetN(),i1,tf3->GetParameter(1));
				tge[i2][1][1]->SetPointError(tge[i2][1][1]->GetN()-1,0,tf3->GetParError(1));
				tge[i2][1][2]->SetPoint(tge[i2][1][2]->GetN(),i1,tf3->GetParameter(1));
				tge[i2][1][2]->SetPointError(tge[i2][1][2]->GetN()-1,0,tf3->GetParError(1));
				cout<<" "<<i2<<" "<<tf3->GetParameter(1)<<" "<<tf3->GetParError(1)<<endl;
				
				sprintf(hname,"%d_%d",Run2List[i1],i2);
				hh->SetName(hname);hh->SetTitle(hname);
				outroot->cd();hh->Write();
			}
			inroot->Close();
		}
	}
	for(int i1=0;i1<NR3L;i1++)
	{
// 		if(Run3List[i1]==131) continue;
		cout<<i1<<" "<<Run3List[i1];
		if(Run3List[i1]<7000)
		{
			sprintf(hname,"%s/Histos/PMTCalibration_%d.root",AnalysisFilePath,Run3List[i1]);
			inroot=new TFile(hname);
			for(int i2=0;i2<3;i2++)
			{
				if(i2==1) continue;
				sprintf(hname,"SignalWindow2_%d",i2);
				inroot->GetObject(hname,hh);
				
				maxbin1=hh->GetMaximumBin();
				hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(maxbin1),30);
				minbin1=hh->GetMinimumBin();
				hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(minbin1),60);
				maxbin2=hh->GetMaximumBin();
				hh->GetXaxis()->SetRangeUser(-50,500);
				
				for(int ik1=maxbin1;ik1>=1;ik1--)
				{
					if((hh->GetBinContent(maxbin1)/hh->GetBinContent(ik1))>100){minbin=ik1;break;}
				}
				for(int ik1=maxbin2;ik1<=hh->GetNbinsX();ik1++)
				{
					if((hh->GetBinContent(maxbin2)/hh->GetBinContent(ik1))>5){maxbin=ik1;break;}
				}
				tf2->SetParameter(0,hh->GetBinContent(maxbin1));
				tf2->SetParameter(1,hh->GetBinCenter(maxbin1));
				tf2->SetParameter(2,5);
				tf2->SetParameter(3,hh->GetBinContent(maxbin2));
				tf2->SetParameter(4,hh->GetBinCenter(maxbin2));
				tf2->SetParameter(5,20);
				tf2->SetParLimits(5,0,50);
				hh->Fit(tf2,"q","q",hh->GetBinCenter(minbin),hh->GetBinCenter(maxbin));
				tge[i2][2][0]->SetPoint(tge[i2][2][0]->GetN(),i1,tf2->GetParameter(4)-tf2->GetParameter(1));
				tge[i2][2][0]->SetPointError(tge[i2][2][0]->GetN()-1,0,sqrt(pow(tf2->GetParError(4),2)+pow(tf2->GetParError(1),2)));
				tge[i2][2][2]->SetPoint(tge[i2][2][2]->GetN(),i1,tf2->GetParameter(4)-tf2->GetParameter(1));
				tge[i2][2][2]->SetPointError(tge[i2][2][2]->GetN()-1,0,sqrt(pow(tf2->GetParError(4),2)+pow(tf2->GetParError(1),2)));
				cout<<" "<<i2<<" "<<(tf2->GetParameter(4)-tf2->GetParameter(1))<<" "<<sqrt(pow(tf2->GetParError(4),2)+pow(tf2->GetParError(1),2))<<endl;
				
				sprintf(hname,"%d_%d",Run3List[i1],i2);
				hh->SetName(hname);hh->SetTitle(hname);
				outroot->cd();hh->Write();
			}
			inroot->Close();
			
			
// 			for(int i2=0;i2<3;i2++)
// 			{
// 				if(i2==1) continue;
// 				sprintf(hname,"SignalWindow2_%d",i2);
// 				inroot->GetObject(hname,hh);
// 				tf1=(TF1*)hh->GetFunction("GGG");
// 				tge[i2][2][0]->SetPoint(tge[i2][2][0]->GetN(),i1,tf1->GetParameter(4));
// 				tge[i2][2][0]->SetPointError(tge[i2][2][0]->GetN()-1,0,tf1->GetParError(4));
// 				cout<<" "<<i2<<" "<<tf1->GetParameter(4)<<" "<<tf1->GetParError(4)<<endl;
// 				sprintf(hname,"%d_%d",Run3List[i1],i2);
// 				hh->SetName(hname);hh->SetTitle(hname);
// 				outroot->cd();hh->Write();
// 			}
// 			inroot->Close();
		}
		else
		{
			sprintf(hname,"%s/Histos/PMTCalibration3_%d.root",AnalysisFilePath,Run3List[i1]);
			inroot=new TFile(hname);
			for(int i2=0;i2<3;i2++)
			{
				if(i2==1) continue;
				sprintf(hname,"Integral_w_%d_6",i2);
				inroot->GetObject(hname,hh);
				
				if(i2==0)
				{
					hh->GetXaxis()->SetRangeUser(0,20);
					maxbin1=hh->GetMaximumBin();
					hh->GetXaxis()->SetRangeUser(25,80);
					maxbin2=hh->GetMaximumBin();
					hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(maxbin1),hh->GetBinCenter(maxbin2));
					minbin1=hh->GetMinimumBin();
					hh->GetXaxis()->SetRangeUser(0,200);
					if(minbin1==maxbin1)
					{
						minbin1-=5;
					}
					tf3->SetParameter(0,hh->GetBinContent(maxbin2));
					tf3->SetParameter(1,hh->GetBinCenter(maxbin2));
					tf3->SetParameter(2,20);
					xmin=hh->GetBinCenter(minbin1);
					xmax=2*hh->GetBinCenter(maxbin2)-hh->GetBinCenter(minbin1);
					hh->Fit(tf3,"q","q",xmin,xmax);
					xmin=tf3->GetParameter(1)-1.5*tf3->GetParameter(2);
					xmax=tf3->GetParameter(1)+1.5*tf3->GetParameter(2);
					if(xmin<hh->GetBinCenter(minbin1)) xmin=hh->GetBinCenter(minbin1);
					hh->Fit(tf3,"q","q",xmin,xmax);
				}
				else
				{
					hh->Fit(tf3,"q","q",0,50);
					xmin=tf3->GetParameter(1)-1.5*tf3->GetParameter(2);
					xmax=tf3->GetParameter(1)+1.5*tf3->GetParameter(2);
					hh->Fit(tf3,"q","q",xmin,xmax);
// 					
// 					if(Run2List[i1]==7241 || Run2List[i1]==7271 || Run2List[i1]==7273)
// 					{
// 						hh->GetXaxis()->SetRangeUser(0,15);
// 						maxbin1=hh->GetMaximumBin();
// 						hh->GetXaxis()->SetRangeUser(15,50);
// 						maxbin2=hh->GetMaximumBin();
// 						hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(maxbin1),hh->GetBinCenter(maxbin2));
// 						minbin1=hh->GetMinimumBin();
// 						hh->GetXaxis()->SetRangeUser(0,200);
// 						tf3->SetParameter(0,hh->GetBinContent(maxbin2));
// 						tf3->SetParameter(1,hh->GetBinCenter(maxbin2));
// 						tf3->SetParameter(2,20);
// 						xmin=hh->GetBinCenter(minbin1);
// 						xmax=2*hh->GetBinCenter(maxbin2)-hh->GetBinCenter(minbin1);
// 						hh->Fit(tf3,"q","q",xmin,xmax);
// 						xmin=tf3->GetParameter(1)-1.5*tf3->GetParameter(2);
// 						xmax=tf3->GetParameter(1)+1.5*tf3->GetParameter(2);
// 						if(xmin<hh->GetBinCenter(minbin1)) xmin=hh->GetBinCenter(minbin1);
// 						hh->Fit(tf3,"q","q",xmin,xmax);
// 					}
				}
				tge[i2][2][1]->SetPoint(tge[i2][2][1]->GetN(),i1,tf3->GetParameter(1));
				tge[i2][2][1]->SetPointError(tge[i2][2][1]->GetN()-1,0,tf3->GetParError(1));
				tge[i2][2][2]->SetPoint(tge[i2][2][2]->GetN(),i1,tf3->GetParameter(1));
				tge[i2][2][2]->SetPointError(tge[i2][2][2]->GetN()-1,0,tf3->GetParError(1));
				cout<<" "<<i2<<" "<<tf3->GetParameter(1)<<" "<<tf3->GetParError(1)<<endl;
				
				sprintf(hname,"%d_%d",Run3List[i1],i2);
				hh->SetName(hname);hh->SetTitle(hname);
				outroot->cd();hh->Write();
			}
			inroot->Close();
			
			
			
// 			for(int i2=0;i2<3;i2++)
// 			{
// 				if(i2==1) continue;
// 				sprintf(hname,"Integral_w_%d_6",i2);
// 				inroot->GetObject(hname,hh);
// 				tf2->SetParameter(0,500);
// 				tf2->SetParameter(1,8);
// 				tf2->SetParameter(2,2);
// 				tf2->SetParameter(3,1000);
// 				tf2->SetParameter(4,30);
// 				tf2->SetParameter(5,10);
// 				hh->Fit(tf2,"q","q",0.,fitulim[i2]);
// 				tge[i2][2][1]->SetPoint(tge[i2][2][1]->GetN(),i1,tf2->GetParameter(4));
// 				tge[i2][2][1]->SetPointError(tge[i2][2][1]->GetN()-1,0,tf2->GetParError(4));
// 				cout<<" "<<i2<<" "<<tf2->GetParameter(4)<<" "<<tf2->GetParError(4)<<endl;
// 				sprintf(hname,"%d_%d",Run3List[i1],i2);
// 				hh->SetName(hname);hh->SetTitle(hname);
// 				outroot->cd();hh->Write();
// 			}
// 			inroot->Close();
		}
	}
	outroot->cd();
	for(int i1=0;i1<3;i1++)
	{
		if(i1==1) continue;
		for(int i2=0;i2<3;i2++)
		{
			for(int i3=0;i3<3;i3++)
			{
				tge[i1][i2][i3]->Write();
			}
		}
	}
	outroot->Close();
}

void PMTCalibrationsAll2()
{
	int NR1L=28;
	int NR2L=18;
	int NR3L=21;
	
	int Run1List[28]={50,51,7074,7076,57,59,60,61,7084,70,71,74,75,76,78,80,7110,82,83,84,7126,85,7131,7132,7136,7141,7142,7143};//7087,7081,7082,7096,7099,7100,7123,7128,7129,7104,7107,7121,7140,7135,7133,7134,7144,7137,7138,7139,7092,77,86,79,81,68,7078,7089,7094,7125,7106
	int Run2List[18]={7241,7245,99,100,7254,101,102,7259,103,7264,104,105,7268,106,7271,107,7273,7281};//7252,7279,7287,92,94,96,97,93,95,98,7257,7263,108,7275,109,110,112,113,114,7266,7277 weird ,7256,7285,7283
	int Run3List[21]={116,117,118,119,120,121,7320,122,7322,7323,123,7325,124,125,126,127,128,129,7369,7394,7397};//,130,132,133,134,135,136,137,138,139 rate ,7371,7374,7375 something wrong with LXe ,7312,7314,7316 /// ,7331,7334,7336,7385,7389,7392 ,7318,7327 115,
	
	TGraphErrors* tge[3][3][3];//pmt runperiod method
	for(int i1=0;i1<3;i1++)
	{
		if(i1==1) continue;
		for(int i2=0;i2<3;i2++)
		{
			for(int i3=0;i3<2;i3++)
			{
				sprintf(hname,"PMT_%s_Run%d_M%d",PMTNames[i1].c_str(),i2+1,i3+1);
				tge[i1][i2][i3]=new TGraphErrors();
				tge[i1][i2][i3]->SetName(hname);tge[i1][i2][i3]->SetTitle(hname);
				tge[i1][i2][i3]->SetMarkerStyle(20+i3);
				tge[i1][i2][i3]->SetMarkerColor(1+i3*3+i1);
				tge[i1][i2][i3]->GetXaxis()->SetTitle("Run ID");tge[i1][i2][i3]->GetXaxis()->CenterTitle();
				tge[i1][i2][i3]->GetYaxis()->SetTitle("Single Photoelectron Charge (arbitrary units)");tge[i1][i2][i3]->GetYaxis()->CenterTitle();
			}
			sprintf(hname,"PMT_%s_Run%d_All",PMTNames[i1].c_str(),i2+1);
			tge[i1][i2][2]=new TGraphErrors();
			tge[i1][i2][2]->SetName(hname);tge[i1][i2][2]->SetTitle(hname);
			tge[i1][i2][2]->SetMarkerStyle(25);
			tge[i1][i2][2]->SetMarkerColor(4);
			tge[i1][i2][2]->GetXaxis()->SetTitle("Run ID");tge[i1][i2][2]->GetXaxis()->CenterTitle();
			tge[i1][i2][2]->GetYaxis()->SetTitle("Single Photoelectron Charge (arbitrary units)");tge[i1][i2][2]->GetYaxis()->CenterTitle();
		}
	}
	
	TFile* outroot=new TFile("PMTCalibrationsAll.root","recreate");
	
	TH1F* hh;TF1* tf1;
	TF1* tf2=new TF1("GG","gaus(0)+gaus(3)",0.,200.);
	TF1* tf3=new TF1("gaus","gaus",0.,300.);
	TF1* tf4=new TF1("gg","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-(2.*[1]))/[4])^2)",0.,300.);
	TF1* tf5=new TF1("ggg","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-(2.*[1]))/[4])^2)+[5]*exp(-0.5*((x-[6])/[7])^2)",0.,300.);
	TF1* tf6=new TF1("ggg2","gaus(0)+gaus(3)+gaus(6)",0.,300.);
// 	TF1* tf7=new TF1("gl","gaus(0)+landau(3)",0.,300.);
	TF1* tf7=new TF1("gg","gaus(0)+gaus(3)",0.,300.);
	TF1* tf8=new TF1("l","landau",0.,300.);
	TFile* inroot;
	float fitulim[3]={45,0,35};float fitlims[2]={0};
	int maxbin1=0;int maxbin2=0;int minbin1=0;
	int minbin=0;int maxbin=0;
	float xmin=0;float xmax=0;
	float fitchi21=0;float fitchi22=0;
	for(int i1=0;i1<NR1L;i1++)
	{
// 		if(Run1List[i1]==52) continue;
		cout<<i1<<" "<<Run1List[i1];
		if(Run1List[i1]<7000)
		{
			sprintf(hname,"%s/Histos/PMTCalibration_%d.root",AnalysisFilePath,Run1List[i1]);
			inroot=new TFile(hname);
			for(int i2=0;i2<3;i2++)
			{
				if(i2==1) continue;
				sprintf(hname,"SignalWindow2_%d",i2);
				inroot->GetObject(hname,hh);
				
				maxbin1=hh->GetMaximumBin();
				hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(maxbin1),30);
				minbin1=hh->GetMinimumBin();
				hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(minbin1),60);
				maxbin2=hh->GetMaximumBin();
				hh->GetXaxis()->SetRangeUser(-50,500);
				
				for(int ik1=maxbin1;ik1>=1;ik1--)
				{
					if((hh->GetBinContent(maxbin1)/hh->GetBinContent(ik1))>100){minbin=ik1;break;}
				}
				for(int ik1=maxbin2;ik1<=hh->GetNbinsX();ik1++)
				{
					if((hh->GetBinContent(maxbin2)/hh->GetBinContent(ik1))>5){maxbin=ik1;break;}
				}
				tf2->SetParameter(0,hh->GetBinContent(maxbin1));
				tf2->SetParameter(1,hh->GetBinCenter(maxbin1));
				tf2->SetParameter(2,5);
				tf2->SetParameter(3,hh->GetBinContent(maxbin2));
				tf2->SetParameter(4,hh->GetBinCenter(maxbin2));
				tf2->SetParameter(5,20);
				tf2->SetParLimits(5,0,50);
				hh->Fit(tf2,"q","q",hh->GetBinCenter(minbin),hh->GetBinCenter(maxbin));
				tge[i2][0][0]->SetPoint(tge[i2][0][0]->GetN(),i1,tf2->GetParameter(4)-tf2->GetParameter(1));
				tge[i2][0][0]->SetPointError(tge[i2][0][0]->GetN()-1,0,sqrt(pow(tf2->GetParError(4),2)+pow(tf2->GetParError(1),2)));
				tge[i2][0][2]->SetPoint(tge[i2][0][2]->GetN(),i1,tf2->GetParameter(4)-tf2->GetParameter(1));
				tge[i2][0][2]->SetPointError(tge[i2][0][2]->GetN()-1,0,sqrt(pow(tf2->GetParError(4),2)+pow(tf2->GetParError(1),2)));
// 				cout<<" "<<i2<<" "<<(tf2->GetParameter(4)-tf2->GetParameter(1))<<" "<<sqrt(pow(tf2->GetParError(4),2)+pow(tf2->GetParError(1),2))<<endl;
				
				sprintf(hname,"%d_%d",Run1List[i1],i2);
				hh->SetName(hname);hh->SetTitle(hname);
				outroot->cd();hh->Write();
			}
			inroot->Close();
		}
		else
		{
			sprintf(hname,"%s/Histos/PMTCalibration3_%d.root",AnalysisFilePath,Run1List[i1]);
			inroot=new TFile(hname);
			for(int i2=0;i2<3;i2++)
			{
				if(i2==1) continue;
				sprintf(hname,"Integral_MAgt2_wgt10_%d",i2);
				inroot->GetObject(hname,hh);
				
				if(i2==2 && Run1List[i1]<=7084)
				{
					fitlims[0]=25;fitlims[1]=55;
					if(Run1List[i1]==7084){fitlims[0]=15;fitlims[1]=45;}
					hh->Fit(tf3,"q","q",fitlims[0],fitlims[1]);
					tge[i2][0][1]->SetPoint(tge[i2][0][1]->GetN(),i1,tf3->GetParameter(1));
					tge[i2][0][1]->SetPointError(tge[i2][0][1]->GetN()-1,0,tf3->GetParError(1));
					tge[i2][0][2]->SetPoint(tge[i2][0][2]->GetN(),i1,tf3->GetParameter(1));
					tge[i2][0][2]->SetPointError(tge[i2][0][2]->GetN()-1,0,tf3->GetParError(1));
				}
				else if(i2==2)
				{
					fitlims[0]=0;fitlims[1]=50;
					hh->Fit(tf8,"q","q",fitlims[0],fitlims[1]);
					tge[i2][0][1]->SetPoint(tge[i2][0][1]->GetN(),i1,tf8->GetParameter(1));
					tge[i2][0][1]->SetPointError(tge[i2][0][1]->GetN()-1,0,tf8->GetParError(1));
					tge[i2][0][2]->SetPoint(tge[i2][0][2]->GetN(),i1,tf8->GetParameter(1));
					tge[i2][0][2]->SetPointError(tge[i2][0][2]->GetN()-1,0,tf8->GetParError(1));
				}
				else
				{
					tf4->SetParameter(0,hh->GetBinContent(hh->GetMaximumBin()));
					tf4->SetParameter(1,hh->GetBinCenter(hh->GetMaximumBin()));
					tf4->SetParameter(2,5);
					tf4->SetParameter(3,hh->GetBinContent(hh->GetMaximumBin())/10);
					tf4->SetParameter(4,20);
					
					hh->Fit(tf4,"q","q",0.,100.);
	// 				hh->Fit(tf4,"q","q",tf4->GetParameter(1)-2*tf4->GetParameter(2),2*tf4->GetParameter(1)+2*tf4->GetParameter(4));
					hh->Fit(tf4,"q","q",tf4->GetParameter(1)-2*tf4->GetParameter(2),2*tf4->GetParameter(1)+20.);
	// 				
					tge[i2][0][1]->SetPoint(tge[i2][0][1]->GetN(),i1,tf4->GetParameter(1));
					tge[i2][0][1]->SetPointError(tge[i2][0][1]->GetN()-1,0,tf4->GetParError(1));
					tge[i2][0][2]->SetPoint(tge[i2][0][2]->GetN(),i1,tf4->GetParameter(1));
					tge[i2][0][2]->SetPointError(tge[i2][0][2]->GetN()-1,0,tf4->GetParError(1));
	// 				cout<<" "<<i2<<" "<<tf3->GetParameter(1)<<" "<<tf3->GetParError(1)<<endl;
				}
				
				sprintf(hname,"%d_%d",Run1List[i1],i2);
				hh->SetName(hname);hh->SetTitle(hname);
				outroot->cd();hh->Write();
			}
			inroot->Close();
		}
	}
	for(int i1=0;i1<NR2L;i1++)
	{
// 		if(Run2List[i1]==88 || Run2List[i1]==91) continue;
// 		cout<<i1<<" "<<Run2List[i1];
		if(Run2List[i1]<7000)
		{
			sprintf(hname,"%s/Histos/PMTCalibration_%d.root",AnalysisFilePath,Run2List[i1]);
			inroot=new TFile(hname);
			for(int i2=0;i2<3;i2++)
			{
				if(i2==1) continue;
				sprintf(hname,"SignalWindow2_%d",i2);
				inroot->GetObject(hname,hh);
				
				maxbin1=hh->GetMaximumBin();
				hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(maxbin1),30);
				minbin1=hh->GetMinimumBin();
				hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(minbin1),60);
				maxbin2=hh->GetMaximumBin();
				hh->GetXaxis()->SetRangeUser(-50,500);
				
				for(int ik1=maxbin1;ik1>=1;ik1--)
				{
					if((hh->GetBinContent(maxbin1)/hh->GetBinContent(ik1))>100){minbin=ik1;break;}
				}
				for(int ik1=maxbin2;ik1<=hh->GetNbinsX();ik1++)
				{
					if((hh->GetBinContent(maxbin2)/hh->GetBinContent(ik1))>5){maxbin=ik1;break;}
				}
				tf2->SetParameter(0,hh->GetBinContent(maxbin1));
				tf2->SetParameter(1,hh->GetBinCenter(maxbin1));
				tf2->SetParameter(2,5);
				tf2->SetParameter(3,hh->GetBinContent(maxbin2));
				tf2->SetParameter(4,hh->GetBinCenter(maxbin2));
				tf2->SetParameter(5,20);
				tf2->SetParLimits(5,0,50);
				hh->Fit(tf2,"q","q",hh->GetBinCenter(minbin),hh->GetBinCenter(maxbin));
				tge[i2][1][0]->SetPoint(tge[i2][1][0]->GetN(),i1,tf2->GetParameter(4)-tf2->GetParameter(1));
				tge[i2][1][0]->SetPointError(tge[i2][1][0]->GetN()-1,0,sqrt(pow(tf2->GetParError(4),2)+pow(tf2->GetParError(1),2)));
				tge[i2][1][2]->SetPoint(tge[i2][1][2]->GetN(),i1,tf2->GetParameter(4)-tf2->GetParameter(1));
				tge[i2][1][2]->SetPointError(tge[i2][1][2]->GetN()-1,0,sqrt(pow(tf2->GetParError(4),2)+pow(tf2->GetParError(1),2)));
// 				cout<<" "<<i2<<" "<<(tf2->GetParameter(4)-tf2->GetParameter(1))<<" "<<sqrt(pow(tf2->GetParError(4),2)+pow(tf2->GetParError(1),2))<<endl;
				
				sprintf(hname,"%d_%d",Run2List[i1],i2);
				hh->SetName(hname);hh->SetTitle(hname);
				outroot->cd();hh->Write();
			}
			inroot->Close();
		}
		else
		{
			sprintf(hname,"%s/Histos/PMTCalibration3_%d.root",AnalysisFilePath,Run2List[i1]);
			inroot=new TFile(hname);
			for(int i2=0;i2<3;i2++)
			{
				if(i2==1) continue;
				sprintf(hname,"Integral_MAgt2_wgt10_%d",i2);
				inroot->GetObject(hname,hh);
				
				tf4->SetParameter(0,hh->GetBinContent(hh->GetMaximumBin()));
				tf4->SetParameter(1,hh->GetBinCenter(hh->GetMaximumBin()));
				tf4->SetParameter(2,5);
				tf4->SetParameter(3,hh->GetBinContent(hh->GetMaximumBin())/10);
				tf4->SetParameter(4,20);
				
				hh->Fit(tf4,"q","q",0.,100.);
// 				hh->Fit(tf4,"q","q",tf4->GetParameter(1)-2*tf4->GetParameter(2),2*tf4->GetParameter(1)+2*tf4->GetParameter(4));
				hh->Fit(tf4,"q","q",tf4->GetParameter(1)-2*tf4->GetParameter(2),2*tf4->GetParameter(1)+20.);
// 				
// 				if(tf4->GetChisquare()/tf4->GetNDF()>10)
// 				{
// 					maxbin1=hh->GetMaximumBin();
// 					hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(maxbin1),hh->GetBinCenter(maxbin1)+20.);
// 					minbin1=hh->GetMinimumBin();
// 					hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(minbin1),60);
// 					maxbin2=hh->GetMaximumBin();
// 					hh->GetXaxis()->SetRangeUser(-50,500);
// 					
// 					hh->Fit(tf4,"q","q",hh->GetBinCenter(minbin1),100.);
// // 					hh->Fit(tf4,"q","q",hh->GetBinCenter(minbin1),2*tf5->GetParameter(1)+2*tf5->GetParameter(4));
// 					hh->Fit(tf4,"q","q",hh->GetBinCenter(minbin1),2*tf4->GetParameter(1)+20.);
// 				}
				tge[i2][1][1]->SetPoint(tge[i2][1][1]->GetN(),i1,tf4->GetParameter(1));
				tge[i2][1][1]->SetPointError(tge[i2][1][1]->GetN()-1,0,tf4->GetParError(1));
				tge[i2][1][2]->SetPoint(tge[i2][1][2]->GetN(),i1,tf4->GetParameter(1));
				tge[i2][1][2]->SetPointError(tge[i2][1][2]->GetN()-1,0,tf4->GetParError(1));
// 				cout<<" "<<i2<<" "<<tf3->GetParameter(1)<<" "<<tf3->GetParError(1)<<endl;
				
				sprintf(hname,"%d_%d",Run2List[i1],i2);
				hh->SetName(hname);hh->SetTitle(hname);
				outroot->cd();hh->Write();
			}
			inroot->Close();
		}
	}
	for(int i1=0;i1<NR3L;i1++)
	{
// 		if(Run3List[i1]==131) continue;
		cout<<i1<<" "<<Run3List[i1];
		if(Run3List[i1]<7000)
		{
			sprintf(hname,"%s/Histos/PMTCalibration_%d.root",AnalysisFilePath,Run3List[i1]);
			inroot=new TFile(hname);
			for(int i2=0;i2<3;i2++)
			{
				if(i2==1) continue;
				sprintf(hname,"SignalWindow2_%d",i2);
				inroot->GetObject(hname,hh);
				
				maxbin1=hh->GetMaximumBin();
				hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(maxbin1),30);
				minbin1=hh->GetMinimumBin();
				hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(minbin1),60);
				maxbin2=hh->GetMaximumBin();
				hh->GetXaxis()->SetRangeUser(-50,500);
				
				for(int ik1=maxbin1;ik1>=1;ik1--)
				{
					if((hh->GetBinContent(maxbin1)/hh->GetBinContent(ik1))>100){minbin=ik1;break;}
				}
				for(int ik1=maxbin2;ik1<=hh->GetNbinsX();ik1++)
				{
					if((hh->GetBinContent(maxbin2)/hh->GetBinContent(ik1))>5){maxbin=ik1;break;}
				}
				tf2->SetParameter(0,hh->GetBinContent(maxbin1));
				tf2->SetParameter(1,hh->GetBinCenter(maxbin1));
				tf2->SetParameter(2,5);
				tf2->SetParameter(3,hh->GetBinContent(maxbin2));
				tf2->SetParameter(4,hh->GetBinCenter(maxbin2));
				tf2->SetParameter(5,20);
				tf2->SetParLimits(5,0,50);
				hh->Fit(tf2,"q","q",hh->GetBinCenter(minbin),hh->GetBinCenter(maxbin));
				tge[i2][2][0]->SetPoint(tge[i2][2][0]->GetN(),i1,tf2->GetParameter(4)-tf2->GetParameter(1));
				tge[i2][2][0]->SetPointError(tge[i2][2][0]->GetN()-1,0,sqrt(pow(tf2->GetParError(4),2)+pow(tf2->GetParError(1),2)));
				tge[i2][2][2]->SetPoint(tge[i2][2][2]->GetN(),i1,tf2->GetParameter(4)-tf2->GetParameter(1));
				tge[i2][2][2]->SetPointError(tge[i2][2][2]->GetN()-1,0,sqrt(pow(tf2->GetParError(4),2)+pow(tf2->GetParError(1),2)));
// 				cout<<" "<<i2<<" "<<(tf2->GetParameter(4)-tf2->GetParameter(1))<<" "<<sqrt(pow(tf2->GetParError(4),2)+pow(tf2->GetParError(1),2))<<endl;
				
				sprintf(hname,"%d_%d",Run3List[i1],i2);
				hh->SetName(hname);hh->SetTitle(hname);
				outroot->cd();hh->Write();
			}
			inroot->Close();
			
		}
		else
		{
			sprintf(hname,"%s/Histos/PMTCalibration3_%d.root",AnalysisFilePath,Run3List[i1]);
			inroot=new TFile(hname);
			for(int i2=0;i2<3;i2++)
			{
				if(i2==1) continue;
				sprintf(hname,"Integral_MAgt2_wgt10_%d",i2);
				inroot->GetObject(hname,hh);
				
				tf4->SetParameter(0,hh->GetBinContent(hh->GetMaximumBin()));
				tf4->SetParameter(1,hh->GetBinCenter(hh->GetMaximumBin()));
				tf4->SetParameter(2,5);
				tf4->SetParameter(3,hh->GetBinContent(hh->GetMaximumBin())/10);
				tf4->SetParameter(4,20);
				
				hh->Fit(tf4,"q","q",0.,100.);
// 				hh->Fit(tf4,"q","q",tf4->GetParameter(1)-2*tf4->GetParameter(2),2*tf4->GetParameter(1)+2*tf4->GetParameter(4));
				hh->Fit(tf4,"q","q",tf4->GetParameter(1)-2*tf4->GetParameter(2),2*tf4->GetParameter(1)+20.);
// 				
// 				if(tf4->GetChisquare()/tf4->GetNDF()>10)
// 				{
// 					maxbin1=hh->GetMaximumBin();
// 					hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(maxbin1),hh->GetBinCenter(maxbin1)+20.);
// 					minbin1=hh->GetMinimumBin();
// 					hh->GetXaxis()->SetRangeUser(hh->GetBinCenter(minbin1),60);
// 					maxbin2=hh->GetMaximumBin();
// 					hh->GetXaxis()->SetRangeUser(-50,500);
// 					
// 					hh->Fit(tf4,"q","q",hh->GetBinCenter(minbin1),100.);
// // 					hh->Fit(tf4,"q","q",hh->GetBinCenter(minbin1),2*tf5->GetParameter(1)+2*tf5->GetParameter(4));
// 					hh->Fit(tf4,"q","q",hh->GetBinCenter(minbin1),2*tf4->GetParameter(1)+20.);
// 				}
				tge[i2][2][1]->SetPoint(tge[i2][2][1]->GetN(),i1,tf4->GetParameter(1));
				tge[i2][2][1]->SetPointError(tge[i2][2][1]->GetN()-1,0,tf4->GetParError(1));
				tge[i2][2][2]->SetPoint(tge[i2][2][2]->GetN(),i1,tf4->GetParameter(1));
				tge[i2][2][2]->SetPointError(tge[i2][2][2]->GetN()-1,0,tf4->GetParError(1));
// 				cout<<" "<<i2<<" "<<tf3->GetParameter(1)<<" "<<tf3->GetParError(1)<<endl;
				
				sprintf(hname,"%d_%d",Run3List[i1],i2);
				hh->SetName(hname);hh->SetTitle(hname);
				outroot->cd();hh->Write();
			}
			inroot->Close();
		}
	}
	outroot->cd();
	for(int i1=0;i1<3;i1++)
	{
		if(i1==1) continue;
		for(int i2=0;i2<3;i2++)
		{
			for(int i3=0;i3<3;i3++)
			{
				tge[i1][i2][i3]->Write();
			}
		}
	}
	outroot->Close();
}





















