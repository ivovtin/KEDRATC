#include <Riostream.h>
#include <sstream>
#include <TROOT.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraph.h>
#include <TTree.h>
#include <TF1.h>
#include <TMath.h>
#include "TFile.h"
#include "TVirtualPad.h"
#include <iomanip>
#pragma hdrstop
#include<stdio.h>
#include<stdlib.h>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TLine.h>
#include <TEventList.h>
#include <TProfile.h>
#include <vector>
#include <TChain.h>
#include "TMinuit.h"
#include "TRandom3.h"
#include <math.h>
#include "TVirtualFitter.h"
#include <algorithm> 
#include <TMultiGraph.h>
#include <TGraphErrors.h>

using namespace std;
string progname;
/*
int Usage(string status)
{
        cout<<"Usage: "<<progname<<"\t"<<"Range ATC->1,2...   Counter First cnt->0-160     End cnt->0-160     Ntrh_{ph.e.}     Data_Region"<<endl;
        exit(0);
}
*/
void eff_ineff_single_cnt()
{

//int main(int argc, char* argv[])
//{  
//    progname=argv[0];
    int reg_aer;
    int first_cnt;
    int end_cnt;
    float npe_trh;
    int region; 
/*    if( argc>1 )
    {
     reg_aer=atoi(argv[1]);
     first_cnt=atoi(argv[2]);  
     end_cnt=atoi(argv[3]);
     npe_trh=atof(argv[4]);
     region=atoi(argv[5]);
    }
    else
    { 
      Usage(progname);
    }
*/
/*      reg_aer=8;
     first_cnt=0;  
     end_cnt=160;
     npe_trh=3.5;
     region=1;
*/
     reg_aer=4;
     first_cnt=20;  
     end_cnt=21;
     npe_trh=0.3;
     region=5;

   TString dirout="/home/ovtin/development/work_KEDR/KEDRATC/eff_ineff/output/";
    
  struct data                                  //структруа с данными АЧС
    {
      int i, t, ntrk, triggered,zero,fitted,estimated,neightrig,wlshit,nearwls,aerogel_region,aerogel_region0,
	 aerogel_region5,aerogel_region20,active_region,active_region0,active_region5,active_region20,test,
	single_aerogel_region,single_aerogel_region0,single_aerogel_region5,single_aerogel_region20,single_active_region,
	single_active_region0,single_active_region5,single_active_region20,single_test,
	in_aerogel_region,in_aerogel_region0,in_aerogel_region5,in_aerogel_region20,in_active_region,in_active_region0,in_active_region5,in_active_region20,in_test,
	   out_aerogel_region,out_aerogel_region0,out_aerogel_region5,out_aerogel_region20,out_active_region,out_active_region0,out_active_region5,out_active_region20,out_test;
	float amp, rtime,time,chi2,npe,npen,tlen,pathwls,rin,phiin,zin,rout,phiout,zout,rwls,phiwls,zwls,neighnpe,Rin_gl,
	Phiin_gl,Zin_gl,Rout_gl,Phiout_gl,Zout_gl;
    };
 data cnt;

  struct data2                                 
  { 
  int t,q,ip,nvec,nvecxy,nvecz,nhits,nhitsxy,nhitsz,nhitsvd;
  float p,pt,theta,phi,chi2,rc,xc,yc,zc,za,ph0,ph1,ph2,x0,x1,x2,y0,y1,y2,z0,z1,z2,vx,vy,vz;
  int emc_ncls, atc_ncnt;
  };
  data2 t;
  
  struct data3                                 
    { int nhits, dchits, namps, ntimes;
      float time[2], length[2], beta[2], phi[2];
      int type[2];
    };
  data3 tof;

  struct data4                                 
    {
      int ncls,ncls_trk,nlkr,ncsi,nstrcls,nstrtrk;
      float energy,elkr,ecsi;
    };
  data4 emc;

  struct data5                                 
    {
      int c,lkr,csi,ncells;
      float e,x,y,z,vx,vy,vz,theta,phi,rho,dtheta,dphi,drho,thetastr,phistr;
      int qlty,ncellsbad,str_ncls,str_ntrk,dc_ntrk,emc_ncls;
    };
  data5 c0;

  struct data6                                 
    {
      int c,lkr,csi,ncells;
      float e,x,y,z,vx,vy,vz,theta,phi,rho,dtheta,dphi,drho,thetastr,phistr;
      int qlty,ncellsbad,str_ncls,str_ntrk,dc_ntrk,emc_ncls;
    };
  data6 c1;

  struct data7                                 
    {
      int nhits,dcmuhits,octant,layer;
      int status;
    };
  data7 mu;

  struct data8                                 
    {   
    int event, evdaq; //succesive event and DAQ event numbers
    int run;     //guess what it is
    int quality; //event quality number, set by user
    };
  data8 ev;

   struct data9                                 
    {
      	int ntrk, nip, nbeam;
	float x, y, z;
	float sig_x, sig_y, sig_z;
	float theta2t, phi2t;
    };
  data9 vrt;

 char branchname[161];
 char branchname1[161];
 char branchname2[161];
 char branchname3[161];
 char branchname4[161];
 char branchname5[161];
 char branchname6[161];
 char branchname7[161];
 char branchname8[161];
/*
   TH1F *h1[160]=new TH1F("h1","Cosmic muon momentum",200,0,3000);
   TH1F *h2[160]=new TH1F("h2","<<Kaons>>",500,0,50);
   TH1F *h3[160]=new TH1F("h3","<<#pi>>",500,0,50);
   TH1F *h4[160]=new TH1F("h4","E/p",200,0,2);
   TH1F *h5[160]=new TH1F("h5","Cosmic muon energy",200,0,1000);
   TProfile *pr0[160]=new TProfile("pr0","Momentum dependence of the	amplitude for #mu",200,0,3000,0,50);
   TProfile *pr1[160]=new TProfile("pr1","Momentum dependence of the	amplitude for K and #pi",200,0,3000,0,50);
   TProfile *pr2[160]=new TProfile("pr2","Momentum dependence of the	for amplitude K and #pi",200,0,3000,0,50);
*/
/*
   TH1F *h1[160];
   TH1F *h2[160];
   TH1F *h3[160];
   TH1F *h4[160];
   TH1F *h5[160];
   TProfile *pr0[160];
   TProfile *pr1[160];
   TProfile *pr2[160];
*/
 TChain *tt1 = new TChain("tt");                            
  if(region==1){
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_1.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_2.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_3.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_4.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_5.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_6.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_7.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_8.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_9.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_9_1.root");
    //}
  //if(region==2){
    //2range
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_10.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_11.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_12.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_13.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_14.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_15.root");
    //}
  //if(region==3){
    //3 range
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_16.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_17.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_18.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_19.root");
    }
  if(region==4){
    //4 range
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_1.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_2.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_3.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_4.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_5.root");
    }
  if(region==5){
    //5 range
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_6.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_7.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_1.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_2.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_3.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_4.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_5.root");
  }
  if(region==6){
    //6 range
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_6.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_7.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_1.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_2.root");
    }
  if(region==7){
    //7 range
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_4.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_5.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_1.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_2.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_1.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_2.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_3.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_4.root");
    }
  if(region==8){
    //8 range
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_1.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_2.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_1.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_2.root");
    }
  if(region==9){
    //9 range
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_1.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_2.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_3.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_4.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_5.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_1.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_2.root");
    }
  if(region==10){
    //10 range
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_1.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_2.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_3.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_4.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_5.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_1.root");
    }
  if(region==11){
    //11 range
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_2.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_3.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_4.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_1.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_2.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_3.root");
    }
  if(region==12){
    //12 range
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_1.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_2.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_3.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_1.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_2.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_3.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_4.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_5.root");
    }
  if(region==13){
    //13 range
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_6.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_7.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_1.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_2.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_3.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_5.root");
    }
  if(region==14){
    //14 range
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_1.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_2.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_3.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_4.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_1.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_2.root");
    }
  if(region==15){
    //15 range
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_3.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_4.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_5.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_6.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_1.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_2.root");
    tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_3.root");
    }

 TChain *tt2 = new TChain("tt");                            
 if(region==1){
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_3.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_4.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_5.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_6.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_7.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_8.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_9.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_9_1.root");
    //}
    //if(region==2){
    //2range
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_10.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_11.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_12.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_13.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_14.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_15.root");
    //}
    //if(region==3){
    //3 range
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_16.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_17.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_18.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_19.root");
    }
  if(region==4){
    //4 range
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_3.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_4.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_5.root");
    }
  if(region==5){
    //5 range
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_6.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_7.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_3.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_4.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_5.root");
    }
  if(region==6){
    //6 range
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_6.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_7.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_2.root");
    }
  if(region==7){
    //7 range
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_4.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_5.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_3.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_4.root");
    }
  if(region==8){
    //8 range
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_2.root");
    }
  if(region==9){
    //9 range
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_3.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_4.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_5.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_2.root");
    }
  if(region==10){
    //10 range
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_3.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_4.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_5.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_1.root");
    }
  if(region==11){
    //11 range
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_3.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_4.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_3.root");
    }
  if(region==12){
    //12 range
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_3.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_3.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_4.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_5.root");
    }
  if(region==13){
    //13 range
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_6.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_7.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_3.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_5.root");
    }
  if(region==14){
    //14 range
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_3.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_4.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_2.root");
    }
  if(region==15){
    //15 range
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_3.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_4.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_5.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_6.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_3.root");
    }

  TChain *tt3=new TChain("emct");
  if(region==1){
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_3.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_4.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_5.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_6.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_7.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_8.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_9.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_9_1.root");
    //}
    //if(region==2){
    //2range
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_10.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_11.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_12.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_13.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_14.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_15.root");
    //}
    //if(region==3){
    //3 range
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_16.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_17.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_18.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_19.root");
    }
  if(region==4){
    //4 range
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_3.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_4.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_5.root");
    }
  if(region==5){
    //5 range
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_6.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_7.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_3.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_4.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_5.root");
    }
  if(region==6){
    //6 range
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_6.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_7.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_2.root");
    }
  if(region==7){
    //7 range
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_4.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_5.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_3.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_4.root");
    }
  if(region==8){
    //8 range
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_2.root");
    }
  if(region==9){
    //9 range
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_3.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_4.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_5.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_2.root");
    }
  if(region==10){
    //10 range
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_3.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_4.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_5.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_1.root");
  }
  if(region==11){
    //11 range
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_3.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_4.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_3.root");
    }
  if(region==12){
    //12 range
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_3.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_3.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_4.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_5.root");
    }
  if(region==13){
    //13 range
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_6.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_7.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_3.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_5.root");
    }
  if(region==14){
    //14 range
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_3.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_4.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_2.root");
    }
  if(region==15){
    //15 range
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_3.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_4.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_5.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_6.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_3.root");
    }

	tt1->SetBranchStatus("*",0);
	tt1->SetBranchStatus("t",1);
	tt1->SetBranchStatus("tof",1);
	tt1->SetBranchStatus("emc",1);
	tt1->SetBranchStatus("mu",1);
	tt1->SetBranchStatus("ev",1);
	tt1->SetBranchStatus("vrt",1);

	tt1->SetBranchAddress("t",&t);              //СЏв”ђСЏв”‚СЏв”ЊРїв•џРїв•«Рїв•џРїв•ЎРїв•©Рїв•¦Рїв•ЎРїв•џРїв•ЈРїв•Є 
	tt1->SetBranchAddress("tof",&tof); 
	tt1->SetBranchAddress("emc",&emc);  
	tt1->SetBranchAddress("mu",&mu); 
	tt1->SetBranchAddress("ev",&ev); 
	tt1->SetBranchAddress("vrt",&vrt); 

	tt3->SetBranchStatus("c0",1);
	tt3->SetBranchStatus("c1",1);

	tt3->SetBranchAddress("c0",&c0);  
	tt3->SetBranchAddress("c1",&c1);        
	
	float P, Pkaon, Ppion, mkaon=493.66, mpion=139.57, mmuon=105.65;
	float eff1, mu_eff1;
	float eff2, mu_eff2;
	float mu_pion, mu_kaon;

	float eff_teoretic, mu_eff;

	//float npe_trh=0.1;
	//float npe_trh=1.7;
	const int kk=20;
   
	float eff_mu[kk];
	float err_eff_mu[kk];
	double num_npenotzero_mu[kk];
	double num_npetotal_mu[kk];
	float p_all[kk];
	float err_p_all[kk];

	float eff_pion[kk];
	float not_eff_pion[kk];
	float err_eff_pion[kk];
	double num_npenotzero_pion[kk];
	double num_npetotal_pion[kk];

	float eff_kaon[kk];
	float not_eff_kaon[kk];
	float err_eff_kaon[kk];
	double num_npenotzero_kaon[kk];
	double num_npetotal_kaon[kk];
	float Ksigma[kk];

   	for(int kk3=0; kk3<=kk-1; kk3++)
   	{
     		eff_mu[kk3]=0;
     		err_eff_mu[kk3]=0;
     		num_npenotzero_mu[kk3]=0;
     		num_npetotal_mu[kk3]=0;

     		eff_pion[kk3]=0;
     		not_eff_pion[kk3]=0;
     		err_eff_pion[kk3]=0;
     		num_npenotzero_pion[kk3]=0;
     		num_npetotal_pion[kk3]=0;

     		eff_kaon[kk3]=0;
     		not_eff_kaon[kk3]=0;
     		err_eff_kaon[kk3]=0;
     		num_npenotzero_kaon[kk3]=0;
     		num_npetotal_kaon[kk3]=0;

     		Ksigma[kk3]=0;
   	}

for (int counter=first_cnt; counter<end_cnt; counter++)
{ 

   	TFile *fout=0;
   	TString fname;
  	fname=TString::Format(dirout+"cnt_single_%d_eff_ineff.root",counter).Data();
   	fout=new TFile(fname,"RECREATE");
   	cout<<fname<<endl; 
/*
   	TCanvas *cc1 = new TCanvas("cc1",TString::Format("cnt_single_%d_eff_ineff.root",counter).Data(),1600,1200);
    	//gStyle->SetOptStat(0);
   	cc1->Divide(4,3);
*/
        sprintf(branchname4,"h1_Cosmic muon momentum_cnt%d",counter);
        sprintf(branchname5,"h2_<<Kaons>>_cnt%d",counter);
        sprintf(branchname6,"h3_<<#pi>>_cnt%d",counter);
        sprintf(branchname7,"h4_E/p_cnt%d",counter);
        sprintf(branchname8,"h5_Cosmic muon energy_cnt%d",counter);
   	TH1F *h1=new TH1F(branchname4,branchname4,200,0,3000);
   	TH1F *h2=new TH1F(branchname5,branchname5,500,0,50);
   	TH1F *h3=new TH1F(branchname6,branchname6,500,0,50);
   	TH1F *h4=new TH1F(branchname7,branchname7,200,0,2);
   	TH1F *h5=new TH1F(branchname8,branchname8,200,0,1000);
        sprintf(branchname1,"pr0_#mu_cnt%d",counter);
        sprintf(branchname2,"pr1_K and #pi_cnt%d",counter);
        sprintf(branchname3,"pr2_K and #pi_cnt%d",counter);
  	TProfile *pr0=new TProfile(branchname1,branchname1,200,0,3000,0,50);
  	TProfile *pr1=new TProfile(branchname2,branchname2,200,0,3000,0,50);
   	TProfile *pr2=new TProfile(branchname3,branchname3,200,0,3000,0,50);

   	int total_entries=0;
       
  for (int n=0; n<9; n++)
  {
	 sprintf(branchname,"cnt%d",n);       
	 tt2->SetBranchStatus("*",0);		
	 tt2->SetBranchStatus(branchname,1);
	 tt2->SetBranchAddress(branchname,&cnt);       //получаем пересечение    
  	 printf("Counter %d\r",n); fflush(stdout); 
    
	 Int_t nentr=tt1->GetEntries();               //число событий в дереве
	 cout<<"crossing of the counter =>"<<n<<"  Nentries="<<nentr<<endl;

    //for(int k=0; k<nentr; k++)                               //цикл по всем событиям в дереве
    for(int k=0; k<4000000; k++)                               //цикл по всем событиям в дереве
    {
      tt1->GetEntry(k);
      tt2->GetEntry(k);
    
      //if( cnt.i>=first_cnt && cnt.i<end_cnt )  //цикл по числу счетчиков
      if( cnt.i==counter )  
      {	     
        	int ATCregion;
        	if (reg_aer==1) {ATCregion=cnt.single_aerogel_region0; }
        	if (reg_aer==2) {ATCregion=cnt.aerogel_region0;  }
       	 	if (reg_aer==3) {ATCregion=cnt.single_aerogel_region5; }
      		if (reg_aer==4) {ATCregion=cnt.aerogel_region5;  }
        	if (reg_aer==5) {ATCregion=cnt.single_aerogel_region; }
        	if (reg_aer==6) {ATCregion=cnt.aerogel_region;  }
        	if (reg_aer==7) {ATCregion=cnt.single_aerogel_region20; }
        	if (reg_aer==8) {ATCregion=cnt.aerogel_region20;  }
        	if (reg_aer==9) {ATCregion=cnt.single_active_region0;  }

        	if( t.p<=1000000 && ATCregion!=0 && cnt.wlshit!=1 && t.nvec>4 && t.chi2<5 && sqrt(pow(t.x0,2)+pow(t.y0,2)+pow(t.z0,2))<40 && t.emc_ncls>1 && tof.nhits>=1 && mu.octant>=1 && mu.status>=1 && mu.nhits>1 && mu.layer>=1 && emc.energy>50 && emc.ncls_trk>=1 ){
		//if( t.p<=1000000 && ATCregion!=0 && cnt.wlshit!=1 && t.nvec>4 && t.chi2<150 && sqrt(pow(t.x0,2)+pow(t.y0,2)+pow(t.z0,2))<40 && t.emc_ncls>1 && tof.nhits>=1 && mu.octant>=1 && mu.status>=1 && mu.nhits>1 && mu.layer>=1 && emc.energy>50 && emc.ncls_trk>=1 ){
            	//if ( t.p<=1000000 && ATCregion!=0 && cnt.wlshit!=1 && t.nvec>4 && t.chi2<150 && sqrt(pow(t.x0,2)+pow(t.y0,2)+pow(t.z0,2))<40 && t.emc_ncls>1 && tof.nhits>=1 && mu.octant>=1 && mu.status>=1 && mu.nhits>1 && mu.layer>=1 && emc.energy>50 && emc.ncls_trk>=1 && (cnt.Zout_gl+cnt.Zout_gl)/2>30){
            	//if (t.p<=3000 && cnt.aerogel_region0!=0 && cnt.wlshit!=1 && t.nvec>4 && t.chi2<100 && sqrt(pow(t.x0,2)+pow(t.y0,2)+pow(t.z0,2))<20 && t.emc_ncls>1 && cnt.npe>0.5){
    			  	tt3->GetEntry(k);
					  //if( c0.e>50 && c1.e>50 ) 
            			//{ 
	             		//cout<<"Run="<<ev.run<<"\t"<<"Event="<<ev.event<<"\t"<<"Evdaq="<<ev.evdaq<<"\t"<<endl;
              			Pkaon=t.p*(mkaon/mmuon);
              			Ppion=t.p*(mpion/mmuon);
              			//cout<<t.p<<"\t"<<Pkaon<<"\t"<<Ppion<<endl;
              			//cout<<cnt.Zout_gl<<endl;

              			h1->Fill(t.p);
              			if(t.p>200 && t.p<300){
              				h2->Fill(cnt.npe);
              				//  cout<<"Run="<<ev.run<<"\t"<<"Event="<<ev.event<<"\t"<<"Evdaq="<<ev.evdaq<<"\t"<<endl;
              				//  cout<<"cnt="<<n<<"\t"<<"npe="<<cnt.npe<<"\t"<<"P="<<t.p<<"\t"<<endl;
              			}
              			if(t.p>700 && t.p<1100){
                			h3->Fill(cnt.npe);
              			}
              			h4->Fill(emc.energy/t.p);
              			h5->Fill(emc.energy);
              			pr0->Fill(t.p,cnt.npe);
              			pr1->Fill(Ppion,cnt.npe);
              			pr2->Fill(Pkaon,cnt.npe);

              			total_entries=++total_entries;
    
           						for(int ii=0; ii<=kk-1; ii++)
     							{
     							    p_all[ii]=100*ii+50;
      							  err_p_all[ii]=50;
      							  //cout<<"p_all[ii]"<<p_all[ii]<<endl;
        							if(t.p>100*ii&&t.p<100+100*ii){ 
        								if(cnt.npe>npe_trh){num_npenotzero_mu[ii]=++num_npenotzero_mu[ii];}
                							  num_npetotal_mu[ii]=++num_npetotal_mu[ii];
                  							  eff_mu[ii]=num_npenotzero_mu[ii]/num_npetotal_mu[ii];
                  							  err_eff_mu[ii]=sqrt(num_npenotzero_mu[ii])/num_npetotal_mu[ii];
                  							  //cout<<ii<<"\t"<<"num_npetotal_mu[ii]="<<num_npetotal_mu[ii]<<"\t"<<"eff_mu[ii]="<<eff_mu[ii]<<endl;
        							}
        							if(Ppion>100*ii&&Ppion<100+100*ii){ 
        								if(cnt.npe>npe_trh){num_npenotzero_pion[ii]=++num_npenotzero_pion[ii];}
                							    num_npetotal_pion[ii]=++num_npetotal_pion[ii];
                  							  eff_pion[ii]=num_npenotzero_pion[ii]/num_npetotal_pion[ii];
                  							  err_eff_pion[ii]=sqrt(num_npenotzero_pion[ii])/num_npetotal_pion[ii];
                  							  not_eff_pion[ii]=1-eff_pion[ii];
                                  //if(eff_pion[ii]<-1 || eff_pion[ii]>1)
                                  //{
                                  //  cout<<ii<<"\t"<<"num_npenotzero_pion[ii]="<<num_npenotzero_pion[ii]<<"\t"<<"num_npetotal_pion[ii]="<<num_npetotal_pion[ii]<<"\t"<<"eff_pion[ii]="<<eff_pion[ii]<<"\t"<<"err_eff_pion[ii]="<<err_eff_pion[ii]<<endl;
                                 // }
        							}
							        if(Pkaon>100*ii&&Pkaon<100+100*ii){ 
 								       	if(cnt.npe>npe_trh){num_npenotzero_kaon[ii]=++num_npenotzero_kaon[ii];}
						                	    num_npetotal_kaon[ii]=++num_npetotal_kaon[ii];
						                  	  eff_kaon[ii]=num_npenotzero_kaon[ii]/num_npetotal_kaon[ii];
						                  	  err_eff_kaon[ii]=sqrt(num_npenotzero_kaon[ii])/num_npetotal_kaon[ii];
						                  	  not_eff_kaon[ii]=1-eff_kaon[ii];
								      }
							         Ksigma[ii]=sqrt(2)*pow(TMath::Erf(-1+2*(1-eff_kaon[ii])),-1)+pow(TMath::Erf(-1+2*eff_pion[ii]),-1);
     					    		} 
          //}
    				}
  			}
		}
	}

//	cc1->cd(1);
	h1->SetXTitle("P, MeV/c");
	h1->Draw("hist");
	TLine *l1=new TLine(200,0,200,4500);
	l1->SetLineColor(kRed);
	l1->Draw();
	TLine *l2=new TLine(300,0,300,4500);
	l2->SetLineColor(kRed);
	l2->Draw();
	TLine *l3=new TLine(700,0,700,4500);
	l3->SetLineColor(kBlue);
	l3->Draw();
	TLine *l4=new TLine(1100,0,1100,4500);
	l4->SetLineColor(kBlue);
	l4->Draw();

//	cc1->cd(2);
	h5->SetXTitle("E, MeV");
	h5->Draw("hist");

//	cc1->cd(3);
	h4->Draw("hist");

//	cc1->cd(4);
	TF1* myfit0=new TF1("myfit0","[0]+[1]*(x^2-[2]^2)/(x^2)*(TMath::Erf(x-[2])+1)/2",200.,3000.);
	myfit0->SetLineColor(kBlue);
	myfit0->SetParameter(0,pr0->GetMinimum());
	myfit0->SetParameter(1,pr0->GetMaximum());
	myfit0->SetParameter(2,200);
	myfit0->SetParNames("#mu_{0}","#mu_{max}","P_{thr}");
	pr0->Fit("myfit0","","",200,3000);
	cout<<"myfit0->GetParameter(0)="<<myfit0->GetParameter(0)<<"\t"<<"myfit0->GetParameter(1)"<<myfit0->GetParameter(1)<<endl;
	pr0->SetXTitle("P, MeV/c");
	pr0->SetYTitle("N_{ph.e.}");
	gStyle->SetOptFit(1011);
	//gStyle->SetOptStat(1011);
	pr0->Draw("prof");

//	cc1->cd(5);
	gPad->SetLogy();
	h2->SetLineColor(kRed);
	h2->Draw("hist");

//	cc1->cd(6);
	h3->Draw("hist");

//	cc1->cd(7);
	TF1* myfit=new TF1("myfit1","[0]+[1]*(x^2-[2]^2)/(x^2)*(TMath::Erf(x-[2])+1)/2",200.,3000.);
	myfit->SetLineColor(kBlue);
	myfit->SetParameter(0,pr1->GetMinimum());
	myfit->SetParameter(1,pr1->GetMaximum());
	myfit->SetParameter(2,200);
	myfit->SetParNames("#mu_{0}","#mu_{max}","P_{thr}");
	pr1->Fit("myfit1","","",200,3000);
	pr1->SetXTitle("P, MeV/c");
	pr1->SetYTitle("N_{ph.e.}");

	TF1* myfit2=new TF1("myfit2","[0]+[1]*(x^2-[2]^2)/(x^2)*(TMath::Erf(x-[2])+1)/2",520.,3000.);
	myfit2->SetLineColor(kRed);
	myfit2->SetParameter(0,pr2->GetMinimum());
	myfit2->SetParameter(1,pr2->GetMaximum());
	myfit2->SetParameter(2,200);
	myfit2->SetParNames("#mu_{0}","#mu_{max}","P_{thr}");
	pr2->Fit("myfit2","","",520,3000);
	pr2->SetLineColor(kRed);
	float Pthr=myfit->GetParameter(2);
	float ref_index;
	ref_index=sqrt(1+pow((105/Pthr),2));
	float mu_0=myfit->GetParameter(0);
	float mu_max=myfit->GetParameter(1);
	//cout<<"mu_ef="<<mu_ef<<"\t"<<"eff="<<eff<<"\t"<<"ineff="<<1-eff<<endl;
	pr1->Draw("prof");
	pr2->Draw("same");

//	cc1->cd(8);
	TMultiGraph *mg1 = new TMultiGraph();
	//mg2->SetTitle("Efficiency&Momentum");
	TGraphErrors* gr1=new TGraphErrors(kk,p_all,eff_mu,err_p_all,err_eff_mu);
	gr1->SetMarkerStyle(20);
	gr1->SetMarkerColor(3);
	gr1->SetLineWidth(2);
	gr1->SetLineColor(3);
	gr1->SetTitle("Efficiency&Momentum");
	gr1->GetXaxis()->SetTitle("P, MeV/c");
	gr1->GetYaxis()->SetTitle("Efficiency");
	//gr1->Draw("ap");

	TGraphErrors* gr2=new TGraphErrors(kk,p_all,eff_pion,err_p_all,err_eff_pion);
	gr2->SetMarkerStyle(20);
	gr2->SetMarkerColor(4);
	gr2->SetLineWidth(2);
	gr2->SetLineColor(4);
	gr2->SetTitle("Efficiency&Momentum");
	gr2->GetXaxis()->SetTitle("P, MeV/c");
	gr2->GetYaxis()->SetTitle("Efficiency");
	//gr2->Draw("same");

	TGraphErrors* gr3=new TGraphErrors(kk,p_all,eff_kaon,err_p_all,err_eff_kaon);
	gr3->SetMarkerStyle(20);
	gr3->SetMarkerColor(2);
	gr3->SetLineWidth(2);
	gr3->SetLineColor(2);
	gr3->SetTitle("Efficiency&Momentum");
 
	mg1->Add(gr1);
	mg1->Add(gr2);
	mg1->Add(gr3);
	mg1->SetTitle("Efficiency&Momentum; P, MeV/c; Efficiency");
	// draw the multigraph
	mg1->Draw("ap");

//	cc1->cd(9);
	TMultiGraph *mg2 = new TMultiGraph();
	TGraphErrors* gr4=new TGraphErrors(kk,p_all,not_eff_pion,err_p_all,0);
	gr4->SetMarkerStyle(20);
	gr4->SetMarkerColor(4);
	gr4->SetLineWidth(2);
	gr4->SetLineColor(4);

	TGraphErrors* gr5=new TGraphErrors(kk,p_all,not_eff_kaon,err_p_all,0);
	gr5->SetMarkerStyle(20);
	gr5->SetMarkerColor(2);
	gr5->SetLineWidth(2);
	gr5->SetLineColor(2);
 
	mg2->Add(gr4);
	mg2->Add(gr5);
	mg2->SetTitle("K efficiency and #pi Misidentification; P, MeV/c; K efficiency and #pi Misidentification");
	// draw the multigraph
	mg2->Draw("ap");
	TLine *l5=new TLine(650,0,650,1);
	l5->SetLineColor(kRed);
	l5->Draw();
	TLine *l6=new TLine(1500,0,1500,1);
	l6->SetLineColor(kRed);
	l6->Draw();

//	cc1->cd(10);
	TMultiGraph *mg3 = new TMultiGraph();
	TGraphErrors* gr6=new TGraphErrors(kk,p_all,Ksigma,err_p_all,0);
	gr6->SetMarkerStyle(20);
	gr6->SetMarkerColor(3);
	gr6->SetLineWidth(2);
	gr6->SetLineColor(3);
	gr6->SetTitle("#sigma");
	mg3->Add(gr6);
	mg3->SetTitle("#sigma&Momentum; P, MeV/c; #sigma");
	mg3->Draw("ap");

//	cc1->cd();
//	cc1->SaveAs(TString::Format("cnt_single_%d_eff_ineff_canvas.root",counter).Data());

	//cout<<"Cnt="<<counter<<"\t"<<"Ksigma[6]="<<Ksigma[6]<<"\t"<<"Ksigma[7]="<<Ksigma[7]<<"\t"<<"Ksigma[8]="<<Ksigma[8]<<"\t"<<"Ksigma[9]="<<Ksigma[9]<<"\t"<<"Ksigma[10]="<<Ksigma[10]<<"\t"<<"Ksigma[11]="<<Ksigma[11]<<"\t"<<"Ksigma[12]="<<Ksigma[12]<<"\t"<<"Ksigma[13]="<<Ksigma[13]<<"\t"<<"Ksigma[14]="<<Ksigma[14]<<"\t"<<"average="<<average_Ksigma<<endl;
	cout<<"Cnt="<<counter<<"\t"<<"Nthr_{ph.e.}="<<npe_trh<<"\t"<<"eff_mu[11]="<<eff_mu[11]<<"\t"<<"eff_pion[11]="<<eff_pion[11]<<"\t"<<"eff_kaon[11]="<<eff_kaon[11]<<"\t"<<"sigma[11]="<<Ksigma[11]<<"\t"<<"total_entries="<<total_entries<<endl;

	for(int kk4=5; kk4<=15; kk4++)
	{
	ofstream of(TString::Format(dirout+"cnt%d-%d_single_npetrh%f_eff_ineff_p%d.dat",first_cnt,end_cnt,npe_trh,kk4).Data(),ios::app);         
	of<<counter<<"\t"<<npe_trh<<"\t"<<myfit0->GetParameter(1)<<"\t"<<eff_mu[kk4]<<"\t"<<eff_pion[kk4]<<"\t"<<eff_kaon[kk4]<<"\t"<<not_eff_pion[kk4]<<"\t"<<not_eff_kaon[kk4]<<"\t"<<Ksigma[kk4]<<"\t"<<total_entries<<endl;
	of.close();
	}
    
   	pr0->Write("N_{ph.e.}&Momentum");
   	mg1->Write("Efficiency&Momentum");
    	mg2->Write("K efficiency and #pi Misidentification");
    	mg3->Write("#sigma&Momentum");

   	fout->Write();                                                      //ÐÉÛÅÍ ÄÁÎÎÙÅ × ÆÁÊÌ É ÚÁËÒÙ×ÁÅÍ ÅÇÏ
   	fout->Close();
   
   	for(int kk1=0; kk1<=kk-1; kk1++)
   	{
	   eff_mu[kk1]=0;
	   err_eff_mu[kk1]=0;
	   num_npenotzero_mu[kk1]=0;
	   num_npetotal_mu[kk1]=0;

	   eff_pion[kk1]=0;
	   not_eff_pion[kk1]=0;
	   err_eff_pion[kk1]=0;
	   num_npenotzero_pion[kk1]=0;
	   num_npetotal_pion[kk1]=0;

	   eff_kaon[kk1]=0;
	   not_eff_kaon[kk1]=0;
	   err_eff_kaon[kk1]=0;
	   num_npenotzero_kaon[kk1]=0;
	   num_npetotal_kaon[kk1]=0;

	   Ksigma[kk1]=0; 
	}
}
}
