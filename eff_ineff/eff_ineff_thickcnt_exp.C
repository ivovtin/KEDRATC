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
//#include <math.h>
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

void eff_ineff_thickcnt_exp()
{
//int main(int argc, char* argv[])
//{  
//    progname=argv[0];
    int reg_aer;
    int counter;
    int end_cnt;
    float npe_trh;
    int region; 
/*    if( argc>1 )
    {
     reg_aer=atoi(argv[1]);
     counter=atoi(argv[2]);
     end_cnt=atoi(argv[3]);  
     npe_trh=atof(argv[4]);
     region=atoi(argv[5]);
    }
    else
    { 
    Usage(progname);
    }
*/
    bool sim=0;  
    reg_aer=4;
    counter=20;
    end_cnt=140;  
    npe_trh=0.7;
    region=11;  

    TString dirout="/home/ovtin/development/work_KEDR/eff_ineff/";

    int first_cnt=counter;

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
  { 
    int nhits, dchits, namps, ntimes;
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

TChain *tt1=new TChain("tt");
if(sim!=1){
//1range
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
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_19.root");
}
//2range
if(region==2){
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_5.root");
}
//3range
if(region==3){
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_4.root");
}
//4range
if(region==4){
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_2.root");
}
//5range
if(region==5){
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_1.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_2.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_3.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_4.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_5.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_1.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_2.root");
}
//6range
if(region==6){
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_1.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_2.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_3.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_4.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_5.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_1.root");
}
//7range
if(region==7){
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_2.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_3.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_4.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_1.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_2.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_3.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_1.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_2.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_3.root");
}
//8range
if(region==8){

tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_1.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_2.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_3.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_4.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_5.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_6.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_7.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_1.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_2.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_3.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_5.root");
}
//9range
if(region==9){
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_1.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_2.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_3.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_4.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_1.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_2.root");
}
//10range
if(region==10){
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_3.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_4.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_5.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_6.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_1.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_2.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_3.root");
}
//11range
if(region==11){
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_1.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_2.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_3.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_4.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_5.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_6.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_7.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_1.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_2.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_3.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_4.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_5.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_6.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_7.root");
}
//12range
if(region==12){
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_1.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_2.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_3.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_4.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_5.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_6.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_7.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_1.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_2.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_3.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_4.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_5.root");
}
//13range
if(region==13){
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_1.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_2.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_3.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_4.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_5.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_1.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_2.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_3.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_4.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_5.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_6.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_1.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_2.root");
tt1->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_3.root");
}
}
else{
  if(region==1){
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_19.root");
  }
  if(region==2){
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_19.root");
  }
  if(region==3){
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_19.root");
  }
  if(region==4){
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_10.root");  
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_19.root");
  }
  if(region==5){
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_19.root");
  }
  if(region==6){
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_19.root");
  }
  if(region==7){
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_19.root");
  }
  if(region==8){
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_19.root");
  }  
  if(region==9){
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_19.root");
  }  
  if(region==10){
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_19.root");
  }  
  if(region==11){
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_19.root");
  }  
  if(region==12){
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_19.root");
  }  
  if(region==13){
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_19.root");
  }  
  if(region==14){
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_19.root");
  }  
  if(region==15){  
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_19.root");
  }
  if(region==16){  
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_19.root");
  }
  if(region==17){  
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_19.root");
  }
  if(region==18){  
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_19.root");
  }
  if(region==19){  
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_19.root");
  }
  if(region==20){  
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_19.root");
  }
  if(region==21){  
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_19.root");
  }
  if(region==22){  
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_19.root");
  }
  if(region==23){  
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_19.root");
  }
  if(region==24){  
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_19.root");
  }
  if(region==25){  
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_1.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_2.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_3.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_4.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_5.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_6.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_7.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_8.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_9.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_10.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_11.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_12.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_13.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_14.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_15.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_16.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_17.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_18.root");
  tt1->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_19.root");
  }             
}

TChain *tt2=new TChain("tt");
if(sim!=1){
//1range
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
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_19.root");
}
//2range
if(region==2){
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_5.root");
}
//3range
if(region==3){
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_4.root");
}
//4range
if(region==4){
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_2.root");
}
//5range
if(region==5){
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_1.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_2.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_3.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_4.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_5.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_1.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_2.root");
}
//6range
if(region==6){
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_1.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_2.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_3.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_4.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_5.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_1.root");
}
//7range
if(region==7){
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_2.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_3.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_4.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_1.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_2.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_3.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_1.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_2.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_3.root");
}
//8range
if(region==8){

tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_1.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_2.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_3.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_4.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_5.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_6.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_7.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_1.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_2.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_3.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_5.root");
}
//9range
if(region==9){
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_1.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_2.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_3.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_4.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_1.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_2.root");
}
//10range
if(region==10){
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_3.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_4.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_5.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_6.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_1.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_2.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_3.root");
}
//11range
if(region==11){
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_1.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_2.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_3.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_4.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_5.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_6.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_7.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_1.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_2.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_3.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_4.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_5.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_6.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_7.root");
}
//12range
if(region==12){
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_1.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_2.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_3.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_4.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_5.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_6.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_7.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_1.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_2.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_3.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_4.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_5.root");
}
//13range
if(region==13){
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_1.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_2.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_3.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_4.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_5.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_1.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_2.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_3.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_4.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_5.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_6.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_1.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_2.root");
tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_3.root");
}
}
else{
    if(region==1){
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_19.root");
  }
  if(region==2){
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_19.root");
  }
  if(region==3){
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_19.root");
  }
  if(region==4){
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_10.root");  
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_19.root");
  }
  if(region==5){
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_19.root");
  }
  if(region==6){
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_19.root");
  }
  if(region==7){
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_19.root");
  }
  if(region==8){
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_19.root");
  }  
  if(region==9){
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_19.root");
  }  
  if(region==10){
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_19.root");
  }  
  if(region==11){
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_19.root");
  }  
  if(region==12){
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_19.root");
  }  
  if(region==13){
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_19.root");
  }  
  if(region==14){
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_19.root");
  }  
  if(region==15){  
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_19.root");
  }
  if(region==16){  
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_19.root");
  }
  if(region==17){  
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_19.root");
  }
  if(region==18){  
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_19.root");
  }
  if(region==19){  
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_19.root");
  }
  if(region==20){  
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_19.root");
  }
  if(region==21){  
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_19.root");
  }
  if(region==22){  
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_19.root");
  }
  if(region==23){  
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_19.root");
  }
  if(region==24){  
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_19.root");
  }
  if(region==25){  
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_8.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_9.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_15.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_19.root");
  }               
}

TChain *tt3=new TChain("tt");
if(sim!=1){
  //1range
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
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_19.root");
}
//2range
if(region==2){
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_5.root");
}
//3range
if(region==3){
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_4.root");
}
//4range
if(region==4){
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_2.root");
}
//5range
if(region==5){
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_1.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_2.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_3.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_4.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_5.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_1.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_2.root");
}
//6range
if(region==6){
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_1.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_2.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_3.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_4.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_5.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_1.root");
}
//7range
if(region==7){
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_2.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_3.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_4.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_1.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_2.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_3.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_1.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_2.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_3.root");
}
//8range
if(region==8){

tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_1.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_2.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_3.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_4.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_5.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_6.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_7.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_1.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_2.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_3.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_5.root");
}
//9range
if(region==9){
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_1.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_2.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_3.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_4.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_1.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_2.root");
}
//10range
if(region==10){
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_3.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_4.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_5.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_6.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_1.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_2.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_3.root");
}
//11range
if(region==11){
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_1.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_2.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_3.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_4.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_5.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_6.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_7.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_1.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_2.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_3.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_4.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_5.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_6.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_7.root");
}
//12range
if(region==12){
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_1.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_2.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_3.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_4.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_5.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_6.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_7.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_1.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_2.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_3.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_4.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_5.root");
}
//13range
if(region==13){
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_1.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_2.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_3.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_4.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_5.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_1.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_2.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_3.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_4.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_5.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_6.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_1.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_2.root");
tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_3.root");
}
}
else{
    if(region==1){
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_1range_19.root");
  }
  if(region==2){
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_2range_19.root");
  }
  if(region==3){
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_3range_19.root");
  }
  if(region==4){
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_10.root");  
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_4range_19.root");
  }
  if(region==5){
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_5range_19.root");
  }
  if(region==6){
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_6range_19.root");
  }
  if(region==7){
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_7range_19.root");
  }
  if(region==8){
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_8range_19.root");
  }  
  if(region==9){
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_9range_19.root");
  }  
  if(region==10){
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_10range_19.root");
  }  
  if(region==11){
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_11range_19.root");
  }  
  if(region==12){
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_12range_19.root");
  }  
  if(region==13){
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_13range_19.root");
  }  
  if(region==14){
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_14range_19.root");
  }  
  if(region==15){  
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_15range_19.root");
  }
  if(region==16){  
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_16range_19.root");
  }
  if(region==17){  
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_17range_19.root");
  }
  if(region==18){  
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_18range_19.root");
  }
  if(region==19){  
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_19range_19.root");
  }
  if(region==20){  
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_20range_19.root");
  }
  if(region==21){  
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_21range_19.root");
  }
  if(region==22){  
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_22range_19.root");
  }
  if(region==23){  
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_23range_19.root");
  }
  if(region==24){  
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_24range_19.root");
  }
  if(region==25){  
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_8.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_9.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_15.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/eff_cnt/sim_cosm/cosm_runs_25range_19.root");
  }         
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
 
float P, Pkaon, Ppion, mkaon=493.66, mpion=139.57, mmuon=105.65;

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

 TFile *fout=0;
 TString fname;
 fname=TString::Format("cnt_thick_%d_%d_eff_ineff_allcnt_%d.root",first_cnt,end_cnt,sim).Data();
 fout=new TFile(fname,"RECREATE");
 cout<<fname<<endl; 

   TH1F *h1=new TH1F("h1","Cosmic muon momentum",200,0,3000);
   TH1F *h2=new TH1F("h2","<<Kaons>>",500,0,50);
   TH1F *h3=new TH1F("h3","<<#pi>>",500,0,50);
   TH1F *h4=new TH1F("h4","E/p",200,0,2);
   TH1F *h5=new TH1F("h5","Cosmic muon energy",200,0,1000);
   TH1F* h6=new TH1F("h6","N_{ph.e.}",1000,0,100);
   
   char name[0];
   if(sim==1)sprintf(name,"Simulation #mu");
   if(sim!=1)sprintf(name,"Experiment #mu");
   TProfile *pr0=new TProfile(name,name,200,0,3000,0,50);

   char name1[0];
   if(sim==1)sprintf(name1,"Simulation #pi");
   if(sim!=1)sprintf(name1,"Experiment #pi");  
   TProfile *pr1=new TProfile(name1,name1,200,0,3000,0,50);

   char name2[0];
   if(sim==1)sprintf(name2,"Simulation K");
   if(sim!=1)sprintf(name2,"Experiment K");
   TProfile *pr2=new TProfile(name2,name2,200,0,3000,0,50);

   TCanvas *cc1 = new TCanvas("cc1",TString::Format("cnt_thick_%d_%d_eff_ineff_all_%d.root",first_cnt,end_cnt,sim).Data(),1600,1200);
   gStyle->SetOptStat(1);
   gStyle->SetOptFit(1);
   gROOT->SetStyle("Plain");
   cc1->Divide(4,3);

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

   float num_npenotzero_cosmmu, num_npetotal_cosmmu, eff_cosmmu, err_eff_cosmmu;
   int total_entries=0;

   int jj=0;

   Int_t nentr=tt1->GetEntries();               //число событий в дереве

//for(int k=0; k<100000; k++)                               
//for(int k=0; k<3000000; k++)                               
for(int k=0; k<nentr; k++)                               //цикл по всем событиям в дереве
{
	tt1->GetEntry(k);
        if( t.p<=1000000 && t.nvec>4 && t.chi2<20 && sqrt(pow(t.x0,2)+pow(t.y0,2)+pow(t.z0,2))<40 && t.emc_ncls>1 && tof.nhits>=1 && mu.octant>=1 && mu.status>=1 && mu.nhits>1 && mu.layer>=1 && emc.energy>50 && t.emc_ncls<=4 && tof.nhits<=8)
//        if( t.p<=1000000 && t.nvec>4 && t.chi2<20 && sqrt(pow(t.x0,2)+pow(t.y0,2)+pow(t.z0,2))<40 && t.emc_ncls>1 && tof.nhits>=1 && mu.octant>=1 && mu.status>=1 && mu.nhits>1 && mu.layer>=1 && emc.energy>50 && emc.ncls_trk>=1 && t.emc_ncls<=4 && tof.nhits<=8)
//        if( t.p<=1000000 && t.nvec>4 && t.chi2<5 && sqrt(pow(t.x0,2)+pow(t.y0,2)+pow(t.z0,2))<40 && t.emc_ncls>1 && tof.nhits>=1 && mu.octant>=1 && mu.status>=1 && mu.nhits>1 && mu.layer>=1 && emc.energy>50 && emc.ncls_trk>=1 && t.emc_ncls<=3 && tof.nhits<=4)
  {
   	int j=0;
   	int counter=-2;
   	int counter2, counter3;
   	float sum_npe=0;
	for (int n=0; n<9; n++)  //цикл по числу пересеченных счетчиков на событии
	{
		tt2->SetBranchStatus("*",0);
		sprintf(branchname,"cnt%d",n);       
		tt2->SetBranchStatus(branchname,1);
		tt2->SetBranchAddress(branchname,&cnt);       //получаем пересечение    
    		//printf("Counter %d\r",n); fflush(stdout);       	
		tt2->GetEntry(k);
		//cout<<"crossing of the counter =>"<<n<<"  Nentries="<<nentr<<endl;
	  	//cout<<"=================================================================="<<endl;
    		//cout<<"cnt.i="<<cnt.i<<"\t"<<"cnt.npe="<<cnt.npe<<endl; 

  		//for ( counter; counter<end_cnt; counter++)
  		//{

		cout<<k<<"\t"<<n<<"\t"<<counter<<"\t"<<"Counter="<<cnt.i<<endl;
	  
	  	if(cnt.i>-1){counter=cnt.i;}
	  
	  	counter2=counter+80;	
	  	counter3=counter+1;
      	  	if(counter==39){
       	  	counter2=counter+80;
          	counter3=counter-19;
          	}
          	if(counter==59){
          	counter2=counter+80;
          	counter3=counter-19;
          	}
          	if(counter==19){
          	counter2=counter+80;
          	counter3=counter-19;
          	}
          	if(counter==79){
          	counter2=counter+80;
          	counter3=counter-19;
          	}
     	  
     	  	if(counter>=80){
     	  	counter2=counter-79;	
	  	counter3=counter+1;
          		if(counter==99){
          		counter2=counter-19;
          		counter3=counter-99;
          		}
          		if(counter==119){
          		counter2=counter-19;
          		counter3=counter-99;
          		}
          		if(counter==139){
          		counter2=counter-19;
          		counter3=counter-99;
          		}
          		if(counter==159){
          		counter2=counter-19;
          		counter3=counter-99;
          		}	  
     	  	}
      		//cout<<"crossing of the counter =>"<<n<<"  Nentries="<<nentr<<"\t"<<"Cnt1="<<counter<<"\t"<<"Cnt2="<<counter2<<"\t"<<"Cnt3="<<counter3<<"\t"<<endl;  

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
 
  	  if( (cnt.i==counter || cnt.i==counter2 || cnt.i==counter3) && ATCregion!=0 && cnt.wlshit!=1 && cnt.nearwls!=1 )  
		{	
			j++;   //считаем число перечеченных счетчиков и затем используем это условие
			sum_npe+=cnt.npe;
       			//cout<<k<<"\t"<<"cnt.i="<<cnt.i<<"\t"<<"cnt.npe="<<cnt.npe<<endl; 
	     		//cout<<k<<"\t"<<"Cnt="<<cnt.i<<"\t"<<"sum_npe="<<sum_npe<<endl;
       			//cout<<"j="<<j<<endl;     
       			//cout<<"k="<<k<<"\t"<<"ncross="<<n<<"\t"<<"cnt.i="<<cnt.i<<"\t"<<"cnt.npe="<<cnt.npe<<"\t"<<"t.t="<<t.t<<"\t"<<"t.x0="<<t.x0<<"\t"<<"t.y0="<<t.y0<<"\t"<<"t.z0="<<t.z0<<endl;  
       			//cout<<"Cnt1="<<counter<<"\t"<<"Cnt2="<<counter2<<"\t"<<"Cnt3="<<counter3<<"\t"<<endl;
 
     			if( j>1 )
       			{
            				printf("k= %d\r",k); fflush(stdout);        
 					//cout<<"j>1"<<"\t"<<k<<"\t"<<"Cnt="<<cnt.i<<"\t"<<"sum_npe="<<sum_npe<<endl;
            				//cout<<"Run="<<ev.run<<"\t"<<"Event="<<ev.event<<"\t"<<"Evdaq="<<ev.evdaq<<"\t"<<endl;
  					Pkaon=t.p*(mkaon/mmuon);
  					Ppion=t.p*(mpion/mmuon);

  					h1->Fill(t.p);
  					if(t.p>200 && t.p<300)
            				{
  						h2->Fill(sum_npe);
				         	/*
						if( sum_npe>5 )
						{
					 		cout<<"Run="<<ev.run<<"\t"<<"Event="<<ev.event<<"\t"<<"Evdaq="<<ev.evdaq<<"\t"<<endl;
					 		cout<<"cnt="<<cnt.i<<"\t"<<"npe="<<cnt.npe<<"\t"<<"P="<<t.p<<"\t"<<"E/p="<<emc.energy/t.p<<"\t"<<"tof.nhits="<<tof.nhits<<endl;
					 	}
						*/
  					}
  					if(t.p>700 && t.p<1100)
            				{
  						h3->Fill(sum_npe);
				         	/*
						if( sum_npe==0 )
						{
					 		cout<<"Run="<<ev.run<<"\t"<<"Event="<<ev.event<<"\t"<<"Evdaq="<<ev.evdaq<<"\t"<<endl;
					 		cout<<"cnt="<<cnt.i<<"\t"<<"npe="<<cnt.npe<<"\t"<<"P="<<t.p<<"\t"<<"E/p="<<emc.energy/t.p<<"\t"<<"tof.nhits="<<tof.nhits<<endl;
					 	}
						*/
  					}
  					h4->Fill(emc.energy/t.p);
  					h5->Fill(emc.energy);
  					pr0->Fill(t.p,sum_npe);
                                        h6->Fill(sum_npe);
  					//cout<<k<<"\t"<<"cnt.i="<<cnt.i<<"\t"<<"cnt.npe="<<cnt.npe<<endl; 
					//cout<<k<<"\t"<<"Cnt="<<cnt.i<<"\t"<<"sum_npe="<<sum_npe<<endl;
            				if(sum_npe>0&&jj<10)
            				{  
              					cout<<t.p<<"\t"<<sum_npe<<endl; 
              					cout<<"Nentries="<<nentr<<endl;
              					jj++;
            				}
					/*if(cnt.npe>5 && t.p<100){
  					cout<<"Run="<<ev.run<<"\t"<<"Event="<<ev.event<<"\t"<<"Evdaq="<<ev.evdaq<<"\t"<<endl;
  					cout<<"cnt="<<n<<"\t"<<"npe="<<cnt.npe<<"\t"<<"P="<<t.p<<"\t"<<endl;
  					}
					*/
  					pr1->Fill(Ppion,sum_npe);
  					pr2->Fill(Pkaon,sum_npe);

             				total_entries=++total_entries;

              				if(t.p>1200)
              				{
              					if(sum_npe>npe_trh){num_npenotzero_cosmmu=++num_npenotzero_cosmmu;}
	      					num_npetotal_cosmmu=++num_npetotal_cosmmu;
    	      				}


					for(int ii=0; ii<=kk-1; ii++)
    				     	{
     						p_all[ii]=100*ii+50;
      				     	     	err_p_all[ii]=50;
        			     		if(t.p>100*ii&&t.p<100+100*ii)
                   				{ 
        	    					if(sum_npe>npe_trh){num_npenotzero_mu[ii]=++num_npenotzero_mu[ii];}
                    					num_npetotal_mu[ii]=++num_npetotal_mu[ii];
                    					eff_mu[ii]=num_npenotzero_mu[ii]/num_npetotal_mu[ii];
                    					err_eff_mu[ii]=sqrt(num_npenotzero_mu[ii])/num_npetotal_mu[ii];
        			     		}
        			     		if(Ppion>100*ii&&Ppion<100+100*ii)
                   				{ 
        	    					if(sum_npe>npe_trh){num_npenotzero_pion[ii]=++num_npenotzero_pion[ii];}
              	    					num_npetotal_pion[ii]=++num_npetotal_pion[ii];
              	    					eff_pion[ii]=num_npenotzero_pion[ii]/num_npetotal_pion[ii];
                    					err_eff_pion[ii]=sqrt(num_npenotzero_pion[ii])/num_npetotal_pion[ii];
                    					not_eff_pion[ii]=1-eff_pion[ii];
        			     		}
                   				if(Pkaon>100*ii&&Pkaon<100+100*ii)
                   				{ 
 							if(sum_npe>npe_trh){num_npenotzero_kaon[ii]=++num_npenotzero_kaon[ii];}
						        num_npetotal_kaon[ii]=++num_npetotal_kaon[ii];
						        eff_kaon[ii]=num_npenotzero_kaon[ii]/num_npetotal_kaon[ii];
						        err_eff_kaon[ii]=sqrt(num_npenotzero_kaon[ii])/num_npetotal_kaon[ii];
						        not_eff_kaon[ii]=1-eff_kaon[ii];
						}
                   				Ksigma[ii]=pow(TMath::Erf(-1+2*(1-eff_kaon[ii])),-1)+pow(TMath::Erf(-1+2*eff_pion[ii]),-1);
					}
       			}
      		}
      		if(n>0 && j>1)
      		{
       			break;
      		}           
       	    }
	}
   }
cc1->cd(1);
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

cc1->cd(2);
h5->SetXTitle("E, MeV");
h5->Draw("hist");

cc1->cd(3);
h4->Draw("hist");

cc1->cd(4);
TF1* myfit0=new TF1("myfit0","[0]+[1]*(x^2-[2]^2)/(x^2)*(TMath::Erf(x-[2])+1)/2",200.,3000.);
myfit0->SetLineColor(kGreen);
myfit0->SetParameter(0,pr0->GetMinimum());
myfit0->SetParameter(1,pr0->GetMaximum());
myfit0->SetParameter(2,200);
myfit0->SetParNames("#mu_{0}","#mu_{max}","P_{thr}");
pr0->Fit("myfit0","","",200,3000);
cout<<"myfit0->GetParameter(0)="<<myfit0->GetParameter(0)<<"\t"<<"myfit0->GetParameter(1)"<<myfit0->GetParameter(1)<<endl;
pr0->SetLineColor(kGreen);
pr0->SetMarkerStyle(21);
if(sim==1){pr0->SetMarkerStyle(24);}
pr0->SetMarkerColor(kGreen); 
pr0->SetMarkerSize(0.8); 
pr0->SetXTitle("P, MeV/c");
pr0->SetYTitle("N_{ph.e.}");
gStyle->SetOptFit(1111);
//gStyle->SetOptStat(1011);
pr0->Draw("prof");

cc1->cd(5);
gPad->SetLogy();
h2->SetLineColor(kRed);
h2->Draw("hist");

cc1->cd(6);
h3->Draw("hist");

cc1->cd(7);
TF1* myfit=new TF1("myfit1","[0]+[1]*(x^2-[2]^2)/(x^2)*(TMath::Erf(x-[2])+1)/2",200.,3000.);
myfit->SetLineColor(kBlue);
myfit->SetParameter(0,pr1->GetMinimum());
myfit->SetParameter(1,pr1->GetMaximum());
myfit->SetParameter(2,200);
myfit->SetParNames("#mu_{0}","#mu_{max}","P_{thr}");
pr1->Fit("myfit1","","",200,3000);
pr1->SetLineColor(kBlue);
pr1->SetMarkerStyle(21);
if(sim==1){pr1->SetMarkerStyle(24);}
pr1->SetMarkerColor(kBlue); 
pr1->SetMarkerSize(0.8); 
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
pr2->SetMarkerStyle(21);
if(sim==1){pr2->SetMarkerStyle(24);}
pr2->SetMarkerColor(kRed); 
pr2->SetMarkerSize(0.8); 
float Pthr=myfit->GetParameter(2);
float ref_index;
ref_index=sqrt(1+pow((105/Pthr),2));
float mu_0=myfit->GetParameter(0);
float mu_max=myfit->GetParameter(1);
//cout<<"mu_ef="<<mu_ef<<"\t"<<"eff="<<eff<<"\t"<<"ineff="<<1-eff<<endl;
pr1->Draw("prof");
pr2->Draw("same");

cc1->cd(8);
TMultiGraph *mg1 = new TMultiGraph();
//mg2->SetTitle("Efficiency&Momentum");
TGraphErrors* gr1=new TGraphErrors(kk,p_all,eff_mu,err_p_all,err_eff_mu);
gr1->SetMarkerStyle(21);
if(sim==1){gr1->SetMarkerStyle(24);}
gr1->SetMarkerSize(0.8); 
gr1->SetMarkerColor(3);
gr1->SetLineWidth(2);
gr1->SetLineColor(3);
gr1->SetTitle("Experiment #mu");
gr1->SetName("Experiment #mu");
if(sim==1)
{
gr1->SetTitle("Simulation #mu");
gr1->SetName("Simulation #mu");
}
gr1->GetXaxis()->SetTitle("P, MeV/c");
gr1->GetYaxis()->SetTitle("Efficiency");
//gr1->Draw("ap");

TGraphErrors* gr2=new TGraphErrors(kk,p_all,eff_pion,err_p_all,err_eff_pion);
gr2->SetMarkerStyle(21);
if(sim==1){gr2->SetMarkerStyle(24);}
gr2->SetMarkerSize(0.8); 
gr2->SetMarkerColor(4);
gr2->SetLineWidth(2);
gr2->SetLineColor(4);
gr2->SetTitle("Experiment #pi");
gr2->SetName("Experiment #pi");
if(sim==1)
{
gr2->SetTitle("Simulation #pi");
gr2->SetName("Simulation #pi");
}
gr2->GetXaxis()->SetTitle("P, MeV/c");
gr2->GetYaxis()->SetTitle("Efficiency");
//gr2->Draw("same");

TGraphErrors* gr3=new TGraphErrors(kk,p_all,eff_kaon,err_p_all,err_eff_kaon);
gr3->SetMarkerStyle(21);
if(sim==1){gr3->SetMarkerStyle(24);}
gr3->SetMarkerSize(0.8); 
gr3->SetMarkerColor(2);
gr3->SetLineWidth(2);
gr3->SetLineColor(2);
gr3->SetTitle("Experiment K");
gr3->SetName("Experiment K");
if(sim==1)
{
gr3->SetTitle("Simulation K");
gr3->SetName("Simulation K");
} 
mg1->Add(gr1);
mg1->Add(gr2);
mg1->Add(gr3);
mg1->SetTitle("Efficiency&Momentum; P, MeV/c; Efficiency");
// draw the multigraph
mg1->Draw("ap");

cc1->cd(9);
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

cc1->cd(10);
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

cc1->cd();
cc1->SaveAs(TString::Format("cnt_thick_%d_%d_eff_ineff_allcnt_canvas_%d.root",first_cnt,end_cnt,sim).Data());

cout<<"Cnt="<<counter<<"\t"<<"Nthr_{ph.e.}="<<npe_trh<<"\t"<<"eff_mu[11]="<<eff_mu[11]<<"\t"<<"eff_pion[11]="<<eff_pion[11]<<"\t"<<"eff_kaon[11]="<<eff_kaon[11]<<"\t"<<"sigma[11]="<<Ksigma[11]<<"\t"<<"total_entries="<<total_entries<<endl;

for(int kk4=5; kk4<=15; kk4++)
{
ofstream of(TString::Format("cnt_thick_%d_%d_npetrh%f_eff_ineff_p%d_allcnt_%d.dat",first_cnt,end_cnt,npe_trh,kk4,sim).Data(),ios::app);         
of<<counter<<"\t"<<npe_trh<<"\t"<<myfit0->GetParameter(1)<<"\t"<<eff_mu[kk4]<<"\t"<<eff_pion[kk4]<<"\t"<<eff_kaon[kk4]<<"\t"<<not_eff_pion[kk4]<<"\t"<<not_eff_kaon[kk4]<<"\t"<<Ksigma[kk4]<<"\t"<<total_entries<<endl;
of.close();
}

eff_cosmmu=num_npenotzero_cosmmu/num_npetotal_cosmmu;
err_eff_cosmmu=sqrt(num_npenotzero_cosmmu)/num_npetotal_cosmmu;
cout<<"Cosmic muons with P>1.2GeV/c ==>"<<"\t"<<"npe_trh="<<npe_trh<<"\t"<<"Eff="<<eff_cosmmu<<"\t"<<"Efferr="<<err_eff_cosmmu<<endl;
ofstream of2(TString::Format("thick_eff_%d_%d_npetrh%f_eff_ineff_allcnt_%d.dat",first_cnt,end_cnt+1,npe_trh,sim).Data(),ios::app);         
of2<<"Cosmic muons with P>1.2GeV/c ==>"<<"\t"<<"npe_trh="<<npe_trh<<"\t"<<"Eff="<<eff_cosmmu<<"\t"<<"Efferr="<<err_eff_cosmmu<<endl;
of2.close();

    pr0->Write("N_{ph.e.}&Momentum_mu");
    pr1->Write("N_{ph.e.}&Momentum_pi");
    pr2->Write("N_{ph.e.}&Momentum_K");
    gr1->Write("eff_mu");
    gr2->Write("eff_pi");
    gr3->Write("eff_K");
    mg1->Write("K efficiency and #pi Misidentification");
    mg2->Write("Efficiency&Momentum");
    h6->Write("Amplitude");

//}
//fout->Write();                                                     
//fout->Close();
}
