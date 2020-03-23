//ЧЩЮЙУМЕОЙЕ ЛПЬЖЖЙГЙЕОФБ УЧЕФПУВПТБ ДМС 5-ФЙ ФЙРПЧ УЮЕФЮЙЛПЧ Й УРТБЧБ Й УМЕЧБ ПФДЕМШОП
//СЂРёСЃРѕРІР°РЅРёРµ С‚СЂРµС…РјРµСЂРЅРѕР№ РіРёСЃС‚РѕРіСЂР°РјРјС‹ x,y,z - РѕС‚СЂРёСЃРѕРІРєР° СЃС‡РµС‚С‡РёРєР° РІ С‚СЂРµС…РјРµСЂРЅРѕРј РІРёРґРµ РІ Р»РѕРєР°Р»СЊРЅС‹С… РєРѕРѕСЂРґРёРЅР°С‚Р°С…
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
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TLine.h>
#include <TEventList.h>
#include <TProfile.h>
#include <vector>
#include <TChain.h>
#include <string>

Double_t poissonf(Double_t*x,Double_t*par)                                         
{                                                                              
  return par[0]*TMath::Poisson(x[0],par[1]);
} 

using namespace std;
string progname;

int Usage(string status)
{
   cout<<"Usage: "<<progname<<"\t"<<"Range data->1,2...  Number area->(AY=1,AS=2,AR=3,AQ=4,EX=5)  Data->(format: 13082017)"<<endl;
   exit(0);
}

int main(int argc, char* argv[])
{  
    progname=argv[0];
    int region;
    int area;
    int day_data;
    
    if( argc>1 )
    {
    region=atoi(argv[1]);
    area=atoi(argv[2]);
    day_data=atoi(argv[3]);
    if(area>5){ Usage(progname); return 0;}  
    }
    else
    { 
      Usage(progname);
    }
    

    struct data                                  //СЏв”‚СЏв”ЊСЏв”ЂСЏв”ђРїв•ЁСЏв”ЊСЏв”ЂСЏв”ђРїв•џ СЏв”‚ Рїв•ўРїв•џРїв•«Рїв•«СЏв–ЂРїв•ЄРїв•¦ Рїв–‘Рїв•–Рїв•‘
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
    { int t,q,ip,nvec,nvecxy,nvecz,nhits,nhitsxy,nhitsz,nhitsvd;
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
      int event,evdaq,run,quality;
    };
data8 ev;

char branchname[161];
char branchname2[1];
char branchname3[1];
char branchname4[1];
char branchname5[1];
char branchname6[1];
char branchname7[1];
char branchname8[1];

float LC, total_LC, average_LC, er_ev,disp, err_LC, phi1, phi2, z1, z2, x1, x2, y1, y2;
int i, j, Ev, AS=1, AQ=1, AR=1, AY=1, EX=1;
int numbin;
int ev_all=0;

TChain *tt=new TChain("tt");
//1range
if(region==1){
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_1.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_2.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_3.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_4.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_5.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_6.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_7.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_8.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_9.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_9_1.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_10.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_11.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_12.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_13.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_14.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_15.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_16.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_17.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_18.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_19.root");
}
//2range
if(region==2){
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_1.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_2.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_3.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_4.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_5.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_6.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_7.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_1.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_2.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_3.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_4.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_5.root");
}
//3range
if(region==3){
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_6.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_7.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_1.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_2.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_4.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_5.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan15_1.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan15_2.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_1.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_2.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_3.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_4.root");
}
//4range
if(region==4){
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr15_1.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr15_2.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may15_1.root");
  tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may15_2.root");
}
//5range
if(region==5){
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_1.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_2.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_3.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_4.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_5.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jul15_1.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jul15_2.root");
}
//6range
if(region==6){
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_1.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_2.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_3.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_4.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_5.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_1.root");
}
//7range
if(region==7){
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_2.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_3.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_4.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec15_1.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec15_2.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec15_3.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan16_1.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan16_2.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan16_3.root");
}
//8range
if(region==8){

tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_1.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_2.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_3.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_4.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_5.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_6.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_7.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_1.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_2.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_3.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_5.root");
}
//9range
if(region==9){
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_1.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_2.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_3.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_4.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_1.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_2.root");
}
//10range
if(region==10){
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_3.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_4.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_5.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_6.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun16_1.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun16_2.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun16_3.root");
}
//11range
if(region==11){
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march17_1.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march17_2.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march17_3.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr17_1.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr17_2.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr17_3.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may17_1.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may17_2.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may17_3.root");
}
//12range
if(region==12){
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_1.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_2.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_3.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_4.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_5.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_6.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_7.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_1.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_2.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_3.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_4.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_5.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_6.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_7.root");
}
//13range
if(region==13){
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_1.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_2.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_3.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_4.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_5.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_6.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_7.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_1.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_2.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_3.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_4.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_5.root");
}
//14range
if(region==14){
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_1.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_2.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_3.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_4.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_5.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_1.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_2.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_3.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_4.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_5.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_6.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun18_1.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun18_2.root");
tt->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun18_3.root");
}


TChain *tt2=new TChain("tt");
//1range
if(region==1){
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_1.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_2.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_3.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_4.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_5.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_6.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_7.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_8.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_9.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_9_1.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_10.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_11.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_12.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_13.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_14.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_15.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_16.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_17.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_18.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_19.root");
}
//2range
if(region==2){
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_1.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_2.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_3.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_4.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_5.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_6.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_7.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_1.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_2.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_3.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_4.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_5.root");
}
//3range
if(region==3){
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_6.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_7.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_1.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_2.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_4.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_5.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan15_1.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan15_2.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_1.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_2.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_3.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_4.root");
}
//4range
if(region==4){
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr15_1.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr15_2.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may15_1.root");
  tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may15_2.root");
}
//5range
if(region==5){
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_1.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_2.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_3.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_4.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_5.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jul15_1.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jul15_2.root");
}
//6range
if(region==6){
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_1.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_2.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_3.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_4.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_5.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_1.root");
}
//7range
if(region==7){
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_2.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_3.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_4.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec15_1.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec15_2.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec15_3.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan16_1.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan16_2.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan16_3.root");
}
//8range
if(region==8){

tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_1.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_2.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_3.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_4.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_5.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_6.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_7.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_1.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_2.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_3.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_5.root");
}
//9range
if(region==9){
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_1.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_2.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_3.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_4.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_1.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_2.root");
}
//10range
if(region==10){
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_3.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_4.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_5.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_6.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun16_1.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun16_2.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun16_3.root");
}
//11range
if(region==11){
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march17_1.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march17_2.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march17_3.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr17_1.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr17_2.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr17_3.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may17_1.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may17_2.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may17_3.root");
}
//12range
if(region==12){
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_1.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_2.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_3.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_4.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_5.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_6.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_7.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_1.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_2.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_3.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_4.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_5.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_6.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_7.root");
}
//13range
if(region==13){
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_1.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_2.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_3.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_4.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_5.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_6.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_7.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_1.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_2.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_3.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_4.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_5.root");
}
//14range
if(region==14){
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_1.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_2.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_3.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_4.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_5.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_1.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_2.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_3.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_4.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_5.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_6.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun18_1.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun18_2.root");
tt2->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun18_3.root");
}

TChain *tt3=new TChain("tt");
//1range
if(region==1){
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_1.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_2.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_3.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_4.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_5.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_6.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_7.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_8.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_9.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_9_1.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_10.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_11.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_12.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_13.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_14.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_15.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_16.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_17.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_18.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_19.root");
}
//2range
if(region==2){
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_1.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_2.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_3.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_4.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_5.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_6.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_7.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_1.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_2.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_3.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_4.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_5.root");
}
//3range
if(region==3){
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_6.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_7.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_1.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_2.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_4.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_5.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan15_1.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan15_2.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_1.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_2.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_3.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_4.root");
}
//4range
if(region==4){
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr15_1.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr15_2.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may15_1.root");
  tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may15_2.root");
}
//5range
if(region==5){
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_1.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_2.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_3.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_4.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_5.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jul15_1.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jul15_2.root");
}
//6range
if(region==6){
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_1.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_2.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_3.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_4.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_5.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_1.root");
}
//7range
if(region==7){
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_2.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_3.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_4.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec15_1.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec15_2.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec15_3.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan16_1.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan16_2.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan16_3.root");
}
//8range
if(region==8){

tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_1.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_2.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_3.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_4.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_5.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_6.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_7.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_1.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_2.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_3.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_5.root");
}
//9range
if(region==9){
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_1.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_2.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_3.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_4.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_1.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_2.root");
}
//10range
if(region==10){
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_3.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_4.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_5.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_6.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun16_1.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun16_2.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun16_3.root");
}
//11range
if(region==11){
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march17_1.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march17_2.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march17_3.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr17_1.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr17_2.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr17_3.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may17_1.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may17_2.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may17_3.root");
}
//12range
if(region==12){
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_1.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_2.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_3.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_4.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_5.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_6.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_7.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_1.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_2.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_3.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_4.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_5.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_6.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_7.root");
}
//13range
if(region==13){
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_1.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_2.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_3.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_4.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_5.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_6.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_7.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_1.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_2.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_3.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_4.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_5.root");
}
//14range
if(region==14){
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_1.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_2.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_3.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_4.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_5.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_1.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_2.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_3.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_4.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_5.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_6.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun18_1.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun18_2.root");
tt3->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun18_3.root");
}

TChain *tt4=new TChain("tt");
//1range
if(region==1){
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_1.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_2.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_3.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_4.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_5.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_6.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_7.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_8.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_9.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_9_1.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_10.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_11.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_12.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_13.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_14.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_15.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_16.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_17.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_18.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_19.root");
}
//2range
if(region==2){
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_1.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_2.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_3.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_4.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_5.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_6.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_7.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_1.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_2.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_3.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_4.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_5.root");
}
//3range
if(region==3){
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_6.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_7.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_1.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_2.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_4.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_5.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan15_1.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan15_2.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_1.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_2.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_3.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_4.root");
}
//4range
if(region==4){
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr15_1.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr15_2.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may15_1.root");
  tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may15_2.root");
}
//5range
if(region==5){
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_1.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_2.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_3.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_4.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_5.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jul15_1.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jul15_2.root");
}
//6range
if(region==6){
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_1.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_2.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_3.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_4.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_5.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_1.root");
}
//7range
if(region==7){
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_2.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_3.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_4.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec15_1.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec15_2.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec15_3.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan16_1.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan16_2.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan16_3.root");
}
//8range
if(region==8){

tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_1.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_2.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_3.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_4.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_5.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_6.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_7.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_1.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_2.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_3.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_5.root");
}
//9range
if(region==9){
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_1.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_2.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_3.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_4.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_1.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_2.root");
}
//10range
if(region==10){
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_3.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_4.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_5.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_6.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun16_1.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun16_2.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun16_3.root");
}
//11range
if(region==11){
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march17_1.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march17_2.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march17_3.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr17_1.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr17_2.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr17_3.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may17_1.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may17_2.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may17_3.root");
}
//12range
if(region==12){
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_1.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_2.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_3.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_4.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_5.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_6.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_7.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_1.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_2.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_3.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_4.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_5.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_6.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_7.root");
}
//13range
if(region==13){
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_1.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_2.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_3.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_4.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_5.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_6.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_7.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_1.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_2.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_3.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_4.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_5.root");
}
//14range
if(region==14){
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_1.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_2.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_3.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_4.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_5.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_1.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_2.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_3.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_4.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_5.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_6.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun18_1.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun18_2.root");
tt4->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun18_3.root");
}


TChain *tt5=new TChain("emct");
//1range
if(region==1){
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_1.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_2.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_3.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_4.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_5.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_6.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_7.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_8.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_9.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_9_1.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_10.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_11.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_12.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_13.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_14.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_15.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_16.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_17.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_18.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_19.root");
}
//2range
if(region==2){
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_1.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_2.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_3.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_4.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_5.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_6.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_7.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_1.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_2.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_3.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_4.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_5.root");
}
//3range
if(region==3){
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_6.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_7.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_1.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_2.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_4.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_5.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan15_1.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan15_2.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_1.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_2.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_3.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_4.root");
}
//4range
if(region==4){
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr15_1.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr15_2.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may15_1.root");
  tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may15_2.root");
}
//5range
if(region==5){
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_1.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_2.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_3.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_4.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_5.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jul15_1.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jul15_2.root");
}
//6range
if(region==6){
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_1.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_2.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_3.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_4.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_5.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_1.root");
}
//7range
if(region==7){
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_2.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_3.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_4.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec15_1.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec15_2.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec15_3.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan16_1.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan16_2.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan16_3.root");
}
//8range
if(region==8){

tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_1.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_2.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_3.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_4.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_5.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_6.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_7.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_1.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_2.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_3.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_5.root");
}
//9range
if(region==9){
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_1.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_2.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_3.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_4.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_1.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_2.root");
}
//10range
if(region==10){
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_3.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_4.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_5.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_6.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun16_1.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun16_2.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun16_3.root");
}
//11range
if(region==11){
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march17_1.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march17_2.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march17_3.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr17_1.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr17_2.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr17_3.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may17_1.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may17_2.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may17_3.root");
}
//12range
if(region==12){
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_1.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_2.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_3.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_4.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_5.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_6.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_7.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_1.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_2.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_3.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_4.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_5.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_6.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_7.root");
}
//13range
if(region==13){
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_1.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_2.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_3.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_4.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_5.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_6.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_7.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_1.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_2.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_3.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_4.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_5.root");
}
//14range
if(region==14){
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_1.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_2.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_3.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_4.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_5.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_1.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_2.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_3.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_4.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_5.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_6.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun18_1.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun18_2.root");
tt5->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun18_3.root");
}

TChain *tt6=new TChain("emct");
//1range
if(region==1){
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_1.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_2.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_3.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_4.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_5.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_6.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_7.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_8.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_9.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_9_1.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_10.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_11.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_12.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_13.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_14.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_15.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_16.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_17.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_18.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_19.root");
}
//2range
if(region==2){
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_1.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_2.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_3.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_4.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_5.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_6.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_7.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_1.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_2.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_3.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_4.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_5.root");
}
//3range
if(region==3){
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_6.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_7.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_1.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_2.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_4.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_5.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan15_1.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan15_2.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_1.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_2.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_3.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_4.root");
}
//4range
if(region==4){
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr15_1.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr15_2.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may15_1.root");
  tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may15_2.root");
}
//5range
if(region==5){
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_1.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_2.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_3.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_4.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_5.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jul15_1.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jul15_2.root");
}
//6range
if(region==6){
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_1.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_2.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_3.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_4.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_5.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_1.root");
}
//7range
if(region==7){
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_2.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_3.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_4.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec15_1.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec15_2.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec15_3.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan16_1.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan16_2.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan16_3.root");
}
//8range
if(region==8){

tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_1.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_2.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_3.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_4.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_5.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_6.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_7.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_1.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_2.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_3.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_5.root");
}
//9range
if(region==9){
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_1.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_2.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_3.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_4.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_1.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_2.root");
}
//10range
if(region==10){
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_3.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_4.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_5.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_6.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun16_1.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun16_2.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun16_3.root");
}
//11range
if(region==11){
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march17_1.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march17_2.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march17_3.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr17_1.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr17_2.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr17_3.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may17_1.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may17_2.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may17_3.root");
}
//12range
if(region==12){
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_1.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_2.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_3.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_4.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_5.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_6.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_7.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_1.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_2.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_3.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_4.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_5.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_6.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_7.root");
}
//13range
if(region==13){
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_1.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_2.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_3.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_4.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_5.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_6.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_7.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_1.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_2.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_3.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_4.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_5.root");
}
//14range
if(region==14){
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_1.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_2.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_3.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_4.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_5.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_1.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_2.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_3.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_4.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_5.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_6.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun18_1.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun18_2.root");
tt6->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun18_3.root");
}


TChain *tt7=new TChain("tt");
//1range
//1range
if(region==1){
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_1.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_2.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_3.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_4.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_5.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_6.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_7.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_8.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_9.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_9_1.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_10.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_11.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_12.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_13.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_14.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_15.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_16.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_17.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_18.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun14_19.root");
}
//2range
if(region==2){
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_1.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_2.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_3.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_4.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_5.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_6.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct14_7.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_1.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_2.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_3.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_4.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_5.root");
}
//3range
if(region==3){
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_6.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov14_7.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_1.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_2.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_4.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec14_5.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan15_1.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan15_2.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_1.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_2.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_3.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb15_4.root");
}
//4range
if(region==4){
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr15_1.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr15_2.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may15_1.root");
  tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may15_2.root");
}
//5range
if(region==5){
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_1.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_2.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_3.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_4.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun15_5.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jul15_1.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jul15_2.root");
}
//6range
if(region==6){
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_1.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_2.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_3.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_4.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_oct15_5.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_1.root");
}
//7range
if(region==7){
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_2.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_3.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_nov15_4.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec15_1.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec15_2.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec15_3.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan16_1.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan16_2.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan16_3.root");
}
//8range
if(region==8){

tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_1.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_2.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_3.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_4.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_5.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_6.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb16_7.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_1.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_2.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_3.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march16_5.root");
}
//9range
if(region==9){
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_1.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_2.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_3.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr16_4.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_1.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_2.root");
}
//10range
if(region==10){
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_3.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_4.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_5.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may16_6.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun16_1.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun16_2.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun16_3.root");
}
//11range
if(region==11){
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march17_1.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march17_2.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march17_3.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr17_1.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr17_2.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr17_3.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may17_1.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may17_2.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may17_3.root");
}
//12range
if(region==12){
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_1.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_2.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_3.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_4.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_5.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_6.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_dec17_7.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_1.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_2.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_3.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_4.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_5.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_6.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jan18_7.root");
}
//13range
if(region==13){
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_1.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_2.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_3.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_4.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_5.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_6.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_feb18_7.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_1.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_2.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_3.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_4.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_march18_5.root");
}
//14range
if(region==14){
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_1.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_2.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_3.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_4.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_apr18_5.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_1.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_2.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_3.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_4.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_5.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_may18_6.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun18_1.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun18_2.root");
tt7->Add("/home/ovtin/development/atc_npe/data_cosm/cosm_runs_jun18_3.root");
}

TCanvas *c11 = new TCanvas();
c11->cd();
TH2F *h31 = new TH2F("h2", "h2 title", 50, -50, 50, 40, -0.18,0.18);

	Int_t nentr=tt->GetEntries();               //СЏв”¤Рїв•¦СЏв”‚Рїв•©Рїв•¬ СЏв”‚Рїв•¬Рїв• СЏв–ЂСЏв”ЊРїв•¦Рїв•§ Рїв•Ў Рїв•ўРїв•ЈСЏв”ЂРїв•ЈРїв•ЎРїв•Ј
	cout<<"Nentries="<<nentr<<endl;

	sprintf(branchname2,"t");
	tt2->SetBranchStatus("*",0);
	tt2->SetBranchStatus(branchname2,1);
	tt2->SetBranchAddress(branchname2,&t);              //СЏв”ђСЏв”‚СЏв”ЊРїв•џРїв•«Рїв•џРїв•ЎРїв•©Рїв•¦Рїв•ЎРїв•џРїв•ЈРїв•Є Рїв•џРїв•ўСЏв”ЂРїв•ЈСЏв”‚СЏв”‚ Рїв•ЎРїв•ЈСЏв”ЊРїв•ЁРїв•¦

	sprintf(branchname3,"tof");
	tt3->SetBranchStatus("*",0);
	tt3->SetBranchStatus(branchname3,1);
	tt3->SetBranchAddress(branchname3,&tof); 

	sprintf(branchname4,"emc");
	tt4->SetBranchStatus("*",0);
	tt4->SetBranchStatus(branchname4,1);
	tt4->SetBranchAddress(branchname4,&emc);  

	sprintf(branchname5,"c0");
	tt5->SetBranchStatus("*",0);
	tt5->SetBranchStatus(branchname5,1);
	tt5->SetBranchAddress(branchname5,&c0);  

	sprintf(branchname6,"c1");
	tt6->SetBranchStatus("*",0);
	tt6->SetBranchStatus(branchname6,1);
	tt6->SetBranchAddress(branchname6,&c1); 

	sprintf(branchname7,"mu");
	tt7->SetBranchStatus("*",0);
	tt7->SetBranchStatus(branchname7,1);
	tt7->SetBranchAddress(branchname7,&mu); 

        cout<<"\t"<<"Region="<<region<<"\t"<<"Area="<<area<<"\t"<<"Data="<<setfill('0')<<setw(8)<<day_data<<endl;

   
   int name_area, size_area;
   int cnt_begin, cnt_end;
   float phi_begin, phi_end, z_begin, z_end, phi_area1, phi_area2, z_area1, z_area2;
   int ATCregion;
   float betta;

   string NAME_str;

   if (area==1){ cnt_begin=20; cnt_end=40; size_area=59; NAME_str=TString::Format("data_AY_%08d_%drange.dat",day_data,region).Data();}  //20-40cnt - AY
   if (area==2){ cnt_begin=40; cnt_end=60; size_area=69; NAME_str=TString::Format("data_AS_%08d_%drange.dat",day_data,region).Data();}  //40-60cnt - AS
   if (area==3){ cnt_begin=100; cnt_end=120; size_area=69; NAME_str=TString::Format("data_AR_%08d_%drange.dat",day_data,region).Data();}  //100-120cnt - AR
   if (area==4){ cnt_begin=120; cnt_end=140; size_area=59; NAME_str=TString::Format("data_AQ_%08d_%drange.dat",day_data,region).Data();}  //120-140cnt - AQ
   if (area==5){ cnt_begin=0; cnt_end=160; size_area=29; NAME_str=TString::Format("data_EX_%08d_%drange.dat",day_data,region).Data();}  //EX

for(name_area=1; name_area<size_area; name_area++)     
{
    LC=0, total_LC=0, average_LC=0, er_ev=0, i=0, numbin=0, j=0, Ev=0, disp=0, err_LC=0;
    
    ofstream of(TString::Format("data_cnt_Klight_area_all_%dregion_%darea.dat",region,area).Data(),ios_base::out);   //очистка файла для последующей записи     
    of.close();

//======================================================BILS,BILD (40-60)============================================================================================
if(area==2){
    if (name_area==1) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-82.709)*3.14/180.; z1=-23.128824+0.40-4.37139; z2=-23.128824+0.40+4.37139; }
    if (name_area==2) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-82.709)*3.14/180.; z1=-23.128824+1*8.74277+0.40-4.37139; z2=-23.128824+1*8.74277+0.40+4.37139;}   
    if (name_area==3) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-82.709)*3.14/180.; z1=-23.128824+2*8.74277+0.40-4.37139; z2=-23.128824+2*8.74277+0.40+4.37139;}  
    if (name_area==4) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-82.709)*3.14/180.; z1=-23.128824+3*8.74277+0.40-4.37139; z2=-23.128824+3*8.74277+0.40+4.37139;}  
    if (name_area==5) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-82.709)*3.14/180.; z1=-23.128824+4*8.74277+0.40-4.37139; z2=-23.128824+4*8.74277+0.40+4.37139;}   
    if (name_area==6) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-82.709)*3.14/180.; z1=-23.128824+5*8.74277+0.40-4.37139; z2=-23.128824+5*8.74277+0.40+4.37139;} 
    if (name_area==7) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-85.6045)*3.14/180.; z1=-23.128824+0.40+0*8.74277-4.37139; z2=-23.128824+0.40+0*8.74277+4.37139;} 
    if (name_area==8) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-85.6045)*3.14/180.; z1=-23.128824+0.40+1*8.74277-4.37139; z2=-23.128824+0.40+1*8.74277+4.37139;}
    if (name_area==9) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-85.6045)*3.14/180.; z1=-23.128824+0.40+2*8.74277-4.37139; z2=-23.128824+0.40+2*8.74277+4.37139;} 
    if (name_area==10) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-85.6045)*3.14/180.; z1=-23.128824+0.40+3*8.74277-4.37139; z2=-23.128824+0.40+3*8.74277+4.37139;}
    if (name_area==11) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-85.6045)*3.14/180.; z1=-23.128824+0.40+4*8.74277-4.37139; z2=-23.128824+0.40+4*8.74277+4.37139;} 
    if (name_area==12) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-85.6045)*3.14/180.; z1=-23.128824+0.40+5*8.74277-4.37139; z2=-23.128824+0.40+5*8.74277+4.37139;} 
    if (name_area==13) {phi1=(90.-85.6045)*3.14/180.; phi2=(90.-87.800)*3.14/180.; z1=-23.128824+0.40+0*8.74277-4.37139; z2=-23.128824+0.40+0*8.74277+4.37139;} 
    if (name_area==14) {phi1=(90.-85.6045)*3.14/180.; phi2=(90.-87.800)*3.14/180.; z1=-23.128824+0.40+1*8.74277-4.37139; z2=-23.128824+0.40+1*8.74277+4.37139;} 
    if (name_area==15) {phi1=(90.-85.6045)*3.14/180.; phi2=(90.-87.800)*3.14/180.; z1=-23.128824+0.40+2*8.74277-4.37139; z2=-23.128824+0.40+2*8.74277+4.37139;} 
    if (name_area==16) {phi1=(90.-85.6045)*3.14/180.; phi2=(90.-87.800)*3.14/180.; z1=-23.128824+0.40+3*8.74277-4.37139; z2=-23.128824+0.40+3*8.74277+4.37139;}  
    if (name_area==17) {phi1=(90.-85.6045)*3.14/180.; phi2=(90.-87.800)*3.14/180.; z1=-23.128824+0.40+4*8.74277-4.37139; z2=-23.128824+0.40+4*8.74277+4.37139;} 
    if (name_area==18) {phi1=(90.-85.6045)*3.14/180.; phi2=(90.-87.800)*3.14/180.; z1=-23.128824+0.40+5*8.74277-4.37139; z2=-23.128824+0.40+5*8.74277+4.37139;} 
    if (name_area==19) {phi1=(90.-87.800)*3.14/180.; phi2=(90.-89.100)*3.14/180.; z1=-23.128824+0.40+0*8.74277-4.37139;  z2=-23.128824+0.40+0*8.74277+4.37139;}
    if (name_area==20) {phi1=(90.-87.800)*3.14/180.; phi2=(90.-89.100)*3.14/180.; z1=-23.128824+0.40+1*8.74277-4.37139;  z2=-23.128824+0.40+1*8.74277+4.37139;}  
    if (name_area==21) {phi1=(90.-87.800)*3.14/180.; phi2=(90.-89.100)*3.14/180.; z1=-23.128824+0.40+2*8.74277-4.37139;  z2=-23.128824+0.40+2*8.74277+4.37139;} 
    if (name_area==22) {phi1=(90.-87.800)*3.14/180.; phi2=(90.-89.100)*3.14/180.; z1=-23.128824+0.40+3*8.74277-4.37139;  z2=-23.128824+0.40+3*8.74277+4.37139;} 
    if (name_area==23) {phi1=(90.-87.800)*3.14/180.; phi2=(90.-89.100)*3.14/180.; z1=-23.128824+0.40+4*8.74277-4.37139;  z2=-23.128824+0.40+4*8.74277+4.37139;}
    if (name_area==24) {phi1=(90.-87.800)*3.14/180.; phi2=(90.-89.100)*3.14/180.; z1=-23.128824+0.40+5*8.74277-4.37139;  z2=-23.128824+0.40+5*8.74277+4.37139;} 
    if (name_area==25) {phi1=(90.-89.100)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-23.128824+0.40+0*8.74277-4.37139;  z2=-23.128824+0.40+0*8.74277+4.37139;} 
    if (name_area==26) {phi1=(90.-89.100)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-23.128824+0.40+1*8.74277-4.37139;  z2=-23.128824+0.40+1*8.74277+4.37139;} 
    if (name_area==27) {phi1=(90.-89.100)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-23.128824+0.40+2*8.74277-4.37139;  z2=-23.128824+0.40+2*8.74277+4.37139;} 
    if (name_area==28) {phi1=(90.-89.100)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-23.128824+0.40+3*8.74277-4.37139;  z2=-23.128824+0.40+3*8.74277+4.37139;} 
    if (name_area==29) {phi1=(90.-89.100)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-23.128824+0.40+4*8.74277-4.37139;  z2=-23.128824+0.40+4*8.74277+4.37139;} 
    if (name_area==30) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-90.900)*3.14/180.; z1=-23.128824+0.40+0*8.74277-4.37139;  z2=-23.128824+0.40+0*8.74277+4.37139;} 
    if (name_area==31) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-90.900)*3.14/180.; z1=-23.128824+0.40+1*8.74277-4.37139;  z2=-23.128824+0.40+1*8.74277+4.37139;} 
    if (name_area==32) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-90.900)*3.14/180.; z1=-23.128824+0.40+2*8.74277-4.37139;  z2=-23.128824+0.40+2*8.74277+4.37139;}  
    if (name_area==33) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-90.900)*3.14/180.; z1=-23.128824+0.40+3*8.74277-4.37139;  z2=-23.128824+0.40+3*8.74277+4.37139;} 
    if (name_area==34) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-90.900)*3.14/180.; z1=-23.128824+0.40+4*8.74277-4.37139;  z2=-23.128824+0.40+4*8.74277+4.37139;} 
    if (name_area==35) {phi1=(90.-90.900)*3.14/180.; phi2=(90.-92.200)*3.14/180.; z1=-23.128824+0.40+0*8.74277-4.37139;  z2=-23.128824+0.40+0*8.74277+4.37139;}  
    if (name_area==36) {phi1=(90.-90.900)*3.14/180.; phi2=(90.-92.200)*3.14/180.; z1=-23.128824+0.40+1*8.74277-4.37139;  z2=-23.128824+0.40+1*8.74277+4.37139;}  
    if (name_area==37) {phi1=(90.-90.900)*3.14/180.; phi2=(90.-92.200)*3.14/180.; z1=-23.128824+0.40+2*8.74277-4.37139;  z2=-23.128824+0.40+2*8.74277+4.37139;} 
    if (name_area==38) {phi1=(90.-90.900)*3.14/180.; phi2=(90.-92.200)*3.14/180.; z1=-23.128824+0.40+3*8.74277-4.37139;  z2=-23.128824+0.40+3*8.74277+4.37139;} 
    if (name_area==39) {phi1=(90.-90.900)*3.14/180.; phi2=(90.-92.200)*3.14/180.; z1=-23.128824+0.40+4*8.74277-4.37139;  z2=-23.128824+0.40+4*8.74277+4.37139;} 
    if (name_area==40) {phi1=(90.-90.900)*3.14/180.; phi2=(90.-92.200)*3.14/180.; z1=-23.128824+0.40+5*8.74277-4.37139;  z2=-23.128824+0.40+5*8.74277+4.37139;} 
    if (name_area==41) {phi1=(90.-92.200)*3.14/180.; phi2=(90.-94.39549)*3.14/180.; z1=-23.128824+0.40+0*8.74277-4.37139; z2=-23.128824+0.40+0*8.74277+4.37139;} 
    if (name_area==42) {phi1=(90.-92.200)*3.14/180.; phi2=(90.-94.39549)*3.14/180.; z1=-23.128824+0.40+1*8.74277-4.37139; z2=-23.128824+0.40+1*8.74277+4.37139;} 
    if (name_area==43) {phi1=(90.-92.200)*3.14/180.; phi2=(90.-94.39549)*3.14/180.; z1=-23.128824+0.40+2*8.74277-4.37139; z2=-23.128824+0.40+2*8.74277+4.37139;}
    if (name_area==44) {phi1=(90.-92.200)*3.14/180.; phi2=(90.-94.39549)*3.14/180.; z1=-23.128824+0.40+3*8.74277-4.37139; z2=-23.128824+0.40+3*8.74277+4.37139;} 
    if (name_area==45) {phi1=(90.-92.200)*3.14/180.; phi2=(90.-94.39549)*3.14/180.; z1=-23.128824+0.40+4*8.74277-4.37139; z2=-23.128824+0.40+4*8.74277+4.37139;} 
    if (name_area==46) {phi1=(90.-92.200)*3.14/180.; phi2=(90.-94.39549)*3.14/180.; z1=-23.128824+0.40+5*8.74277-4.37139; z2=-23.128824+0.40+5*8.74277+4.37139;} 
    if (name_area==47) {phi1=(90.-94.39549)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-23.128824+0.40+0*8.74277-4.37139; z2=-23.128824+0.40+0*8.74277+4.37139;} 
    if (name_area==48) {phi1=(90.-94.39549)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-23.128824+0.40+1*8.74277-4.37139; z2=-23.128824+0.40+1*8.74277+4.37139;} 
    if (name_area==49) {phi1=(90.-94.39549)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-23.128824+0.40+2*8.74277-4.37139; z2=-23.128824+0.40+2*8.74277+4.37139;} 
    if (name_area==50) {phi1=(90.-94.39549)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-23.128824+0.40+3*8.74277-4.37139; z2=-23.128824+0.40+3*8.74277+4.37139;} 
    if (name_area==51) {phi1=(90.-94.39549)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-23.128824+0.40+4*8.74277-4.37139; z2=-23.128824+0.40+4*8.74277+4.37139;}   
    if (name_area==52) {phi1=(90.-94.39549)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-23.128824+0.40+5*8.74277-4.37139; z2=-23.128824+0.40+5*8.74277+4.37139;}  
    if (name_area==53) {phi1=(90.-97.291)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-23.128824+0.40+0*8.74277-4.37139; z2=-23.128824+0.40+0*8.74277+4.37139;}   
    if (name_area==54) {phi1=(90.-97.291)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-23.128824+0.40+1*8.74277-4.37139; z2=-23.128824+0.40+1*8.74277+4.37139;}  
    if (name_area==55) {phi1=(90.-97.291)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-23.128824+0.40+2*8.74277-4.37139; z2=-23.128824+0.40+2*8.74277+4.37139;} 
    if (name_area==56) {phi1=(90.-97.291)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-23.128824+0.40+3*8.74277-4.37139; z2=-23.128824+0.40+3*8.74277+4.37139;}  
    if (name_area==57) {phi1=(90.-97.291)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-23.128824+0.40+4*8.74277-4.37139; z2=-23.128824+0.40+4*8.74277+4.37139;} 
    if (name_area==58) {phi1=(90.-97.291)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-23.128824+0.40+5*8.74277-4.37139; z2=-23.128824+0.40+5*8.74277+4.37139;} 
    if (name_area==59) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-82.709)*3.14/180.; z1=-23.128824+0.40+5*8.74277+4.37139; z2=-23.128824+0.40+5*8.74277+4.37139+4.654-0.55;}
    if (name_area==60) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-86.212)*3.14/180.; z1=-23.128824+0.40+5*8.74277+4.37139; z2=-23.128824+0.40+5*8.74277+4.37139+(4.654-0.55)*2/3;}
    if (name_area==61) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-86.212)*3.14/180.; z1=-23.128824+0.40+5*8.74277+4.37139+(4.654-0.55)*2/3; z2=-23.128824+0.40+5*8.74277+4.37139+4.654-0.55;} 
    if (name_area==62) {phi1=(90.-93.788)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-23.128824+0.40+5*8.74277+4.37139; z2=-23.128824+0.40+5*8.74277+4.37139+(4.654-0.55)*2/3;}   
    if (name_area==63) {phi1=(90.-93.788)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-23.128824+0.40+5*8.74277+4.37139+(4.654-0.55)*2/3; z2=-23.128824+0.40+5*8.74277+4.37139+4.654-0.55;}  
    if (name_area==64) {phi1=(90.-97.291)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-23.128824+0.40+5*8.74277+4.37139; z2=-23.128824+0.40+5*8.74277+4.37139+4.654-0.55;}         
    if (name_area==65) {phi1=(90.-89.100)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-23.128824+0.40+4*8.74277+4.37139; z2=-23.128824+0.40+4*8.74277+4.37139+(8.74277-3.0);}     
    if (name_area==66) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-90.900)*3.14/180.; z1=-23.128824+0.40+4*8.74277+4.37139; z2=-23.128824+0.40+4*8.74277+4.37139+(8.74277-3.0);}     
    if (name_area==67) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-23.128824+0.40-4.37139-(9.1229/4+0.12); z2=-23.128824+0.40-4.37139;}                       
    if (name_area==68) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-23.128824+0.40-4.37139-(9.1229/4+0.12); z2=-23.128824+0.40-4.37139;}                       //AS68
}
//================================================================================================================================================================================


//======================================================BOSS,BOSD (121-140)=======================================================================================================
if(area==4){
//   cout<<"check area"<<area<<endl;
    if (name_area==1) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-82.709)*3.14/180.; z1=-17.42032+0.38-4.0220909;  z2=-17.42032+0.38+4.0220909;}     //name_area1
    if (name_area==2) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-82.709)*3.14/180.; z1=-17.42032+0.38-4.0220909+1*8.044182; z2=-17.42032+0.38+4.0220909+1*8.044182;}    //name_area5
    if (name_area==3) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-82.709)*3.14/180.; z1=-17.42032+0.38-4.0220909+2*8.044182; z2=-17.42032+0.38+4.0220909+2*8.044182;}    //name_area5
    if (name_area==4) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-82.709)*3.14/180.; z1=-17.42032+0.38-4.0220909+3*8.044182; z2=-17.42032+0.38+4.0220909+3*8.044182;}    //name_area5
    if (name_area==5) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-82.709)*3.14/180.; z1=-17.42032+0.38-4.0220909+4*8.044182; z2=-17.42032+0.38+4.0220909+4*8.044182;}    //name_area5
    if (name_area==6) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-85.6045)*3.14/180.; z1=-17.42032+0.38-4.0220909+0*8.044182; z2=-17.42032+0.38+4.0220909+0*8.044182;}     //name_area10
    if (name_area==7) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-85.6045)*3.14/180.; z1=-17.42032+0.38-4.0220909+1*8.044182; z2=-17.42032+0.38+4.0220909+1*8.044182;}     //name_area10
    if (name_area==8) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-85.6045)*3.14/180.; z1=-17.42032+0.38-4.0220909+2*8.044182; z2=-17.42032+0.38+4.0220909+2*8.044182;}     //name_area10
    if (name_area==9) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-85.6045)*3.14/180.; z1=-17.42032+0.38-4.0220909+3*8.044182; z2=-17.42032+0.38+4.0220909+3*8.044182;}     //name_area10
    if (name_area==10) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-85.6045)*3.14/180.; z1=-17.42032+0.38-4.0220909+4*8.044182; z2=-17.42032+0.38+4.0220909+4*8.044182;}     //name_area10
    if (name_area==11) {phi1=(90.-85.6045)*3.14/180.; phi2=(90.-87.800)*3.14/180.; z1=-17.42032+0.38-4.0220909+0*8.044182; z2=-17.42032+0.38+4.0220909+0*8.044182;}     //name_area15
    if (name_area==12) {phi1=(90.-85.6045)*3.14/180.; phi2=(90.-87.800)*3.14/180.; z1=-17.42032+0.38-4.0220909+1*8.044182; z2=-17.42032+0.38+4.0220909+1*8.044182;}     //name_area15
    if (name_area==13) {phi1=(90.-85.6045)*3.14/180.; phi2=(90.-87.800)*3.14/180.; z1=-17.42032+0.38-4.0220909+2*8.044182; z2=-17.42032+0.38+4.0220909+2*8.044182;}     //name_area15
    if (name_area==14) {phi1=(90.-85.6045)*3.14/180.; phi2=(90.-87.800)*3.14/180.; z1=-17.42032+0.38-4.0220909+3*8.044182; z2=-17.42032+0.38+4.0220909+3*8.044182;}     //name_area15
    if (name_area==15) {phi1=(90.-85.6045)*3.14/180.; phi2=(90.-87.800)*3.14/180.; z1=-17.42032+0.38-4.0220909+4*8.044182; z2=-17.42032+0.38+4.0220909+4*8.044182;}     //name_area15
    if (name_area==16) {phi1=(90.-87.800)*3.14/180.; phi2=(90.-89.100)*3.14/180.; z1=-17.42032+0.38-4.0220909+0*8.044182; z2=-17.42032+0.38+4.0220909+0*8.044182;}     //name_area20
    if (name_area==17) {phi1=(90.-87.800)*3.14/180.; phi2=(90.-89.100)*3.14/180.; z1=-17.42032+0.38-4.0220909+1*8.044182; z2=-17.42032+0.38+4.0220909+1*8.044182;}     //name_area20
    if (name_area==18) {phi1=(90.-87.800)*3.14/180.; phi2=(90.-89.100)*3.14/180.; z1=-17.42032+0.38-4.0220909+2*8.044182; z2=-17.42032+0.38+4.0220909+2*8.044182;}     //name_area20
    if (name_area==19) {phi1=(90.-87.800)*3.14/180.; phi2=(90.-89.100)*3.14/180.; z1=-17.42032+0.38-4.0220909+3*8.044182; z2=-17.42032+0.38+4.0220909+3*8.044182;}     //name_area20
    if (name_area==20) {phi1=(90.-87.800)*3.14/180.; phi2=(90.-89.100)*3.14/180.; z1=-17.42032+0.38-4.0220909+4*8.044182; z2=-17.42032+0.38+4.0220909+4*8.044182;}     //name_area20
    if (name_area==21) {phi1=(90.-89.100)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-17.42032+0.38-4.0220909+0*8.044182; z2=-17.42032+0.38+4.0220909+0*8.044182;}     //name_area24
    if (name_area==22) {phi1=(90.-89.100)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-17.42032+0.38-4.0220909+1*8.044182; z2=-17.42032+0.38+4.0220909+1*8.044182;}     //name_area24
    if (name_area==23) {phi1=(90.-89.100)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-17.42032+0.38-4.0220909+2*8.044182; z2=-17.42032+0.38+4.0220909+2*8.044182;}     //name_area24
    if (name_area==24) {phi1=(90.-89.100)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-17.42032+0.38-4.0220909+3*8.044182; z2=-17.42032+0.38+4.0220909+3*8.044182;}     //name_area24
    if (name_area==25) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-90.900)*3.14/180.; z1=-17.42032+0.38-4.0220909+0*8.044182; z2=-17.42032+0.38+4.0220909+0*8.044182;}     //name_area28
    if (name_area==26) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-90.900)*3.14/180.; z1=-17.42032+0.38-4.0220909+1*8.044182; z2=-17.42032+0.38+4.0220909+1*8.044182;}     //name_area28
    if (name_area==27) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-90.900)*3.14/180.; z1=-17.42032+0.38-4.0220909+2*8.044182; z2=-17.42032+0.38+4.0220909+2*8.044182;}     //name_area28
    if (name_area==28) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-90.900)*3.14/180.; z1=-17.42032+0.38-4.0220909+3*8.044182; z2=-17.42032+0.38+4.0220909+3*8.044182;}     //name_area28
    if (name_area==29) {phi1=(90.-90.900)*3.14/180.; phi2=(90.-92.200)*3.14/180.; z1=-17.42032+0.38-4.0220909+0*8.044182; z2=-17.42032+0.38+4.0220909+0*8.044182;}     //name_area33
    if (name_area==30) {phi1=(90.-90.900)*3.14/180.; phi2=(90.-92.200)*3.14/180.; z1=-17.42032+0.38-4.0220909+1*8.044182; z2=-17.42032+0.38+4.0220909+1*8.044182;}     //name_area33
    if (name_area==31) {phi1=(90.-90.900)*3.14/180.; phi2=(90.-92.200)*3.14/180.; z1=-17.42032+0.38-4.0220909+2*8.044182; z2=-17.42032+0.38+4.0220909+2*8.044182;}     //name_area33
    if (name_area==32) {phi1=(90.-90.900)*3.14/180.; phi2=(90.-92.200)*3.14/180.; z1=-17.42032+0.38-4.0220909+3*8.044182; z2=-17.42032+0.38+4.0220909+3*8.044182;}     //name_area33
    if (name_area==33) {phi1=(90.-90.900)*3.14/180.; phi2=(90.-92.200)*3.14/180.; z1=-17.42032+0.38-4.0220909+4*8.044182; z2=-17.42032+0.38+4.0220909+4*8.044182;}     //name_area33
    if (name_area==34) {phi1=(90.-92.200)*3.14/180.; phi2=(90.-94.39549)*3.14/180.; z1=-17.42032+0.38-4.0220909+0*8.044182; z2=-17.42032+0.38+4.0220909+0*8.044182;}    //name_area38
    if (name_area==35) {phi1=(90.-92.200)*3.14/180.; phi2=(90.-94.39549)*3.14/180.; z1=-17.42032+0.38-4.0220909+1*8.044182; z2=-17.42032+0.38+4.0220909+1*8.044182;}    //name_area38
    if (name_area==36) {phi1=(90.-92.200)*3.14/180.; phi2=(90.-94.39549)*3.14/180.; z1=-17.42032+0.38-4.0220909+2*8.044182; z2=-17.42032+0.38+4.0220909+2*8.044182;}    //name_area38
    if (name_area==37) {phi1=(90.-92.200)*3.14/180.; phi2=(90.-94.39549)*3.14/180.; z1=-17.42032+0.38-4.0220909+3*8.044182; z2=-17.42032+0.38+4.0220909+3*8.044182;}    //name_area38
    if (name_area==38) {phi1=(90.-92.200)*3.14/180.; phi2=(90.-94.39549)*3.14/180.; z1=-17.42032+0.38-4.0220909+4*8.044182; z2=-17.42032+0.38+4.0220909+4*8.044182;}    //name_area38
    if (name_area==39) {phi1=(90.-94.39549)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-17.42032+0.38-4.0220909+0*8.044182; z2=-17.42032+0.38+4.0220909+0*8.044182;}   //name_area43
    if (name_area==40) {phi1=(90.-94.39549)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-17.42032+0.38-4.0220909+1*8.044182; z2=-17.42032+0.38+4.0220909+1*8.044182;}   //name_area43
    if (name_area==41) {phi1=(90.-94.39549)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-17.42032+0.38-4.0220909+2*8.044182; z2=-17.42032+0.38+4.0220909+2*8.044182;}   //name_area43
    if (name_area==42) {phi1=(90.-94.39549)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-17.42032+0.38-4.0220909+3*8.044182; z2=-17.42032+0.38+4.0220909+3*8.044182;}   //name_area43
    if (name_area==43) {phi1=(90.-94.39549)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-17.42032+0.38-4.0220909+4*8.044182; z2=-17.42032+0.38+4.0220909+4*8.044182;}   //name_area43
    if (name_area==44) {phi1=(90.-97.291)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-17.42032+0.38-4.0220909+0*8.044182; z2=-17.42032+0.38+4.0220909+0*8.044182;}     //name_area48
    if (name_area==45) {phi1=(90.-97.291)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-17.42032+0.38-4.0220909+1*8.044182; z2=-17.42032+0.38+4.0220909+1*8.044182;}     //name_area48
    if (name_area==46) {phi1=(90.-97.291)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-17.42032+0.38-4.0220909+2*8.044182; z2=-17.42032+0.38+4.0220909+2*8.044182;}     //name_area48
    if (name_area==47) {phi1=(90.-97.291)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-17.42032+0.38-4.0220909+3*8.044182; z2=-17.42032+0.38+4.0220909+3*8.044182;}     //name_area48
    if (name_area==48) {phi1=(90.-97.291)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-17.42032+0.38-4.0220909+4*8.044182; z2=-17.42032+0.38+4.0220909+4*8.044182;}     //name_area48
    if (name_area==49) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-82.709)*3.14/180.; z1=-17.42032+0.38+4.0220909+4*8.044182; z2=-17.42032+0.38+4.0220909+4*8.044182+4.654-0.55;}  //name_area49
    if (name_area==50) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-86.212)*3.14/180.; z1=-17.42032+0.38+4.0220909+4*8.044182; z2=-17.42032+0.38+4.0220909+4*8.044182+(4.654-0.55)*2/3;}  //name_area50
    if (name_area==51) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-86.212)*3.14/180.; z1=-17.42032+0.38+4.0220909+4*8.044182+(4.654-0.55)*2/3; z2=-17.42032+0.38+4.0220909+4*8.044182+4.654-0.55;} //name_area51
    if (name_area==52) {phi1=(90.-93.788)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-17.42032+0.38+4.0220909+4*8.044182; z2=-17.42032+0.38+4.0220909+4*8.044182+(4.654-0.55)*2/3;}  //name_area52
    if (name_area==53) {phi1=(90.-93.788)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-17.42032+0.38+4.0220909+4*8.044182+(4.654-0.55)*2/3; z2=-17.42032+0.38+4.0220909+4*8.044182+4.654-0.55;}  //name_area53
    if (name_area==54) {phi1=(90.-97.291)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-17.42032+0.38+4.0220909+4*8.044182; z2=-17.42032+0.38+4.0220909+4*8.044182+4.654-0.55;}  //name_area54
    if (name_area==55) {phi1=(90.-89.100)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-17.42032+0.38+4.0220909+3*8.044182; z2=-17.42032+0.38+4.0220909+3*8.044182+(8.044182-3.0);}  //name_area55
    if (name_area==56) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-90.900)*3.14/180.; z1=-17.42032+0.38+4.0220909+3*8.044182; z2=-17.42032+0.38+4.0220909+3*8.044182+(8.044182-3.0);}     //name_area56
    if (name_area==57) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-17.42032+0.38-4.0220909-(8.46756/4+0.12); z2=-17.42032+0.38-4.0220909;}                       //name_area57
    if (name_area==58) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-17.42032+0.38-4.0220909-(8.46756/4+0.12); z2=-17.42032+0.38-4.0220909;}                       //AQ58
}
//======================================================BOLS,BOLD (101-120)===========================================================================
if(area==3){
  if (name_area==1) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-83.567)*3.14/180.; z1=-21.6842479-0.295-3.611752; z2=-21.6842479-0.295+3.611752;}    //name_area1
  if (name_area==2) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-82.709)*3.14/180.; z1=-21.6842479-0.295-0.800+8.742779-4.3713895; z2=-21.6842479-0.295-0.800+8.742779+4.3713895;}   //name_area2
  if (name_area==3) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-82.709)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+2*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+2*8.742779;}     //name_area6
  if (name_area==4) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-82.709)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+3*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+3*8.742779;}     //name_area6
  if (name_area==5) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-82.709)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+4*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+4*8.742779;}     //name_area6
  if (name_area==6) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-82.709)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+5*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+5*8.742779;}     //name_area6
  if (name_area==7) {phi1=(90.-83.567)*3.14/180.; phi2=(90.-85.6045)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895; z2=-21.6842479-0.295-0.800+4.3713895;}     //name_area7
  if (name_area==8) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-85.6045)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+1*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+1*8.742779;}    //name_area12
  if (name_area==9) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-85.6045)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+2*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+2*8.742779;}    //name_area12
  if (name_area==10) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-85.6045)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+3*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+3*8.742779;}    //name_area12
  if (name_area==11) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-85.6045)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+4*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+4*8.742779;}    //name_area12
  if (name_area==12) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-85.6045)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+5*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+5*8.742779;}    //name_area12
  if (name_area==13) {phi1=(90.-85.6045)*3.14/180.; phi2=(90.-87.800)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+0*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+0*8.742779;}    //name_area18
  if (name_area==14) {phi1=(90.-85.6045)*3.14/180.; phi2=(90.-87.800)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+1*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+1*8.742779;}    //name_area18
  if (name_area==15) {phi1=(90.-85.6045)*3.14/180.; phi2=(90.-87.800)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+2*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+2*8.742779;}    //name_area18
  if (name_area==16) {phi1=(90.-85.6045)*3.14/180.; phi2=(90.-87.800)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+3*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+3*8.742779;}    //name_area18
  if (name_area==17) {phi1=(90.-85.6045)*3.14/180.; phi2=(90.-87.800)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+4*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+4*8.742779;}    //name_area18
  if (name_area==18) {phi1=(90.-85.6045)*3.14/180.; phi2=(90.-87.800)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+5*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+5*8.742779;}    //name_area18
  if (name_area==19) {phi1=(90.-87.800)*3.14/180.; phi2=(90.-89.10000)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+0*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+0*8.742779;}     //name_area24
  if (name_area==20) {phi1=(90.-87.800)*3.14/180.; phi2=(90.-89.10000)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+1*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+1*8.742779;}     //name_area24
  if (name_area==21) {phi1=(90.-87.800)*3.14/180.; phi2=(90.-89.10000)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+2*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+2*8.742779;}     //name_area24
  if (name_area==22) {phi1=(90.-87.800)*3.14/180.; phi2=(90.-89.10000)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+3*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+3*8.742779;}     //name_area24
  if (name_area==23) {phi1=(90.-87.800)*3.14/180.; phi2=(90.-89.10000)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+4*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+4*8.742779;}     //name_area24
  if (name_area==24) {phi1=(90.-87.800)*3.14/180.; phi2=(90.-89.10000)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+5*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+5*8.742779;}     //name_area24
  if (name_area==25) {phi1=(90.-89.100)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+0*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+0*8.742779;}     //name_area29
  if (name_area==26) {phi1=(90.-89.100)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+1*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+1*8.742779;}     //name_area29
  if (name_area==27) {phi1=(90.-89.100)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+2*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+2*8.742779;}     //name_area29
  if (name_area==28) {phi1=(90.-89.100)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+3*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+3*8.742779;}     //name_area29
  if (name_area==29) {phi1=(90.-89.100)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+4*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+4*8.742779;}     //name_area29
  if (name_area==30) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-90.900)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+0*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+0*8.742779;}     //name_area34
  if (name_area==31) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-90.900)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+1*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+1*8.742779;}     //name_area34
  if (name_area==32) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-90.900)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+2*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+2*8.742779;}     //name_area34
  if (name_area==33) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-90.900)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+3*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+3*8.742779;}     //name_area34
  if (name_area==34) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-90.900)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+4*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+4*8.742779;}     //name_area34
  if (name_area==35) {phi1=(90.-90.900)*3.14/180.; phi2=(90.-92.200)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+0*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+0*8.742779;}     //name_area40
  if (name_area==36) {phi1=(90.-90.900)*3.14/180.; phi2=(90.-92.200)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+1*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+1*8.742779;}     //name_area40
  if (name_area==37) {phi1=(90.-90.900)*3.14/180.; phi2=(90.-92.200)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+2*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+2*8.742779;}     //name_area40
  if (name_area==38) {phi1=(90.-90.900)*3.14/180.; phi2=(90.-92.200)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+3*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+3*8.742779;}     //name_area40
  if (name_area==39) {phi1=(90.-90.900)*3.14/180.; phi2=(90.-92.200)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+4*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+4*8.742779;}     //name_area40
  if (name_area==40) {phi1=(90.-90.900)*3.14/180.; phi2=(90.-92.200)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+5*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+5*8.742779;}     //name_area40
  if (name_area==41) {phi1=(90.-92.200)*3.14/180.; phi2=(90.-94.39549)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+0*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+0*8.742779;}   //name_area46
  if (name_area==42) {phi1=(90.-92.200)*3.14/180.; phi2=(90.-94.39549)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+1*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+1*8.742779;}   //name_area46
  if (name_area==43) {phi1=(90.-92.200)*3.14/180.; phi2=(90.-94.39549)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+2*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+2*8.742779;}   //name_area46
  if (name_area==44) {phi1=(90.-92.200)*3.14/180.; phi2=(90.-94.39549)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+3*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+3*8.742779;}   //name_area46
  if (name_area==45) {phi1=(90.-92.200)*3.14/180.; phi2=(90.-94.39549)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+4*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+4*8.742779;}   //name_area46
  if (name_area==46) {phi1=(90.-92.200)*3.14/180.; phi2=(90.-94.39549)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+5*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+5*8.742779;}   //name_area46
  if (name_area==47) {phi1=(90.-94.39549)*3.14/180.; phi2=(90.-96.433)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895; z2=-21.6842479-0.295-0.800+4.3713895;}                      //name_area47
  if (name_area==48) {phi1=(90.-94.39549)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+1*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+1*8.742779;}     //name_area52
  if (name_area==49) {phi1=(90.-94.39549)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+2*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+2*8.742779;}     //name_area52
  if (name_area==50) {phi1=(90.-94.39549)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+3*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+3*8.742779;}     //name_area52
  if (name_area==51) {phi1=(90.-94.39549)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+4*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+4*8.742779;}     //name_area52
  if (name_area==52) {phi1=(90.-94.39549)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+5*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+5*8.742779;}     //name_area52
  if (name_area==53) {phi1=(90.-96.433)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-21.6842479-0.295-3.611752; z2=-21.6842479-0.295+3.611752;}                                     //name_area53
  if (name_area==54) {phi1=(90.-97.291)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+1*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+1*8.742779;}    //name_area58
  if (name_area==55) {phi1=(90.-97.291)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+2*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+2*8.742779;}    //name_area58
  if (name_area==56) {phi1=(90.-97.291)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+3*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+3*8.742779;}    //name_area58
  if (name_area==57) {phi1=(90.-97.291)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+4*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+4*8.742779;}    //name_area58
  if (name_area==58) {phi1=(90.-97.291)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895+5*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+5*8.742779;}    //name_area58
  if (name_area==59) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-82.709)*3.14/180.; z1=-21.6842479-0.295-0.800+4.3713895+5*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+5*8.742779+(4.654-0.55);}     //name_area59
  if (name_area==60) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-86.212)*3.14/180.; z1=-21.6842479-0.295-0.800+4.3713895+5*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+5*8.742779+(4.654-0.55)*2/3;} //name_area60
  if (name_area==61) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-86.212)*3.14/180.; z1=-21.6842479-0.295-0.800+4.3713895+5*8.742779+(4.654-0.55)*2/3; z2=-21.6842479-0.295-0.800+4.3713895+5*8.742779+4.654-0.55;} //name_area61
  if (name_area==62) {phi1=(90.-93.788)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-21.6842479-0.295-0.800+4.3713895+5*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+5*8.742779+(4.654-0.55)*2/3;}   //name_area62
  if (name_area==63) {phi1=(90.-93.788)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-21.6842479-0.295-0.800+4.3713895+5*8.742779+(4.654-0.55)*2/3; z2=-21.6842479-0.295-0.800+4.3713895+5*8.742779+(4.654-0.55);}    //name_area63
  if (name_area==64) {phi1=(90.-97.291)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-21.6842479-0.295-0.800+4.3713895+5*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+5*8.742779+(4.654-0.55);}    //name_area64
  if (name_area==65) {phi1=(90.-89.100)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-21.6842479-0.295-0.800+4.3713895+4*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+4*8.742779+(8.742779-3.0);}  //name_area65
  if (name_area==66) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-90.900)*3.14/180.; z1=-21.6842479-0.295-0.800+4.3713895+4*8.742779; z2=-21.6842479-0.295-0.800+4.3713895+4*8.742779+(8.742779-3.0);}     //name_area66
  if (name_area==67) {phi1=(90.-83.567)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895-(9.1229/4+0.12); z2=-21.6842479-0.295-0.800-4.3713895;}    //name_area67
  if (name_area==68) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-96.433)*3.14/180.; z1=-21.6842479-0.295-0.800-4.3713895-(9.1229/4+0.12); z2=-21.6842479-0.295-0.800-4.3713895;}     //AR68
}
//=======================================================================================================================================================

//======================================================BISS,BISD (21-40)===========================================================================
if(area==1){
     if (name_area==1) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-83.567)*3.14/180.; z1=-16.790464+0.60-3.1805359; z2=-16.790464+0.60+3.1805359;}                                  
     if (name_area==2) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-82.709)*3.14/180.; z1=-16.790464+0.60-0.85+8.044182-4.0220909; z2=-16.790464+0.60-0.85+8.044182+4.0220909;}      
     if (name_area==3) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-82.709)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+2*8.044182; z2=-16.790464+0.60-0.85+4.0220909+2*8.044182;}    
     if (name_area==4) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-82.709)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+3*8.044182; z2=-16.790464+0.60-0.85+4.0220909+3*8.044182;}    
     if (name_area==5) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-82.709)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+4*8.044182; z2=-16.790464+0.60-0.85+4.0220909+4*8.044182;}    
     if (name_area==6) {phi1=(90.-83.567)*3.14/180.; phi2=(90.-85.6045)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+0*8.044182; z2=-16.790464+0.60-0.85+4.0220909+0*8.044182;}  
     if (name_area==7) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-85.6045)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+1*8.044182; z2=-16.790464+0.60-0.85+4.0220909+1*8.044182;}   
     if (name_area==8) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-85.6045)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+2*8.044182; z2=-16.790464+0.60-0.85+4.0220909+2*8.044182;}   
     if (name_area==9) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-85.6045)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+3*8.044182; z2=-16.790464+0.60-0.85+4.0220909+3*8.044182;}   
     if (name_area==10) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-85.6045)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+4*8.044182; z2=-16.790464+0.60-0.85+4.0220909+4*8.044182;}   
     if (name_area==11) {phi1=(90.-85.6045)*3.14/180.; phi2=(90.-87.800)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+0*8.044182; z2=-16.790464+0.60-0.85+4.0220909+0*8.044182;}  
     if (name_area==12) {phi1=(90.-85.6045)*3.14/180.; phi2=(90.-87.800)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+1*8.044182; z2=-16.790464+0.60-0.85+4.0220909+1*8.044182;}   
     if (name_area==13) {phi1=(90.-85.6045)*3.14/180.; phi2=(90.-87.800)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+2*8.044182; z2=-16.790464+0.60-0.85+4.0220909+2*8.044182;}   
     if (name_area==14) {phi1=(90.-85.6045)*3.14/180.; phi2=(90.-87.800)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+3*8.044182; z2=-16.790464+0.60-0.85+4.0220909+3*8.044182;}   
     if (name_area==15) {phi1=(90.-85.6045)*3.14/180.; phi2=(90.-87.800)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+4*8.044182; z2=-16.790464+0.60-0.85+4.0220909+4*8.044182;}   
     if (name_area==16) {phi1=(90.-87.800)*3.14/180.; phi2=(90.-89.100)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+0*8.044182;  z2=-16.790464+0.60-0.85+4.0220909+0*8.044182;}    
     if (name_area==17) {phi1=(90.-87.800)*3.14/180.; phi2=(90.-89.100)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+1*8.044182;  z2=-16.790464+0.60-0.85+4.0220909+1*8.044182;}     
     if (name_area==18) {phi1=(90.-87.800)*3.14/180.; phi2=(90.-89.100)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+2*8.044182;  z2=-16.790464+0.60-0.85+4.0220909+2*8.044182;}     
     if (name_area==19) {phi1=(90.-87.800)*3.14/180.; phi2=(90.-89.100)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+3*8.044182;  z2=-16.790464+0.60-0.85+4.0220909+3*8.044182;}     
     if (name_area==20) {phi1=(90.-87.800)*3.14/180.; phi2=(90.-89.100)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+4*8.044182;  z2=-16.790464+0.60-0.85+4.0220909+4*8.044182;}     
     if (name_area==21) {phi1=(90.-89.100)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+0*8.044182;  z2=-16.790464+0.60-0.85+4.0220909+0*8.044182;}     
     if (name_area==22) {phi1=(90.-89.100)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+1*8.044182;  z2=-16.790464+0.60-0.85+4.0220909+1*8.044182;}     
     if (name_area==23) {phi1=(90.-89.100)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+2*8.044182;  z2=-16.790464+0.60-0.85+4.0220909+2*8.044182;}     
     if (name_area==24) {phi1=(90.-89.100)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+3*8.044182;  z2=-16.790464+0.60-0.85+4.0220909+3*8.044182;}     
     if (name_area==25) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-90.900)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+0*8.044182;  z2=-16.790464+0.60-0.85+4.0220909+0*8.044182;}     
     if (name_area==26) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-90.900)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+1*8.044182;  z2=-16.790464+0.60-0.85+4.0220909+1*8.044182;}     
     if (name_area==27) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-90.900)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+2*8.044182;  z2=-16.790464+0.60-0.85+4.0220909+2*8.044182;}    
     if (name_area==28) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-90.900)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+3*8.044182;  z2=-16.790464+0.60-0.85+4.0220909+3*8.044182;}     
     if (name_area==29) {phi1=(90.-90.900)*3.14/180.; phi2=(90.-92.200)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+0*8.044182;  z2=-16.790464+0.60-0.85+4.0220909+0*8.044182;}     
     if (name_area==30) {phi1=(90.-90.900)*3.14/180.; phi2=(90.-92.200)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+1*8.044182;  z2=-16.790464+0.60-0.85+4.0220909+1*8.044182;}     
     if (name_area==31) {phi1=(90.-90.900)*3.14/180.; phi2=(90.-92.200)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+2*8.044182;  z2=-16.790464+0.60-0.85+4.0220909+2*8.044182;}     
     if (name_area==32) {phi1=(90.-90.900)*3.14/180.; phi2=(90.-92.200)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+3*8.044182;  z2=-16.790464+0.60-0.85+4.0220909+3*8.044182;}     
     if (name_area==33) {phi1=(90.-90.900)*3.14/180.; phi2=(90.-92.200)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+4*8.044182;  z2=-16.790464+0.60-0.85+4.0220909+4*8.044182;}   
     if (name_area==34) {phi1=(90.-92.200)*3.14/180.; phi2=(90.-94.39549)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+0*8.044182; z2=-16.790464+0.60-0.85+4.0220909+0*8.044182;}     
     if (name_area==35) {phi1=(90.-92.200)*3.14/180.; phi2=(90.-94.39549)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+1*8.044182; z2=-16.790464+0.60-0.85+4.0220909+1*8.044182;}     
     if (name_area==36) {phi1=(90.-92.200)*3.14/180.; phi2=(90.-94.39549)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+2*8.044182; z2=-16.790464+0.60-0.85+4.0220909+2*8.044182;}    
     if (name_area==37) {phi1=(90.-92.200)*3.14/180.; phi2=(90.-94.39549)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+3*8.044182; z2=-16.790464+0.60-0.85+4.0220909+3*8.044182;}     
     if (name_area==38) {phi1=(90.-92.200)*3.14/180.; phi2=(90.-94.39549)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+4*8.044182; z2=-16.790464+0.60-0.85+4.0220909+4*8.044182;}     
     if (name_area==39) {phi1=(90.-94.39549)*3.14/180.; phi2=(90.-96.433)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909; z2=-16.790464+0.60-0.85+4.0220909;}                          
     if (name_area==40) {phi1=(90.-94.39549)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+1*8.044182; z2=-16.790464+0.60-0.85+4.0220909+1*8.044182;}     
     if (name_area==41) {phi1=(90.-94.39549)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+2*8.044182; z2=-16.790464+0.60-0.85+4.0220909+2*8.044182;}     
     if (name_area==42) {phi1=(90.-94.39549)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+3*8.044182; z2=-16.790464+0.60-0.85+4.0220909+3*8.044182;}     
     if (name_area==43) {phi1=(90.-94.39549)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+4*8.044182; z2=-16.790464+0.60-0.85+4.0220909+4*8.044182;}    
     if (name_area==44) {phi1=(90.-96.433)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-16.790464+0.60-3.1805359; z2=-16.790464+0.60+3.1805359;}                                     
     if (name_area==45) {phi1=(90.-97.291)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+1*8.044182; z2=-16.790464+0.60-0.85+4.0220909+1*8.044182;}     
     if (name_area==46) {phi1=(90.-97.291)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+2*8.044182; z2=-16.790464+0.60-0.85+4.0220909+2*8.044182;}     
     if (name_area==47) {phi1=(90.-97.291)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+3*8.044182; z2=-16.790464+0.60-0.85+4.0220909+3*8.044182;}     
     if (name_area==48) {phi1=(90.-97.291)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909+4*8.044182; z2=-16.790464+0.60-0.85+4.0220909+4*8.044182;}     
     if (name_area==49) {phi1=(90.-81.607)*3.14/180.; phi2=(90.-82.709)*3.14/180.; z1=-16.790464+0.60-0.85+4.0220909+4*8.044182; z2=-16.790464+0.60-0.85+4.0220909+4*8.044182+(4.654-0.55);}  
     if (name_area==50) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-86.212)*3.14/180.; z1=-16.790464+0.60-0.85+4.0220909+4*8.044182; z2=-16.790464+0.60-0.85+4.0220909+4*8.044182+(4.654-0.55)*2/3;}    
     if (name_area==51) {phi1=(90.-82.709)*3.14/180.; phi2=(90.-86.212)*3.14/180.; z1=-16.790464+0.60-0.85+4.0220909+4*8.044182+(4.654-0.55)*2/3; z2=-16.790464+0.60-0.85+4.0220909+4*8.044182+(4.654-0.55);}    
     if (name_area==52) {phi1=(90.-93.788)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-16.790464+0.60-0.85+4.0220909+4*8.044182; z2=-16.790464+0.60-0.85+4.0220909+4*8.044182+(4.654-0.55)*2/3;}    
     if (name_area==53) {phi1=(90.-93.788)*3.14/180.; phi2=(90.-97.291)*3.14/180.; z1=-16.790464+0.60-0.85+4.0220909+4*8.044182+(4.654-0.55)*2/3; z2=-16.790464+0.60-0.85+4.0220909+4*8.044182+(4.654-0.55);}   
     if (name_area==54) {phi1=(90.-97.291)*3.14/180.; phi2=(90.-98.393)*3.14/180.; z1=-16.790464+0.60-0.85+4.0220909+4*8.044182; z2=-16.790464+0.60-0.85+4.0220909+4*8.044182+(4.654-0.55);}     
     if (name_area==55) {phi1=(90.- 89.100)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-16.790464+0.60-0.85+4.0220909+3*8.044182; z2=-16.790464+0.60-0.85+4.0220909+3*8.044182+(8.044182-3.0);}  
     if (name_area==56) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-90.900)*3.14/180.; z1=-16.790464+0.60-0.85+4.0220909+3*8.044182; z2=-16.790464+0.60-0.85+4.0220909+3*8.044182+(8.044182-3.0);}  
     if (name_area==57) {phi1=(90.-83.567)*3.14/180.; phi2=(90.-89.800)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909-(8.46756/4+0.12); z2=-16.790464+0.60-0.85-4.0220909;}     
     if (name_area==58) {phi1=(90.-90.200)*3.14/180.; phi2=(90.-96.433)*3.14/180.; z1=-16.790464+0.60-0.85-4.0220909-(8.46756/4+0.12); z2=-16.790464+0.60-0.85-4.0220909;}     
}
//===============================================================================================================================================================================================

//======================================================ECDI, ECSI========================================================================================
if(area==5){
   if (name_area==1) {x2=10; x1=0.495-0.090; y2=21.9759+3.5; y1=21.9759-3.5; }  
   if (name_area==2) {x2=1.5749+2.85*0.5-0.090; x1=1.5749-2.85*0.5-0.090; y2=29.3132+3.83747; y1=29.3132-3.83747;}    
   if (name_area==3) {x2=1.5749+2.85*0.5-0.090; x1=1.5749-2.85*0.5-0.090; y2=36.98814+3.83747; y1=36.98814-3.83747;}  
   if (name_area==4) {x2=1.5749+2.85*0.5-0.090; x1=1.5749-2.85*0.5-0.090; y2=44.66308+3.83747; y1=44.66308-3.83747;}  
   if (name_area==5) {x2=1.5749+2.85*0.5-0.090; x1=1.5749-2.85*0.5-0.090; y2=52.33802+3.83747; y1=52.33802-3.83747;}  
   if (name_area==6) {x2=10; x1=1.5749+2.85*0.5-0.090; y2=33.97589+8.5;  y1=33.97589-8.5;}   
   if (name_area==7) {x2=4.3824+2.765*0.5-0.090; x1=4.3824-2.765*0.5-0.090; y2=45.9758+3.5; y1=45.9758-3.5;}     
   if (name_area==8) {x2=4.3824+2.765*0.5-0.090; x1=4.3824-2.765*0.5-0.090; y2=54.47579+5.0; y1=54.47579-5.0;} 
   if (name_area==9) {x2=15; x1=4.3824+2.765*0.5-0.090; y2=50.97588+8.5; y1=50.97587-8.5;}   
   if (name_area==10) {x2=15; x1=4.3824+2.765*0.5-0.090; y2=64.4758+5.0; y1=64.4758-5.0;}   
   if (name_area==11) {x2=4.3824+2.765*0.5-0.090; x1=4.3824+2.765*0.5-3.35-0.090; y2=67.8758+1.6; y1=67.8758-1.6;}       
   if (name_area==12) {x2=4.3824+2.765*0.5-0.090; x1=4.3824+2.765*0.5-3.55-0.090; y2=64.87579+1.4; y1=64.87579-1.4;}      
   if (name_area==13) {x2=4.3824+2.765*0.5-0.090; x1=4.3824+2.765*0.5-4.95-0.090; y2=61.525789+1.95; y1=61.525789-1.95;}  
   if (name_area==14) {x2=4.3824-2.765*0.5-0.090; x1=4.3824-2.765*0.5-2.15-0.090; y2=57.925789+1.65; y1=57.925789-1.65;}  
   if (name_area==15) {x2=-0.495-0.090; x1=-10; y2=21.9759+3.5; y1=21.9759-3.5;}  //area29 <-area1
   if (name_area==16) {x2=-1.5749+2.85*0.5-0.090; x1=-1.5749-2.85*0.5-0.090; y2=29.3132+3.83747; y1=29.3132-3.83747;}      
   if (name_area==17) {x2=-1.5749+2.85*0.5-0.090; x1=-1.5749-2.85*0.5-0.090; y2=36.98814+3.83747; y1=36.98814-3.83747;}   
   if (name_area==18) {x2=-1.5749+2.85*0.5-0.090; x1=-1.5749-2.85*0.5-0.090; y2=44.66308+3.83747; y1=44.66308-3.83747;}   
   if (name_area==19) {x2=-1.5749+2.85*0.5-0.090; x1=-1.5749-2.85*0.5-0.090; y2=52.33802+3.83747; y1=52.33802-3.83747;}   
   if (name_area==20) {x2=-1.5749-2.85*0.5-0.090; x1=-15; y2=33.97589+8.5;  y1=33.97589-8.5;} 
   if (name_area==21) {x2=-4.3824+2.765*0.5-0.090; x1=-4.3824-2.765*0.5-0.090; y2=45.9758+3.5; y1=45.9758-3.5;}    
   if (name_area==22) {x2=-4.3824+2.765*0.5-0.090; x1=-4.3824-2.765*0.5-0.090; y2=54.47579+5.0; y1=54.47579-5.0;}  
   if (name_area==23) {x2=-4.3824-2.765*0.5-0.090; x1=-15; y2=50.97588+8.5; y1=50.97587-8.5;}   
   if (name_area==24) {x2=-4.3824-2.765*0.5-0.090; x1=-15; y2=64.4758+5.0; y1=64.4758-5.0;}     
   if (name_area==25) {x1=-4.3824-2.765*0.5-0.090; x2=-4.3824-2.765*0.5+3.35-0.090; y2=67.8758+1.6;  y1=67.8758-1.6;}       
   if (name_area==26) {x1=-4.3824-2.765*0.5-0.090; x2=-4.3824-2.765*0.5+3.55-0.090; y2=64.87579+1.4; y1=64.87579-1.4;}      
   if (name_area==27) {x1=-4.3824-2.765*0.5-0.090; x2=-4.3824-2.765*0.5+4.95-0.090; y2=61.525789+1.95; y1=61.525789-1.95;}   
   if (name_area==28) {x1=-4.3824+2.765*0.5-0.090; x2=-4.3824+2.765*0.5+2.15-0.090; y2=57.925789+1.65; y1=57.925789-1.65;}   // ex55 <-ex27
}
//=======================================================================================================================================================

//for (int n=cnt_begin; n<cnt_end; n++)
//{
for (int n=0; n<9; n++)
{         
  //if( area==5 && (n==20||n==100)){n=n+40;}            //0-20 60-100  140-160
// if(n<20 || n>=60 && n<100 || n>=140)
//{
	sprintf(branchname,"cnt%d",n);
	tt->SetBranchStatus("*",0);
	tt->SetBranchStatus(branchname,1);
	tt->SetBranchAddress(branchname,&cnt);              //СЏв”ђСЏв”‚СЏв”ЊРїв•џРїв•«Рїв•џРїв•ЎРїв•©Рїв•¦Рїв•ЎРїв•џРїв•ЈРїв•Є Рїв•џРїв•ўСЏв”ЂРїв•ЈСЏв”‚СЏв”‚ Рїв•ЎРїв•ЈСЏв”ЊРїв•ЁРїв•¦

	cout<<"Counter="<<n<<endl;

//TProfile* pr=new TProfile(branchname,branchname,1200,-600,600,0,200);        //z

for(int k=0; k<nentr; k++)                               //СЏв”њРїв•¦Рїв•ЁРїв•© РїВ©Рїв•¬ Рїв•ЎСЏв”‚Рїв•ЈРїв•Є СЏв”‚Рїв•¬Рїв• СЏв–ЂСЏв”ЊРїв•¦СЏв–ђРїв•Є Рїв•Ў Рїв•ўРїв•ЈСЏв”ЂРїв•ЈРїв•ЎРїв•Ј	
{
    tt->GetEntry(k);
    
         if( area==5 && (cnt.i>=20&&cnt.i<=60)){cnt_begin=60;  cnt_end=100;}            //0-20 60-100  140-160
         if( area==5 && (cnt.i>=140)){cnt_begin=140;  cnt_end=160;}            //0-20 60-100  140-160

    	 if( cnt.i>=cnt_begin && cnt.i<cnt_end )  //цикл по числу счетчиков
    {     
 
 if(area==1){ phi_begin=(cnt.phiin+cnt.phiout)/2; phi_end=(cnt.phiin+cnt.phiout)/2; z_begin=(cnt.zin+cnt.zout)/2; z_end=(cnt.zin+cnt.zout)/2; phi_area1=phi1-0.0015; phi_area2=phi2-0.0015; z_area1=z1; z_area2=z2; ATCregion=cnt.aerogel_region0; }
 if(area==2){ phi_begin=(cnt.phiin+cnt.phiout)/2; phi_end=(cnt.phiin+cnt.phiout)/2; z_begin=(cnt.zin+cnt.zout)/2; z_end=(cnt.zin+cnt.zout)/2; phi_area1=phi1+0.0055; phi_area2=phi2+0.00055; z_area1=z1; z_area2=z2; ATCregion=cnt.aerogel_region0;}
 if(area==3){ phi_begin=(cnt.phiin+cnt.phiout)/2; phi_end=(cnt.phiin+cnt.phiout)/2; z_begin=(cnt.zin+cnt.zout)/2; z_end=(cnt.zin+cnt.zout)/2; phi_area1=phi1-0.002785; phi_area2=phi2-0.002785; z_area1=z1; z_area2=z2; ATCregion=cnt.aerogel_region0;}
 if(area==4){ phi_begin=(cnt.phiin+cnt.phiout)/2; phi_end=(cnt.phiin+cnt.phiout)/2; z_begin=(cnt.zin+cnt.zout)/2; z_end=(cnt.zin+cnt.zout)/2; phi_area1=phi1+0.001874; phi_area2=phi2+0.001874; z_area1=z1; z_area2=z2; ATCregion=cnt.aerogel_region0;}
 //cout<<"check area"<<area<<endl; }
 if(area==5){ phi_begin=(cnt.rin*sin(cnt.phiin)+cnt.rout*sin(cnt.phiout))/2; phi_end=(cnt.rin*sin(cnt.phiin)+cnt.rout*sin(cnt.phiout))/2; z_begin=(cnt.rin*cos(cnt.phiin)+cnt.rout*cos(cnt.phiout))/2; z_end=(cnt.rin*cos(cnt.phiin)+cnt.rout*cos(cnt.phiout))/2; phi_area1=x2; phi_area2=x1; z_area1=y1; z_area2=y2; ATCregion=cnt.single_aerogel_region0;} 
  
//if( cnt.phiin!=0&&cnt.phiout!=0 && ATCregion!=0 && cnt.wlshit!=1 && cnt.nearwls!=1 && cnt.tlen>2.0 && phi_begin<phi_area1 && phi_end>phi_area2 && z_begin>z_area1 && z_end<z_area2 ) { 
if( ATCregion!=0 && cnt.wlshit!=1 && cnt.nearwls!=1 && phi_begin<phi_area1 && phi_end>phi_area2 && z_begin>z_area1 && z_end<z_area2 ) { 

	 tt2->GetEntry(k);
  //if (t.nvec>2 && t.chi2<150 && sqrt(pow(t.x0,2)+pow(t.y0,2)+pow(t.z0,2))<40 && t.emc_ncls>1 && t.p>1000 ) {
  //if (t.nvec>4 && t.chi2<150 && sqrt(pow(t.x0,2)+pow(t.y0,2)+pow(t.z0,2))<40 && t.emc_ncls>1 ) {
  //if (t.nvec>4 && t.chi2<20 && sqrt(pow(t.x0,2)+pow(t.y0,2)+pow(t.z0,2))<40 && t.emc_ncls>1 ) {
  if (t.nvec>4 && t.chi2<50 && t.emc_ncls>1 ) {            //*****************************************************************
  //if (t.nvec>4 && t.chi2<50 && t.emc_ncls>1 && t.p>1100 ) {
	    tt3->GetEntry(k);
	if (tof.nhits>=2) {    //for endcap
	    tt7->GetEntry(k);  
	    if (mu.status>=1 && mu.nhits>5 && mu.layer>=1) {
	    tt4->GetEntry(k);
	if (emc.energy>100 && emc.ncls_trk>=2 && emc.ncls>=2) {
	    tt5->GetEntry(k);
	    tt6->GetEntry(k);
	if( c0.e>50 && c1.e>50 ) { 

         betta=sqrt(1-1/pow(((emc.energy+105.65)/105.65),2));
 
 //1,2,3,5,13,15,19,23,26,37,39,40,57,60,63,64,66,93,109,113,125,129,152
    if( cnt.i==0 || cnt.i==1 ||  cnt.i==2 || cnt.i==4 || cnt.i==12 || cnt.i==14 || cnt.i==18 || cnt.i==22 || cnt.i==25 || cnt.i==36 || cnt.i==38 || cnt.i==39 || cnt.i==56 || cnt.i==59 || cnt.i==62 ||  cnt.i==63 ||  cnt.i==65 || cnt.i==92 ||  cnt.i==108 || cnt.i==112 || cnt.i==124 || cnt.i==128 ||  cnt.i==151) {                                            //номера счетчиков с плотным аэрогелем
       LC=(cnt.npe/cnt.tlen)/(1.-1./pow(1.05576*betta,2)); 
	}                
       else{ 
    //LC=(cnt.npe/cnt.tlen)/(1.-1./pow(1.05,2));
    LC=(cnt.npe/cnt.tlen)/(1.-1./pow(1.04752*betta,2));
	   } 
	   i++;
	   ev_all++;    
            //cout<<LC<<"\t"<<cnt.npe<<"\t"<<t.p<<"\t"<<emc.energy<<"\t"<<betta<<"\t"<<TMath::C()<<"\t"<<cnt.tlen<<endl;        
            ofstream of(TString::Format("data_cnt_Klight_area_all_%dregion_%darea.dat",region,area).Data(),ios::app);  //запись LC в текстовый файл для последующего фитирования    
 	    of<<i<<"\t"<<LC<<endl;
            of.close();
//cout<<i<<"\t"<<"j="<<j<<"\t"<<"LC="<<LC<<"\t"<<"Total_LC="<<total_LC<<"\t"<<"npe="<<cnt.npe<<"\t"<<"wlshit="<<cnt.wlshit<<"\t"<<"nearwls="<<cnt.nearwls<<"\t"<<"Average_LC="<<average_LC<<endl;
	
	h31->Fill((cnt.zin+cnt.zout)/2,(cnt.phiin+cnt.phiout)/2);
	
	}
	}
	}
	}
	}
    }
  }   //if(n<20 || n>=60 && n<100 || n>=140)
}
// }  //if cnt   
}  //cnt

h31->SetXTitle("Zin/out, cm");
h31->SetYTitle("Phiin/out, rad");
h31->SetLineColor(3);          //red
h31->Draw("colz");
//c11->SaveAs(TString::Format("test_biss_26072017_3range.root").Data());

//   float Weight;
//   if(area==1 && region==12 ){Weight=i/392142.;}   //AY       
//   if(area==2 && region==12 ){Weight=i/640371.;}   //AS       
//   if(area==3 && region==12 ){Weight=i/579737.;}   //AR      
//   if(area==4 && region==12 ){Weight=i/408414.;}   //AQ      
//   if(area==5 && region==12 ){Weight=i/2389739.;}   //EX       
//cout<<"Weight="<<Weight<<endl;
//   cout<<"i="<<i<<"\t"<<std::setprecision(6)<<"ev_all="<<ev_all<<"\t"<<std::setprecision(6)<<"Weight="<<Weight<<endl;
//   ofstream of3(TString::Format("data_weight_%dregion_%darea.dat",region,area).Data(),ios::app);  //запись LC в текстовый файл для последующего фитирования    
//   of3<<Weight<<"\t"<<LC<<endl;
//   of3.close();

  //========фитирование распределения Klight(LC)===========================================================================
   float x0[i],y0[i];
   ostringstream o_str;
   const string DIR_REF = TString::Format("data_cnt_Klight_area_all_%dregion_%darea.dat",region,area).Data();              

   if(i<100){ numbin=3*i; }   
     else{
     	numbin=i;             //число событий в области    
   }
   
   TH1F *h1=new TH1F("h1","The distribution Klihgt",numbin,0,100);
   
   string f_in_ref;
   o_str << DIR_REF;
   f_in_ref = o_str.str(); o_str.str("");

   FILE *pF;
   pF=fopen(f_in_ref.c_str(), "r");
   for (int k=0; k<i; k++)
   {
      fscanf(pF," %f %f ",&x0[k],&y0[k]);
      h1->Fill(y0[k]);
     //cout<<x0[k]<<"\t"<<y0[k]<<endl;
   }
   fclose(pF);

   h1->Draw("hist");

//   TF1* myfit1=new TF1("myfit1","gaus");
//   h1->Fit("myfit1","","",0.1,25);

   TF1 myfit1("myfit1",poissonf,0,10,2); // x in [0;10], 2parameters                  
   
   myfit1.SetParName(0,"Const");                                                
   myfit1.SetParName(1,"#mu");    
   
   myfit1.SetParameter(0,0.01);                                              
   myfit1.SetParameter(1,0.1);    
  
   h1->Fit("myfit1","","",1,25);    


   //average_LC=Weight*myfit1.GetParameter(1);                 //средний коэффициент светосбора в области - профитированный
   //err_LC=Weight*myfit1.GetParError(1); 
   average_LC=myfit1.GetParameter(1);                 //средний коэффициент светосбора в области - профитированный
   err_LC=myfit1.GetParError(1); 
   cout<<"Mean="<<myfit1.GetParameter(1)<<"\t"<<"average_LC="<<average_LC<<"\t"<<"err_LC="<<myfit1.GetParError(1)<<endl;
   //h1->delete();
   //myfit1->delete();
    
  //====================================================================================

//ofstream of2(TString::Format("data_%d_13082017_2range_test.dat",area).Data(),ios::app);   
ofstream of2(NAME_str.c_str(),ios::app);   
of2<<name_area<<"\t"<<average_LC<<"\t"<<err_LC<<endl;                 
of2.close();

}  //area

} 

