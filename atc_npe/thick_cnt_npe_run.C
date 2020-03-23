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
#include <stdio.h>
#include <stdlib.h>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TLine.h>
#include <TEventList.h>
#include <TProfile.h>
#include <vector>
#include <TChain.h>
#include <string>
#include <iostream>
#include <fstream>
#include <time.h>
#include "TMinuit.h"
#include "TRandom3.h"
#include <math.h>
#include "TVirtualFitter.h"
#include <algorithm> 
#include <TGraphErrors.h>

using namespace std;
string progname;
//void npe_run(int counter)
//void npe_run()
int main(int argc, char* argv[])
{
  struct data                                  //СЏв”‚СЏв”ЊСЏв”ЂСЏв”ђРїв•ЁСЏв”ЊСЏв”ЂСЏв”ђРїв•џ СЏв”‚ Рїв•ўРїв•џРїв•«Рїв•«СЏв–ЂРїв•ЄРїв•¦ Рїв–‘Рїв•–Рїв•
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
  { UInt_t nhits, dchits, namps, ntimes;
    float time[2], length[2], beta[2], phi[2];
    UInt_t type[2];
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
    UInt_t nhits,dcmuhits,octant,layer;
    UInt_t status;
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

  TChain *tt=new TChain("tt");
  //1range
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_1.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_2.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_3.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_4.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_5.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_6.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_7.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_8.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_9.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_9_1.root");

  //2range
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_10.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_11.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_12.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_13.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_14.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_15.root");

  //3 range
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_16.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_17.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_18.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_19.root");

  //4 range
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_1.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_2.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_3.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_4.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_5.root");

  //5 range
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_6.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_7.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_1.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_2.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_3.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_4.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_5.root");

  //6 range
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_6.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_7.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_1.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_2.root");

  //7 range
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_4.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_5.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_1.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_2.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_1.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_2.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_3.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_4.root");

  //8 range
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_1.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_2.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_1.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_2.root");

  //9 range
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_1.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_2.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_3.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_4.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_5.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_1.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_2.root");
  
  //10 range
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_1.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_2.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_3.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_4.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_5.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_1.root");
  
  //11 range
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_2.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_3.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_4.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_1.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_2.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_3.root");

  //12 range
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_1.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_2.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_3.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_1.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_2.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_3.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_4.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_5.root");

  //13 range
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_6.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_7.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_1.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_2.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_3.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_5.root");

  //14 range
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_1.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_2.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_3.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_4.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_1.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_2.root");

  //15 range
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_3.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_4.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_5.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_6.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_1.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_2.root");
  tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_3.root");


  TChain *tt2=new TChain("emct");
  //1range
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

  //2range
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_10.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_11.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_12.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_13.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_14.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_15.root");

  //3 range
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_16.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_17.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_18.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_19.root");

  //4 range
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_5.root");

  //5 range
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_5.root");

  //6 range
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_2.root");

  //7 range
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_4.root");

  //8 range
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_2.root");

  //9 range
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_2.root");
  
  //10 range
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_1.root");
  
  //11 range
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_3.root");

  //12 range
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_5.root");

  //13 range
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_7.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_5.root");

  //14 range
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_2.root");

  //15 range
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_3.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_4.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_5.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_6.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_1.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_2.root");
  tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_3.root");

  TChain *tt3=new TChain("tt");
  //1range
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

  //2range
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_10.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_11.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_12.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_13.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_14.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_15.root");

  //3 range
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_16.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_17.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_18.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_19.root");

  //4 range
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_5.root");

  //5 range
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_5.root");

  //6 range
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_2.root");

  //7 range
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_4.root");

  //8 range
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_2.root");

  //9 range
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_2.root");
  
  //10 range
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_1.root");
  
  //11 range
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_3.root");

  //12 range
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_5.root");

  //13 range
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_7.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_5.root");

  //14 range
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_2.root");

  //15 range
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_3.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_4.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_5.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_6.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_1.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_2.root");
  tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_3.root");


  TCanvas *c11 = new TCanvas();
  c11->cd();
  gROOT->SetStyle("Plain");
  gStyle->SetTimeOffset(0);
  //gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  gStyle->SetGridStyle(1);
  gStyle->SetGridWidth(1);
  gStyle->SetGridColor(kGray);
  //TH2F *h31 = new TH2F("h2", "h2 title", 40, -26, 28, 40, -0.18,0.18);

	Int_t nentr=tt->GetEntries();               //СЏв”¤Рїв•¦СЏв”‚Рїв•©Рїв•¬ СЏв”‚Рїв•¬Рїв• СЏв–ЂСЏв”ЊРїв•¦Рїв•§ Рїв•Ў Рїв•ўРїв•ЈСЏв”ЂРїв•ЈРїв•ЎРїв•Ј
	cout<<"Nentries="<<nentr<<endl;

	sprintf(branchname2,"t");
	sprintf(branchname3,"tof");
	sprintf(branchname4,"emc");
	sprintf(branchname7,"mu");
	sprintf(branchname8,"ev");

	tt->SetBranchStatus("*",0);
  /*tt->SetBranchStatus(branchname2,1);
	tt->SetBranchStatus(branchname3,1);
	tt->SetBranchStatus(branchname4,1);
	tt->SetBranchStatus(branchname7,1);
	tt->SetBranchStatus(branchname8,1);
  */
	tt->SetBranchStatus("t",1);
	tt->SetBranchStatus("tof",1);
	tt->SetBranchStatus("emc",1);
	tt->SetBranchStatus("mu",1);
	tt->SetBranchStatus("ev",1);
  /*
	tt->SetBranchAddress(branchname2,&t);              //СЏв”ђСЏв”‚СЏв”ЊРїв•џРїв•«Рїв•џРїв•ЎРїв•©Рїв•¦Рїв•ЎРїв•џРїв•ЈРїв•Є 
	tt->SetBranchAddress(branchname3,&tof); 
	tt->SetBranchAddress(branchname4,&emc);  
	tt->SetBranchAddress(branchname7,&mu); 
	tt->SetBranchAddress(branchname8,&ev); 
  */
	tt->SetBranchAddress("t",&t);              //СЏв”ђСЏв”‚СЏв”ЊРїв•џРїв•«Рїв•џРїв•ЎРїв•©Рїв•¦Рїв•ЎРїв•џРїв•ЈРїв•Є 
	tt->SetBranchAddress("tof",&tof); 
	tt->SetBranchAddress("emc",&emc);  
	tt->SetBranchAddress("mu",&mu); 
	tt->SetBranchAddress("ev",&ev); 

	sprintf(branchname5,"c0");
	sprintf(branchname6,"c1");
	tt2->SetBranchStatus("*",0);
  /*	
	tt2->SetBranchStatus(branchname5,1);
	tt2->SetBranchStatus(branchname6,1);

	tt2->SetBranchAddress(branchname5,&c0);  
	tt2->SetBranchAddress(branchname6,&c1); 
  */
	tt2->SetBranchStatus("c0",1);
	tt2->SetBranchStatus("c1",1);

	tt2->SetBranchAddress("c0",&c0);  
	tt2->SetBranchAddress("c1",&c1); 

  int next_run=0;
  int ttime;
  float day=31.;
  float kx1,kx2,kx3,kx4;
  float eff;
  double num_npezero=0;
  double num_npenotzero=0;
  double num_npetotal=0;      
  
  time_t rawtime;
  time ( &rawtime ); // текущая дата в секундах
  //std::cout << "time="<<  ctime (&rawtime); 
  std::cout << "time="<<  rawtime <<"\n"<<endl; 

  float Npe_ref=0.0;
  int jj=0;
  float ref_sum_npe=0.0;

  TFile *fout=0;
  //fout = new TFile(TString::Format("npe_time_80-160cnt_active_region20.root").Data(),"RECREATE");
  fout = new TFile(TString::Format("thick_npe_time_0-80cnt_active_region20.root").Data(),"RECREATE");
  //fout = new TFile(TString::Format("thick_npe_time_0-80cnt_aerogel_region20.root").Data(),"RECREATE");
  
  //TProfile* pr=new TProfile(branchname,branchname,150,19635,20000,0,200);        //z
  //TProfile* pr=new TProfile("npe&time","npe&time",200,0,25,0,200);        //z
  TProfile* pr=new TProfile("npe&time","npe&time",1000,1400778000,1465603200,0,200);        //z
  //TProfile* pr=new TProfile(branchname,branchname,50,19600,20850,0,200);        //z
  //TProfile* pr=new TProfile(branchname,branchname,1250,19600,20850,0,200);        //z
  //TProfile* pr2=new TProfile("eff&time","eff&time",200,0,25,0,100);        //z
  TProfile* pr2=new TProfile("eff&time","eff&time",1000,1400778000,1465603200,0,100);        //z
  TProfile* pr3=new TProfile("npe/nperef&time","npe/nperef&time",1000,1400778000,1465603200,0,3.4);        //z
  std::vector<double> vec1;
  std::vector<double> vec_time;
  
    for(int k=0; k<nentr; k++)                               
    //for(int k=0; k<600000; k++)                               
    {
      tt->GetEntry(k);

      if( t.p>1200 && t.chi2<20 && sqrt(pow(t.x0,2)+pow(t.y0,2)+pow(t.z0,2))<40 && t.nvec>4 && mu.status>=1 && mu.nhits>5 && mu.layer>=1 && tof.nhits>=1 && emc.energy>100 && emc.ncls_trk>=2 && emc.ncls>=2 && t.emc_ncls>=2 && t.emc_ncls<=4 && tof.nhits<=8 ) 
    {
      int j=0;
      float sum_npe=0;

  	int counter=-2;
  	//int end_cnt=25;

      //for ( counter; counter<end_cnt; counter++)
      //{ 
      int counter2, counter3;
            
      for (int n=0; n<7; n++)
      {        
          //printf("k= %d\r",k); fflush(stdout);        		 
    	  tt3->SetBranchStatus("*",0);
	  sprintf(branchname,"cnt%d",n);
	  tt3->SetBranchStatus(branchname,1);
	  tt3->SetBranchAddress(branchname,&cnt);              
    	  //tt3->SetBranchStatus("cnt%d",n,1);
    	  //tt3->SetBranchAddress("cnt%d",n,&cnt); 
    	  tt3->GetEntry(k);              

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
 
      	  //if( (cnt.i==counter || cnt.i==counter2 || cnt.i==counter3) && cnt.aerogel_region20!=0 && cnt.wlshit!=1 && cnt.nearwls!=1 ) 
      	  if( (cnt.i==counter || cnt.i==counter2 || cnt.i==counter3) && cnt.active_region20!=0 && cnt.wlshit==1 && cnt.nearwls==1 ) 
          {         
	       //tt2->GetEntry(k);
	       //if( c0.e>50 && c1.e>50 )
               //{	
                  j++;   //считаем число перечеченных счетчиков и затем используем это условие
		  sum_npe+=cnt.npe;
               
                  if( j>1 )
       		  {
       		  printf("k= %d\r",k); fflush(stdout); 
       		  
                  if( next_run!=ev.run )
                  {	
                     //cout<<ev.run<<"\t"<<endl;
                     //of<<ev.run<<"\t"<<endl;
	             ifstream f_in("list_of_run_time.dat");
	             string s;	
	             while( getline(f_in, s) )
                     {
    		        istringstream iss(s);
    		        iss >> kx1 >> kx2 >> kx3 >> kx4;
    		        if (kx1==ev.run) 
                        {
                       	  struct tm tm;
			  time_t ts = 0;
  			  memset(&tm, 0, sizeof(tm));

    			  //strptime("1998-04-11", "%Y-%m-%d", &tm);
    			  tm.tm_year   = kx4 - 1900;
    			  tm.tm_mon    = kx2 - 1;    //months since January - [0,11]
    			  tm.tm_mday   = kx3;
    			  ts = mktime(&tm);

    			  printf("%d \n", (int)ts); //unix time-stamp
    			  printf("%s \n", ctime(&ts)); //human readable date
  
                          ttime=ts;     
                          /*  
    		          if((kx2==6 || kx2==9 || kx2==11) && kx4==2014){ day=30.;} 
    		          if(kx2==2 && 2015){ day=28.;} 
    		          if((kx2==4 || kx2==6 || kx2==9 || kx2==11) && (kx4==2015 || 2016)){ day=30.;} 
    		          if(kx2==2 && 2016){ day=29.;}   		        
    		       	  ttime=kx2+kx3/day-5.74194;
    		       	  if(kx4==2015){ttime=ttime+12.;}
    		       	  if(kx4==2016){ttime=ttime+24.;}
                          */    		       		
			  //cout<<kx1<<"\t"<<kx2<<"\t"<<kx3<<"\t"<<kx4<<"\t"<<ttime<<endl;
    		          break; 
                       }
	            }    
	            f_in.close();
            
            	    eff=0;
            	    num_npezero=0;
            	    num_npenotzero=0;
            	    num_npetotal=0;      
		 }
          	 next_run=ev.run;
          
          	 pr->Fill(ttime,sum_npe);
          	 //pr->Fill(ev.run,sum_npe);

          	 //if(sum_npe==0)
          	 //{
            	 //	num_npezero=++num_npezero;
          	 //} 
          	 if(sum_npe>0.3){num_npenotzero=++num_npenotzero;}
          	 num_npetotal=++num_npetotal;
          	 eff=num_npenotzero/num_npetotal;
          	 //cout<<fixed<<setprecision(10)<<eff<<"\t"<<num_npetotal<<"\t"<<num_npezero<<"\t"<<num_npenotzero<<"\t"<<ttime<<endl;
          	 pr2->Fill(ttime,eff*100);
          	 //sum_npe+=cnt.npe;
          	 if(sum_npe>0&&jj<20000)
          	 {
          	     //if(jj<20000){
          	     jj++;
          	     ref_sum_npe+=sum_npe;
          	     Npe_ref=ref_sum_npe/jj;         
          	     //cout<<jj<<"\t"<<"cnt.npe="<<cnt.npe<<"\t"<<"Npe_ref="<<Npe_ref<<"\t"<<cnt.npe/Npe_ref<<"\t"<<"ref_sum_npe="<<ref_sum_npe<<endl;
          	} 
          	pr3->Fill(ttime,sum_npe/Npe_ref);	      
	
		}	      
	      //}
	    }
	    if(n>0 && j>1)
      	    {
       		break;
      	    }
      	   }         
    //}
   }  
  }
 //}    
  TF1* myfit=new TF1("myfit","[0]*exp(-x/[1])+[2]",1400778000,1465603200);

  float npe0=pr->GetBinContent(1);
  cout<<"npe0="<<npe0<<endl;

  myfit->SetParameter(1,npe0);   

  myfit->SetParLimits(0,npe0,15);                        
  myfit->SetParLimits(1,0.,10000.);                        
  myfit->SetParLimits(2,0.,10000.);                         

  myfit->SetLineColor(kBlue);

  //myfit->SetParameter(0,"tau",1.,0.1,0,0);
  //myfit->SetParameter(1,"off",0.5,0.1,0,0);

  pr->Fit("myfit","","",1400778000,1465603200);
  //cout<<"pr->GetMinimum()="<<pr->GetMinimum()<<"\t"<<"pr->GetMaximum()"<<pr->GetMaximum()<<endl;
  pr->SetMarkerStyle(20);
  pr->SetMarkerSize(0.5);
  pr->SetMarkerColor(2);
  pr->SetLineWidth(2);
  pr->SetLineColor(2);
  pr->GetXaxis()->SetTimeDisplay(1);
  pr->GetXaxis()->SetTimeFormat("%d/%m/%y%F1970-01-01 00:00:00");
  pr->GetXaxis()->SetTitleSize(0.05);
  pr->GetXaxis()->SetTitleOffset(1.0);
  pr->GetXaxis()->SetLabelSize(0.05);
  pr->GetXaxis()->SetNdivisions(205);
  pr->GetYaxis()->SetTitleSize(0.05);
  pr->GetYaxis()->SetTitleOffset(1.00);
  pr->GetYaxis()->SetLabelSize(0.05);
  pr->GetYaxis()->SetNdivisions(205);
  pr->GetYaxis()->SetDecimals();
  pr->GetXaxis()->SetTitle("Time(dd/mm/yy)");
  pr->GetYaxis()->SetTitle("N_{ph.e.}");
  pr->SetTitle("");
 //
  pr2->SetMarkerStyle(20);
  pr2->SetMarkerSize(0.5);
  pr2->SetMarkerColor(4);
  pr2->SetLineWidth(2);
  pr2->SetLineColor(4);
  pr2->GetXaxis()->SetTimeDisplay(1);
  pr2->GetXaxis()->SetTimeFormat("%d/%m/%y%F1970-01-01 00:00:00");
  pr2->GetXaxis()->SetTitleSize(0.05);
  pr2->GetXaxis()->SetTitleOffset(1.0);
  pr2->GetXaxis()->SetLabelSize(0.05);
  pr2->GetXaxis()->SetNdivisions(205);
  pr2->GetYaxis()->SetTitleSize(0.05);
  pr2->GetYaxis()->SetTitleOffset(1.00);
  pr2->GetYaxis()->SetLabelSize(0.05);
  pr2->GetYaxis()->SetNdivisions(205);
  pr2->GetYaxis()->SetDecimals();
  pr2->GetXaxis()->SetTitle("Time(dd/mm/yy)");
  pr2->GetYaxis()->SetTitle("Efficiency, %");
  pr2->SetTitle("");
  
  pr3->SetMarkerStyle(20);
  pr3->SetMarkerSize(0.5);
  pr3->SetMarkerColor(2);
  pr3->SetLineWidth(2);
  pr3->SetLineColor(2);
  pr3->GetXaxis()->SetTimeDisplay(1);
  pr3->GetXaxis()->SetTimeFormat("%d/%m/%y%F1970-01-01 00:00:00");
  pr3->GetXaxis()->SetTitleSize(0.05);
  pr3->GetXaxis()->SetTitleOffset(1.0);
  pr3->GetXaxis()->SetLabelSize(0.05);
  pr3->GetXaxis()->SetNdivisions(205);
  pr3->GetYaxis()->SetTitleSize(0.05);
  pr3->GetYaxis()->SetTitleOffset(1.00);
  pr3->GetYaxis()->SetLabelSize(0.05);
  pr3->GetYaxis()->SetNdivisions(205);
  pr3->GetYaxis()->SetDecimals();
  pr3->GetXaxis()->SetTitle("Time(dd/mm/yy)");
  pr3->GetYaxis()->SetTitle("N_{ph.e.}/Nref_{ph.e.}");
  pr3->SetTitle("");

  //pr->SetLineColor(5);         //blue
  pr->Draw("prof");
/*  
  TGraph *gr1;
  gr1=new TGraph(vec_time.size(),&vec_time[0],&vec1[0]);
  gr1->SetTitle("gr1");
  gr1->SetName("gr1");
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(0.5);
  gr1->SetMarkerColor(2);
  gr1->SetLineWidth(2);
  gr1->SetLineColor(2);
  gr1->GetXaxis()->SetTimeDisplay(1);
  gr1->GetXaxis()->SetTimeFormat("%d/%m/%y%F1970-01-01 00:00:00");
  gr1->GetXaxis()->SetTitleSize(0.05);
  gr1->GetXaxis()->SetTitleOffset(0.95);
  gr1->GetXaxis()->SetLabelSize(0.05);
  gr1->GetXaxis()->SetNdivisions(205);
  gr1->GetYaxis()->SetTitleSize(0.05);
  gr1->GetYaxis()->SetTitleOffset(1.00);
  gr1->GetYaxis()->SetLabelSize(0.05);
  gr1->GetYaxis()->SetDecimals();
  gr1->GetYaxis()->SetNdivisions(205);
  gr1->GetXaxis()->SetTitle("Time(dd/mm/yy)");
  gr1->GetYaxis()->SetTitle("N_{ph.e.}/Nref_{ph.e.}");
  gr1->SetTitle("");
  gr1->Write();
*/
  //of.close();
  fout->Write();                                                //ÐÉÛÅÍ ÄÁÎÎÙÅ × ÆÁÊÌ É ÚÁËÒÙ×ÁÅÍ ÅÇÏ
  fout->Close();
} 


