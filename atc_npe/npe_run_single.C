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
/*
int Usage(string status)
{
       cout<<"Usage: "<<progname<<"\t"<<"First cnt  End cnt"<<endl;
       exit(0);
}
*/
//void npe_run(int counter)
//void npe_run_single()
int main(int argc, char* argv[])
{

  // progname=argv[0];

   int first_cnt=0;
   int end_cnt=160;
   //int first_cnt=20;
   //int end_cnt=21;

/*    if( argc>1 )
    {
      first_cnt=atoi(argv[1]);  
      end_cnt=atoi(argv[2]);
      if(first_cnt<0 || first_cnt>160){ Usage(progname); return 0;}  
      if(end_cnt<0 || end_cnt>160){ Usage(progname); return 0;}
    }
    else
    { 
      Usage(progname);
    }
*/    
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
  char name0[161];
  char name1[161];
  char name2[161];
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
    //16range 
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march17_1.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march17_2.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march17_3.root");
    //17range 
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr17_1.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr17_2.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr17_3.root");
    //18range
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may17_1.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may17_2.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may17_3.root");
    //19range
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_1.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_2.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_3.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_4.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_5.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_6.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_7.root");
    //20range
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_1.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_2.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_3.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_4.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_5.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_6.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_7.root");
    //21range
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_1.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_2.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_3.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_4.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_5.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_6.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_7.root");
    //22range
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_1.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_2.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_3.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_4.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_5.root");
    //23range
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_1.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_2.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_3.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_4.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_5.root");
    //24range
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_1.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_2.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_3.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_4.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_5.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_6.root");
    //25range
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_1.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_2.root");
    tt->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_3.root");
    
    TChain *tt2=new TChain("tt");
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
    //16range 
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march17_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march17_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march17_3.root");
    //17range 
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr17_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr17_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr17_3.root");
    //18range
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may17_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may17_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may17_3.root");
    //19range
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_3.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_4.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_5.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_6.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_7.root");
    //20range
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_3.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_4.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_5.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_6.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_7.root");
    //21range
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_3.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_4.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_5.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_6.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_7.root");
    //22range
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_3.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_4.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_5.root");
    //23range
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_3.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_4.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_5.root");
    //24range
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_3.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_4.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_5.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_6.root");
    //25range
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_1.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_2.root");
    tt2->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_3.root");


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
    //16range 
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march17_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march17_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march17_3.root");
    //17range 
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr17_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr17_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr17_3.root");
    //18range
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may17_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may17_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may17_3.root");
    //19range
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_3.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_4.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_5.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_6.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_7.root");
    //20range
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_3.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_4.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_5.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_6.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_7.root");
    //21range
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_3.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_4.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_5.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_6.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_7.root");
    //22range
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_3.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_4.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_5.root");
    //23range
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_3.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_4.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_5.root");
    //24range
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_3.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_4.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_5.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_6.root");
    //25range
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_1.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_2.root");
    tt3->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_3.root");
    
    TChain *tt4=new TChain("tt");
    //1range
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_3.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_4.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_5.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_6.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_7.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_8.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_9.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_9_1.root");

    //2range
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_10.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_11.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_12.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_13.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_14.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_15.root");

    //3 range
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_16.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_17.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_18.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_19.root");

    //4 range
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_3.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_4.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_5.root");

    //5 range
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_6.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_7.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_3.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_4.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_5.root");

    //6 range
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_6.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_7.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_2.root");

    //7 range
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_4.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_5.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_3.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_4.root");

    //8 range
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_2.root");

     //9 range
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_3.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_4.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_5.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_2.root");
  
    //10 range
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_3.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_4.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_5.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_1.root");
  
    //11 range
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_3.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_4.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_3.root");

    //12 range
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_3.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_3.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_4.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_5.root");

    //13 range
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_6.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_7.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_3.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_5.root");

    //14 range
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_3.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_4.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_2.root");

    //15 range
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_3.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_4.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_5.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_6.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_3.root");
    //16range 
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march17_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march17_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march17_3.root");
    //17range 
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr17_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr17_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr17_3.root");
    //18range
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may17_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may17_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may17_3.root");
    //19range
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_3.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_4.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_5.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_6.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_7.root");
    //20range
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_3.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_4.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_5.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_6.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_7.root");
    //21range
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_3.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_4.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_5.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_6.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_7.root");
    //22range
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_3.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_4.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_5.root");
    //23range
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_3.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_4.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_5.root");
    //24range
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_3.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_4.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_5.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_6.root");
    //25range
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_1.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_2.root");
    tt4->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_3.root");
    
    TChain *tt5=new TChain("emct");
    //1range
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_3.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_4.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_5.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_6.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_7.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_8.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_9.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_9_1.root");

    //2range
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_10.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_11.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_12.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_13.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_14.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_15.root");

    //3 range
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_16.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_17.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_18.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_19.root");

    //4 range
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_3.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_4.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_5.root");

    //5 range
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_6.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_7.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_3.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_4.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_5.root");

    //6 range
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_6.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_7.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_2.root");

    //7 range
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_4.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_5.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_3.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_4.root");

    //8 range
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_2.root");

     //9 range
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_3.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_4.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_5.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_2.root");
  
    //10 range
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_3.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_4.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_5.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_1.root");
  
    //11 range
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_3.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_4.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_3.root");

    //12 range
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_3.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_3.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_4.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_5.root");

    //13 range
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_6.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_7.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_3.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_5.root");

    //14 range
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_3.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_4.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_2.root");

    //15 range
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_3.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_4.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_5.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_6.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_3.root");
    //16range 
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march17_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march17_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march17_3.root");
    //17range 
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr17_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr17_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr17_3.root");
    //18range
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may17_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may17_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may17_3.root");
    //19range
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_3.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_4.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_5.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_6.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_7.root");
    //20range
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_3.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_4.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_5.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_6.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_7.root");
    //21range
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_3.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_4.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_5.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_6.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_7.root");
    //22range
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_3.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_4.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_5.root");
    //23range
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_3.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_4.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_5.root");
    //24range
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_3.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_4.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_5.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_6.root");
    //25range
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_1.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_2.root");
    tt5->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_3.root");
    
    TChain *tt6=new TChain("emct");
    //1range
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_3.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_4.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_5.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_6.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_7.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_8.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_9.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_9_1.root");

    //2range
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_10.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_11.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_12.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_13.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_14.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_15.root");

    //3 range
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_16.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_17.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_18.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_19.root");

    //4 range
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_3.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_4.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_5.root");

    //5 range
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_6.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_7.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_3.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_4.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_5.root");

    //6 range
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_6.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_7.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_2.root");

    //7 range
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_4.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_5.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_3.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_4.root");

    //8 range
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_2.root");

     //9 range
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_3.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_4.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_5.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_2.root");
  
    //10 range
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_3.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_4.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_5.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_1.root");
  
    //11 range
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_3.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_4.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_3.root");

    //12 range
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_3.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_3.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_4.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_5.root");

    //13 range
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_6.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_7.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_3.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_5.root");

    //14 range
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_3.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_4.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_2.root");

    //15 range
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_3.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_4.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_5.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_6.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_3.root");
    //16range 
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march17_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march17_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march17_3.root");
    //17range 
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr17_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr17_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr17_3.root");
    //18range
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may17_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may17_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may17_3.root");
    //19range
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_3.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_4.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_5.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_6.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_7.root");
    //20range
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_3.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_4.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_5.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_6.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_7.root");
    //21range
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_3.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_4.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_5.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_6.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_7.root");
    //22range
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_3.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_4.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_5.root");
    //23range
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_3.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_4.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_5.root");
    //24range
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_3.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_4.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_5.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_6.root");
    //25range
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_1.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_2.root");
    tt6->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_3.root");
    
    TChain *tt7=new TChain("tt");
    //1range
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_3.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_4.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_5.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_6.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_7.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_8.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_9.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_9_1.root");

    //2range
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_10.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_11.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_12.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_13.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_14.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_15.root");

    //3 range
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_16.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_17.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_18.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_19.root");

    //4 range
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_3.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_4.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_5.root");

    //5 range
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_6.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_7.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_3.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_4.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_5.root");

    //6 range
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_6.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_7.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_2.root");

    //7 range
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_4.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_5.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_3.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_4.root");

    //8 range
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_2.root");

     //9 range
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_3.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_4.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_5.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_2.root");
  
    //10 range
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_3.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_4.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_5.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_1.root");
  
    //11 range
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_3.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_4.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_3.root");

    //12 range
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_3.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_3.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_4.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_5.root");

    //13 range
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_6.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_7.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_3.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_5.root");

    //14 range
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_3.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_4.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_2.root");

    //15 range
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_3.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_4.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_5.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_6.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_3.root");
    //16range 
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march17_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march17_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march17_3.root");
    //17range 
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr17_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr17_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr17_3.root");
    //18range
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may17_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may17_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may17_3.root");
    //19range
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_3.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_4.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_5.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_6.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_7.root");
    //20range
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_3.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_4.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_5.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_6.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_7.root");
    //21range
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_3.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_4.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_5.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_6.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_7.root");
    //22range
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_3.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_4.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_5.root");
    //23range
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_3.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_4.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_5.root");
    //24range
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_3.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_4.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_5.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_6.root");
    //25range
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_1.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_2.root");
    tt7->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_3.root");
    
    TChain *tt8=new TChain("tt");
    //1range
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_3.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_4.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_5.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_6.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_7.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_8.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_9.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_9_1.root");

    //2range
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_10.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_11.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_12.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_13.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_14.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_15.root");

    //3 range
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_16.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_17.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_18.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun14_19.root");

    //4 range
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_3.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_4.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_5.root");

    //5 range
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_6.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct14_7.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_3.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_4.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_5.root");

    //6 range
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_6.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov14_7.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_2.root");

    //7 range
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_4.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec14_5.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan15_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_3.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb15_4.root");

    //8 range
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr15_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may15_2.root");

     //9 range
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_3.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_4.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun15_5.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jul15_2.root");
  
    //10 range
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_3.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_4.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_oct15_5.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_1.root");
  
    //11 range
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_3.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_nov15_4.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec15_3.root");

    //12 range
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan16_3.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_3.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_4.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_5.root");

    //13 range
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_6.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb16_7.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_3.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march16_5.root");

    //14 range
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_3.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr16_4.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_2.root");

    //15 range
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_3.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_4.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_5.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may16_6.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun16_3.root");
    //16range 
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march17_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march17_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march17_3.root");
    //17range 
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr17_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr17_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr17_3.root");
    //18range
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may17_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may17_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may17_3.root");
    //19range
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_3.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_4.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_5.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_6.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_dec17_7.root");
    //20range
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_3.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_4.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_5.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_6.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jan18_7.root");
    //21range
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_3.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_4.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_5.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_6.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_feb18_7.root");
    //22range
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_3.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_4.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_march18_5.root");
    //23range
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_3.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_4.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_apr18_5.root");
    //24range
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_3.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_4.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_5.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_may18_6.root");
    //25range
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_1.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_2.root");
    tt8->Add("/home/ovtin/development/work_KEDR/data_cosm/cosm_runs_jun18_3.root");

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

	sprintf(branchname8,"ev");
	tt8->SetBranchStatus("*",0);
	tt8->SetBranchStatus(branchname8,1);
	tt8->SetBranchAddress(branchname8,&ev); 

  TFile *fout=0;
  //fout = new TFile(TString::Format("npe_time_0-80cnt_active_region20.root").Data(),"RECREATE");
  //fout = new TFile(TString::Format("npe_time_0-80cnt_aerogel_region20.root").Data(),"RECREATE"); 
  // fout = new TFile(TString::Format("npe_time_0-80cnt_aerogel_region.root").Data(),"RECREATE");
  //fout = new TFile(TString::Format("npe_time_0-80cnt_aerogel_region5.root").Data(),"RECREATE");
  //fout = new TFile(TString::Format("npe_time_0-80cnt_aerogel_region0.root").Data(),"RECREATE");
  
  fout = new TFile(TString::Format("/home/ovtin/development/work_KEDR/KEDRATC/atc_npe/npe_time_aerogel.root").Data(),"RECREATE");
  //fout = new TFile(TString::Format("/home/ovtin/development/atc_npe/npe_time_shifter.root").Data(),"RECREATE");
  //fout = new TFile(TString::Format("npe_time_80-160cnt_active_region0.root").Data(),"RECREATE");
  //fout = new TFile(TString::Format("npe_time_80-160cnt_active_region20.root").Data(),"RECREATE");
  //fout = new TFile(TString::Format("npe_time_80-160cnt_aerogel_region0.root").Data(),"RECREATE");
  //fout = new TFile(TString::Format("npe_time_80-160cnt_aerogel_region5.root").Data(),"RECREATE");
  //fout = new TFile(TString::Format("npe_time_20-60cnt_aerogel_region.root").Data(),"RECREATE");
  //fout = new TFile(TString::Format("npe_time_80-160cnt_aerogel_region20_test.root").Data(),"RECREATE");
  
  //ofstream of("list_of_run.dat",ios::app);
  //for (int n=counter; n<counter+1; n++)
  
  ofstream Result(TString::Format("/home/ovtin/development/work_KEDR/KEDRATC/atc_npe/result_cnt_npe_degr_aerogel.dat"),ios_base::out); 
  //ofstream Result(TString::Format("/home/ovtin/development/atc_npe/result_cnt_npe_degr_shifter.dat"),ios_base::out); 


  for (int j=first_cnt; j<end_cnt; j++)
 {
	sprintf(name0,"Cnt%d_Nphe",j);       
	sprintf(name1,"Cnt%d_Eff",j);       
	sprintf(name2,"Cnt%d_Nphenorm",j);       
	//sprintf(branchname10,"N_{ph.e.}&#beta#gamma cnt%d ",j);       
	//TProfile* pr=new TProfile(branchname,branchname,100,0,momentum,0,20); 
	//TProfile* pr=new TProfile(branchname,branchname,300,0,momentum,0,20); 
	//TProfile* pr2=new TProfile(branchname10,branchname10,100,0,15,0,30); 
	cout<<"Counter="<<j<<endl;
	
	int next_run=0;
  	int ttime;
  	float day=31.;
  	float kx1,kx2,kx3,kx4,kx5,kx6,kx7;
  	float eff;
  	double num_npezero=0;
  	double num_npenotzero=0;
  	double num_npetotal=0;      
  
  	time_t rawtime;
  	time ( &rawtime ); // текущая дата в секундах
  	//std::cout << "time="<<  ctime (&rawtime); 
  	std::cout << "time="<<  rawtime <<"\n"<<endl; 

  	float Npe_ref=1.0;
  	int jj=0;
  	float sum_npe=0.0;

  	//TProfile* pr=new TProfile(branchname,branchname,150,19635,20000,0,200);        //z
  	//TProfile* pr=new TProfile("npe&time","npe&time",200,0,25,0,200);        //z
  	TProfile* pr=new TProfile(name0,name0,3000,1400778000,1529296837,0,200);        //z
  	//TProfile* pr=new TProfile("Shifter - 20 mm","npe&time",1000,1400778000,1465603200,0,200);        //z
  	//TProfile* pr=new TProfile(branchname,branchname,50,19600,20850,0,200);        //z
  	//TProfile* pr=new TProfile(branchname,branchname,1250,19600,20850,0,200);        //z
  	//TProfile* pr2=new TProfile("eff&time","eff&time",200,0,25,0,100);        //z
  	TProfile* pr2=new TProfile(name1,name1,3000,1400778000,1529296837,0,100);        //z
  	TProfile* pr3=new TProfile(name2,name2,3000,1400778000,1529296837,0,3.0);        //z
  	//std::vector<double> vec1;
  	//std::vector<double> vec_time;
  	
        ofstream myFile(TString::Format("/home/ovtin/development/work_KEDR/KEDRATC/atc_npe/data_aerogel_cnt_%d.dat",j),ios_base::out); 
        //ofstream myFile(TString::Format("/home/ovtin/development/atc_npe/data_shifter_cnt_%d.dat",j),ios_base::out); 

  //for (int n=0; n<1; n++)
  for (int n=0; n<7; n++)
  //for (int n=0; n<9; n++)
  {        
    //if(n<20 || n>=60 && n<100 || n>=140)
    //{
	  sprintf(branchname,"cnt%d",n);
    	  tt->SetBranchStatus("*",0);
	  tt->SetBranchStatus(branchname,1);
	  tt->SetBranchAddress(branchname,&cnt);              

    	  //tt->SetBranchStatus("cnt%d",n,1);
    	  //tt->SetBranchAddress("cnt%d",n,&cnt);              	
	  cout<<"Сrossing="<<n<<endl;

    for(int k=0; k<nentr; k++)                               
    //for(int k=0; k<1000000; k++)                               
    {
        tt->GetEntry(k);
        
	//if( cnt.i>=first_cnt && cnt.i<end_cnt )  //цикл по числу счетчиков
	//if( (cnt.i>=0 && cnt.i<20) || (cnt.i>=60 && cnt.i<80) )  //цикл по числу счетчиков
	//if( (cnt.i>=80 && cnt.i<100) || (cnt.i>=140 && cnt.i<160) )  //цикл по числу счетчиков
        if( cnt.i==j && cnt.aerogel_region!=0 && cnt.wlshit!=1 && cnt.nearwls!=1 ){    //aerogel
        //if( cnt.i==j && cnt.active_region0!=0 && cnt.wlshit!=0 && cnt.nearwls!=0 ){  //shifter
	 tt2->GetEntry(k);      
 	//if( t.p>1200 && t.chi2<20 && sqrt(pow(t.x0,2)+pow(t.y0,2)+pow(t.z0,2))<40 && t.nvec>4 && t.emc_ncls>=2 && t.emc_ncls<=4 )
 	//if( t.p>1200 && t.chi2<120 && sqrt(pow(t.x0,2)+pow(t.y0,2)+pow(t.z0,2))<80 && t.nvec>4 && t.emc_ncls>=2 && t.emc_ncls<=4 )
 	if( t.p>1200 && t.chi2<150 && t.nvec>4 && t.emc_ncls>=2 && t.emc_ncls<=4 )
	{
         tt3->GetEntry(k);
  	if( tof.nhits>=1 && tof.nhits<=8 ) 
        {
         tt7->GetEntry(k);      
	 //cout<<"cnt.npe="<<cnt.npe<<endl;
        if( mu.status>=1 && mu.nhits>5 && mu.layer>=1 ) 
        {  
         tt4->GetEntry(k);
	if( emc.energy>100 && emc.ncls_trk>=2 && emc.ncls>=2 ) 
        {  
               tt5->GetEntry(k);
    	       tt6->GetEntry(k);
               tt8->GetEntry(k);
	      if( c0.e>50 && c1.e>50 )
        {	
          if(n>0){cout<<"ev.run="<<ev.run<<endl;} 
          if( next_run!=ev.run )
          {	
            //cout<<ev.run<<"\t"<<endl;
            //of<<ev.run<<"\t"<<endl;
	         ifstream f_in("/home/ovtin/development/work_KEDR/KEDRATC/atc_npe/list_of_run_time.dat");
	         string s;	
	          while( getline(f_in, s) )
            {
    		      istringstream iss(s);
    		      iss >> kx1 >> kx2 >> kx3 >> kx4 >>kx5 >>kx6 >>kx7;
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
          
          	pr->Fill(ttime,cnt.npe);
          	//pr->Fill(ev.run,cnt.npe);
          	myFile<< ttime <<"\t"<< cnt.npe <<endl;

          	if(cnt.npe==0)
          	{
            		num_npezero=++num_npezero;
          	} 
          	if(cnt.npe>0){num_npenotzero=++num_npenotzero;}
          	num_npetotal=++num_npetotal;
          	eff=num_npenotzero/num_npetotal;
          	//cout<<fixed<<setprecision(10)<<eff<<"\t"<<num_npetotal<<"\t"<<num_npezero<<"\t"<<num_npenotzero<<"\t"<<ttime<<endl;
          	pr2->Fill(ttime,eff*100);
          	//sum_npe+=cnt.npe;
          	if(cnt.npe>0&&jj<20000){
          	//if(jj<20000){
          		jj++;
          		sum_npe+=cnt.npe;
          		Npe_ref=sum_npe/jj;         
          	//cout<<jj<<"\t"<<"cnt.npe="<<cnt.npe<<"\t"<<"Npe_ref="<<Npe_ref<<"\t"<<cnt.npe/Npe_ref<<"\t"<<"sum_npe="<<sum_npe<<endl;
          	} 
          	//cout<<jj<<"\t"<<"cnt.npe="<<cnt.npe<<"\t"<<"Npe_ref="<<Npe_ref<<"\t"<<cnt.npe/Npe_ref<<endl;
          	pr3->Fill(ttime,cnt.npe/Npe_ref);	      
//		vec1.push_back(cnt.npe/Npe_ref);
//                vec_time.push_back(ttime);
//          	pr3->Fill(ttime,cnt.npe);	      
	      }
	    }
	   }
	  }
	}
      }
    } 
  }
  myFile.close(); 
      // }    //*************
  
  //float npe_begin=(pr->GetBinContent(1)+pr->GetBinContent(2)+pr->GetBinContent(3)+pr->GetBinContent(4)+pr->GetBinContent(5)+pr->GetBinContent(6)+pr->GetBinContent(7)+pr->GetBinContent(8)+pr->GetBinContent(9)+pr->GetBinContent(10))/10;
  float sum_npe_begin=0;
  int nn=0;
  for(int ii=0; ii<200; ii++)
  {
     //cout<<"pr->GetBinContent(ii)="<<pr->GetBinContent(ii)<<endl;      
     if(pr->GetBinContent(ii)>0)
       {
        sum_npe_begin+=pr->GetBinContent(ii);
        nn++;
       } 
  }
  float npe_begin=0;
  if(nn>0){npe_begin=sum_npe_begin/nn;}
  cout<<"npe_begin="<<npe_begin<<endl;  
  
  int bin=pr->GetNbinsX();  
  //float npe_end=(pr->GetBinContent(bin-10)+pr->GetBinContent(bin-9)+pr->GetBinContent(bin-8)+pr->GetBinContent(bin-7)+pr->GetBinContent(bin-6)+pr->GetBinContent(bin-5)+pr->GetBinContent(bin-4)+pr->GetBinContent(bin-3)+pr->GetBinContent(bin-2)+pr->GetBinContent(bin-1))/10;
  float sum_npe_end=0;
  nn=0;
  for(int ii=0; ii<200; ii++)
  {
     //cout<<"pr->GetBinContent(bin-ii)="<<pr->GetBinContent(bin-ii)<<endl; 
       if(pr->GetBinContent(bin-ii)>0)
       {     
        sum_npe_end+=pr->GetBinContent(bin-ii); 
        nn++;
       } 
  }
  float npe_end=0;
  if(nn>0){npe_end=sum_npe_end/nn;}

  cout<<"npe_end="<<npe_end<<endl;  
  
  float degr=0;
  if(npe_begin>0){degr=(1-(npe_end/npe_begin))*100;}
  
  /*
  TF1* myfit=new TF1("myfit","[0]*exp(-(x-1400778000)/[1])",1400778000,1529296837);
  myfit->SetParLimits(0,0,10000);
  myfit->SetParLimits(1,500000,1e12);
  myfit->SetLineColor(kBlue);
  pr->Fit("myfit","","",1400778000,1529296837);
  */
  //cout<<"pr->GetMinimum()="<<pr->GetMinimum()<<"\t"<<"pr->GetMaximum()"<<pr->GetMaximum()<<endl;
  pr->SetMarkerStyle(20);
  pr->SetMarkerSize(0.5);
  pr->SetMarkerColor(4);
  pr->SetLineWidth(2);
  pr->SetLineColor(4);
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
  
  
  //float npe_begin_n=(pr3->GetBinContent(1)+pr3->GetBinContent(2)+pr3->GetBinContent(3)+pr3->GetBinContent(4)+pr3->GetBinContent(5)+pr3->GetBinContent(6)+pr3->GetBinContent(7)+pr3->GetBinContent(8)+pr3->GetBinContent(9)+pr3->GetBinContent(10))/10;
  //cout<<"npe_begin_n="<<npe_begin_n<<endl;  
  //float npe_end_n=(pr3->GetBinContent(bin-10)+pr3->GetBinContent(bin-9)+pr3->GetBinContent(bin-8)+pr3->GetBinContent(bin-7)+pr3->GetBinContent(bin-6)+pr3->GetBinContent(bin-5)+pr3->GetBinContent(bin-4)+pr3->GetBinContent(bin-3)+pr3->GetBinContent(bin-2)+pr3->GetBinContent(bin-1))/10;  
  float sum_npe_begin_n=0;
  nn=0;
  for(int ii=0; ii<200; ii++)
  {
     //cout<<"pr3->GetBinContent(ii)="<<pr3->GetBinContent(ii)<<endl;  
       if(pr3->GetBinContent(ii)>0)
     {     
     sum_npe_begin_n+=pr3->GetBinContent(ii);
     nn++;
     } 
  }
  float npe_begin_n=0;
  if(nn>0){npe_begin_n=sum_npe_begin_n/nn;}
  cout<<"npe_begin_n="<<npe_begin_n<<endl;  
  
  bin=pr3->GetNbinsX();  
  //float npe_end=(pr->GetBinContent(bin-10)+pr->GetBinContent(bin-9)+pr->GetBinContent(bin-8)+pr->GetBinContent(bin-7)+pr->GetBinContent(bin-6)+pr->GetBinContent(bin-5)+pr->GetBinContent(bin-4)+pr->GetBinContent(bin-3)+pr->GetBinContent(bin-2)+pr->GetBinContent(bin-1))/10;
  float sum_npe_end_n=0;
  nn=0;
  for(int ii=0; ii<200; ii++)
  {
     //cout<<"pr3->GetBinContent(bin-ii)="<<pr3->GetBinContent(bin-ii)<<endl;      
     if(pr3->GetBinContent(bin-ii)>0)
     {     
     sum_npe_end_n+=pr3->GetBinContent(bin-ii); 
     nn++;
     }
  }
  float npe_end_n=0;
  if(nn>0){npe_end_n=sum_npe_end_n/nn;}
  cout<<"npe_end_n="<<npe_end_n<<endl;  

  float degr_n=0;
  if(npe_begin_n>0){degr_n=(1-(npe_end_n/npe_begin_n))*100;}
  
  //Result<< j <<"\t"<< degr <<"\t"<< degr_n <<endl;
  Result<< j <<"\t"<< degr << endl;
  
  /*
  TF1* myfit3=new TF1("myfit3","[0]*exp(-(x-1400778000)/[1])",1400778000,1529296837);
  myfit3->SetParLimits(0,0,10000);
  myfit3->SetParLimits(1,500000,1e12);
  myfit3->SetLineColor(kBlue);
  pr3->Fit("myfit3","","",1400778000,1529296837);
  */

  pr3->SetMarkerStyle(20);
  pr3->SetMarkerSize(0.5);
  pr3->SetMarkerColor(4);
  pr3->SetLineWidth(2);
  pr3->SetLineColor(4);
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
  //pr->Draw("prof");
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
  }
  fout->Write();                                                //ÐÉÛÅÍ ÄÁÎÎÙÅ × ÆÁÊÌ É ÚÁËÒÙ×ÁÅÍ ÅÇÏ
  fout->Close();
} 


