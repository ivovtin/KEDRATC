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

//void npetrh_eff_ineff(int cnt_i)
void npetrh_eff_ineff_thick()
{
  int cnt_i=0;

  int cnt;
  float npetrh;
  float npe;
  float eff_mu, eff_pion, eff_kaon, not_eff_pion, not_eff_kaon, Ksigma;
  int total_entries;
 
  TFile *fout=0;
  TString fname;
  fname=TString::Format("cnt_thick_misindentification_allcnt.root").Data();
  fout=new TFile(fname,"RECREATE");
  cout<<fname<<endl;

   //TCanvas *cc1 = new TCanvas("cc1","cnt_thick_misidentification_eff_ineff_all.root",1600,1200);
   //gStyle->SetOptStat(0);
   //gROOT->SetStyle("Plain");
 
  vector<float> eff_P,ineff_P,eff_K,ineff_K,trh;
  vector<float> eff_P_1st,ineff_P_1st,eff_K_1st,ineff_K_1st,trh_1st;
  vector<float> eff_P_2nd,ineff_P_2nd,eff_K_2nd,ineff_K_2nd,trh_2nd;

  FILE *pF;
   //pF=fopen("thick_npetrh0.0/cnt_thick_20_60_npetrh0.000000_eff_ineff_p10.dat", "r");
   pF=fopen("thick_npetrh0.0/cnt_thick_0_80_npetrh0.000000_eff_ineff_p10_allcnt.dat", "r");
   //for (int k=1; k<41; k++)
   for (int k=1; k<2; k++)
  {
    fscanf(pF," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF); 

   FILE *pF01;
   pF01=fopen("thick_npetrh0.3/cnt_thick_0_80_npetrh0.300000_eff_ineff_p10_allcnt.dat", "r");
   //for (int k=1; k<41; k++)
   for (int k=1; k<2; k++)
  {
    fscanf(pF01," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF01); 

  FILE *pF1;
   //pF1=fopen("thick_npetrh0.5/cnt_thick_20_60_npetrh0.500000_eff_ineff_p10.dat", "r");
   pF1=fopen("thick_npetrh0.5/cnt_thick_0_80_npetrh0.500000_eff_ineff_p10_allcnt.dat", "r");
   for (int k1=1; k1<2; k1++)
  {
    fscanf(pF1," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF1); 

  FILE *pF2;
   //pF2=fopen("thick_npetrh1.0/cnt_thick_20_60_npetrh1.000000_eff_ineff_p10.dat", "r");
   pF2=fopen("thick_npetrh1.0/cnt_thick_0_80_npetrh1.000000_eff_ineff_p10_allcnt.dat", "r");
   for (int k2=1; k2<2; k2++)
  {
    fscanf(pF2," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF2); 

  FILE *pF3;
   //pF3=fopen("thick_npetrh2.0/cnt_thick_20_60_npetrh2.000000_eff_ineff_p10.dat", "r");
   pF3=fopen("thick_npetrh2.0/cnt_thick_0_80_npetrh2.000000_eff_ineff_p10_allcnt.dat", "r");
   for (int k3=1; k3<2; k3++)
  {
    fscanf(pF3," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF3);
    
  FILE *pF4;
   //pF4=fopen("thick_npetrh2.5/cnt_thick_20_60_npetrh2.500000_eff_ineff_p10.dat", "r");
   pF4=fopen("thick_npetrh2.5/cnt_thick_0_80_npetrh2.500000_eff_ineff_p10_allcnt.dat", "r");
   for (int k4=1; k4<2; k4++)
  {
    fscanf(pF4," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF4); 

  FILE *pF5;
   //pF5=fopen("thick_npetrh3.0/cnt_thick_20_60_npetrh3.000000_eff_ineff_p10.dat", "r");
   pF5=fopen("thick_npetrh3.0/cnt_thick_0_80_npetrh3.000000_eff_ineff_p10_allcnt.dat", "r");
   for (int k5=1; k5<2; k5++)
  {
    fscanf(pF5," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF5); 

  FILE *pF6;
   //pF6=fopen("thick_npetrh3.5/cnt_thick_20_60_npetrh3.500000_eff_ineff_p10.dat", "r");
   pF6=fopen("thick_npetrh3.5/cnt_thick_0_80_npetrh3.500000_eff_ineff_p10_allcnt.dat", "r");
   for (int k6=1; k6<2; k6++)
  {
    fscanf(pF6," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF6); 


   FILE *pF7;
   //pF7=fopen("thick_npetrh1.5/cnt_thick_20_60_npetrh1.500000_eff_ineff_p10.dat", "r");
   pF7=fopen("thick_npetrh1.5/cnt_thick_0_80_npetrh1.500000_eff_ineff_p10_allcnt.dat", "r");
   for (int k7=1; k7<2; k7++)
  {
    fscanf(pF7," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF7); 

   
    FILE *pF48;
   pF48=fopen("thick_npetrh4.0/cnt_thick_0_80_npetrh4.000000_eff_ineff_p10_allcnt.dat", "r");
   for (int k48=1; k48<2; k48++)
  {
    fscanf(pF48," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF48); 

     FILE *pF49;
   pF49=fopen("thick_npetrh4.5/cnt_thick_0_80_npetrh4.500000_eff_ineff_p10_allcnt.dat", "r");
   for (int k49=1; k49<2; k49++)
  {
    fscanf(pF49," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF49);


//==========================================================================================================
FILE *pF8;
   //pF8=fopen("thick_npetrh0.0/cnt_thick_20_60_npetrh0.000000_eff_ineff_p11.dat", "r");
   pF8=fopen("thick_npetrh0.0/cnt_thick_0_80_npetrh0.000000_eff_ineff_p11_allcnt.dat", "r");
   for (int k8=1; k8<2; k8++)
  {
    fscanf(pF8," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF8); 

    FILE *pF02;
   pF02=fopen("thick_npetrh0.3/cnt_thick_0_80_npetrh0.300000_eff_ineff_p11_allcnt.dat", "r");
   //for (int k=1; k<41; k++)
   for (int k=1; k<2; k++)
  {
    fscanf(pF02," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF02); 


  FILE *pF9;
   //pF9=fopen("thick_npetrh0.5/cnt_thick_20_60_npetrh0.500000_eff_ineff_p11.dat", "r");
   pF9=fopen("thick_npetrh0.5/cnt_thick_0_80_npetrh0.500000_eff_ineff_p11_allcnt.dat", "r");
   for (int k9=1; k9<2; k9++)
  {
    fscanf(pF9," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF9); 

  FILE *pF10;
   //pF10=fopen("thick_npetrh1.0/cnt_thick_20_60_npetrh1.000000_eff_ineff_p11.dat", "r");
   pF10=fopen("thick_npetrh1.0/cnt_thick_0_80_npetrh1.000000_eff_ineff_p11_allcnt.dat", "r");
   for (int k10=1; k10<2; k10++)
  {
    fscanf(pF10," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF10); 

  FILE *pF11;
   //pF11=fopen("thick_npetrh2.0/cnt_thick_20_60_npetrh2.000000_eff_ineff_p11.dat", "r");
   pF11=fopen("thick_npetrh2.0/cnt_thick_0_80_npetrh2.000000_eff_ineff_p11_allcnt.dat", "r");
   for (int k11=1; k11<2; k11++)
  {
    fscanf(pF11," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF11);
    
  FILE *pF12;
   //pF12=fopen("thick_npetrh2.5/cnt_thick_20_60_npetrh2.500000_eff_ineff_p11.dat", "r");
   pF12=fopen("thick_npetrh2.5/cnt_thick_0_80_npetrh2.500000_eff_ineff_p11_allcnt.dat", "r");
   for (int k12=1; k12<2; k12++)
  {
    fscanf(pF12," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF12); 

  FILE *pF13;
   //pF13=fopen("thick_npetrh3.0/cnt_thick_20_60_npetrh3.000000_eff_ineff_p11.dat", "r");
   pF13=fopen("thick_npetrh3.0/cnt_thick_0_80_npetrh3.000000_eff_ineff_p11_allcnt.dat", "r");
   for (int k13=1; k13<2; k13++)
  {
    fscanf(pF13," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
    if(cnt<80)
    {
    ineff_P_1st.push_back(not_eff_pion);
    eff_P_1st.push_back(eff_pion);
    ineff_K_1st.push_back(not_eff_kaon);
    eff_K_1st.push_back(eff_kaon);
    trh_1st.push_back(npetrh);
    }  
    if(cnt>80)
    {
    ineff_P_2nd.push_back(not_eff_pion);
    eff_P_2nd.push_back(eff_pion);
    ineff_K_2nd.push_back(not_eff_kaon);
    eff_K_2nd.push_back(eff_kaon);
    trh_2nd.push_back(npetrh);
    }  
   }
   fclose(pF13); 

  FILE *pF14;
   //pF14=fopen("thick_npetrh3.5/cnt_thick_20_60_npetrh3.500000_eff_ineff_p11.dat", "r");
   pF14=fopen("thick_npetrh3.5/cnt_thick_0_80_npetrh3.500000_eff_ineff_p11_allcnt.dat", "r");
   for (int k14=1; k14<2; k14++)
  {
    fscanf(pF14," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
    if(cnt<80)
    {
    ineff_P_1st.push_back(not_eff_pion);
    eff_P_1st.push_back(eff_pion);
    ineff_K_1st.push_back(not_eff_kaon);
    eff_K_1st.push_back(eff_kaon);
    trh_1st.push_back(npetrh);
    }  
    if(cnt>80)
    {
    ineff_P_2nd.push_back(not_eff_pion);
    eff_P_2nd.push_back(eff_pion);
    ineff_K_2nd.push_back(not_eff_kaon);
    eff_K_2nd.push_back(eff_kaon);
    trh_2nd.push_back(npetrh);
    }  
   }
   fclose(pF14); 


   FILE *pF15;
   //pF15=fopen("thick_npetrh1.5/cnt_thick_20_60_npetrh1.500000_eff_ineff_p11.dat", "r");
   pF15=fopen("thick_npetrh1.5/cnt_thick_0_80_npetrh1.500000_eff_ineff_p11_allcnt.dat", "r");
   for (int k15=1; k15<2; k15++)
  {
    fscanf(pF15," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
    if(cnt<80)
    {
    ineff_P_1st.push_back(not_eff_pion);
    eff_P_1st.push_back(eff_pion);
    ineff_K_1st.push_back(not_eff_kaon);
    eff_K_1st.push_back(eff_kaon);
    trh_1st.push_back(npetrh);
    }  
    if(cnt>80)
    {
    ineff_P_2nd.push_back(not_eff_pion);
    eff_P_2nd.push_back(eff_pion);
    ineff_K_2nd.push_back(not_eff_kaon);
    eff_K_2nd.push_back(eff_kaon);
    trh_2nd.push_back(npetrh);
    }  
   }
   fclose(pF15); 

    FILE *pF46;
   pF46=fopen("thick_npetrh4.0/cnt_thick_0_80_npetrh4.000000_eff_ineff_p11_allcnt.dat", "r");
   for (int k46=1; k46<2; k46++)
  {
    fscanf(pF46," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
    if(cnt<80)
    {
    ineff_P_1st.push_back(not_eff_pion);
    eff_P_1st.push_back(eff_pion);
    ineff_K_1st.push_back(not_eff_kaon);
    eff_K_1st.push_back(eff_kaon);
    trh_1st.push_back(npetrh);
    }  
    if(cnt>80)
    {
    ineff_P_2nd.push_back(not_eff_pion);
    eff_P_2nd.push_back(eff_pion);
    ineff_K_2nd.push_back(not_eff_kaon);
    eff_K_2nd.push_back(eff_kaon);
    trh_2nd.push_back(npetrh);
    }  
   }
   fclose(pF46); 

     FILE *pF47;
   pF47=fopen("thick_npetrh4.5/cnt_thick_0_80_npetrh4.500000_eff_ineff_p11_allcnt.dat", "r");
   for (int k47=1; k47<2; k47++)
  {
    fscanf(pF47," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
    if(cnt<80)
    {
    ineff_P_1st.push_back(not_eff_pion);
    eff_P_1st.push_back(eff_pion);
    ineff_K_1st.push_back(not_eff_kaon);
    eff_K_1st.push_back(eff_kaon);
    trh_1st.push_back(npetrh);
    }  
    if(cnt>80)
    {
    ineff_P_2nd.push_back(not_eff_pion);
    eff_P_2nd.push_back(eff_pion);
    ineff_K_2nd.push_back(not_eff_kaon);
    eff_K_2nd.push_back(eff_kaon);
    trh_2nd.push_back(npetrh);
    }  
   }
   fclose(pF47);


//==========================================================================================================
FILE *pF16;
   //pF16=fopen("thick_npetrh0.0/cnt_thick_20_60_npetrh0.000000_eff_ineff_p12.dat", "r");
   pF16=fopen("thick_npetrh0.0/cnt_thick_0_80_npetrh0.000000_eff_ineff_p12_allcnt.dat", "r");
   for (int k16=1; k16<2; k16++)
  {
    fscanf(pF16," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
    if(cnt<80)
    {
    ineff_P_1st.push_back(not_eff_pion);
    eff_P_1st.push_back(eff_pion);
    ineff_K_1st.push_back(not_eff_kaon);
    eff_K_1st.push_back(eff_kaon);
    trh_1st.push_back(npetrh);
    }  
    if(cnt>80)
    {
    ineff_P_2nd.push_back(not_eff_pion);
    eff_P_2nd.push_back(eff_pion);
    ineff_K_2nd.push_back(not_eff_kaon);
    eff_K_2nd.push_back(eff_kaon);
    trh_2nd.push_back(npetrh);
    }  
   }
   fclose(pF16); 

    FILE *pF03;
   pF03=fopen("thick_npetrh0.3/cnt_thick_0_80_npetrh0.300000_eff_ineff_p12_allcnt.dat", "r");
   //for (int k=1; k<41; k++)
   for (int k=1; k<2; k++)
  {
    fscanf(pF03," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
    if(cnt<80)
    {
    ineff_P_1st.push_back(not_eff_pion);
    eff_P_1st.push_back(eff_pion);
    ineff_K_1st.push_back(not_eff_kaon);
    eff_K_1st.push_back(eff_kaon);
    trh_1st.push_back(npetrh);
    }  
    if(cnt>80)
    {
    ineff_P_2nd.push_back(not_eff_pion);
    eff_P_2nd.push_back(eff_pion);
    ineff_K_2nd.push_back(not_eff_kaon);
    eff_K_2nd.push_back(eff_kaon);
    trh_2nd.push_back(npetrh);
    }  
   }
   fclose(pF03); 


  FILE *pF17;
   //pF17=fopen("thick_npetrh0.5/cnt_thick_20_60_npetrh0.500000_eff_ineff_p12.dat", "r");
   pF17=fopen("thick_npetrh0.5/cnt_thick_0_80_npetrh0.500000_eff_ineff_p12_allcnt.dat", "r");
   for (int k17=1; k17<2; k17++)
  {
    fscanf(pF17," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
    if(cnt<80)
    {
    ineff_P_1st.push_back(not_eff_pion);
    eff_P_1st.push_back(eff_pion);
    ineff_K_1st.push_back(not_eff_kaon);
    eff_K_1st.push_back(eff_kaon);
    trh_1st.push_back(npetrh);
    }  
    if(cnt>80)
    {
    ineff_P_2nd.push_back(not_eff_pion);
    eff_P_2nd.push_back(eff_pion);
    ineff_K_2nd.push_back(not_eff_kaon);
    eff_K_2nd.push_back(eff_kaon);
    trh_2nd.push_back(npetrh);
    }  
   }
   fclose(pF17); 

  FILE *pF18;
   //pF18=fopen("thick_npetrh1.0/cnt_thick_20_60_npetrh1.000000_eff_ineff_p12.dat", "r");
   pF18=fopen("thick_npetrh1.0/cnt_thick_0_80_npetrh1.000000_eff_ineff_p12_allcnt.dat", "r");
   for (int k18=1; k18<2; k18++)
  {
    fscanf(pF18," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
    if(cnt<80)
    {
    ineff_P_1st.push_back(not_eff_pion);
    eff_P_1st.push_back(eff_pion);
    ineff_K_1st.push_back(not_eff_kaon);
    eff_K_1st.push_back(eff_kaon);
    trh_1st.push_back(npetrh);
    }  
    if(cnt>80)
    {
    ineff_P_2nd.push_back(not_eff_pion);
    eff_P_2nd.push_back(eff_pion);
    ineff_K_2nd.push_back(not_eff_kaon);
    eff_K_2nd.push_back(eff_kaon);
    trh_2nd.push_back(npetrh);
    }  
   }
   fclose(pF18); 

  FILE *pF19;
   //pF19=fopen("thick_npetrh2.0/cnt_thick_20_60_npetrh2.000000_eff_ineff_p12.dat", "r");
   pF19=fopen("thick_npetrh2.0/cnt_thick_0_80_npetrh2.000000_eff_ineff_p12_allcnt.dat", "r");
   for (int k19=1; k19<2; k19++)
  {
    fscanf(pF19," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
    if(cnt<80)
    {
    ineff_P_1st.push_back(not_eff_pion);
    eff_P_1st.push_back(eff_pion);
    ineff_K_1st.push_back(not_eff_kaon);
    eff_K_1st.push_back(eff_kaon);
    trh_1st.push_back(npetrh);
    }  
    if(cnt>80)
    {
    ineff_P_2nd.push_back(not_eff_pion);
    eff_P_2nd.push_back(eff_pion);
    ineff_K_2nd.push_back(not_eff_kaon);
    eff_K_2nd.push_back(eff_kaon);
    trh_2nd.push_back(npetrh);
    }  
   }
   fclose(pF19);
    
  FILE *pF20;
   //pF20=fopen("thick_npetrh2.5/cnt_thick_20_60_npetrh2.500000_eff_ineff_p12.dat", "r");
   pF20=fopen("thick_npetrh2.5/cnt_thick_0_80_npetrh2.500000_eff_ineff_p12_allcnt.dat", "r");
   for (int k20=1; k20<2; k20++)
  {
    fscanf(pF20," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
    if(cnt<80)
    {
    ineff_P_1st.push_back(not_eff_pion);
    eff_P_1st.push_back(eff_pion);
    ineff_K_1st.push_back(not_eff_kaon);
    eff_K_1st.push_back(eff_kaon);
    trh_1st.push_back(npetrh);
    }  
    if(cnt>80)
    {
    ineff_P_2nd.push_back(not_eff_pion);
    eff_P_2nd.push_back(eff_pion);
    ineff_K_2nd.push_back(not_eff_kaon);
    eff_K_2nd.push_back(eff_kaon);
    trh_2nd.push_back(npetrh);
    }  
   }
   fclose(pF20); 

  FILE *pF21;
   //pF21=fopen("thick_npetrh3.0/cnt_thick_20_60_npetrh3.000000_eff_ineff_p12.dat", "r");
   pF21=fopen("thick_npetrh3.0/cnt_thick_0_80_npetrh3.000000_eff_ineff_p12_allcnt.dat", "r");
   for (int k21=1; k21<2; k21++)
  {
    fscanf(pF21," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
    if(cnt<80)
    {
    ineff_P_1st.push_back(not_eff_pion);
    eff_P_1st.push_back(eff_pion);
    ineff_K_1st.push_back(not_eff_kaon);
    eff_K_1st.push_back(eff_kaon);
    trh_1st.push_back(npetrh);
    }  
    if(cnt>80)
    {
    ineff_P_2nd.push_back(not_eff_pion);
    eff_P_2nd.push_back(eff_pion);
    ineff_K_2nd.push_back(not_eff_kaon);
    eff_K_2nd.push_back(eff_kaon);
    trh_2nd.push_back(npetrh);
    }  
   }
   fclose(pF21); 

  FILE *pF22;
   //pF22=fopen("thick_npetrh3.5/cnt_thick_20_60_npetrh3.500000_eff_ineff_p12.dat", "r");
   pF22=fopen("thick_npetrh3.5/cnt_thick_0_80_npetrh3.500000_eff_ineff_p12_allcnt.dat", "r");
   for (int k22=1; k22<2; k22++)
  {
    fscanf(pF22," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF22); 


   FILE *pF23;
   //pF23=fopen("thick_npetrh1.5/cnt_thick_20_60_npetrh1.500000_eff_ineff_p12.dat", "r");
   pF23=fopen("thick_npetrh1.5/cnt_thick_0_80_npetrh1.500000_eff_ineff_p12_allcnt.dat", "r");
   for (int k23=1; k23<2; k23++)
  {
    fscanf(pF23," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF23); 

   FILE *pF44;
   pF44=fopen("thick_npetrh4.0/cnt_thick_0_80_npetrh4.000000_eff_ineff_p12_allcnt.dat", "r");
   for (int k44=1; k44<2; k44++)
  {
    fscanf(pF44," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF44); 

     FILE *pF45;
   pF45=fopen("thick_npetrh4.5/cnt_thick_0_80_npetrh4.500000_eff_ineff_p12_allcnt.dat", "r");
   for (int k45=1; k45<2; k45++)
  {
    fscanf(pF45," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF45);

//==========================================================================================================
FILE *pF24;
   //pF24=fopen("thick_npetrh0.0/cnt_thick_20_60_npetrh0.000000_eff_ineff_p13.dat", "r");
   pF24=fopen("thick_npetrh0.0/cnt_thick_0_80_npetrh0.000000_eff_ineff_p13_allcnt.dat", "r");
   for (int k24=1; k24<2; k24++)
  {
    fscanf(pF24," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF24); 

   
    FILE *pF04;
   pF04=fopen("thick_npetrh0.3/cnt_thick_0_80_npetrh0.300000_eff_ineff_p13_allcnt.dat", "r");
   //for (int k=1; k<41; k++)
   for (int k=1; k<2; k++)
  {
    fscanf(pF04," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF04); 


  FILE *pF25;
   //pF25=fopen("thick_npetrh0.5/cnt_thick_20_60_npetrh0.500000_eff_ineff_p13.dat", "r");
   pF25=fopen("thick_npetrh0.5/cnt_thick_0_80_npetrh0.500000_eff_ineff_p13_allcnt.dat", "r");
   for (int k25=1; k25<2; k25++)
  {
    fscanf(pF25," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF25); 

  FILE *pF26;
   //pF26=fopen("thick_npetrh1.0/cnt_thick_20_60_npetrh1.000000_eff_ineff_p13.dat", "r");
   pF26=fopen("thick_npetrh1.0/cnt_thick_0_80_npetrh1.000000_eff_ineff_p13_allcnt.dat", "r");
   for (int k26=1; k26<2; k26++)
  {
    fscanf(pF26," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF26); 

  FILE *pF27;
   //pF27=fopen("thick_npetrh2.0/cnt_thick_20_60_npetrh2.000000_eff_ineff_p13.dat", "r");
   pF27=fopen("thick_npetrh2.0/cnt_thick_0_80_npetrh2.000000_eff_ineff_p13_allcnt.dat", "r");
   for (int k27=1; k27<2; k27++)
  {
    fscanf(pF27," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF27);
    
  FILE *pF28;
   //pF28=fopen("thick_npetrh2.5/cnt_thick_20_60_npetrh2.500000_eff_ineff_p13.dat", "r");
   pF28=fopen("thick_npetrh2.5/cnt_thick_0_80_npetrh2.500000_eff_ineff_p13_allcnt.dat", "r");
   for (int k28=1; k28<2; k28++)
  {
    fscanf(pF28," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF28); 

  FILE *pF29;
   //pF29=fopen("thick_npetrh3.0/cnt_thick_20_60_npetrh3.000000_eff_ineff_p13.dat", "r");
   pF29=fopen("thick_npetrh3.0/cnt_thick_0_80_npetrh3.000000_eff_ineff_p13_allcnt.dat", "r");
   for (int k29=1; k29<2; k29++)
  {
    fscanf(pF29," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF29); 

  FILE *pF30;
   //pF30=fopen("thick_npetrh3.5/cnt_thick_20_60_npetrh3.500000_eff_ineff_p13.dat", "r");
   pF30=fopen("thick_npetrh3.5/cnt_thick_0_80_npetrh3.500000_eff_ineff_p13_allcnt.dat", "r");
   for (int k30=1; k30<2; k30++)
  {
    fscanf(pF30," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF30); 


   FILE *pF31;
   //pF31=fopen("thick_npetrh1.5/cnt_thick_20_60_npetrh1.500000_eff_ineff_p13.dat", "r");
   pF31=fopen("thick_npetrh1.5/cnt_thick_0_80_npetrh1.500000_eff_ineff_p13_allcnt.dat", "r");
   for (int k31=1; k31<2; k31++)
  {
    fscanf(pF31," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF31); 

     FILE *pF42;
   pF42=fopen("thick_npetrh4.0/cnt_thick_0_80_npetrh4.000000_eff_ineff_p13_allcnt.dat", "r");
   for (int k42=1; k42<2; k42++)
  {
    fscanf(pF42," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF42); 

     FILE *pF43;
   pF43=fopen("thick_npetrh4.5/cnt_thick_0_80_npetrh4.500000_eff_ineff_p13_allcnt.dat", "r");
   for (int k43=1; k43<2; k43++)
  {
    fscanf(pF43," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF43);



//==========================================================================================================
FILE *pF32;
   //pF32=fopen("thick_npetrh0.0/cnt_thick_20_60_npetrh0.000000_eff_ineff_p14.dat", "r");
   pF32=fopen("thick_npetrh0.0/cnt_thick_0_80_npetrh0.000000_eff_ineff_p14_allcnt.dat", "r");
   for (int k32=1; k32<2; k32++)
  {
    fscanf(pF32," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF32); 

       FILE *pF05;
   pF05=fopen("thick_npetrh0.3/cnt_thick_0_80_npetrh0.300000_eff_ineff_p14_allcnt.dat", "r");
   //for (int k=1; k<41; k++)
   for (int k=1; k<2; k++)
  {
    fscanf(pF05," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF05); 

  FILE *pF33;
   //pF33=fopen("thick_npetrh0.5/cnt_thick_20_60_npetrh0.500000_eff_ineff_p14.dat", "r");
   pF33=fopen("thick_npetrh0.5/cnt_thick_0_80_npetrh0.500000_eff_ineff_p14_allcnt.dat", "r");
   for (int k33=1; k33<2; k33++)
  {
    fscanf(pF33," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF33); 

  FILE *pF34;
   //pF34=fopen("thick_npetrh1.0/cnt_thick_20_60_npetrh1.000000_eff_ineff_p14.dat", "r");
   pF34=fopen("thick_npetrh1.0/cnt_thick_0_80_npetrh1.000000_eff_ineff_p14_allcnt.dat", "r");
   for (int k34=1; k34<2; k34++)
  {
    fscanf(pF34," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF34); 

  FILE *pF35;
   //pF35=fopen("thick_npetrh2.0/cnt_thick_20_60_npetrh2.000000_eff_ineff_p14.dat", "r");
   pF35=fopen("thick_npetrh2.0/cnt_thick_0_80_npetrh2.000000_eff_ineff_p14_allcnt.dat", "r");
   for (int k35=1; k35<2; k35++)
  {
    fscanf(pF35," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF35);
    
  FILE *pF36;
   //pF36=fopen("thick_npetrh2.5/cnt_thick_20_60_npetrh2.500000_eff_ineff_p14.dat", "r");
   pF36=fopen("thick_npetrh2.5/cnt_thick_0_80_npetrh2.500000_eff_ineff_p14_allcnt.dat", "r");
   for (int k36=1; k36<2; k36++)
  {
    fscanf(pF36," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF36); 

  FILE *pF37;
   //pF37=fopen("thick_npetrh3.0/cnt_thick_20_60_npetrh3.000000_eff_ineff_p14.dat", "r");
   pF37=fopen("thick_npetrh3.0/cnt_thick_0_80_npetrh3.000000_eff_ineff_p14_allcnt.dat", "r");
   for (int k37=1; k37<2; k37++)
  {
    fscanf(pF37," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF37); 

  FILE *pF38;
   //pF38=fopen("thick_npetrh3.5/cnt_thick_20_60_npetrh3.500000_eff_ineff_p14.dat", "r");
   pF38=fopen("thick_npetrh3.5/cnt_thick_0_80_npetrh3.500000_eff_ineff_p14_allcnt.dat", "r");
   for (int k38=1; k38<2; k38++)
  {
    fscanf(pF38," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF38); 


   FILE *pF39;
   //pF39=fopen("thick_npetrh1.5/cnt_thick_20_60_npetrh1.500000_eff_ineff_p14.dat", "r");
   pF39=fopen("thick_npetrh1.5/cnt_thick_0_80_npetrh1.500000_eff_ineff_p14_allcnt.dat", "r");
   for (int k39=1; k39<2; k39++)
  {
    fscanf(pF39," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF39); 

  FILE *pF40;
   pF40=fopen("thick_npetrh4.0/cnt_thick_0_80_npetrh4.000000_eff_ineff_p14_allcnt.dat", "r");
   for (int k40=1; k40<2; k40++)
  {
    fscanf(pF40," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF40); 

     FILE *pF41;
   pF41=fopen("thick_npetrh4.5/cnt_thick_0_80_npetrh4.500000_eff_ineff_p14_allcnt.dat", "r");
   for (int k41=1; k41<2; k41++)
  {
    fscanf(pF41," %i %f %f %f %f %f %f %f %f %i ",&cnt,&npetrh,&npe,&eff_mu,&eff_pion,&eff_kaon,&not_eff_pion,&not_eff_kaon,&Ksigma,&total_entries);
    if(cnt==cnt_i)
    {
    cout<<"cnt="<<cnt<<"\t"<<"eff_pion="<<eff_pion<<endl;
    ineff_P.push_back(not_eff_pion);
    eff_P.push_back(eff_pion);
    ineff_K.push_back(not_eff_kaon);
    eff_K.push_back(eff_kaon);
    trh.push_back(npetrh);
    }
   }
   fclose(pF41);



//===========================================================================================================
/*
  TFile *fout=0;
  TString fname;
  fname=TString::Format("cnt_thick_misindentification_allcnt.root").Data();
  fout=new TFile(fname,"RECREATE");
  cout<<fname<<endl;
*/

   TProfile *pr0=new TProfile("pr0","Misidentification",100,0,5,0,5);
   TProfile *pr1=new TProfile("pr1","Misidentification",100,0,5,0,5);

   gROOT->SetStyle("Plain");
  
 //  TCanvas *cc1 = new TCanvas("cc1","cnt_thick_misidentification_eff_ineff_all.root",1600,1200);

   for( int i=0; i<eff_P.size(); i++ )
   {
    pr0->Fill(trh.at(i),eff_K.at(i));
    pr1->Fill(trh.at(i),ineff_P.at(i));
    //cout<<"ineff_p="<<1-eff_P.at(i)<<"\t"<<"eff_K="<<eff_K.at(i)<<"\t"<<"ineff="<<ineff_K.at(i)<<endl;
   }

   pr0->SetMarkerStyle(20);
   pr0->SetMarkerColor(3);
   pr0->SetLineWidth(2);
   pr0->SetLineColor(3);
   pr0->SetTitle("Misidentification; Threshold, N_{ph.e.}; Misidentification");
   pr0->GetYaxis()->SetRangeUser(0,0.35);
   //pr0->SetTitle("Misidentification for 1st layer ATC; Threshold, N_{ph.e.}; Misidentification");
  
   pr1->SetMarkerStyle(20);
   pr1->SetMarkerColor(4);
   pr1->SetLineWidth(2);
   pr1->SetLineColor(4);
   //pr1->SetTitle("Misidentification for 1st layer ATC; Threshold, N_{ph.e.}; Misidentification");
  
   //pr0->Draw("prof");
   //pr1->Draw("same");
   pr0->Print("all");
   pr1->Print("all");
    
   //======================================================================================
   TMultiGraph *mg1 = new TMultiGraph();
   //TGraphErrors* gr1=new TGraphErrors(eff_P.size(), &trh[0], &eff_K[0]);
   TGraphErrors *gr1=new TGraphErrors("for_gr1.dat");   
   gr1->SetMarkerStyle(20);
   gr1->SetMarkerColor(2);
   gr1->SetMarkerSize(1.4);
   gr1->SetLineWidth(2);
   gr1->SetLineColor(2);
   gr1->SetName("gr1");
   gr1->SetTitle("gr1");
  
   //TGraphErrors* gr2=new TGraphErrors(eff_P.size(), &trh[0], &ineff_P[0]);
   TGraphErrors *gr2=new TGraphErrors("for_gr2.dat");   
   gr2->SetMarkerStyle(25);
   gr2->SetMarkerColor(2);
   gr2->SetMarkerSize(1.4);
   gr2->SetLineWidth(2);
   gr2->SetLineColor(2);
   gr2->SetName("gr2");
   gr2->SetTitle("gr2");
  
   TGraphErrors *gr3=new TGraphErrors("for_gr1_OR.dat");   
   gr3->SetMarkerStyle(25);
   gr3->SetMarkerColor(4);
   gr3->SetMarkerSize(1.4);
   gr3->SetLineWidth(2);
   gr3->SetLineColor(4);
   gr3->SetName("gr3");
   gr3->SetTitle("gr3");

   TGraphErrors *gr4=new TGraphErrors("for_gr2_OR.dat");   
   gr4->SetMarkerStyle(20);
   gr4->SetMarkerColor(4);
   gr4->SetMarkerSize(1.4);
   gr4->SetLineWidth(2);
   gr4->SetLineColor(4);
   gr4->SetName("gr4");
   gr4->SetTitle("gr4");

   TGraphErrors *gr5=new TGraphErrors("for_gr2_AND.dat");   
   gr5->SetMarkerStyle(25);
   gr5->SetMarkerColor(3);
   gr5->SetMarkerSize(1.4);
   gr5->SetLineWidth(2);
   gr5->SetLineColor(3);
   gr5->SetName("gr5");
   gr5->SetTitle("gr5");

   TGraphErrors *gr6=new TGraphErrors("for_gr1_AND.dat");   
   gr6->SetMarkerStyle(20);
   gr6->SetMarkerColor(3);
   gr6->SetMarkerSize(1.4);
   gr6->SetLineWidth(2);
   gr6->SetLineColor(3);
   gr6->SetName("gr6");
   gr6->SetTitle("gr6");
  
   mg1->Add(gr1);
   mg1->Add(gr2);
   mg1->Add(gr3);
   mg1->Add(gr4);
   mg1->Add(gr5);
   mg1->Add(gr6);
   //mg1->GetYaxis()->SetRangeUser(0,0.35); 
   //mg1->GetYaxis()->SetRangeUser(0,1.0); 
   mg1->GetXaxis()->SetTitleOffset(1.20);
   mg1->GetXaxis()->SetLabelSize(0.05);
   mg1->GetYaxis()->SetTitleSize(0.05);
   mg1->GetYaxis()->SetTitleOffset(1.05);
   mg1->GetYaxis()->SetLabelSize(0.05);
   mg1->GetYaxis()->SetTitle("Misidentification");
   mg1->GetXaxis()->SetTitle("Threshold, N_{ph.e.}");
   mg1->Draw("apl");

   auto legend = new TLegend(0.1,0.75,0.40,0.9);
   //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   //legend->AddEntry(h1,"Histogram filled with random numbers","f");
   //legend->AddEntry("f1","Function abs(#frac{sin(x)}{x})","l");
   legend->AddEntry("gr1","Thick - K","lep");
   legend->AddEntry("gr2","Thick - #pi ","lep");
   legend->AddEntry("gr4","OR - K","lep");
   legend->AddEntry("gr3","OR - #pi","lep");
   legend->AddEntry("gr6","AND - K","lep");
   legend->AddEntry("gr5","AND - #pi","lep");
   //legend->AddEntry("pr2","2^{nd} layer","lep");
   legend->Draw("same");

/*
pr0->Write();                                                     
pr1->Write();                                                     
pr2->Write();                                                     
pr3->Write();                                                     
legend->Write();                                                     
mg1->Write();                                                     
//fout->Write();                                                     
fout->Close();
*/

//cc1->cd();
//cc1->SaveAs("cnt_thick_misidentification_allcnt_canvas.root");

}







