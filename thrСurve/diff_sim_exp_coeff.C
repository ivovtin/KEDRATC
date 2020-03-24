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
#include <TStyle.h>
#include <TLine.h>
#include <TEventList.h>
#include <TProfile.h>
#include <vector>

int region;

void change_Kcnt()
{
const int n=161;
float x1[n],x2[n],x3[n],x4[n];
float y1[n],y2[n];

float ampl_pr[n],eff_pr[n];
float sigma_A[n], sigma_eff[n];
std::vector< float > cnt;
std::vector< float > Kcnt1;

std::vector< float > cnt2;
std::vector< float > Kcnt2;

ofstream of(TString::Format("change_data_Kcnt_all_31122018_%drange_barrel.dat",region).Data(),ios_base::out);   //очистка файла для последующей записи  
of.close();

ostringstream o_str;

//const string DIR_REF = TString::Format("/home/ovtin/development/ATC/cosmic/cosm_18012016/result_diff_exp_sim_31122018_%drange.dat",region).Data();
const string DIR_REF = TString::Format("result_diff_exp_sim_31122018_%drange.dat",region).Data();
string f_in_ref;
o_str << DIR_REF;
f_in_ref = o_str.str(); o_str.str("");

FILE *pF;
pF=fopen(f_in_ref.c_str(), "r");
for (int k=1; k<n; k++)
 {
     fscanf(pF," %f %f %f %f ",&x1[k],&x2[k],&x3[k],&x4[k]);

       cnt.push_back(x1[k]);    
       Kcnt1.push_back(x4[k]);
// cout<<x1[k]<<"\t"<<x2[k]<<"\t"<<x3[k]<<"\t"<<x4[k]<<endl;
 }
fclose(pF);

 string DIR_CAL;

if(region<10)
{
     DIR_CAL=TString::Format("result/data_Kcnt_all_0%drange.dat",region).Data();
}
else{
     DIR_CAL=TString::Format("result/data_Kcnt_all_%drange.dat",region).Data();
}

string f_in_cal;
o_str << DIR_CAL;
f_in_cal = o_str.str(); o_str.str("");

FILE *pF2;
pF2=fopen(f_in_cal.c_str(), "r");
for (int j=1; j<n; j++)
{
     fscanf(pF2," %f %f ",&y1[j],&y2[j]);

       cnt2.push_back(y1[j]);   
       Kcnt2.push_back(y2[j]);       
//       cout<<y1[j]<<"\t"<<y2[j]<<"\t"<<endl;  
 }
fclose(pF2);

       cout<<"counter"<<"\t"<<"Kcnt1"<<"\t"<<"Kcnt2"<<"\t"<<"Kcnt1*Kcnt2"<<"\t"<<endl;
  float result;
  
for( int i=0; i<cnt.size(); i++ )
{
	ofstream myFile(TString::Format("change_data_Kcnt_all_31122018_%drange_barrel.dat",region).Data(),ios::app); 
      //if( (i>=20 && i<60) || (i>=100 && i<140) )  //цикл по числу счетчиков
      //if(Kcnt2.at(i)>0){
     // if( (i<20 || (i>=60 && i<100) || i>=140 ) )  //цикл по числу счетчиков
     // {
	if(Kcnt1.at(i)<0 || Kcnt1.at(i)>100 || ((i<20 || (i>=60 && i<100) || i>=140 )) )
	{
	result=Kcnt2.at(i);
	}
	else{
	result=Kcnt1.at(i)*Kcnt2.at(i);
	//result=Kcnt2.at(i)/Kcnt1.at(i);
	}
        cout<<cnt.at(i)<<"\t"<<Kcnt1.at(i)<<"\t"<<"\t"<<Kcnt2.at(i)<<"\t"<<"\t"<<result<<endl;

     // }
     // else{
     //    result=Kcnt2.at(i);	
     // }
	myFile<<cnt.at(i)<<"\t"<<result<<endl;
	myFile.close(); 
      
      //}
}

}



void diff_sim_exp_coeff_1()
{
const int n=161;
float max_exp[n], x1[n],x2[n],x3[n],x4[n],x5[n],x6[n],x7[n];
float x8[n];//,x9[n];
float max_sim[n], y1[n],y2[n],y3[n],y4[n],y5[n],y6[n],y7[n];
float y8[n];//,y9[n];

float ampl_pr[n],eff_pr[n];
float sigma_A[n], sigma_eff[n];
std::vector< float > cnt;
std::vector< float > ev_exp;
std::vector< float > npe_max;

std::vector< float > cnt2;
std::vector< double > ev_sim;
std::vector< float > npe_max2;

ofstream of(TString::Format("result_diff_exp_sim_31122018_%drange.dat",region).Data(),ios_base::out);   //очистка файла для последующей записи  
of.close();

ostringstream o_str;

//const string DIR_REF = TString::Format("result_eff_cnt_single_exp_15032018_%drange_2.dat",region).Data();
const string DIR_REF = TString::Format("result_eff_cnt_single_exp_15032018_%drange_4.dat",region).Data();
string f_in_ref;
o_str << DIR_REF;
f_in_ref = o_str.str(); o_str.str("");

FILE *pF;
pF=fopen(f_in_ref.c_str(), "r");
for (int k=1; k<n; k++)
 {
     fscanf(pF," %f %f %f %f %f %f %f %f",&x1[k],&x2[k],&x3[k],&x4[k],&x5[k],&x6[k],&x7[k],&x8[k]);

       cnt.push_back(x1[k]);
       max_exp[k]=x4[k]+x3[k]*((pow(2.4,2)-pow(x2[k]*0.001,2))/pow(2.4,2));
       npe_max.push_back(max_exp[k]);
       ev_exp.push_back(x8[k]);

   cout<<x1[k]<<"\t"<<x2[k]<<"\t"<<x3[k]<<"\t"<<x4[k]<<"\t"<<x5[k]<<"\t"<<x6[k]<<"\t"<<x7[k]<<"\t"<<x8[k]<<endl;
 }
fclose(pF);


//const string DIR_CAL = TString::Format("result_eff_cnt_single_sim_15032018_%drange_2.dat",region).Data();
const string DIR_CAL = TString::Format("result_eff_cnt_single_sim_15032018_%drange_4.dat",region).Data();
string f_in_cal;
o_str << DIR_CAL;
f_in_cal = o_str.str(); o_str.str("");

FILE *pF2;
pF2=fopen(f_in_cal.c_str(), "r");
for (int j=1; j<n; j++)
{
     fscanf(pF2," %f %f %f %f %f %f %f %f",&y1[j],&y2[j],&y3[j],&y4[j],&y5[j],&y6[j],&y7[j],&y8[j]);

       cnt2.push_back(y1[j]);
       max_sim[j]=y4[j]+y3[j]*((pow(2.4,2)-pow(y2[j]*0.001,2))/pow(2.4,2));   
       npe_max2.push_back(max_sim[j]);
       ev_sim.push_back(y8[j]);
       
       //cout<<y1[j]<<"\t"<<y2[j]<<"\t"<<y3[j]<<"\t"<<y4[j]<<"\t"<<y5[j]<<"\t"<<y6[j]<<"\t"<<y7[j]<<"\t"<<endl;  
 }
fclose(pF2);

       cout<<"counter"<<"\t"<<"Npe_max_exp"<<"\t"<<"Npe_max_sim"<<"\t"<<"Npe_max_exp/Npe_max_sim"<<"\t"<<endl;
float npe_res;
for( int i=0; i<cnt.size(); i++ )
{
   //    if(npe_max2.at(i)>0){
        cout<<cnt.at(i)<<"\t"<<npe_max.at(i)<<"\t"<<"\t"<<npe_max2.at(i)<<"\t"<<"\t"<<npe_max.at(i)/npe_max2.at(i)<<endl;
	ofstream myFile(TString::Format("result_diff_exp_sim_31122018_%drange.dat",region).Data(),ios::app);
	npe_res=npe_max.at(i)/npe_max2.at(i); 
	if( npe_max2.at(i)>15.0 || ev_sim.at(i)<1100 || ev_exp.at(i)<1100 ){npe_res=1.0;}
	//if( npe_max2.at(i)>15.0 ){npe_res=1.0;}
	myFile<<cnt.at(i)<<"\t"<<npe_max.at(i)<<"\t"<<npe_max2.at(i)<<"\t"<<npe_res<<endl;
	myFile.close(); 
    //   }
}

change_Kcnt();

}


void diff_sim_exp_coeff()
{
 for(region=1; region<=13; region++)
 //for(region=12; region<=13; region++)
 {
   cout<<"----------------------------------"<<"START REGION="<<region<<endl;
   diff_sim_exp_coeff_1();
   cout<<"----------------------------------"<<"END REGION="<<region<<endl;
 }
}


