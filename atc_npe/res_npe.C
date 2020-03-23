void res_npe()
{       
        int Ncnt=160;
        //int Ncnt=91;
        int cnt;
        int kx1;
        float kx2,kx3,kx4,kx5,kx6;
        float sum_degr1=0,sum_degr2=0;
        float aer_degr1=0,aer_degr2=0;
        float sh_degr1=0,sh_degr2=0;
        int n1=0,n2=0;
        int Nevaer=0,Nevsh=0;      
        float aer[160],sh[160];  
        for(int cnt=0; cnt<160; cnt++)
        {
         aer[cnt]=0; sh[cnt]=0;
        }
       
        TFile *fout=0;
        fout = new TFile(TString::Format("npe_time_aerogel_shifter.root").Data(),"RECREATE");
        char name0[161];
        char name1[161];
	
        for( int cnt=0; cnt<Ncnt; cnt++ )
        {
                sprintf(name0,"aer_cnt%d_Nphe",cnt);
                TProfile* pr=new TProfile(name0,name0,1000,1400778000,1529296837,0,100); 
                ifstream f_inaer(TString::Format("/home/ovtin/development/work_KEDR/KEDRATC/atc_npe/data_aerogel_cnt_%d.dat",cnt).Data());
                string saer;
                while( getline(f_inaer, saer) )
                {
     		    istringstream issaer(saer);
    		    issaer >> kx3 >> kx4;
                    Nevaer++;
                    pr->Fill(kx3,kx4);
	        }
                f_inaer.close();
                cout<<Nevaer<<endl;

                TF1* myfit=new TF1("myfit","[0]*exp(-(x-1400778000)/[1])",1400778000,1529296837);
  		myfit->SetParLimits(0,0,10000);
  		myfit->SetParLimits(1,500000,1e12);
  		myfit->SetLineColor(kRed);
  		pr->Fit("myfit","","",1400778000,1529296837);
 		pr->SetMarkerStyle(20);
  		pr->SetMarkerSize(0.5);
  		pr->SetMarkerColor(2);  //Red
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
                //pr->Draw("prof");

                aer[cnt]=myfit->GetParameter(0)*exp(-(1529296837-1400778000)/myfit->GetParameter(1))/myfit->GetParameter(0);

                cout<<cnt<<"\t"<<aer[cnt]<<endl; 
                
    		if( cnt<80 )
                //if( cnt>=20 && cnt<60 && myfit->GetChisquare()/myfit->GetNDF()<20 )
                {
                 n1++;
                 sum_degr1+=aer[cnt];
                }  
    		if( cnt>=80 )
                //if( cnt>=100 && cnt<140 && myfit->GetChisquare()/myfit->GetNDF()<20 )
                {
                 n2++;
                 sum_degr2+=aer[cnt];
                }  
                Nevaer=0;
        }
        if(cnt<80)aer_degr1=(1-sum_degr1/n1)*100;
        if(cnt>=80)aer_degr2=(1-sum_degr2/n2)*100;
        //cout<<sum_degr1<<"\t"<<n1<<endl;                  
        //cout<<sum_degr2<<"\t"<<n2<<endl;                  
        
        sum_degr1=0; 
        sum_degr2=0;
        n1=0;
        n2=0;

	for( int cnt=0; cnt<Ncnt; cnt++ )
        {
                sprintf(name1,"sh_cnt%d_Nphe",cnt);
                TProfile* pr1=new TProfile(name1,name1,1000,1400778000,1529296837,0,100); 
                ifstream f_insh(TString::Format("/home/ovtin/development/work_KEDR/KEDRATC/atc_npe/data_shifter_cnt_%d.dat",cnt).Data());
                string ssh;
                while( getline(f_insh, ssh) )
		{ 
   		    istringstream isssh(ssh);
    		    isssh >> kx5 >> kx6;
                    Nevsh++;
                    pr1->Fill(kx5,kx6);
		}    
	        f_insh.close();
                cout<<Nevsh<<endl;

                TF1* myfit1=new TF1("myfit1","[0]*exp(-(x-1400778000)/[1])",1400778000,1529296837);
  		myfit1->SetParLimits(0,0,10000);
  		myfit1->SetParLimits(1,500000,1e12);
  		myfit1->SetLineColor(kBlue);
  		pr1->Fit("myfit1","","",1400778000,1529296837);
 		pr1->SetMarkerStyle(20);
  		pr1->SetMarkerSize(0.5);
  		pr1->SetMarkerColor(4);  //Blue
  		pr1->SetLineWidth(2);
  		pr1->SetLineColor(4);
  		pr1->GetXaxis()->SetTimeDisplay(1);
  		pr1->GetXaxis()->SetTimeFormat("%d/%m/%y%F1970-01-01 00:00:00");
  		pr1->GetXaxis()->SetTitleSize(0.05);
  		pr1->GetXaxis()->SetTitleOffset(1.0);
  		pr1->GetXaxis()->SetLabelSize(0.05);
  		pr1->GetXaxis()->SetNdivisions(205);
  		pr1->GetYaxis()->SetTitleSize(0.05);
  		pr1->GetYaxis()->SetTitleOffset(1.00);
  		pr1->GetYaxis()->SetLabelSize(0.05);
  		pr1->GetYaxis()->SetNdivisions(205);
  		pr1->GetYaxis()->SetDecimals();
  		pr1->GetXaxis()->SetTitle("Time(dd/mm/yy)");
  		pr1->GetYaxis()->SetTitle("N_{ph.e.}");
  		pr1->SetTitle("");

                sh[cnt]=myfit1->GetParameter(0)*exp(-(1529296837-1400778000)/myfit1->GetParameter(1))/myfit1->GetParameter(0);
    	
                if(cnt<80)
		//if( cnt>=20 && cnt<60 && myfit1->GetChisquare()/myfit1->GetNDF()<20 )
                {
                sum_degr1+=sh[cnt];
                n1++;
                }  
                if(cnt>=80)
    		//if( cnt>=100 && cnt<140 && myfit1->GetChisquare()/myfit1->GetNDF()<20 )
                {
                sum_degr2+=sh[cnt];
                n2++;
                }                              
                Nevsh=0;
        }
        if(cnt<80)sh_degr1=(1-sum_degr1/n1)*100;
        if(cnt>=80)sh_degr2=(1-sum_degr2/n2)*100;
        cout<<n1<<"\t"<<n2<<"\t"<<"Aerogel + Shifter -->"<<"\t"<<"Degr1="<<aer_degr1<<"\t"<<"Degr2="<<aer_degr2<<endl;
        cout<<n1<<"\t"<<n2<<"\t"<<"Shifter (QE of PMT) -->"<<"\t"<<"Degr1="<<sh_degr1<<"\t"<<"Degr2="<<sh_degr2<<endl;
        cout<<n1<<"\t"<<n2<<"\t"<<"Aerogel -->"<<"\t"<<"Degr1="<<aer_degr1-sh_degr1<<"\t"<<"Degr2="<<aer_degr2-sh_degr2<<endl;
        for(int cnt=0; cnt<Ncnt; cnt++)
        {
        cout<<cnt<<"\t"<<(1-aer[cnt])*100<<"\t"<<(1-sh[cnt])*100<<"\t"<<(1-aer[cnt])*100-(1-sh[cnt])*100<<endl;
        }  
        fout->Write(); 
        fout->Close();
}


















