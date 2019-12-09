using namespace std;

#include "TFile.h" 
#include "TMath.h"
#include <cmath>
#include "iostream"
#include "fstream"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"

void FitBeamOffset(){
    
    
    double OffsetInitialx=-10;
    double OffsetFinalx=0;
    double OffsetInitialy=-5;
    double OffsetFinaly=5;
    double interval=(OffsetFinalx-OffsetInitialx)/20;
    
    int TD=80;
    float TT=49.99;
    
    TString matrix = Form("_%.0f_%.0f_%.0f_%.0f_",OffsetInitialx,OffsetFinalx,OffsetInitialy,OffsetFinaly);
    
	TFile *f = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/carbon/C_pedestal_TRIUMF_DL_BeamOffset" + matrix + Form("TargetDistance%imm_TargetThickness%.2fum_NewSectorGeometry.root",TD,TT),"READ");
	TF1* fit_func1 = new TF1 ( "fit_func1","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))",-10,10);
    TF1* fit_func3 = new TF1 ( "fit_func3","[0]*TMath::Exp(-pow((x-[1]),2)/(2*[2]*[2]))+[3]*TMath::Exp(-pow((x-[4]),2)/(2*[2]*[2]))+[5]*TMath::Exp(-pow((x-[6]),2)/(2*[2]*[2]))+[7]*TMath::Exp([8]*x)",-10,10);
	
    TFile *f_out = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/carbon/C_pedestal_TRIUMF_DL_BeamOffset" + matrix + Form("TargetDistance%imm_TargetThickness%.2fum_NewSectorGeometry_withFit.root",TD,TT),"RECREATE");
    
	ofstream qvalfitpar;
	//qvalfitpar.open ( "/home/jerome/12Be_exp/Analysis/CarbonGain/TES_QValues.txt" ); //change the path and the name of the file acordingly
	qvalfitpar.open ( "/home/jerome/12Be_exp/Analysis/BeamOffset/carbon/GS_QValues_BeamOffset" + matrix + Form("TargetDistance%imm_TargetThickness%.2fum_NewSectorGeometry.txt",TD,TT) ); //change the path and the name of the file acordingly
	//qvalfitpar << " x      y    	 Amplitude 	 Mean 	 SD  	Deviation\n";
	
    int OffsetLengthx=1+(OffsetFinalx-OffsetInitialx)/interval;
    int OffsetLengthy=1+(OffsetFinaly-OffsetInitialy)/interval;
    int Length=OffsetLengthx*OffsetLengthy;

	TCanvas *c1 = new TCanvas ( "c1" );
    c1->Divide ( OffsetLengthy,OffsetLengthx);
    
    f_out->cd();
    
    double meanGS[Length], meanSES[Length], meanTES[Length], sd[Length], devMeanGS[Length], devMeanSES[Length], devMeanTES[Length], fom[Length];
	
	for(int i=0;i<Length;i++){
		TH1D *hQ = (TH1D*)f->Get(Form("hQval_%d",i));
		c1->cd(i+1);
		//hQ->Draw ();
		//hQ->Rebin(2);
        
        
        
		//GS
        
		fit_func1->SetParLimits ( 2,0.05,0.5);
		fit_func1->SetParLimits ( 0,1,25);
        fit_func1->SetParLimits ( 1,2.1,3.2);
		hQ->Fit ( "fit_func1","QE","",1,4);
        
        hQ->GetXaxis()->SetRangeUser(1.5,4.5);
        string name2 = Form ( "hQval_%i",i );
        hQ->Write(name2.c_str());
        
        meanGS[i]=fit_func1->GetParameter ( 1 );
        devMeanGS[i]=fabs(2.722-meanGS[i]);//fabs
        sd[i]=fit_func1->GetParameter ( 2 );
        fom[i]=sqrt(fabs(sd[i]*devMeanGS[i]));
        
        qvalfitpar << i << "	 " << OffsetInitialx+((i-i%OffsetLengthy)/OffsetLengthy)*interval << "	 " << OffsetInitialy+(i%OffsetLengthy)*interval << "	 " <<  sd[i] <<  "	 " <<  devMeanGS[i]  << "	 " << fom[i] << endl;
        
        
        
        //TES
		/*fit_func1->SetParLimits ( 2,0.1,0.5);
		if(i>1){
			fit_func1->SetParLimits ( 0,25,80);
		}
		hQ->Fit ( "fit_func1","","",-3,1);
		
		hQ->GetXaxis()->SetRangeUser(-3,1);
		deviation[i]=abs(fit_func1->GetParameter ( 1 )-C_TESQval);*/
        
        
        
        
        /*
        // 3-Gaussian Fit for Carbon
        //fit_func3->SetParLimits ( 0,25,200);
        fit_func3->SetParLimits ( 1,-1.7,-0.8);
        //fit_func3->SetParLimits ( 3,5,40);
        fit_func3->SetParLimits (4,-0.7,-0.1);//( 4,-0.6,0.4);
        //fit_func3->SetParLimits ( 5,5,50);
        fit_func3->SetParLimits ( 6,2.1,3.2);
        fit_func3->SetParLimits ( 2,0.1,0.5);
        
        hQ->Fit ( "fit_func3","RQE","",-3,4);
    
		//hQ->GetXaxis()->SetRangeUser(1.5,4.5);
        string name2 = Form ( "hQval_%i",i );
        hQ->Write(name2.c_str());
        
        meanGS[i]=fit_func3->GetParameter ( 6 );
        meanSES[i]=fit_func3->GetParameter ( 4 );
        meanTES[i]=fit_func3->GetParameter ( 1 );
        sd[i]=fit_func3->GetParameter ( 2 );
        devMeanGS[i]=fabs(2.722-meanGS[i]);//fabs
        devMeanSES[i]=fabs(-0.36-meanSES[i]);
        devMeanTES[i]=fabs(-1.132-meanTES[i]);
        //fom[i]=sqrt((devMeanTES[i]*devMeanTES[i])+(devMeanGS[i]*devMeanGS[i]));
        fom[i]=sqrt(sd[i]*devMeanTES[i]*devMeanSES[i]*devMeanGS[i]);
		
		
        
        
		
		//qvalfitpar << OffsetInitialx+((i-i%OffsetLengthx)/OffsetLengthx)*interval << "	 " << OffsetInitialy-(i%OffsetLengthy)*interval << "	 " << fit_func1->GetParameter ( 0 ) << " 	" <<  fit_func1->GetParameter ( 1 ) << "	 " <<  fit_func1->GetParameter ( 2 ) << endl;
		qvalfitpar << i << "	 " << OffsetInitialx+((i-i%OffsetLengthy)/OffsetLengthy)*interval << "	 " << OffsetInitialy+(i%OffsetLengthy)*interval << "	 " <<  sd[i] <<  "	 " <<  devMeanGS[i]  <<  "	 " << devMeanSES[i]  <<  "	 " << devMeanTES[i]   <<  "	 " << fom[i] << endl;
        */
	}
	f_out->Close();
}
