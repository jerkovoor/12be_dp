#include <fstream>
#include <sstream>

void solidAngle() {
    
    TFile *g = new TFile("/home/jerome/12Be_exp/Analysis/BeamOffset/Be_pedestal_TRIUMF_DL_NoOffset_Shift0_0885_IC_Cut_RandomAngle_TargetDistance80_88_Yu_Mod_UnevenBinning_SolidAngle_n4.root","READ");
    TH1D *hAngle_1 = (TH1D*)g->Get("hAngle_1");
    TH1D *hAngle_2 = (TH1D*)g->Get("hAngle_2");
    TH1D *hAngle_3 = (TH1D*)g->Get("hAngle_3");
    TH1D *hAngle_4 = (TH1D*)g->Get("hAngle_4");

    //read in the geometry file to calculate the angles
	ifstream geometry;
	geometry.open ("/home/jerome/12Be_exp/scripts/geometry_s1506.txt"); //open the geometry file; change the path and name accordingly
	if(!geometry.is_open()){
		cout << " No Geometry file found " << endl;
		return -1;
	}

	string read_geometry;
	istringstream iss;
	string name, dummy;
	float Ydr0, Ydr1, Yuz; // Ydz distance from the target, Ydr0 inner radius, Ydz outer radius, same for the Sd

	while ( getline ( geometry,read_geometry ) ){ //in the calibration file named geometry start reading the lines

		if ( read_geometry.find ( "YD_INNER_RADIUS",0 ) !=string::npos ){
			iss.clear();
			iss.str ( read_geometry );
			iss >> name >> dummy >> Ydr0;
		}

		if ( read_geometry.find ( "YD_OUTER_RADIUS",0 ) !=string::npos ){
			iss.clear();
			iss.str ( read_geometry );
			iss >> name >> dummy >> Ydr1;
		}

		if ( read_geometry.find ( "YU_DISTANCE",0 ) !=string::npos ){
			iss.clear();
			iss.str ( read_geometry );
			iss >> name >> dummy >> Yuz;
		}
	}
	//end of the geometry file
    
    double width=(Ydr1-Ydr0)/16;
    
    int n=4;
    int distLength=(16/n)+1, omegaLength=(16/n);
    double dist[distLength], omega[omegaLength];
    double Yubins[17], solidAngleBins[distLength];
    
    for (int i=0; i<17; i++){
		Yubins[i]=180-TMath::RadToDeg()*TMath::ATan((Ydr0+(16-i)*width)/(0-Yuz));
		//cout << Yubins[i] << endl;
	}
    
    for (int i=0;i<distLength;i++){
        dist[i]=Ydr0+n*i*width;
        solidAngleBins[i]=180-TMath::RadToDeg()*TMath::ATan((Ydr0+(distLength-1-i)*n*width)/(0-Yuz));
        //cout << solidAngleBins[i] << endl;
    }
    
    for (int i=0;i<omegaLength;i++){
        omega[i]=42*(M_PI/180)*(0-Yuz)*((1/TMath::Sqrt(pow(dist[i],2)+Yuz*Yuz))-(1/TMath::Sqrt(pow(dist[i+1],2)+Yuz*Yuz)));
        //cout << omega[i] << endl;
    }
    
    double counts1[16], counts2[16], counts3[16], angle1[16], angle2[16], angle3[16], angle4[16], counts4[16], error1[16], error2[16], error3[16], error4[16];
    
    TGraphErrors *gr1 = new TGraphErrors();
    TGraphErrors *gr2 = new TGraphErrors();
    TGraphErrors *gr3 = new TGraphErrors();
    TGraphErrors *gr4 = new TGraphErrors();
    
    for (int i=0;i<omegaLength;i++){
       counts1[i]=hAngle_1->GetBinContent(i);
       angle1[i]=hAngle_1->GetXaxis()->GetBinCenter(i);
       error1[i]=hAngle_1->GetBinError(i);
       gr1->SetPoint(i,angle1[i],counts1[i]/omega[i]);
       gr1->SetPointError(i,0,error1[i]/omega[i]);
       
       counts2[i]=hAngle_2->GetBinContent(i);
       angle2[i]=hAngle_2->GetXaxis()->GetBinCenter(i);
       error2[i]=hAngle_2->GetBinError(i);
       gr2->SetPoint(i,angle2[i],counts2[i]/omega[i]);
       gr2->SetPointError(i,0,error2[i]/omega[i]);
       
       counts3[i]=hAngle_3->GetBinContent(i);
       angle3[i]=hAngle_3->GetXaxis()->GetBinCenter(i);
       error3[i]=hAngle_3->GetBinError(i);
       gr3->SetPoint(i,angle3[i],counts3[i]/omega[i]);
       gr3->SetPointError(i,0,error3[i]/omega[i]);
       
       
       counts4[i]=hAngle_4->GetBinContent(i);
       angle4[i]=hAngle_4->GetXaxis()->GetBinCenter(i);
       error4[i]=hAngle_4->GetBinError(i);
       gr4->SetPoint(i,angle4[i],counts4[i]/omega[i]);
       gr4->SetPointError(i,0,error4[i]/omega[i]);
    }
       
       TCanvas *c1 = new TCanvas ( "c1" ); //create a canvas
       c1->Divide (2,2);
       
       c1->cd(1);
       gr1->Draw("AP*");
       gr1->SetTitle("State 1");
       
       c1->cd(2);
       gr2->Draw("AP*");
       gr2->SetTitle("State 2");
       
       c1->cd(3);
       gr3->Draw("AP*");
       gr3->SetTitle("State 3");
       
       c1->cd(4);
       gr4->Draw("AP*");
       gr4->SetTitle("State 4");
        
}
