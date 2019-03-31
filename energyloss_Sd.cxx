#include <iostream>
#include <fstream>
#include "TMath.h"
using namespace std;

void energyloss_Sd()
{
  double thickness1=61, thickness2=493;
  int entries=535; //change entries to 534 for Boron
  double energy[entries], range1[entries], range2[entries], range3[entries], range4[entries], range5[entries], range6[entries], range7[entries], range8[entries], range9[entries], eff1[24]={0}, eff2[24]={0}, ener_entries1[24]={0}, ener_loss1[24]={0}, ener_entries2[24]={0}, ener_loss2[24]={0};
  //double angle[24]={31.684,32.872,34.029,35.155,36.25,37.324,38.352,39.36,40.339,41.291,42.215,43.114,43.987,44.834,45.658,46.458};

  ifstream test1;
  test1.open("/home/jerome/12Be_exp/scripts/C_in_Si.txt"); //change Al to B when required
    if(test1.is_open())
    {
      for(int i=0;i<entries;i++)
      {
      test1 >> energy[i] >> range1[i] >> range2[i] >> range3[i] >> range4[i] >> range5[i] >> range6[i] >> range7[i] >> range8[i] >> range9[i];
      }
    }
   
  
  
  
  double angl1[24]={0};
    for (int i=0; i<24; i++){angl1[i]=TMath::RadToDeg()*TMath::ATan((11+(i+0.5))/600);
        
     cout << angl1[i] << endl;   
    }

    for (int i=0; i<24; i++){eff1[i]=thickness1/(TMath::Cos(TMath::DegToRad()*angl1[i]));
     //cout << eff1[i] << endl;
    
    }
    
    
  ofstream loss1;
  loss1.open("Sd1_loss.txt"); //change name accordingly
  //loss1 << "Thickness = " << thickness << " micron" << endl;
  //loss1 << "Eff. thickness (microns)    Energy loss (MeV/u)  \n";
    for (int i=0;i<24;i++)
    {
      for (int j=0;j<entries;j++)
      {
	if (range1[j]-0.00025<eff1[i] && eff1[i]<range1[j]+0.00025) //change error bar to 0.00008 for Boron
	{ener_entries1[i]=j;
        ener_loss1[i]=energy[j];
	}
    
      }
      //cout << ener_loss1[i] << endl;

      loss1 << angl1[i] << " " << eff1[i] << " " << ener_loss1[i] << endl;

    }
    
    
    
    double angl2[24]={0};
    for (int i=0; i<24; i++){angl2[i]=TMath::RadToDeg()*TMath::ATan((11+(i+0.5))/690);
        
     cout << angl2[i] << endl;   
    }

    for (int i=0; i<24; i++){eff2[i]=thickness2/(TMath::Cos(TMath::DegToRad()*angl2[i]));
     //cout << eff2[i] << endl;
    
    }
    
    
  ofstream loss2;
  loss2.open("Sd2_loss.txt"); //change name accordingly
  //loss2 << "Thickness = " << thickness << " micron" << endl;
  //loss2 << "Eff. thickness (microns)    Energy loss (MeV/u)  \n";
    for (int i=0;i<24;i++)
    {
      for (int j=0;j<entries;j++)
      {
	if (range1[j]-0.00025<eff2[i] && eff2[i]<range1[j]+0.00025) //change error bar to 0.00008 for Boron
	{ener_entries2[i]=j;
        ener_loss2[i]=energy[j];
	}
    
      }
      //cout << ener_loss2[i] << endl;

      loss2 <<  angl2[i] << " " << eff2[i] << " " << ener_loss2[i] << endl;

    }


//For finding the combined loss

 /* double thicknes6[24], los6[24];
  ifstream bor1;
  bor1.open("/home/jerome/Software/Energy_Loss/Boron1.txt");
    if(bor1.is_open())
    {
      for(int i=0; i<24; i++)
                
      {
        bor1 >> thicknes6[i] >> los6[i];
      }
    }

  double thicknes1[24], los1[24], comb_loss1[24];
  ifstream alum1;
  alum1.open("/home/jerome/Software/Energy_Loss/Aluminium1.txt");
    if(alum1.is_open())
    {
      for(int i=0; i<24; i++)
                
      {
        alum1 >> thicknes1[i] >> los1[i];
        comb_loss1[i]=los1[i]+los6[i];
    
      }
    }*/
}