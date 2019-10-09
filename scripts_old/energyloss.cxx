#include <iostream>
#include <fstream>
#include "TMath.h"
using namespace std;

void energyloss()
{
double thickness=0.1;
int entries=535; //change entries to 534 for Boron
double energy[entries], range1[entries], range2[entries], range3[entries], range4[entries], range5[entries], range6[entries], range7[entries], range8[entries], range9[entries], eff[16]={0}, ener_entries[16]={0}, ener_loss[16]={0};
//double angle[16]={31.684,32.872,34.029,35.155,36.25,37.316,38.352,39.36,40.339,41.291,42.215,43.114,43.987,44.834,45.658,46.458};

ifstream test1;
test1.open("/home/jerome/Software/Energy_Loss/He_in_Al.txt"); //change Al to B when required
if(test1.is_open())
{
  for(int i=0;i<entries;i++)
  {
    test1 >> energy[i] >> range1[i] >> range2[i] >> range3[i] >> range4[i] >> range5[i] >> range6[i] >> range7[i] >> range8[i] >> range9[i];
  }
}
   
    
    
    double angl[16]={0};
    for (int i=0; i<16; i++){angl[i]=TMath::RadToDeg()*TMath::ATan((50+(i+0.5)*(4.94))/85);
        
     //cout << angl[i] << endl;   
    }

    for (int i=0; i<16; i++){eff[i]=thickness/(TMath::Cos(TMath::DegToRad()*angl[i]));
    
    }
    
    
ofstream loss;
  loss.open("Aluminium2.txt"); //change name accordingly
  //loss << "Thickness = " << thickness << " micron" << endl;
  //loss << "Eff. thickness (microns)    Energy loss (MeV/u)  \n";
for (int i=0;i<16;i++)
{for (int j=0;j<entries;j++)
{
    if (range1[j]-0.003<eff[i] && eff[i]<range1[j]+0.003) //change error bar to 0.00008 for Boron
    {ener_entries[i]=j;
        ener_loss[i]=energy[j];
    }
    
}
//cout << ener_entries[i] << endl;

loss << eff[i] << " " << ener_loss[i] << endl;

}


//For finding the combined loss

double thicknes6[16], los6[16];
ifstream bor1;
bor1.open("/home/jerome/Software/Energy_Loss/Boron1.txt");
if(bor1.is_open())
{for(int i=0; i<16; i++)
                
    {
        bor1 >> thicknes6[i] >> los6[i];
    }
}

double thicknes1[16], los1[16], comb_loss1[16];
ifstream alum1;
alum1.open("/home/jerome/Software/Energy_Loss/Aluminium1.txt");
if(alum1.is_open())
{for(int i=0; i<16; i++)
                
    {
        alum1 >> thicknes1[i] >> los1[i];
        comb_loss1[i]=los1[i]+los6[i];
    
    }
}



double thicknes2[16], los2[16], comb_loss2[16];
ifstream alum2;
alum2.open("/home/jerome/Software/Energy_Loss/Aluminium2.txt");
if(alum2.is_open())
{for(int i=0; i<16; i++)
                
    {
        alum2 >> thicknes2[i] >> los2[i];
        comb_loss2[i]=los2[i]+los6[i];
        //cout << comb_loss2[i] << endl;
    }
}

double thicknes3[16], los3[16], comb_loss3[16];
ifstream alum3;
alum3.open("/home/jerome/Software/Energy_Loss/Aluminium3.txt");
if(alum3.is_open())
{for(int i=0; i<16; i++)
                
    {
        alum3 >> thicknes3[i] >> los3[i];
        comb_loss3[i]=los3[i]+los6[i];
    }
}

double thicknes5[16], los5[16], comb_loss5[16];
ifstream alum5;
alum5.open("/home/jerome/Software/Energy_Loss/Aluminium5.txt");
if(alum5.is_open())
{for(int i=0; i<16; i++)
                
    {
        alum5 >> thicknes5[i] >> los5[i];
        comb_loss5[i]=los5[i]+los6[i];
    }
}


ofstream combloss;
combloss.open("Combined_loss");
for(int i=0;i<16;i++)
{combloss << i+1 << " " << comb_loss5[i] << endl;}


TCanvas *c1 = new TCanvas("c1","c1",600, 400);
TMultiGraph *mg = new TMultiGraph();

TGraph *g = new TGraph(16,thicknes1,comb_loss1);
g->SetName("g");
g->SetLineColor(kBlue);
g->SetLineWidth(3);
g->SetTitle("Al 0.1 + B 0.1 micron");

TGraph *h = new TGraph(16,thicknes2,comb_loss2);
h->SetName("h"); 
h->SetLineColor(kGreen);
h->SetLineWidth(3);
h->SetTitle("Al 0.2 + B 0.1 micron");

TGraph *i = new TGraph(16,thicknes3,comb_loss3);
i->SetName("i"); 
i->SetLineColor(kRed);
i->SetLineWidth(3);
i->SetTitle("Al 0.3 + B 0.1 micron");

TGraph *j = new TGraph(16,thicknes5,comb_loss5);
j->SetName("j"); 
j->SetLineColor(kOrange);
j->SetLineWidth(3);
j->SetTitle("Al 0.5 + B 0.1 micron");

mg->Add(g);
mg->Add(h);
mg->Add(i);
mg->Add(j);
mg->Draw("a");
mg->GetXaxis()->SetTitle("Effective Thickness (microns)");
mg->GetYaxis()->SetTitle("Energy loss (MeV/u)");
//mg->GetXaxis()->SetLimits(0,1);
c1->BuildLegend();
}
