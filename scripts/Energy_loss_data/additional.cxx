double a[16], b[16], c[16], loss[16]

ifstream Alfit;
Alfit.open("/home/jerome/12Be_exp/scripts/scripts/Energy_loss_data/Al_Dead_layer/Al_fitparameters.txt");

if(Alfit.is_open())
{
    Alfit.ignore(256,'\n');
    
    for(int i=0;i<16;i++)
    {
      Alfit >> a[i] >> b[i] >> c[i];
    }
    }

    
	   
	   for(int j=0; j<8; j++)
    {if(YuChannel->at(0)==j*16)
     {
    loss1=(-b[0]+TMath::Sqrt((b[0])*b[0]-4*(a[0]-YuEnergy->at(0))*c[0]))/(2*(a[0]-YuEnergy->at(0)));
    loss2=(-e[0]+TMath::Sqrt((e[0])*e[0]-4*(d[0]-(YuEnergy->at(0)+loss1))*f[0]))/(2*(d[0]-(YuEnergy->at(0)+loss1)));
    loss3=(-h[0]+TMath::Sqrt((h[0])*h[0]-4*(g[0]-(YuEnergy->at(0)+loss1+loss2))*k[0]))/(2*(g[0]-(YuEnergy->at(0)+loss1+loss2)));
    loss=loss1+loss2+loss3;
    hYu[YuChannel->at ( 0 )]->Fill ( YuEnergy->at ( 0 )+loss);
    }
    
    
    else if ( YuChannel->at(0)==16*j+1)
    {
    loss1=(-b[1]+TMath::Sqrt((b[1])*b[1]-4*(a[1]-YuEnergy->at(0))*c[1]))/(2*(a[1]-YuEnergy->at(0)));
    loss2=(-e[1]+TMath::Sqrt((e[1])*e[1]-4*(d[1]-(YuEnergy->at(0)+loss1))*f[1]))/(2*(d[1]-(YuEnergy->at(0)+loss1)));
    loss3=(-h[1]+TMath::Sqrt((h[1])*h[1]-4*(g[1]-(YuEnergy->at(0)+loss1+loss2))*k[1]))/(2*(g[1]-(YuEnergy->at(0)+loss1+loss2)));
    loss=loss1+loss2+loss3;
    hYu[YuChannel->at ( 0 )]->Fill ( YuEnergy->at ( 0 )+loss);
    }
    
    else if ( YuChannel->at(0)==16*j+2)
    {
    loss1=(-b[2]+TMath::Sqrt((b[2])*b[2]-4*(a[2]-YuEnergy->at(0))*c[2]))/(2*c[2]);
    loss2=(-e[2]+TMath::Sqrt((e[2])*e[2]-4*(d[2]-(YuEnergy->at(0)+loss1))*f[2]))/(2*(d[2]-(YuEnergy->at(0)+loss1)));
    loss3=(-h[2]+TMath::Sqrt((h[2])*h[2]-4*(g[2]-(YuEnergy->at(0)+loss1+loss2))*k[2]))/(2*(g[2]-(YuEnergy->at(0)+loss1+loss2)));
    loss=loss1+loss2+loss3;
    hYu[YuChannel->at ( 0 )]->Fill ( YuEnergy->at ( 0 )+loss);
    }
    
    else if ( YuChannel->at(0)==16*j+3)
     {
    loss1=(-b[3]+TMath::Sqrt((b[3])*b[3]-4*(a[3]-YuEnergy->at(0))*c[3]))/(2*c[3]);
    loss2=(-e[3]+TMath::Sqrt((e[3])*e[3]-4*(d[3]-(YuEnergy->at(0)+loss1))*f[3]))/(2*(d[3]-(YuEnergy->at(0)+loss1)));
    loss3=(-h[3]+TMath::Sqrt((h[3])*h[3]-4*(g[3]-(YuEnergy->at(0)+loss1+loss2))*k[3]))/(2*(g[3]-(YuEnergy->at(0)+loss1+loss2)));
    loss=loss1+loss2+loss3;
    hYu[YuChannel->at ( 0 )]->Fill ( YuEnergy->at ( 0 )+loss);
    }
    
    else if ( YuChannel->at(0)==16*j+4)
     {
    loss1=(-b[4]+TMath::Sqrt((b[4])*b[4]-4*(a[4]-YuEnergy->at(0))*c[4]))/(2*c[4]);
    loss2=(-e[4]+TMath::Sqrt((e[4])*e[4]-4*(d[4]-(YuEnergy->at(0)+loss1))*f[4]))/(2*(d[4]-(YuEnergy->at(0)+loss1)));
    loss3=(-h[4]+TMath::Sqrt((h[4])*h[4]-4*(g[4]-(YuEnergy->at(0)+loss1+loss2))*k[4]))/(2*(g[4]-(YuEnergy->at(0)+loss1+loss2)));
    loss=loss1+loss2+loss3;
    hYu[YuChannel->at ( 0 )]->Fill ( YuEnergy->at ( 0 )+loss);
    }
    
    else if ( YuChannel->at(0)==16*j+5)
     {
    loss1=(-b[5]+TMath::Sqrt((b[5])*b[5]-4*(a[5]-YuEnergy->at(0))*c[5]))/(2*c[5]);
    loss2=(-e[5]+TMath::Sqrt((e[5])*e[5]-4*(d[5]-(YuEnergy->at(0)+loss1))*f[5]))/(2*(d[5]-(YuEnergy->at(0)+loss1)));
    loss3=(-h[5]+TMath::Sqrt((h[5])*h[5]-4*(g[5]-(YuEnergy->at(0)+loss1+loss2))*k[5]))/(2*(g[5]-(YuEnergy->at(0)+loss1+loss2)));
    loss=loss1+loss2+loss3;
    hYu[YuChannel->at ( 0 )]->Fill ( YuEnergy->at ( 0 )+loss);
    }
    
    else if ( YuChannel->at(0)==16*j+6)
     {
    loss1=(-b[6]+TMath::Sqrt((b[6])*b[6]-4*(a[6]-YuEnergy->at(0))*c[6]))/(2*c[6]);
    loss2=(-e[6]+TMath::Sqrt((e[6])*e[6]-4*(d[6]-(YuEnergy->at(0)+loss1))*f[6]))/(2*(d[6]-(YuEnergy->at(0)+loss1)));
    loss3=(-h[6]+TMath::Sqrt((h[6])*h[6]-4*(g[6]-(YuEnergy->at(0)+loss1+loss2))*k[6]))/(2*(g[6]-(YuEnergy->at(0)+loss1+loss2)));
    loss=loss1+loss2+loss3;
    hYu[YuChannel->at ( 0 )]->Fill ( YuEnergy->at ( 0 )+loss);
    }
    
    else if ( YuChannel->at(0)==16*j+7)
     {
    loss1=(-b[7]+TMath::Sqrt((b[7])*b[7]-4*(a[7]-YuEnergy->at(0))*c[7]))/(2*c[7]);
    loss2=(-e[7]+TMath::Sqrt((e[7])*e[7]-4*(d[7]-(YuEnergy->at(0)+loss1))*f[7]))/(2*(d[7]-(YuEnergy->at(0)+loss1)));
    loss3=(-h[7]+TMath::Sqrt((h[7])*h[7]-4*(g[7]-(YuEnergy->at(0)+loss1+loss2))*k[7]))/(2*(g[7]-(YuEnergy->at(0)+loss1+loss2)));
    loss=loss1+loss2+loss3;
    hYu[YuChannel->at ( 0 )]->Fill ( YuEnergy->at ( 0 )+loss);
    }
    
    else if ( YuChannel->at(0)==16*j+8)
     {
    loss1=(-b[8]+TMath::Sqrt((b[8])*b[8]-4*(a[8]-YuEnergy->at(0))*c[8]))/(2*c[8]);
    loss2=(-e[8]+TMath::Sqrt((e[8])*e[8]-4*(d[8]-(YuEnergy->at(0)+loss1))*f[8]))/(2*(d[8]-(YuEnergy->at(0)+loss1)));
    loss3=(-h[8]+TMath::Sqrt((h[8])*h[8]-4*(g[8]-(YuEnergy->at(0)+loss1+loss2))*k[8]))/(2*(g[8]-(YuEnergy->at(0)+loss1+loss2)));
    loss=loss1+loss2+loss3;
    hYu[YuChannel->at ( 0 )]->Fill ( YuEnergy->at ( 0 )+loss);
    }
    
    else if ( YuChannel->at(0)==16*j+9)
     {
    loss1=(-b[9]+TMath::Sqrt((b[9])*b[9]-4*(a[9]-YuEnergy->at(0))*c[9]))/(2*c[9]);
    loss2=(-e[9]+TMath::Sqrt((e[9])*e[9]-4*(d[9]-(YuEnergy->at(0)+loss1))*f[9]))/(2*(d[9]-(YuEnergy->at(0)+loss1)));
    loss3=(-h[9]+TMath::Sqrt((h[9])*h[9]-4*(g[9]-(YuEnergy->at(0)+loss1+loss2))*k[9]))/(2*(g[9]-(YuEnergy->at(0)+loss1+loss2)));
    loss=loss1+loss2+loss3;
    hYu[YuChannel->at ( 0 )]->Fill ( YuEnergy->at ( 0 )+loss);
    }
    
    else if ( YuChannel->at(0)==16*j+10)
     {
    loss1=(-b[10]+TMath::Sqrt((b[10])*b[10]-4*(a[10]-YuEnergy->at(0))*c[10]))/(2*c[10]);
    loss2=(-e[10]+TMath::Sqrt((e[10])*e[10]-4*(d[10]-(YuEnergy->at(0)+loss1))*f[10]))/(2*(d[10]-(YuEnergy->at(0)+loss1)));
    loss3=(-h[10]+TMath::Sqrt((h[10])*h[10]-4*(g[10]-(YuEnergy->at(0)+loss1+loss2))*k[10]))/(2*(g[10]-(YuEnergy->at(0)+loss1+loss2)));
    loss=loss1+loss2+loss3;
    hYu[YuChannel->at ( 0 )]->Fill ( YuEnergy->at ( 0 )+loss);
    }
    
    else if ( YuChannel->at(0)==16*j+11)
     {
    loss1=(-b[11]+TMath::Sqrt((b[11])*b[11]-4*(a[11]-YuEnergy->at(0))*c[11]))/(2*c[11]);
    loss2=(-e[11]+TMath::Sqrt((e[11])*e[11]-4*(d[11]-(YuEnergy->at(0)+loss1))*f[11]))/(2*(d[11]-(YuEnergy->at(0)+loss1)));
    loss3=(-h[11]+TMath::Sqrt((h[11])*h[11]-4*(g[11]-(YuEnergy->at(0)+loss1+loss2))*k[11]))/(2*(g[11]-(YuEnergy->at(0)+loss1+loss2)));
    loss=loss1+loss2+loss3;
    hYu[YuChannel->at ( 0 )]->Fill ( YuEnergy->at ( 0 )+loss);
    }
    
    else if ( YuChannel->at(0)==16*j+12)
     {
    loss1=(-b[12]+TMath::Sqrt((b[12])*b[12]-4*(a[12]-YuEnergy->at(0))*c[12]))/(2*c[12]);
    loss2=(-e[12]+TMath::Sqrt((e[12])*e[12]-4*(d[12]-(YuEnergy->at(0)+loss1))*f[12]))/(2*(d[12]-(YuEnergy->at(0)+loss1)));
    loss3=(-h[12]+TMath::Sqrt((h[12])*h[12]-4*(g[12]-(YuEnergy->at(0)+loss1+loss2))*k[12]))/(2*(g[12]-(YuEnergy->at(0)+loss1+loss2)));
    loss=loss1+loss2+loss3;
    hYu[YuChannel->at ( 0 )]->Fill ( YuEnergy->at ( 0 )+loss);
    }
    
    else if ( YuChannel->at(0)==16*j+13)
     {
    loss1=(-b[13]+TMath::Sqrt((b[13])*b[13]-4*(a[13]-YuEnergy->at(0))*c[13]))/(2*c[13]);
    loss2=(-e[13]+TMath::Sqrt((e[13])*e[13]-4*(d[13]-(YuEnergy->at(0)+loss1))*f[13]))/(2*(d[13]-(YuEnergy->at(0)+loss1)));
    loss3=(-h[13]+TMath::Sqrt((h[13])*h[13]-4*(g[13]-(YuEnergy->at(0)+loss1+loss2))*k[13]))/(2*(g[13]-(YuEnergy->at(0)+loss1+loss2)));
    loss=loss1+loss2+loss3;
    hYu[YuChannel->at ( 0 )]->Fill ( YuEnergy->at ( 0 )+loss);
    }
    
    else if ( YuChannel->at(0)==16*j+14)
     {
    loss1=(-b[14]+TMath::Sqrt((b[14])*b[14]-4*(a[14]-YuEnergy->at(0))*c[14]))/(2*c[14]);
    loss2=(-e[14]+TMath::Sqrt((e[14])*e[14]-4*(d[14]-(YuEnergy->at(0)+loss1))*f[14]))/(2*(d[14]-(YuEnergy->at(0)+loss1)));
    loss3=(-h[14]+TMath::Sqrt((h[14])*h[14]-4*(g[14]-(YuEnergy->at(0)+loss1+loss2))*k[14]))/(2*(g[14]-(YuEnergy->at(0)+loss1+loss2)));
    loss=loss1+loss2+loss3;
    hYu[YuChannel->at ( 0 )]->Fill ( YuEnergy->at ( 0 )+loss);
    }
    
    else if ( YuChannel->at(0)==16*j+15)
     {
    loss1=(-b[15]+TMath::Sqrt((b[15])*b[15]-4*(a[15]-YuEnergy->at(0))*c[15]))/(2*c[15]);
    loss2=(-e[15]+TMath::Sqrt((e[15])*e[15]-4*(d[15]-(YuEnergy->at(0)+loss1))*f[15]))/(2*(d[15]-(YuEnergy->at(0)+loss1)));
    loss3=(-h[15]+TMath::Sqrt((h[15])*h[15]-4*(g[15]-(YuEnergy->at(0)+loss1+loss2))*k[15]))/(2*(g[15]-(YuEnergy->at(0)+loss1+loss2)));
    loss=loss1+loss2+loss3;
    hYu[YuChannel->at ( 0 )]->Fill ( YuEnergy->at ( 0 )+loss);
    }
    }