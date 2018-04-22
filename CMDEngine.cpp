#include <iostream>
#include <fstream>
#include <math.h>
#include <conio.h>

#include <cstdlib>

//------------------------------------------------------------------------------
//FUNKCJE MATEMATYCZNE
//------------------------------------------------------------------------------

//dzielenie---------------------------------------------------------------------
double dziel (double D1, double D2)
{
double iloraz;
iloraz=D1*pow(D2,-1);
return iloraz;       
}
//pierwiastkowanie--------------------------------------------------------------
double pierw (double D1)
{
double pierwiastek;
pierwiastek=pow(D1,0.5);
return pierwiastek;       
}
//kwadratowanie-----------------------------------------------------------------
double kw (double D1)
{
double kwadrat;
kwadrat=pow(D1,2);
return kwadrat;       
}
//zaokraglanie------------------------------------------------------------------
double przybliz (double D1)
{
double zaokraglenie;
zaokraglenie=int(D1+0.5);
return zaokraglenie;
} 

//------------------------------------------------------------------------------
//KLASA KULKA
//------------------------------------------------------------------------------
class kuleLJ
{
      
//------------------------------------------------------------------------------
              
      public:
             
//liczby opisujace stan czastki------------------------------------------------- 
             double x, y, z, x0, y0, z0, xp, yp, zp;
             double px, py, pz, px0, py0, pz0, ppx, ppy, ppz;
             double fx, fy, fz, fxp, fyp, fzp;
             double m, CKFp, CKF; //calkowity kwadrat sily dzialajacy na i-ta czastke
//liczby opisujace stan calego ukladu niezbedne do posrednich obliczen----------
             double UC, UCp, D2FC, D2FCp, dzeta, dzeta0, dp, s, s0, sp;
//lista sasiadow
             double licznik, znacznik, lista[500];             
};

//------------------------------------------------------------------------------
//FUNKCJE FIZYCZNE DLA UKLADU CZASTEK
//------------------------------------------------------------------------------

//energia potencjalna-----------------------------------------------------------
double UP (double dx, double dy, double dz)
{
double r, N1, N2, a;
double energia;
r=dx*dx+dy*dy+dz*dz;
N1=12;N2=6;a=4; 

              
        energia=a*((1/pow(r,N1/2))-(1/pow(r,N2/2)));	
        //energia=a*((1/pow(r,N1/2))-(1/pow(r,N2/2)))-4*(pow(zasieg,-12)-pow(zasieg,-6));				              
    
return energia;           
}
//plus U (sciana energetyczna)--------------------------------------------------
double plusU (double x, double y, double z)
{
double r, b;
double energia;
r=x*x+y*y+z*z;
b=0.000001;
              
        energia=b*pow(r,8);				              
    
return energia;           
}
//sila X------------------------------------------------------------------------
double FX (double dx, double dy, double dz)
{
double silaX;
double r, N1, N2, a;
r=dx*dx+dy*dy+dz*dz;
N1=12;N2=6;a=4; 			
        
        silaX=a*(((N1*dx)/(pow(r,(N1+2)/2)))-((N2*dx)/(pow(r,(N2+2)/2))));            
    
return silaX;           
}
//plus FX (sciana energetyczna)--------------------------------------------------
double plusFX (double x, double y, double z)
{
double silaX;
double r, b;
r=x*x+y*y+z*z;
b=0.000001; 			
        
        silaX=-16*b*x*pow(r,7);          
    
return silaX;           
}
//sila Y------------------------------------------------------------------------
double FY (double dx, double dy, double dz)
{
double silaY;
double r, N1, N2, a;
r=dx*dx+dy*dy+dz*dz;
N1=12;N2=6;a=4; 
        
        silaY=a*(((N1*dy)/(pow(r,(N1+2)/2)))-((N2*dy)/(pow(r,(N2+2)/2))));          
    
return silaY;           
}
//plus FY (sciana energetyczna)-------------------------------------------------
double plusFY (double x, double y, double z)
{
double silaY;
double r, b;
r=x*x+y*y+z*z;
b=0.000001; 			
        
        silaY=-16*b*y*pow(r,7);          
    
return silaY;           
}
//sila Z------------------------------------------------------------------------
double FZ (double dx, double dy, double dz)
{
double silaZ;
double r, N1, N2, a;
r=dx*dx+dy*dy+dz*dz;
N1=12;N2=6;a=4;             
        
        silaZ=a*(((N1*dz)/(pow(r,(N1+2)/2)))-((N2*dz)/(pow(r,(N2+2)/2))));             
    
return silaZ;           
}
//plus FZ (sciana energetyczna)-------------------------------------------------
double plusFZ (double x, double y, double z)
{
double silaZ;
double r, b;
r=x*x+y*y+z*z;
b=0.000001; 			
        
        silaZ=-16*b*z*pow(r,7);          
    
return silaZ;           
}
//druga pochodna potencjalu-----------------------------------------------------
double D2F (double dx, double dy, double dz)
{
double Dsila;
double r, N1, N2, a;
r=dx*dx+dy*dy+dz*dz;
N1=12;N2=6;a=4; 
                                  
        Dsila=24*(26.*pow(r,-(N1+2)/2)-7.*pow(r,-(N2+2)/2));                
    
return Dsila;           
}
//plus druga pochodna potencjalu (sciana energetyczna)---------------------------
double plusD2F (double x, double y, double z)
{
double Dsila;
double r, b;
r=x*x+y*y+z*z;
b=0.000001; 			
        
        Dsila=240*b*pow(r,7);          
    
return Dsila;           
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//FUNKCJE KLASY KULKA
//------------------------------------------------------------------------------

//ENERGIA POTENCJALNA UKLADU CZASTEK--------------------------------------------
double energiaP(class kuleLJ kulka[], int LiczbaCzastek)
{
double potencjalna;
potencjalna=0;
for(int i=0;i<LiczbaCzastek;i++)
{
potencjalna=potencjalna+kulka[i].UC;
}
potencjalna=0.5*potencjalna;
return potencjalna;        
}
//ENERGIA KINETYCZNA UKLADU CZASTEK---------------------------------------------
double energiaK(class kuleLJ kulka[], int LiczbaCzastek)
{
double kinetyczna; 
kinetyczna=0;
for(int i=0;i<LiczbaCzastek;i++)
{
kinetyczna=kinetyczna+dziel(kw(kulka[i].px)+kw(kulka[i].py)+kw(kulka[i].pz),2*kulka[i].m);
}
return kinetyczna;        
}
//PROBNA ENERGIA KINETYCZNA UKLADU CZASTEK--------------------------------------
double EnKinP(class kuleLJ kulka[], int LiczbaCzastek)
{
double kinP; 
kinP=0;
for(int i=0;i<LiczbaCzastek;i++)
{
kinP=kinP+dziel(kw(kulka[i].ppx)+kw(kulka[i].ppy)+kw(kulka[i].ppz),2*kulka[i].m);
}
return kinP;        
}
//ENERGIA CALKOWITA UKLADU CZASTEK---------------------------------------------
double energiaC(double D1, double D2)
{
double calkowita;
calkowita=0; 
calkowita=D1+D2;
return calkowita;       
}

//WCZYTANIE POZYCJI-------------------------------------------------------------

	   
//POCZATKOWE USTAWIANIE KULEK---------------------------------------------------
void ustaw(class kuleLJ kulka[], int LiczbaCzastek, double Rozstawienie, double Predkosci, double Wychylenia, double PrzyrMasy, double sciana, double temp, int listX)
{
       
//------------------------------------------------------------------------------       
int const nmax=LiczbaCzastek;

/*double pozycjeX[13]={-0.5*Rozstawienie, 0.5*Rozstawienie, 0, 0, 0};             
double pozycjeY[13]={-Rozstawienie*pow(3,0.5)*pow(6,-1), -Rozstawienie*pow(3,0.5)*pow(6,-1), Rozstawienie*pow(3,0.5)*pow(3,-1), 0, 0};
double pozycjeZ[13]={0,0,0,Rozstawienie*pow(3,0.5)*pow(3,-1),-Rozstawienie*pow(3,0.5)*pow(3,-1)}; */

/*double pozycjeX[13]={0.9637,-0.1209,0.3710,-0.2182,0.1485,0.1516,-0.4832,0.9155,0.1030,-0.9276,-0.8096,-0.9128,0.8190};
double pozycjeY[13]={-0.0843,-0.1361,-0.9532,-0.8891,0.9179,0.2265,0.8838,0.8296,0.0761,-0.0875,0.3933,-0.7825,-0.3946};
double pozycjeZ[13]={0.6422,1.1551,0.3621,-0.5750,0.9203,-1.1522,-0.5984,0.0023,-0.0382,-0.6387,0.4164,0.1927,-0.6886};*/

double pozycjeX[nmax];
double pozycjeY[nmax];
double pozycjeZ[nmax];

//fcc---------------------------------------------------------------------------
double DBOX=0, DST=0;
int NUNIT=0, m=0;

DBOX=sciana;

		if (LiczbaCzastek==4) 
		{
            NUNIT=1;
        }
        if (LiczbaCzastek==32) 
		{
            NUNIT=2;
        }
        if (LiczbaCzastek==108) 
		{
            NUNIT=3;
        }
        if (LiczbaCzastek==256) 
		{
            NUNIT=4;
        }
        if (LiczbaCzastek==500) 
		{
            NUNIT=5;
        }

DST=DBOX*pow(NUNIT,-1);

pozycjeX[0]=0;
pozycjeY[0]=0;
pozycjeZ[0]=0;

pozycjeX[1]=0;
pozycjeY[1]=DST*pow(2,-1);
pozycjeZ[1]=DST*pow(2,-1);

pozycjeX[2]=DST*pow(2,-1);
pozycjeY[2]=0;
pozycjeZ[2]=DST*pow(2,-1);

pozycjeX[3]=DST*pow(2,-1);
pozycjeY[3]=DST*pow(2,-1);
pozycjeZ[3]=0;


m=0;
        for (int i=0; i<NUNIT; i++) 
        {
            for (int j=0; j<NUNIT; j++) 
            {
                for (int k=0; k<NUNIT; k++) 
                {
                    for (int ij=0; ij<4; ij++) 
                    {
                        pozycjeX[ij+m]=pozycjeX[ij]+DST*(k);
                        pozycjeY[ij+m]=pozycjeY[ij]+DST*(j);
                        pozycjeZ[ij+m]=pozycjeZ[ij]+DST*(i);
                    }
                    m=m+4;
                }
            }
        }



        for (int i=0; i<LiczbaCzastek; i++) 
        {
            pozycjeX[i]=pozycjeX[i]+DST*pow(4,-1)-DBOX*pow(2,-1);
            pozycjeY[i]=pozycjeY[i]+DST*pow(4,-1)-DBOX*pow(2,-1);
            pozycjeZ[i]=pozycjeZ[i]+DST*pow(4,-1)-DBOX*pow(2,-1);
        }
//------------------------------------------------------------------------------
double predkosciX[13]={-0.116246782551502,0.696058419719532,-0.712062198998200,0.476971307058052,-0.854598618462655,0.242120872225855,0.475877594586946,0.200870528646558,-0.663176047488368,0.199838305877608,0.778514206869891,-0.191770755923049,-0.532396831560660};
double predkosciY[13]={-1.135458085431340,0.563443539901198,0.331586731150320,-0.800846231525055,-0.125239386160960,0.118923582600583,-0.489948579851302,0.422656239240425,1.137557693720790,-0.262357227387297,-0.554545520467662,0.666514184787888,0.127713059422233};
double predkosciZ[13]={0.180716022977857,0.073348683780944,-0.316360325376074,-0.312181780467466,0.914800484436235,-0.478165183462807,0.040050893306352,-0.096641871380843,-0.008166810730054,-0.049109827032349,0.429766056279282,-0.769799662575877,0.391743320244766};

double randX[13]={0.0242, 0.0400, 0.0070, 0.0210, 0.0457, 0.0396, 0.0479, 0.0327, 0.0017, 0.0424, 0.0466, 0.0339, 0.0378};
double randY[13]={0.0017, 0.0219, 0.0190, 0.0382, 0.0397, 0.0093, 0.0244, 0.0222, 0.0323, 0.0354, 0.0377, 0.0138, 0.0339};
double randZ[13]={0.0327, 0.0081, 0.0059, 0.0249, 0.0479, 0.0170, 0.0292, 0.0111, 0.0375, 0.0127, 0.0252, 0.0349, 0.0445};       
//------------------------------------------------------------------------------       

double GRAN=0, alfa=temp;
double pedpx, pedpy, pedpz, MasaCalkowita, DeltaMasy;
pedpx=0; pedpy=0; pedpz=0; MasaCalkowita=0;
double skala=Rozstawienie;       
       for(int i=0;i<LiczbaCzastek;i++)
       { 
	   kulka[i].licznik=listX;	
       kulka[i].m=1+i*PrzyrMasy;  
       
       kulka[i].x0=pozycjeX[i];
       kulka[i].y0=pozycjeY[i];
       kulka[i].z0=pozycjeZ[i];
    
       kulka[i].px0=0;
       kulka[i].py0=0;
       kulka[i].pz0=0;
       
       kulka[i].dzeta0=0;
       kulka[i].s0=1;
       
       GRAN=cos(2*M_PI*(rand()*pow(RAND_MAX,-1)))*pow(-2*log((rand()*pow(RAND_MAX,-1))),(1/2));
       kulka[i].px0=GRAN*alfa;
       GRAN=cos(2*M_PI*(rand()*pow(RAND_MAX,-1)))*pow(-2*log((rand()*pow(RAND_MAX,-1))),(1/2));
       kulka[i].py0=GRAN*alfa;
       GRAN=cos(2*M_PI*(rand()*pow(RAND_MAX,-1)))*pow(-2*log((rand()*pow(RAND_MAX,-1))),(1/2));
       kulka[i].pz0=GRAN*alfa;
    
       pedpx=pedpx+kulka[i].px0;
       pedpy=pedpy+kulka[i].py0;
       pedpz=pedpz+kulka[i].pz0;
       
       MasaCalkowita=MasaCalkowita+kulka[i].m;
       
       
       }
      
       for (int i=0; i<LiczbaCzastek; i++)
       {  
       kulka[i].px0=kulka[i].px0-dziel(pedpx,MasaCalkowita);
       kulka[i].py0=kulka[i].py0-dziel(pedpy,MasaCalkowita);
       kulka[i].pz0=kulka[i].pz0-dziel(pedpz,MasaCalkowita);
       
       kulka[i].x=kulka[i].x0;
       kulka[i].y=kulka[i].y0;
       kulka[i].z=kulka[i].z0;
       kulka[i].px=kulka[i].px0;
       kulka[i].py=kulka[i].py0;
       kulka[i].pz=kulka[i].pz0;
       
       kulka[i].dzeta=kulka[i].dzeta0;
       
       kulka[i].s=kulka[i].s0;
       
       
       }

//------------------------------------------------------------------------------       
      if(LiczbaCzastek==2)
      {
         for(int i=0;i<LiczbaCzastek;i++)
         { 
         kulka[i].y0=0+Wychylenia*randY[i];
         kulka[i].z0=0+Wychylenia*randZ[i]; 
         }                   
      }
//------------------------------------------------------------------------------      
}


/*
//STAN CZASTKI (nowy - wyglada ok)------------------------------------------------------------------
double stan(class kuleLJ kulka[], int LiczbaCzastek, double sciana, double zasieg, double shift)
{
double dx=0, dy=0, dz=0, r=0;   
double EU=0;
double Fx=0, Fy=0;
double Fz=0, F2=0;
double F2calk[LiczbaCzastek], Fcalk[LiczbaCzastek], EUsum[LiczbaCzastek];
double Fxx[LiczbaCzastek], Fyy[LiczbaCzastek], Fzz[LiczbaCzastek];
   

for (int i=0;i<LiczbaCzastek;i++)
{
Fxx[i]=0.0;
Fyy[i]=0.0;
Fzz[i]=0.0;
F2calk[i]=0.0;
Fcalk[i]=0.0; 
EUsum[i]=0.0;

EU=0.0;
Fx=0.0;
Fy=0.0;
Fz=0.0;
F2=0.0;   

    for(int j=0;j<LiczbaCzastek;j++)
    {
		dx=0; dy=0; dz=0; r=0;        
            
		dx=kulka[i].x-kulka[j].x;
		dy=kulka[i].y-kulka[j].y;
		dz=kulka[i].z-kulka[j].z; 
    
		dx=dx-sciana*round(dx*pow(sciana,-1));
		dy=dy-sciana*round(dy*pow(sciana,-1));
		dz=dz-sciana*round(dz*pow(sciana,-1));
    
		r=dx*dx+dy*dy+dz*dz;
		
		if(i!=j)
		{
			if(r<zasieg*zasieg)
			{
			EU=UP(dx,dy,dz)+shift;				
      
			Fx=FX(dx,dy,dz);
			Fy=FY(dx,dy,dz);
			Fz=FZ(dx,dy,dz);
                    
			F2=D2F(dx,dy,dz); 
			
			}   
			else if(r>=zasieg*zasieg)
			{
			EU=0.0;				
					
			Fx=0.0;
			Fy=0.0;
			Fz=0.0;
                    
			F2=0.0; 	  
			 
			}   

			Fxx[i]+=Fx;
			Fyy[i]+=Fy;
			Fzz[i]+=Fz;
    
			F2calk[i]+=F2;
			EUsum[i]+=EU;

			Fcalk[i]+=Fxx[i]*Fxx[i]+Fyy[i]*Fyy[i]+Fzz[i]*Fzz[i];   
  
			kulka[i].fx=Fxx[i];
			kulka[i].fy=Fyy[i];
			kulka[i].fz=Fzz[i];
  
			kulka[i].D2FC=F2calk[i];
			kulka[i].CKF=Fcalk[i];
			kulka[i].UC=EUsum[i]; 

		}


	   
    }   
 
}
}
//------------------------------------------------------------------------------
double stanP(class kuleLJ kulka[], int LiczbaCzastek, double sciana, double zasieg, double shift)
{
double dx=0, dy=0, dz=0, r=0;   
double EU;
double Fx, Fy;
double Fz, F2;
double F2calk[LiczbaCzastek], Fcalk[LiczbaCzastek], EUsum[LiczbaCzastek];
double Fxx[LiczbaCzastek], Fyy[LiczbaCzastek], Fzz[LiczbaCzastek];
   

for (int i=0;i<LiczbaCzastek;i++)
{
Fxx[i]=0.0;
Fyy[i]=0.0;
Fzz[i]=0.0;
F2calk[i]=0.0;
Fcalk[i]=0.0; 
EUsum[i]=0.0;


EU=0.0;
Fx=0.0;
Fy=0.0;
Fz=0.0;
F2=0.0;   

		for(int j=0;j<LiczbaCzastek;j++)
		{
		dx=0; dy=0; dz=0; r=0;       
            
		dx=kulka[i].xp-kulka[j].xp;  
		dy=kulka[i].yp-kulka[j].yp;
		dz=kulka[i].zp-kulka[j].zp;
     
		dx=dx-sciana*round(dx*pow(sciana,-1));
		dy=dy-sciana*round(dy*pow(sciana,-1));
		dz=dz-sciana*round(dz*pow(sciana,-1));
    
		r=dx*dx+dy*dy+dz*dz;
		
		if(i!=j)
		{
    
			if(r<zasieg*zasieg)
			{
			EU=UP(dx,dy,dz)+shift;				
			
			Fx=FX(dx,dy,dz);
			Fy=FY(dx,dy,dz);
			Fz=FZ(dx,dy,dz);
                    
			F2=D2F(dx,dy,dz); 
			}   
			else if(r>=zasieg*zasieg)
			{
			EU=0.0;				
					
			Fx=0.0;
			Fy=0.0;
			Fz=0.0;
                    
			F2=0.0;      
			}
			
			
			    Fxx[i]+=Fx;
    			Fyy[i]+=Fy;
  				Fzz[i]+=Fz;
    
    			F2calk[i]+=F2;
    			EUsum[i]+=EU;
			
			
		}			
		}   
    
Fcalk[i]+=Fxx[i]*Fxx[i]+Fyy[i]*Fyy[i]+Fzz[i]*Fzz[i];   
  
  kulka[i].fxp=Fxx[i];
  kulka[i].fyp=Fyy[i];
  kulka[i].fzp=Fzz[i];
  
  kulka[i].D2FCp=F2calk[i];
  kulka[i].CKFp=Fcalk[i];
  kulka[i].UCp=EUsum[i];


}
    
}
*/








//STAN CZASTKI (nowy - wyglada ok)------------------------------------------------------------------
double stan(class kuleLJ kulka[], int LiczbaCzastek, double sciana, double zasieg, double shift, int listX, double skora)
{
double dx=0, dy=0, dz=0, r=0;   
double EU=0;
double Fx=0, Fy=0;
double Fz=0, F2=0;
double F2calk[LiczbaCzastek], Fcalk[LiczbaCzastek], EUsum[LiczbaCzastek];
double Fxx[LiczbaCzastek], Fyy[LiczbaCzastek], Fzz[LiczbaCzastek];
int NLL;
   
if(kulka[0].licznik==listX)
{
kulka[0].licznik=0;
                    
        for(int i=0;i<LiczbaCzastek;i++)
		{
            NLL=-1;
            kulka[i].znacznik=0;
            for(int j=0;j<LiczbaCzastek;j++)
			{
                            
            dx=0; dy=0; dz=0; r=0;        
            
			dx=kulka[i].x-kulka[j].x;
			dy=kulka[i].y-kulka[j].y;
			dz=kulka[i].z-kulka[j].z; 
    
			dx=dx-sciana*round(dx*pow(sciana,-1));
			dy=dy-sciana*round(dy*pow(sciana,-1));
			dz=dz-sciana*round(dz*pow(sciana,-1));
    
			r=dx*dx+dy*dy+dz*dz;
                            
                if (r<((zasieg+skora)*(zasieg+skora)) && i!=j) 
				{
                    NLL++;
                    kulka[i].lista[NLL]=j;
                }
            }
            kulka[i].znacznik=NLL+1;        
        }
               	
}


int j;   

for (int i=0;i<LiczbaCzastek;i++)
{
Fxx[i]=0.0;
Fyy[i]=0.0;
Fzz[i]=0.0;
F2calk[i]=0.0;
Fcalk[i]=0.0; 
EUsum[i]=0.0;

EU=0.0;
Fx=0.0;
Fy=0.0;
Fz=0.0;
F2=0.0;   
j=0;

		for(int jj=0;jj<kulka[i].znacznik;jj++)
		{
            
        j=kulka[i].lista[jj];

    
		dx=0; dy=0; dz=0; r=0;        
            
		dx=kulka[i].x-kulka[j].x;
		dy=kulka[i].y-kulka[j].y;
		dz=kulka[i].z-kulka[j].z; 
    
		dx=dx-sciana*round(dx*pow(sciana,-1));
		dy=dy-sciana*round(dy*pow(sciana,-1));
		dz=dz-sciana*round(dz*pow(sciana,-1));
    
		r=dx*dx+dy*dy+dz*dz;
		
		if(i!=j)
		{
			if(r<zasieg*zasieg)
			{
			EU=UP(dx,dy,dz)+shift;				
      
			Fx=FX(dx,dy,dz);
			Fy=FY(dx,dy,dz);
			Fz=FZ(dx,dy,dz);
                    
			F2=D2F(dx,dy,dz); 
			
			}   
			else if(r>=zasieg*zasieg)
			{
			EU=0.0;				
					
			Fx=0.0;
			Fy=0.0;
			Fz=0.0;
                    
			F2=0.0; 	  
			 
			}   

			Fxx[i]+=Fx;
			Fyy[i]+=Fy;
			Fzz[i]+=Fz;
    
			F2calk[i]+=F2;
			EUsum[i]+=EU;

			Fcalk[i]+=Fxx[i]*Fxx[i]+Fyy[i]*Fyy[i]+Fzz[i]*Fzz[i];   
  
			kulka[i].fx=Fxx[i];
			kulka[i].fy=Fyy[i];
			kulka[i].fz=Fzz[i];
  
			kulka[i].D2FC=F2calk[i];
			kulka[i].CKF=Fcalk[i];
			kulka[i].UC=EUsum[i]; 

		}


	   
    }   
 
}
kulka[0].licznik=kulka[0].licznik+1;
}
//------------------------------------------------------------------------------
double stanP(class kuleLJ kulka[], int LiczbaCzastek, double sciana, double zasieg, double shift, int listX, double skora)
{
double dx=0, dy=0, dz=0, r=0;   
double EU;
double Fx, Fy;
double Fz, F2;
double F2calk[LiczbaCzastek], Fcalk[LiczbaCzastek], EUsum[LiczbaCzastek];
double Fxx[LiczbaCzastek], Fyy[LiczbaCzastek], Fzz[LiczbaCzastek];
int NLL;
   
if(kulka[0].licznik==listX)
{
kulka[0].licznik=0;
                    
        for(int i=0;i<LiczbaCzastek;i++)
		{
            NLL=-1;
            kulka[i].znacznik=0;
            for(int j=0;j<LiczbaCzastek;j++)
			{
                            
            dx=0; dy=0; dz=0; r=0;       
            
			dx=kulka[i].xp-kulka[j].xp;  
			dy=kulka[i].yp-kulka[j].yp;
			dz=kulka[i].zp-kulka[j].zp;
     
			dx=dx-sciana*round(dx*pow(sciana,-1));
			dy=dy-sciana*round(dy*pow(sciana,-1));
			dz=dz-sciana*round(dz*pow(sciana,-1));
    
			r=dx*dx+dy*dy+dz*dz;
                            
                if (r<((zasieg+skora)*(zasieg+skora)) && i!=j) 
				{
                    NLL++;
                    kulka[i].lista[NLL]=j;
                }
            }
            kulka[i].znacznik=NLL+1;        
        }
               	
}
	
	
	
	
int j;	
	for (int i=0;i<LiczbaCzastek;i++)
	{
	Fxx[i]=0.0;
	Fyy[i]=0.0;
	Fzz[i]=0.0;
	F2calk[i]=0.0;
	Fcalk[i]=0.0; 
	EUsum[i]=0.0;
	
	EU=0.0;
	Fx=0.0;
	Fy=0.0;
	Fz=0.0;
	F2=0.0;   
j=0;
			for(int jj=0;jj<kulka[i].znacznik;jj++)
			{
            
            j=kulka[i].lista[jj];
					
			dx=0; dy=0; dz=0; r=0;       
            
			dx=kulka[i].xp-kulka[j].xp;  
			dy=kulka[i].yp-kulka[j].yp;
			dz=kulka[i].zp-kulka[j].zp;
     
			dx=dx-sciana*round(dx*pow(sciana,-1));
			dy=dy-sciana*round(dy*pow(sciana,-1));
			dz=dz-sciana*round(dz*pow(sciana,-1));
    
			r=dx*dx+dy*dy+dz*dz;
		
			if(i!=j)
			{
    
				if(r<zasieg*zasieg)
				{
				EU=UP(dx,dy,dz)+shift;				
			
				Fx=FX(dx,dy,dz);
				Fy=FY(dx,dy,dz);
				Fz=FZ(dx,dy,dz);
                    
				F2=D2F(dx,dy,dz); 
				}   
				else if(r>=zasieg*zasieg)
				{
				EU=0.0;				
					
				Fx=0.0;
				Fy=0.0;
				Fz=0.0;
                    
				F2=0.0;      
				}
			
			
			    	Fxx[i]+=Fx;
    				Fyy[i]+=Fy;
  					Fzz[i]+=Fz;
    
   	 				F2calk[i]+=F2;
    				EUsum[i]+=EU;
			
			
			}			
			}   
    
	Fcalk[i]+=Fxx[i]*Fxx[i]+Fyy[i]*Fyy[i]+Fzz[i]*Fzz[i];   
  
  	kulka[i].fxp=Fxx[i];
  	kulka[i].fyp=Fyy[i];
  	kulka[i].fzp=Fzz[i];
  
  	kulka[i].D2FCp=F2calk[i];
  	kulka[i].CKFp=Fcalk[i];
  	kulka[i].UCp=EUsum[i];


	}



kulka[0].licznik=kulka[0].licznik+1;
}






//------------------------------------------------------------------------------
//PRZESUWANIE CZASTEK
//------------------------------------------------------------------------------

//ROWNANIA NEWTONA--------------------------------------------------------------
double NewtonRK4(class kuleLJ kulka[], int LiczbaCzastek, double KrokCzasowy, double sciana, double zasieg, double shift, int listX, double skora)
{
       
double k1X[LiczbaCzastek], k2X[LiczbaCzastek], k3X[LiczbaCzastek], k4X[LiczbaCzastek];
double k1Y[LiczbaCzastek], k2Y[LiczbaCzastek], k3Y[LiczbaCzastek], k4Y[LiczbaCzastek];
double k1Z[LiczbaCzastek], k2Z[LiczbaCzastek], k3Z[LiczbaCzastek], k4Z[LiczbaCzastek];
       
double w1X[LiczbaCzastek], w2X[LiczbaCzastek], w3X[LiczbaCzastek], w4X[LiczbaCzastek];
double w1Y[LiczbaCzastek], w2Y[LiczbaCzastek], w3Y[LiczbaCzastek], w4Y[LiczbaCzastek];
double w1Z[LiczbaCzastek], w2Z[LiczbaCzastek], w3Z[LiczbaCzastek], w4Z[LiczbaCzastek];       
       
       
//BLOK 1------------------------------------------------------------------------       
      
       for(int i=0;i<LiczbaCzastek;i++)
       {
       kulka[i].x0=kulka[i].x;    
       kulka[i].y0=kulka[i].y;
       kulka[i].z0=kulka[i].z; 
      
       kulka[i].px0=kulka[i].px;    
       kulka[i].py0=kulka[i].py;
       kulka[i].pz0=kulka[i].pz;
      
       kulka[i].xp=0;    
       kulka[i].yp=0;
       kulka[i].zp=0; 
      
       kulka[i].xp=kulka[i].x0;    
       kulka[i].yp=kulka[i].y0;
       kulka[i].zp=kulka[i].z0; 
      
       kulka[i].ppx=kulka[i].px0;    
       kulka[i].ppy=kulka[i].py0;
       kulka[i].ppz=kulka[i].pz0;
       }
      
       stanP(kulka, LiczbaCzastek, sciana, zasieg, shift, listX, skora);
      
       for(int i=0;i<LiczbaCzastek;i++)
       {            
       w1X[i]=dziel(kulka[i].ppx,kulka[i].m);
       w1Y[i]=dziel(kulka[i].ppy,kulka[i].m);
       w1Z[i]=dziel(kulka[i].ppz,kulka[i].m);
       k1X[i]=kulka[i].fxp;
       k1Y[i]=kulka[i].fyp;
       k1Z[i]=kulka[i].fzp;
       }
       
//BLOK 2------------------------------------------------------------------------       
       
       for(int i=0;i<LiczbaCzastek;i++)
       {
       kulka[i].xp=0;    
       kulka[i].yp=0;
       kulka[i].zp=0;        
               
       kulka[i].xp=kulka[i].x0+0.5*KrokCzasowy*w1X[i];
       kulka[i].yp=kulka[i].y0+0.5*KrokCzasowy*w1Y[i];
       kulka[i].zp=kulka[i].z0+0.5*KrokCzasowy*w1Z[i];
       
       kulka[i].xp=kulka[i].xp-sciana*round(kulka[i].xp*pow(sciana,-1));
       kulka[i].yp=kulka[i].yp-sciana*round(kulka[i].yp*pow(sciana,-1));
       kulka[i].zp=kulka[i].zp-sciana*round(kulka[i].zp*pow(sciana,-1));
       
       kulka[i].ppx=kulka[i].px0+0.5*KrokCzasowy*k1X[i];
       kulka[i].ppy=kulka[i].py0+0.5*KrokCzasowy*k1Y[i];
       kulka[i].ppz=kulka[i].pz0+0.5*KrokCzasowy*k1Z[i];
       }

       stanP(kulka, LiczbaCzastek, sciana, zasieg, shift, listX, skora);
       
       for(int i=0;i<LiczbaCzastek;i++)
       {            
       w2X[i]=dziel(kulka[i].ppx,kulka[i].m);
       w2Y[i]=dziel(kulka[i].ppy,kulka[i].m);
       w2Z[i]=dziel(kulka[i].ppz,kulka[i].m);
       k2X[i]=kulka[i].fxp;
       k2Y[i]=kulka[i].fyp;
       k2Z[i]=kulka[i].fzp;
       }
       
//BLOK 3------------------------------------------------------------------------

       for(int i=0;i<LiczbaCzastek;i++)
       {
       kulka[i].xp=0;    
       kulka[i].yp=0;
       kulka[i].zp=0;        
               
       kulka[i].xp=kulka[i].x0+0.5*KrokCzasowy*w2X[i];
       kulka[i].yp=kulka[i].y0+0.5*KrokCzasowy*w2Y[i];
       kulka[i].zp=kulka[i].z0+0.5*KrokCzasowy*w2Z[i];
       
       kulka[i].xp=kulka[i].xp-sciana*round(kulka[i].xp*pow(sciana,-1));
       kulka[i].yp=kulka[i].yp-sciana*round(kulka[i].yp*pow(sciana,-1));
       kulka[i].zp=kulka[i].zp-sciana*round(kulka[i].zp*pow(sciana,-1));
       
       kulka[i].ppx=kulka[i].px0+0.5*KrokCzasowy*k2X[i];
       kulka[i].ppy=kulka[i].py0+0.5*KrokCzasowy*k2Y[i];
       kulka[i].ppz=kulka[i].pz0+0.5*KrokCzasowy*k2Z[i];
       }
       
       stanP(kulka,LiczbaCzastek, sciana, zasieg, shift, listX, skora);

       for(int i=0;i<LiczbaCzastek;i++)
       {            
       w3X[i]=dziel(kulka[i].ppx,kulka[i].m);
       w3Y[i]=dziel(kulka[i].ppy,kulka[i].m);
       w3Z[i]=dziel(kulka[i].ppz,kulka[i].m);
       k3X[i]=kulka[i].fxp;
       k3Y[i]=kulka[i].fyp;
       k3Z[i]=kulka[i].fzp;
       }
       
//BLOK 4------------------------------------------------------------------------  
       for(int i=0;i<LiczbaCzastek;i++)
       {
       kulka[i].xp=0;    
       kulka[i].yp=0;
       kulka[i].zp=0;        
               
       kulka[i].xp=kulka[i].x0+KrokCzasowy*w3X[i];
       kulka[i].yp=kulka[i].y0+KrokCzasowy*w3Y[i];
       kulka[i].zp=kulka[i].z0+KrokCzasowy*w3Z[i];
       
       kulka[i].xp=kulka[i].xp-sciana*round(kulka[i].xp*pow(sciana,-1));
       kulka[i].yp=kulka[i].yp-sciana*round(kulka[i].yp*pow(sciana,-1));
       kulka[i].zp=kulka[i].zp-sciana*round(kulka[i].zp*pow(sciana,-1));
       
       kulka[i].ppx=kulka[i].px0+KrokCzasowy*k3X[i];
       kulka[i].ppy=kulka[i].py0+KrokCzasowy*k3Y[i];
       kulka[i].ppz=kulka[i].pz0+KrokCzasowy*k3Z[i];
       }
       
       stanP(kulka,LiczbaCzastek, sciana, zasieg, shift, listX, skora);
       
       for(int i=0;i<LiczbaCzastek;i++)
       {            
       w4X[i]=dziel(kulka[i].ppx,kulka[i].m);
       w4Y[i]=dziel(kulka[i].ppy,kulka[i].m);
       w4Z[i]=dziel(kulka[i].ppz,kulka[i].m);
       k4X[i]=kulka[i].fxp;
       k4Y[i]=kulka[i].fyp;
       k4Z[i]=kulka[i].fzp;
       }
       
//NOWE POZYCJE I PREDKOSCI------------------------------------------------------
       
       for(int i=0;i<LiczbaCzastek;i++)
       {
       kulka[i].px=kulka[i].px0+dziel(KrokCzasowy,6)*(k1X[i]+2*k2X[i]+2*k3X[i]+k4X[i]);
       kulka[i].py=kulka[i].py0+dziel(KrokCzasowy,6)*(k1Y[i]+2*k2Y[i]+2*k3Y[i]+k4Y[i]);
       kulka[i].pz=kulka[i].pz0+dziel(KrokCzasowy,6)*(k1Z[i]+2*k2Z[i]+2*k3Z[i]+k4Z[i]);
       
       kulka[i].x=0;    
       kulka[i].y=0;
       kulka[i].z=0;
       
       kulka[i].x=kulka[i].x0+dziel(KrokCzasowy,6)*(w1X[i]+2*w2X[i]+2*w3X[i]+w4X[i]);
       kulka[i].y=kulka[i].y0+dziel(KrokCzasowy,6)*(w1Y[i]+2*w2Y[i]+2*w3Y[i]+w4Y[i]);
       kulka[i].z=kulka[i].z0+dziel(KrokCzasowy,6)*(w1Z[i]+2*w2Z[i]+2*w3Z[i]+w4Z[i]);
       
       kulka[i].x=kulka[i].x-sciana*round(kulka[i].x*pow(sciana,-1));
       kulka[i].y=kulka[i].y-sciana*round(kulka[i].y*pow(sciana,-1));
       kulka[i].z=kulka[i].z-sciana*round(kulka[i].z*pow(sciana,-1));
       
       }
       
}
//ROWNANIA NOSEGO-HOOVERA-------------------------------------------------------

double NoseHooverRK4(class kuleLJ kulka[],int LiczbaCzastek, double Q, double temp, double g, double kB, double KrokCzasowy, double sciana, double zasieg, double shift, int listX, double skora)
{
double PP=0;       
       
double kd1, ks1, kd2, ks2, kd3, ks3, kd4, ks4;
      
double k1X[LiczbaCzastek], k2X[LiczbaCzastek], k3X[LiczbaCzastek], k4X[LiczbaCzastek];
double k1Y[LiczbaCzastek], k2Y[LiczbaCzastek], k3Y[LiczbaCzastek], k4Y[LiczbaCzastek];
double k1Z[LiczbaCzastek], k2Z[LiczbaCzastek], k3Z[LiczbaCzastek], k4Z[LiczbaCzastek];
       
double w1X[LiczbaCzastek], w2X[LiczbaCzastek], w3X[LiczbaCzastek], w4X[LiczbaCzastek];
double w1Y[LiczbaCzastek], w2Y[LiczbaCzastek], w3Y[LiczbaCzastek], w4Y[LiczbaCzastek];
double w1Z[LiczbaCzastek], w2Z[LiczbaCzastek], w3Z[LiczbaCzastek], w4Z[LiczbaCzastek];       
       
       
//BLOK 1------------------------------------------------------------------------       
      
       for(int i=0;i<LiczbaCzastek;i++)
       {
       kulka[i].x0=kulka[i].x;    
       kulka[i].y0=kulka[i].y;
       kulka[i].z0=kulka[i].z; 
      
       kulka[i].px0=kulka[i].px;    
       kulka[i].py0=kulka[i].py;
       kulka[i].pz0=kulka[i].pz;
      
       kulka[i].xp=0;    
       kulka[i].yp=0;
       kulka[i].zp=0; 
       
       kulka[i].xp=kulka[i].x0;    
       kulka[i].yp=kulka[i].y0;
       kulka[i].zp=kulka[i].z0; 
      
       kulka[i].ppx=kulka[i].px0;    
       kulka[i].ppy=kulka[i].py0;
       kulka[i].ppz=kulka[i].pz0;
       
       kulka[i].dzeta0=kulka[i].dzeta;
       kulka[i].s0=kulka[i].s;
       
       kulka[i].dp=kulka[i].dzeta0;
       kulka[i].sp=kulka[i].s0;
       }
       
       
      
       stanP(kulka, LiczbaCzastek, sciana, zasieg, shift, listX, skora);
       PP=EnKinP(kulka, LiczbaCzastek); 
      
       for(int i=0;i<LiczbaCzastek;i++)
       {            
       w1X[i]=dziel(kulka[i].ppx,kulka[i].m);
       w1Y[i]=dziel(kulka[i].ppy,kulka[i].m);
       w1Z[i]=dziel(kulka[i].ppz,kulka[i].m);
       k1X[i]=kulka[i].fxp-kulka[0].dp*kulka[i].ppx;
       k1Y[i]=kulka[i].fyp-kulka[0].dp*kulka[i].ppy;
       k1Z[i]=kulka[i].fzp-kulka[0].dp*kulka[i].ppz;
       }
       
   
   
kd1=dziel(1,Q)*(2*PP-g*kB*temp);
ks1=kulka[0].sp*kulka[0].dp;
//BLOK 2------------------------------------------------------------------------       
       
       for(int i=0;i<LiczbaCzastek;i++)
       {
       kulka[i].xp=0;    
       kulka[i].yp=0;
       kulka[i].zp=0;         
               
       kulka[i].xp=kulka[i].x0+0.5*KrokCzasowy*w1X[i];
       kulka[i].yp=kulka[i].y0+0.5*KrokCzasowy*w1Y[i];
       kulka[i].zp=kulka[i].z0+0.5*KrokCzasowy*w1Z[i];
       
       kulka[i].xp=kulka[i].xp-sciana*round(kulka[i].xp*pow(sciana,-1));
       kulka[i].yp=kulka[i].yp-sciana*round(kulka[i].yp*pow(sciana,-1));
       kulka[i].zp=kulka[i].zp-sciana*round(kulka[i].zp*pow(sciana,-1));
       
       kulka[i].ppx=kulka[i].px0+0.5*KrokCzasowy*k1X[i];
       kulka[i].ppy=kulka[i].py0+0.5*KrokCzasowy*k1Y[i];
       kulka[i].ppz=kulka[i].pz0+0.5*KrokCzasowy*k1Z[i];
       
       kulka[i].dp=kulka[i].dzeta0+0.5*KrokCzasowy*kd1;
       kulka[i].sp=kulka[i].s0+0.5*KrokCzasowy*ks1;
       }
       


       stanP(kulka, LiczbaCzastek, sciana, zasieg, shift, listX, skora);
       PP=EnKinP(kulka, LiczbaCzastek); 
       
       for(int i=0;i<LiczbaCzastek;i++)
       {            
       w2X[i]=dziel(kulka[i].ppx,kulka[i].m);
       w2Y[i]=dziel(kulka[i].ppy,kulka[i].m);
       w2Z[i]=dziel(kulka[i].ppz,kulka[i].m);
       k2X[i]=kulka[i].fxp-kulka[0].dp*kulka[i].ppx;
       k2Y[i]=kulka[i].fyp-kulka[0].dp*kulka[i].ppy;
       k2Z[i]=kulka[i].fzp-kulka[0].dp*kulka[i].ppz;
       }
kd2=dziel(1,Q)*(2*PP-g*kB*temp);
ks2=kulka[0].sp*kulka[0].dp;
       
//BLOK 3------------------------------------------------------------------------

       for(int i=0;i<LiczbaCzastek;i++)
       {
       kulka[i].xp=0;    
       kulka[i].yp=0;
       kulka[i].zp=0;       
               
       kulka[i].xp=kulka[i].x0+0.5*KrokCzasowy*w2X[i];
       kulka[i].yp=kulka[i].y0+0.5*KrokCzasowy*w2Y[i];
       kulka[i].zp=kulka[i].z0+0.5*KrokCzasowy*w2Z[i];
       
       kulka[i].xp=kulka[i].xp-sciana*round(kulka[i].xp*pow(sciana,-1));
       kulka[i].yp=kulka[i].yp-sciana*round(kulka[i].yp*pow(sciana,-1));
       kulka[i].zp=kulka[i].zp-sciana*round(kulka[i].zp*pow(sciana,-1));
       
       kulka[i].ppx=kulka[i].px0+0.5*KrokCzasowy*k2X[i];
       kulka[i].ppy=kulka[i].py0+0.5*KrokCzasowy*k2Y[i];
       kulka[i].ppz=kulka[i].pz0+0.5*KrokCzasowy*k2Z[i];
       
       kulka[i].dp=kulka[i].dzeta0+0.5*KrokCzasowy*kd2;
       kulka[i].sp=kulka[i].s0+0.5*KrokCzasowy*ks2;
       }
       


       stanP(kulka, LiczbaCzastek, sciana, zasieg, shift, listX, skora);
       PP=EnKinP(kulka, LiczbaCzastek);

       for(int i=0;i<LiczbaCzastek;i++)
       {            
       w3X[i]=dziel(kulka[i].ppx,kulka[i].m);
       w3Y[i]=dziel(kulka[i].ppy,kulka[i].m);
       w3Z[i]=dziel(kulka[i].ppz,kulka[i].m);
       k3X[i]=kulka[i].fxp-kulka[0].dp*kulka[i].ppx;
       k3Y[i]=kulka[i].fyp-kulka[0].dp*kulka[i].ppy;
       k3Z[i]=kulka[i].fzp-kulka[0].dp*kulka[i].ppz;
       }
kd3=dziel(1,Q)*(2*PP-g*kB*temp);
ks3=kulka[0].sp*kulka[0].dp;
       
//BLOK 4------------------------------------------------------------------------  
       
       for(int i=0;i<LiczbaCzastek;i++)
       {
       kulka[i].xp=0;    
       kulka[i].yp=0;
       kulka[i].zp=0;         
               
       kulka[i].xp=kulka[i].x0+KrokCzasowy*w3X[i];
       kulka[i].yp=kulka[i].y0+KrokCzasowy*w3Y[i];
       kulka[i].zp=kulka[i].z0+KrokCzasowy*w3Z[i];
       
       kulka[i].xp=kulka[i].xp-sciana*round(kulka[i].xp*pow(sciana,-1));
       kulka[i].yp=kulka[i].yp-sciana*round(kulka[i].yp*pow(sciana,-1));
       kulka[i].zp=kulka[i].zp-sciana*round(kulka[i].zp*pow(sciana,-1));
       
       kulka[i].ppx=kulka[i].px0+KrokCzasowy*k3X[i];
       kulka[i].ppy=kulka[i].py0+KrokCzasowy*k3Y[i];
       kulka[i].ppz=kulka[i].pz0+KrokCzasowy*k3Z[i];
       
       kulka[i].dp=kulka[i].dzeta0+KrokCzasowy*kd3;
       kulka[i].sp=kulka[i].s0+KrokCzasowy*ks3;
       }
       

       stanP(kulka, LiczbaCzastek, sciana, zasieg, shift, listX, skora);
       PP=EnKinP(kulka, LiczbaCzastek);
       
       for(int i=0;i<LiczbaCzastek;i++)
       {            
       w4X[i]=dziel(kulka[i].ppx,kulka[i].m);
       w4Y[i]=dziel(kulka[i].ppy,kulka[i].m);
       w4Z[i]=dziel(kulka[i].ppz,kulka[i].m);
       k4X[i]=kulka[i].fxp-kulka[0].dp*kulka[i].ppx;
       k4Y[i]=kulka[i].fyp-kulka[0].dp*kulka[i].ppy;
       k4Z[i]=kulka[i].fzp-kulka[0].dp*kulka[i].ppz;
       }
kd4=dziel(1,Q)*(2*PP-g*kB*temp);
ks4=kulka[0].sp*kulka[0].dp;
       
//NOWE POZYCJE I PREDKOSCI------------------------------------------------------
       
       for(int i=0;i<LiczbaCzastek;i++)
       {
       kulka[i].px=kulka[i].px0+dziel(KrokCzasowy,6)*(k1X[i]+2*k2X[i]+2*k3X[i]+k4X[i]);
       kulka[i].py=kulka[i].py0+dziel(KrokCzasowy,6)*(k1Y[i]+2*k2Y[i]+2*k3Y[i]+k4Y[i]);
       kulka[i].pz=kulka[i].pz0+dziel(KrokCzasowy,6)*(k1Z[i]+2*k2Z[i]+2*k3Z[i]+k4Z[i]);
       
       kulka[i].x=0;    
       kulka[i].y=0;
       kulka[i].z=0; 
       
       kulka[i].x=kulka[i].x0+dziel(KrokCzasowy,6)*(w1X[i]+2*w2X[i]+2*w3X[i]+w4X[i]);
       kulka[i].y=kulka[i].y0+dziel(KrokCzasowy,6)*(w1Y[i]+2*w2Y[i]+2*w3Y[i]+w4Y[i]);
       kulka[i].z=kulka[i].z0+dziel(KrokCzasowy,6)*(w1Z[i]+2*w2Z[i]+2*w3Z[i]+w4Z[i]);
       
       kulka[i].x=kulka[i].x-sciana*round(kulka[i].x*pow(sciana,-1));
       kulka[i].y=kulka[i].y-sciana*round(kulka[i].y*pow(sciana,-1));
       kulka[i].z=kulka[i].z-sciana*round(kulka[i].z*pow(sciana,-1));
       
       kulka[i].dzeta=kulka[i].dzeta0+dziel(KrokCzasowy,6)*(kd1+2*kd2+2*kd3+kd4);
       kulka[i].s=kulka[i].s0+dziel(KrokCzasowy,6)*(ks1+2*ks2+2*ks3+ks4);
       }
       
       
}
//ROWNANIA TRAVISA-BRAGA--------------------------------------------------------

double TravisBragaRK4(class kuleLJ kulka[],int LiczbaCzastek, double Q, double temp, double kB, double KrokCzasowy, double sciana, double zasieg, double shift, int listX, double skora)
{      
       
double kd1, ks1, kd2, ks2, kd3, ks3, kd4, ks4;
double SumCKF;
      
double k1X[LiczbaCzastek], k2X[LiczbaCzastek], k3X[LiczbaCzastek], k4X[LiczbaCzastek];
double k1Y[LiczbaCzastek], k2Y[LiczbaCzastek], k3Y[LiczbaCzastek], k4Y[LiczbaCzastek];
double k1Z[LiczbaCzastek], k2Z[LiczbaCzastek], k3Z[LiczbaCzastek], k4Z[LiczbaCzastek];
       
double w1X[LiczbaCzastek], w2X[LiczbaCzastek], w3X[LiczbaCzastek], w4X[LiczbaCzastek];
double w1Y[LiczbaCzastek], w2Y[LiczbaCzastek], w3Y[LiczbaCzastek], w4Y[LiczbaCzastek];
double w1Z[LiczbaCzastek], w2Z[LiczbaCzastek], w3Z[LiczbaCzastek], w4Z[LiczbaCzastek];       
       
       
//BLOK 1------------------------------------------------------------------------       
      
       for(int i=0;i<LiczbaCzastek;i++)
       {
       kulka[i].x0=kulka[i].x;    
       kulka[i].y0=kulka[i].y;
       kulka[i].z0=kulka[i].z; 
      
       kulka[i].px0=kulka[i].px;    
       kulka[i].py0=kulka[i].py;
       kulka[i].pz0=kulka[i].pz;
       
       kulka[i].xp=0;    
       kulka[i].yp=0;
       kulka[i].zp=0;
      
       kulka[i].xp=kulka[i].x0;    
       kulka[i].yp=kulka[i].y0;
       kulka[i].zp=kulka[i].z0; 
      
       kulka[i].ppx=kulka[i].px0;    
       kulka[i].ppy=kulka[i].py0;
       kulka[i].ppz=kulka[i].pz0;
       
       kulka[i].dzeta0=kulka[i].dzeta;
       kulka[i].s0=kulka[i].s;
       
       kulka[i].dp=kulka[i].dzeta0;
       kulka[i].sp=kulka[i].s0;
       }
       
       
       stanP(kulka, LiczbaCzastek, sciana, zasieg, shift, listX, skora);
      
      
       SumCKF=0;
       for(int i=0;i<LiczbaCzastek;i++)
       {            
       w1X[i]=dziel(kulka[i].ppx,kulka[i].m)+kulka[0].dp*kulka[i].fxp;
       w1Y[i]=dziel(kulka[i].ppy,kulka[i].m)+kulka[0].dp*kulka[i].fyp;
       w1Z[i]=dziel(kulka[i].ppz,kulka[i].m)+kulka[0].dp*kulka[i].fzp;
       k1X[i]=kulka[i].fxp;
       k1Y[i]=kulka[i].fyp;
       k1Z[i]=kulka[i].fzp;
       SumCKF=SumCKF+kulka[i].CKFp;
       }
       
     
kd1=dziel(1,Q)*(SumCKF-temp*kulka[0].D2FCp);
ks1=kulka[0].sp*kulka[0].dp*kulka[0].D2FCp;

//BLOK 2------------------------------------------------------------------------       
       
       for(int i=0;i<LiczbaCzastek;i++)
       {
       kulka[i].xp=0;    
       kulka[i].yp=0;
       kulka[i].zp=0;        
               
       kulka[i].xp=kulka[i].x0+0.5*KrokCzasowy*w1X[i];
       kulka[i].yp=kulka[i].y0+0.5*KrokCzasowy*w1Y[i];
       kulka[i].zp=kulka[i].z0+0.5*KrokCzasowy*w1Z[i];
       
       kulka[i].xp=kulka[i].xp-sciana*round(kulka[i].xp*pow(sciana,-1));
       kulka[i].yp=kulka[i].yp-sciana*round(kulka[i].yp*pow(sciana,-1));
       kulka[i].zp=kulka[i].zp-sciana*round(kulka[i].zp*pow(sciana,-1));
       
       kulka[i].ppx=kulka[i].px0+0.5*KrokCzasowy*k1X[i];
       kulka[i].ppy=kulka[i].py0+0.5*KrokCzasowy*k1Y[i];
       kulka[i].ppz=kulka[i].pz0+0.5*KrokCzasowy*k1Z[i];
       
       kulka[i].dp=kulka[i].dzeta0+0.5*KrokCzasowy*kd1;
       kulka[i].sp=kulka[i].s0+0.5*KrokCzasowy*ks1;
       }
       

       stanP(kulka, LiczbaCzastek, sciana, zasieg, shift, listX, skora);
       
       
       SumCKF=0;
       for(int i=0;i<LiczbaCzastek;i++)
       {            
       w2X[i]=dziel(kulka[i].ppx,kulka[i].m)+kulka[0].dp*kulka[i].fxp;
       w2Y[i]=dziel(kulka[i].ppy,kulka[i].m)+kulka[0].dp*kulka[i].fyp;
       w2Z[i]=dziel(kulka[i].ppz,kulka[i].m)+kulka[0].dp*kulka[i].fzp;
       k2X[i]=kulka[i].fxp;
       k2Y[i]=kulka[i].fyp;
       k2Z[i]=kulka[i].fzp;
       SumCKF=SumCKF+kulka[i].CKFp;
       }
       
kd2=dziel(1,Q)*(SumCKF-temp*kulka[0].D2FCp);
ks2=kulka[0].sp*kulka[0].dp*kulka[0].D2FCp;       
     
       
//BLOK 3------------------------------------------------------------------------

       for(int i=0;i<LiczbaCzastek;i++)
       {
       kulka[i].xp=0;
       kulka[i].yp=0;
       kulka[i].zp=0;        
               
       kulka[i].xp=kulka[i].x0+0.5*KrokCzasowy*w2X[i];
       kulka[i].yp=kulka[i].y0+0.5*KrokCzasowy*w2Y[i];
       kulka[i].zp=kulka[i].z0+0.5*KrokCzasowy*w2Z[i];
       
       kulka[i].xp=kulka[i].xp-sciana*round(kulka[i].xp*pow(sciana,-1));
       kulka[i].yp=kulka[i].yp-sciana*round(kulka[i].yp*pow(sciana,-1));
       kulka[i].zp=kulka[i].zp-sciana*round(kulka[i].zp*pow(sciana,-1));
       
       kulka[i].ppx=kulka[i].px0+0.5*KrokCzasowy*k2X[i];
       kulka[i].ppy=kulka[i].py0+0.5*KrokCzasowy*k2Y[i];
       kulka[i].ppz=kulka[i].pz0+0.5*KrokCzasowy*k2Z[i];
       
       kulka[i].dp=kulka[i].dzeta0+0.5*KrokCzasowy*kd2;
       kulka[i].sp=kulka[i].s0+0.5*KrokCzasowy*ks2;
       }


       stanP(kulka, LiczbaCzastek, sciana, zasieg, shift, listX, skora);


       SumCKF=0;
       for(int i=0;i<LiczbaCzastek;i++)
       {            
       w3X[i]=dziel(kulka[i].ppx,kulka[i].m)+kulka[0].dp*kulka[i].fxp;
       w3Y[i]=dziel(kulka[i].ppy,kulka[i].m)+kulka[0].dp*kulka[i].fyp;
       w3Z[i]=dziel(kulka[i].ppz,kulka[i].m)+kulka[0].dp*kulka[i].fzp;
       k3X[i]=kulka[i].fxp;
       k3Y[i]=kulka[i].fyp;
       k3Z[i]=kulka[i].fzp;
       SumCKF=SumCKF+kulka[i].CKFp;
       }
       
kd3=dziel(1,Q)*(SumCKF-temp*kulka[0].D2FCp);
ks3=kulka[0].sp*kulka[0].dp*kulka[0].D2FCp;       
       
//BLOK 4------------------------------------------------------------------------  
       
       for(int i=0;i<LiczbaCzastek;i++)
       {
       kulka[i].xp=0;
       kulka[i].yp=0;
       kulka[i].zp=0;        
               
       kulka[i].xp=kulka[i].x0+KrokCzasowy*w3X[i];
       kulka[i].yp=kulka[i].y0+KrokCzasowy*w3Y[i];
       kulka[i].zp=kulka[i].z0+KrokCzasowy*w3Z[i];
       
       kulka[i].xp=kulka[i].xp-sciana*round(kulka[i].xp*pow(sciana,-1));
       kulka[i].yp=kulka[i].yp-sciana*round(kulka[i].yp*pow(sciana,-1));
       kulka[i].zp=kulka[i].zp-sciana*round(kulka[i].zp*pow(sciana,-1));
       
       kulka[i].ppx=kulka[i].px0+KrokCzasowy*k3X[i];
       kulka[i].ppy=kulka[i].py0+KrokCzasowy*k3Y[i];
       kulka[i].ppz=kulka[i].pz0+KrokCzasowy*k3Z[i];
       
       kulka[i].dp=kulka[i].dzeta0+KrokCzasowy*kd3;
       kulka[i].sp=kulka[i].s0+KrokCzasowy*ks3;
       }
       

       stanP(kulka, LiczbaCzastek, sciana, zasieg, shift, listX, skora);
       
       
       SumCKF=0;
       for(int i=0;i<LiczbaCzastek;i++)
       {            
       w4X[i]=dziel(kulka[i].ppx,kulka[i].m)+kulka[0].dp*kulka[i].fxp;
       w4Y[i]=dziel(kulka[i].ppy,kulka[i].m)+kulka[0].dp*kulka[i].fyp;
       w4Z[i]=dziel(kulka[i].ppz,kulka[i].m)+kulka[0].dp*kulka[i].fzp;
       k4X[i]=kulka[i].fxp;
       k4Y[i]=kulka[i].fyp;
       k4Z[i]=kulka[i].fzp;
       SumCKF=SumCKF+kulka[i].CKFp;
       }
       
kd4=dziel(1,Q)*(SumCKF-temp*kulka[0].D2FCp);
ks4=kulka[0].sp*kulka[0].dp*kulka[0].D2FCp; 
       
//NOWE POZYCJE I PREDKOSCI------------------------------------------------------
       
       for(int i=0;i<LiczbaCzastek;i++)
       {
       kulka[i].px=kulka[i].px0+dziel(KrokCzasowy,6)*(k1X[i]+2*k2X[i]+2*k3X[i]+k4X[i]);
       kulka[i].py=kulka[i].py0+dziel(KrokCzasowy,6)*(k1Y[i]+2*k2Y[i]+2*k3Y[i]+k4Y[i]);
       kulka[i].pz=kulka[i].pz0+dziel(KrokCzasowy,6)*(k1Z[i]+2*k2Z[i]+2*k3Z[i]+k4Z[i]);
       
       kulka[i].x=0;
       kulka[i].y=0;
       kulka[i].z=0;
       
       kulka[i].x=kulka[i].x0+dziel(KrokCzasowy,6)*(w1X[i]+2*w2X[i]+2*w3X[i]+w4X[i]);
       kulka[i].y=kulka[i].y0+dziel(KrokCzasowy,6)*(w1Y[i]+2*w2Y[i]+2*w3Y[i]+w4Y[i]);
       kulka[i].z=kulka[i].z0+dziel(KrokCzasowy,6)*(w1Z[i]+2*w2Z[i]+2*w3Z[i]+w4Z[i]);
       
       kulka[i].x=kulka[i].x-sciana*round(kulka[i].x*pow(sciana,-1));
       kulka[i].y=kulka[i].y-sciana*round(kulka[i].y*pow(sciana,-1));
       kulka[i].z=kulka[i].z-sciana*round(kulka[i].z*pow(sciana,-1));
       
       kulka[i].dzeta=kulka[i].dzeta0+dziel(KrokCzasowy,6)*(kd1+2*kd2+2*kd3+kd4);
       kulka[i].s=kulka[i].s0+dziel(KrokCzasowy,6)*(ks1+2*ks2+2*ks3+ks4);
       }       
       
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//DEKLARACJE ZMIENNYCH
//------------------------------------------------------------------------------

int const czastki=108;
int petla=0, skok, odczyt, potencjal, termostat, konsola, listX, dzielnik;

double iteracje, rozstaw, vp, wych, dm, dt, EnP, EnK, EnC;
double rho, Volume, sciana, zasieg, shift, skora;
class kuleLJ kulka[czastki];

double VX, VY, VZ;
double calka, Q, temp, g, kB;


double TEMP, CISNIENIE, TCONF, SUMd2Ur, SUMKWSIL;
double dx, dy, dz, r, PXX, PYY, PZZ;

double SUMTK, SUMTC, SUMP, SUME, TKS, TCS, PS, EPS;

//------------------------------------------------------------------------------      
//FUNKCJA GLOWNA
//------------------------------------------------------------------------------      
using namespace std;	

main()
{

g=(3*czastki)-3; 
kB=1;

temp=2;	
Q=100;
dt=0.00025;
iteracje=1000;
rho=0.5;
Volume=czastki*pow(rho,-1);
sciana=pow(Volume, pow(3,-1));	
skora=0.5;
listX=12;
	
odczyt=0;
skok=100;
konsola=2;//0-nic; 1-dane; 2-postep;
	
rozstaw=1;
vp=0;    
wych=0;
dm=0;  

termostat=1;//0-N; 1-NH; 2-TB;
potencjal=0;//0-LJ; 1-WCA; 2-KZLJ;

	if(potencjal==0)//LJ
	{
	shift=0; zasieg=2.5;	
	}
	
	else if(potencjal==1)//WCA
	{
	zasieg=1.122462048309373; shift=1;
	}
	
	else if(potencjal==2)//KZLJ
	{
	zasieg=2.5; shift=-4*(pow(zasieg,-12)-pow(zasieg,-6));
	}
//------------------------------------------------------------------------------
ofstream pozycje;
ofstream pedy;
ofstream energie;

pozycje.open("pozycjeNHT2R05");
if(!pozycje)return 0;

pedy.open("pedyNHT2R05");
if(!pedy)return 0;

energie.open("energieNHT2R05");
if(!energie)return 0;      

//------------------------------------------------------------------------------




//WCZYTYWANIE DANYCH STARTOWYCH-------------------------------------------------
if(odczyt==1)
{
fstream plik_poloz1("input.txt", ios::in | ios::out );
double pedpx, pedpy, pedpz, MasaCalkowita, DeltaMasy;
pedpx=0; pedpy=0; pedpz=0; MasaCalkowita=0;
int LiczbaCzastek=czastki;
        
		double pomoc1,pomoc2,pomoc3,pomoc4,pomoc5,pomoc6;
		for (int i=0; i<LiczbaCzastek; i++) 
		{
			kulka[i].m=1;
			kulka[i].licznik=listX;
			
			plik_poloz1 >> pomoc1 >> pomoc2 >> pomoc3;
            plik_poloz1 >> pomoc4 >> pomoc5 >> pomoc6;
			plik_poloz1.flush();
			kulka[i].x0=pomoc1;
			kulka[i].y0=pomoc2;
			kulka[i].z0=pomoc3;
            kulka[i].px0=pomoc4;
			kulka[i].py0=pomoc5;
			kulka[i].pz0=pomoc6;
		}
        plik_poloz1 >> pomoc1 >> pomoc2;
        for(int i=0;i<LiczbaCzastek;i++)
        {
		kulka[i].dzeta=pomoc1;
        kulka[i].s=pomoc2;
        
        kulka[i].dzeta0=pomoc1;
        kulka[i].s0=pomoc2;
        
        pedpx=pedpx+kulka[i].px0;
       	pedpy=pedpy+kulka[i].py0;
       	pedpz=pedpz+kulka[i].pz0;
       
       	MasaCalkowita=MasaCalkowita+kulka[i].m;
       
       	kulka[i].znacznik=i;
		}
		plik_poloz1.close();

        // Zerowanie poczatkowego ped
    
		for (int i=0; i<LiczbaCzastek; i++)
       {  
       kulka[i].px0=kulka[i].px0-dziel(pedpx,MasaCalkowita);
       kulka[i].py0=kulka[i].py0-dziel(pedpy,MasaCalkowita);
       kulka[i].pz0=kulka[i].pz0-dziel(pedpz,MasaCalkowita);
       
       kulka[i].x=kulka[i].x0;
       kulka[i].y=kulka[i].y0;
       kulka[i].z=kulka[i].z0;
       kulka[i].px=kulka[i].px0;
       kulka[i].py=kulka[i].py0;
       kulka[i].pz=kulka[i].pz0;
       
       kulka[i].dzeta=kulka[i].dzeta0;
       
       kulka[i].s=kulka[i].s0;
       }
}
else if(odczyt==0)
{
	ustaw(kulka,czastki,rozstaw,vp,wych,dm, sciana, temp, listX);   
	
}
//------------------------------------------------------------------------------      


     

  




//------------------------------------------------------------------------------  
//PETLA GLOWNA     
//------------------------------------------------------------------------------
SUMTK=0, SUMTC=0, SUMP=0, SUME=0; dzielnik=0;
for(int i=0;i<iteracje;i++)
{    
if(i==100000)
{
	dt=0.0005;
	SUMTK=0;
	SUMTC=0;
	SUMP=0;
	SUME=0;
	dzielnik=0;
}
//cout<<kulka[0].licznik<<endl;

stan(kulka, czastki, sciana, zasieg, shift, listX, skora); 


EnP=energiaP(kulka, czastki); 
EnK=energiaK(kulka, czastki); 
EnC=energiaC(EnP, EnK); 



//WIELKOSCI FIZYCZNE------------------------------------------------------------
TEMP=0, CISNIENIE=0, TCONF=0, SUMd2Ur=0, SUMKWSIL=0;
dx=0, dy=0, dz=0, PXX=0, PYY=0, PZZ=0;
TKS=0, TCS=0, PS=0, EPS=0;


TEMP=2*pow(g,-1)*EnK;

for(int k=0;k<czastki;k++)
{
	for(int j=0;j<czastki;j++)
	{
	dx=0; dy=0; dz=0; r=0;  
            
	dx=kulka[k].x-kulka[j].x;
	dy=kulka[k].y-kulka[j].y;
	dz=kulka[k].z-kulka[j].z; 
	
	dx=dx-sciana*round(dx*pow(sciana,-1));
	dy=dy-sciana*round(dy*pow(sciana,-1));
	dz=dz-sciana*round(dz*pow(sciana,-1));
    
	r=dx*dx+dy*dy+dz*dz;
		
	if(k!=j)
	
	if(r<zasieg*zasieg)
	{
	PXX+=FX(dx,dy,dz)*dx;
	PYY+=FY(dx,dy,dz)*dy;
	PZZ+=FZ(dx,dy,dz)*dz;  
	}
	}

SUMd2Ur=SUMd2Ur+kulka[k].D2FC;
SUMKWSIL=SUMKWSIL+kulka[k].CKF;	
}
//cout<<PXX+PYY+PZZ<<endl;
CISNIENIE=dziel(1,(3*Volume))*(0.5*(PXX+PYY+PZZ))+rho*TEMP;
TCONF=dziel(SUMKWSIL,SUMd2Ur); 
TCONF=dziel(TCONF,czastki);

SUMTK=SUMTK+TEMP;
SUMTC=SUMTC+TCONF;
SUMP=SUMP+CISNIENIE;
SUME=SUME+EnP;

TKS=dziel(SUMTK,(dzielnik+1));
TCS=dziel(SUMTC,(dzielnik+1));
PS=dziel(SUMP,(dzielnik+1));
EPS=dziel(SUME,(dzielnik+1));

dzielnik=dzielnik+1;

	if(termostat==0)
	{
		calka=dziel(EnC,czastki);
		NewtonRK4(kulka,czastki, dt, sciana, zasieg, shift, listX, skora);   
	}

	else if(termostat==1)
	{
		calka=(EnC+0.5*Q*kulka[0].dzeta*kulka[0].dzeta+g*kB*temp*log(kulka[0].s))*pow(czastki,-1);
		NoseHooverRK4(kulka, czastki, Q, temp, g, kB, dt, sciana, zasieg, shift, listX, skora);  
	}

	else if(termostat==2)
	{//calka=EnC;
		calka=(EnC+0.5*Q*kulka[0].dzeta*kulka[0].dzeta+kB*temp*log(kulka[0].s))*pow(czastki,-1);
		TravisBragaRK4(kulka, czastki, Q, temp, kB, dt, sciana, zasieg, shift, listX, skora);
	}

if (konsola==1)
{
	cout.precision(16);
	cout<<calka<<"   "<<TKS<<"   "<<EPS*pow(czastki,-1)<<"   "<<endl;     
}
else if (konsola==2)
{
	cout.precision(16);
	cout << "Postep: " << round(((i+1)*pow(iteracje,-1))*100) << "% \r";  
}

energie.precision(16);
energie<<calka<<"   "<<EnP<<"   "<<EnK<<"   "<<TEMP<<"   "<<TCONF<<"   "<<CISNIENIE<<"   "<<EPS<<"   "<<TKS<<"   "<<TCS<<"   "<<PS<<endl;
energie.flush();


//wielkosci termodynamiczne i wartosci srednie


//ZAPIS DO PLIKU

        petla=i-1;
        if(petla%skok==0)
        {
        pozycje.precision(16);

        pozycje<<i+1<<" "; //czas

        for(int k=0;k<czastki;k++){pozycje<<kulka[k].x<<" ";}
        for(int k=0;k<czastki;k++){pozycje<<kulka[k].y<<" ";}
        for(int k=0;k<czastki;k++){pozycje<<kulka[k].z<<" ";}

        pozycje<<" "<<endl;
        pozycje.flush();
        
        }
        


}     
fstream plik_poloz("wyjscie.txt", ios::out );
    if (plik_poloz.good()) {  
            for(int i=0;i<czastki;i++){
                plik_poloz.precision(16);
                plik_poloz << kulka[i].x<< " \t " << kulka[i].y << " \t " << kulka[i].z << "\n";
                plik_poloz << kulka[i].px<< " \t " << kulka[i].py << " \t " << kulka[i].pz << "\n";
                plik_poloz.flush();
            }
            plik_poloz << kulka[0].dzeta<< " \t " << kulka[0].s << "\n";
            //zamkniecie pliku
            plik_poloz.close();
        }


pozycje.close();
energie.close();
//cout<<"koniec"<<endl;
getch();     
}
