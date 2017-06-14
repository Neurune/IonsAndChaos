#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <ctype.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <random>
using namespace std;



int main(){    

  for (unsigned int ig = 0; ig < 1; ig++){ 
    // cout << ig << endl;
    double dt = 0.001;
    double TimeStop = 300000.0;
    double time = 0.0;

    ostringstream Newname;
    Newname << "Datafiles/CheckOnlyAll_Poin_cc"; Newname << ig; Newname << ".txt";
    std::ofstream Newfile (Newname.str().c_str());
    
    ostringstream Trajname;
    Trajname << "Datafiles/CheckOnlyAll_Traj_cc"; Trajname << ig; Trajname << ".txt";
    std::ofstream Trajfile (Trajname.str().c_str());
    
    int click = 1;
    
    ///% Currents
    double I [9];
    double Iext [9];
    double VKrk[10]; 
    double vk[10][4];
    ///% Variables
    double m2, h2, n3, m4, h4, m5, m6, m7, m8, h9, s1, x2, s2, s3, Ca, Kpl;
    m2 = 0.01; h2 = 0.045; n3 = 0.54; m4 = 0.1; h4 = 0.045; 
    m5 = 0.34; m6 = 0.01; m7 = 0.01; m8 = 0.01; h9 = 0.01;
    s1 = 0.01; x2 = 0.01; s2 = 0.01; s3 = 0.01; Ca = 1.0;
    ///% Conductances
    double g [9];
    g[0] = 0.03573; g[1] = 12.2438; g[2] = 2.61868; g[3] = 1.79259; g[4] = 0.0350135; 
    g[5] = 0.0256867; g[6] = 2.34906; g[7] = 0.0717984; g[8] = 0.0166454;

    g[1] = 2.5;
    g[6] = 2.3;

    int cpoin = 0;
    int sw = 0;
    double Poin1, Poin2;

    g[ig] = g[ig]*1.0;

    double Koncen [10];
    Koncen[0] = 3.9*1000;    Koncen[1] = 140*1000;    Koncen[2] = 140*1000;    Koncen[3] = 7*1000;
    Koncen[4] = 140*1000;    Koncen[5] = 7*1000;    Koncen[6] = 1.35*1000;    Koncen[7] = 1.20*1000;
    Koncen[8] = 0.8*1000;     Koncen[9] = 0.8*1000;

    double gext [3];
    gext[0] = 0.513425; gext[1] = 0.00434132; gext[2] = 0.00252916;
    
    ///% Potentials (in mV!!)
    double V [9]; 
    //    double Vr = -45.0;
    double Vr = -45.0 + 0.000001*ig;
    Ca = 1.0;
    V[0] = -60.95; V[1] = 55.0; V[2] = -100.0; V[5] = 120.0;
    V[3] = V[2]; V[4] = V[2]; V[6] = V[2]; V[7] = V[1]; V[8] = V[2];
    double Vext [3]; Vext[0] = 0; Vext[1] = 0; Vext[2] = -70.0;
    
    ///% Other Constants
    double C, A, tau4, K, Ts1, Tx2, Ts2, Ts3, TCa, aCa;
    C = 1.0;  A = 0.02;  tau4 = 15.0; K = 30.0;
    Ts1 = 2.0; Tx2 = 2.0; Ts2 = 100.0; Ts3 = 10.0; TCa = 121.4;
    aCa = 0.5;
    
    ///% New Parameters
    double aKpl, TKpl, Kons, fV;
    aKpl = 0.9;  TKpl = 10000.0;  Kons = 26.7137;
    
    double ccKo, ccKi, ccNao, ccNai, ccClo, ccCli, ccCao, ccCai, ccMgo, ccMgi, ccKoNew, ccCaoNew, ccMgoNew;
    double UpK, UpNa, UpCl, a, b, UpCa, Vk1;
    double k1,k2,k3,k4;
        
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0.0 , 4.00000001);
    
    std::random_device rd2;
    std::mt19937 gen2(rd2());
    std::normal_distribution<> d2(0.0 , 0.4);
    
    int swpea = 0;

    ccKo = Koncen[0]; 
    ccKi = Koncen[1]; /// Kalium
    ccNao = Koncen[2]; 
    ccNai = Koncen[3]; /// Natrium
    ccClo = Koncen[4]; 
    ccCli = Koncen[5]; /// Chlorid
    ccCao = Koncen[6]; ///% Calcium
    ccMgi = Koncen[8];
    ccMgo = Koncen[9];

    int cliswi = 0; double Vmax = -1000.0; double TimeMax = 0; double TimeMax0 = 0;

    ///%% Here we change the concentrations
    ccKoNew = (4.4)*1000; 
    ccKo = (4.4)*1000; 
    ccKi = Koncen[1]; /// Kalium
    ccNao = Koncen[2]; 
    ccNai = Koncen[3]; /// Natrium
    ccClo = Koncen[4]; 
    ccCli = Koncen[5]; /// Chlorid
    ccCao = Koncen[6]; ///% Calcium
    ccCao = 1.20*1000; ///% Calcium
    ccCaoNew = (1.20)*1000;
    g[6] = 2.3*(1.0);
    ccMgoNew = (0.7)*1000;
    ccMgo = (0.7)*1000;
    /////////////////////////////////////////
    
    
    while (time < TimeStop){
      time = time + dt;

      for (unsigned int rk = 0; rk< 4; rk++){
	if (rk == 0){
	  VKrk[0] = Vr;   VKrk[1] = Ca;   VKrk[2] = h2;   VKrk[3] = n3;
	  VKrk[4] = h4;   VKrk[5] = m5;   VKrk[6] = s1;   VKrk[7] = x2;
	  VKrk[8] = s2;           VKrk[9] = s3;
	  
	}
	else if (rk == 1){
	  VKrk[0] = Vr + 0.5*dt*vk[0][0];    VKrk[1] = Ca + 0.5*dt*vk[1][0];    VKrk[2] = h2 + 0.5*dt*vk[2][0];    VKrk[3] = n3 + 0.5*dt*vk[3][0];
	  VKrk[4] = h4 + 0.5*dt*vk[4][0];    VKrk[5] = m5 + 0.5*dt*vk[5][0];    VKrk[6] = s1 + 0.5*dt*vk[6][0];    VKrk[7] = x2 + 0.5*dt*vk[7][0];
	  VKrk[8] = s2 + 0.5*dt*vk[8][0];    VKrk[9] = s3 + 0.5*dt*vk[9][0];
	  
	}
	else if (rk == 2){
	  VKrk[0] = Vr + 0.5*dt*vk[0][1];    VKrk[1] = Ca + 0.5*dt*vk[1][1];    VKrk[2] = h2 + 0.5*dt*vk[2][1];    VKrk[3] = n3 + 0.5*dt*vk[3][1];
	  VKrk[4] = h4 + 0.5*dt*vk[4][1];    VKrk[5] = m5 + 0.5*dt*vk[5][1];    VKrk[6] = s1 + 0.5*dt*vk[6][1];    VKrk[7] = x2 + 0.5*dt*vk[7][1];
	  VKrk[8] = s2 + 0.5*dt*vk[8][1];    VKrk[9] = s3 + 0.5*dt*vk[9][1];
	  
	}
	else if (rk == 3){
	  VKrk[0] = Vr + dt*vk[0][2];        VKrk[1] = Ca + dt*vk[1][2];        VKrk[2] = h2 + dt*vk[2][2];        VKrk[3] = n3 + dt*vk[3][2];
	  VKrk[4] = h4 + dt*vk[4][2];        VKrk[5] = m5 + dt*vk[5][2];        VKrk[6] = s1 + dt*vk[6][2];        VKrk[7] = x2 + dt*vk[7][2];
	  VKrk[8] = s2 + dt*vk[8][2];        VKrk[9] = s3 + dt*vk[9][2];
	}
	ccCai = VKrk[1];
	UpK = 1.00; UpNa = 0.08; UpCl = 0.1; 	

	//////% Current 1
	V[0] = Kons*log( (UpK*ccKo + UpNa*ccNao + UpCl*ccCli) / (UpK*ccKi + UpNa*ccNai + UpCl*ccClo));
	I[0] = g[0]*(VKrk[0] - V[0]);
	
	a = 0.1*(VKrk[0] + 33.)/(1 - exp(-(VKrk[0] + 33.)/10.0));
	b = 4.0*exp(-(VKrk[0] + 53.7)/12.0);
	m2 = a/(a + b);    
	a = 0.07*exp(-(VKrk[0] + 50.)/10.0);
	b = 1.0/(1.0 + exp(-(VKrk[0] + 20.0)/10.0));
	vk[2][rk] = dt*( 4.*a*(1 - VKrk[2]) - b*VKrk[2] ); 
	
	//////% Current 2
	V[1] = Kons*log( (ccNao) / (ccNai)); ///% 
	I[1] = g[1]*(pow(m2 , 3))*VKrk[2]*(VKrk[0]-V[1]);
	
	//////% Current 3
	
	a = 0.01 * (VKrk[0] + 34.0)/(1 - exp(-(VKrk[0] + 34.0)/10.0));
	b = 0.125*exp(-(VKrk[0] + 44.0)/25.0);
	vk[3][rk] = dt*( 4.0*(a*(1.0 - VKrk[3]) - b*VKrk[3]) );
	V[2] = Kons*log( (ccKo) / (ccKi));     ///% ?????
	I[2] = g[2]*(pow(VKrk[3],4) )*(VKrk[0]-V[2]);
	
	//////% Current 4 
	m4 = 1.0/( 1.0 + exp(-(VKrk[0] + 50.0 )/20.0 ) );
	a = 1.0/( 1.0 + exp( (VKrk[0] + 80.0)/ 6.0 )  );
	vk[4][rk] = dt*( (a - VKrk[4])/tau4 );
	V[3] = Kons*log( (ccKo) / (ccKi));
	I[3] = g[3]*(pow(m4 , 3) )*VKrk[4]*(VKrk[0] - V[4]);
	
	//////% Current 5   
	a = 1.0/( 1.0 + exp( -( VKrk[0] + 34.0 )/6.5  )  );
	b = 8.0/( exp(-( VKrk[0] + 55.0)/30.0 ) + exp(( VKrk[0] + 55.0)/30.0 )  );
	vk[5][rk] = dt*( (a - VKrk[5])/b );
	V[4] = Kons*log( (ccKo) / (ccKi));
	I[4] = g[4]*VKrk[5]*(VKrk[0] - V[4]);
	
	//////% Current 6     
	
	m6 = 1.0/( 1.0 + exp(-( VKrk[0] + 20.0)/9.0) );
	V[5] = Kons*log( (ccCaoNew) / (ccCai));
	I[5] = g[5]*(pow(m6 , 2) )*(VKrk[0] - V[5]);	
	vk[1][rk] = dt*( -aCa*(10.0*A*I[5] + Iext[1] ) - VKrk[1]/TCa);	
	//////% Current 7
	
	m7 = 1.0/( 1.0 + pow(K/VKrk[1] , 3.5 ) );
	V[6] = Kons*log( (ccKo) / (ccKi));
	I[6] = g[6]*m7*(VKrk[0] - V[6]);
	
	//////% Current 8     
	m8 = 1.0/(1.0 + exp(-( VKrk[0] + 55.7 )/7.7 )   );
	V[7] = Kons*log( (ccNao) / (ccNai));
	I[7] = g[7]*pow(m8,3)*(VKrk[0]-V[7]);
	//////% Current 9    
	h9 = 1.0/( 1.0 + exp( (VKrk[0] + 75.0)/4.0)  );
	V[8] = Kons*log( (ccKo) / (ccKi));
	I[8] = g[8]*h9*(VKrk[0] - V[8]);
	
	////////////%%  Extrinsic currents  
	
	fV = 1.0/(1.0 + exp(-(VKrk[0] - 20.0)/2.0));
	
	vk[6][rk] = dt*( 3.48*fV - VKrk[6]/Ts1 );
	UpK = 1.0; UpNa = 1.0; 
	Vext[0] = Kons*log( (UpK*ccKo + UpNa*ccNao) / (UpK*ccKi + UpNa*ccNai) );
	Iext[0] = gext[0]*VKrk[6]*(VKrk[0] - Vext[0]);
	
	vk[7][rk] = dt*( 3.48*fV - VKrk[7]/Tx2 );
	vk[8][rk] = dt*( 0.5*VKrk[7]*(1.0 - VKrk[8]) - VKrk[8]/Ts2 );
	UpK = 1.0; UpNa = 1.0; UpCa = 1.0;
	Vext[1] = Kons*log( (UpK*ccKo + UpNa*ccNao + UpCa*ccCao) / (UpK*ccKi + UpNa*ccNai + UpCa*ccCai) );
	Iext[1] = 1.1/(1.0 + ccMgo/(8.0*1000) )*gext[1]*VKrk[8]*(VKrk[0] - Vext[1]);
	
	vk[9][rk] = dt*( fV - VKrk[9]/Ts3 );
	Iext[2] = gext[2]*VKrk[3]*(VKrk[0] - Vext[2]);


	vk[0][rk] = dt*(-10.0*A*(I[0] + I[1] +I[2] +I[3] +I[4] +I[5] +I[6] +I[7] +I[8])/(10.0*A*C) - (Iext[0] + Iext[1] + Iext[2])/(10.0*A*C));	
      }
      
      
      Vr = Vr + 1.0/6.0*(vk[0][0] + 2.0*vk[0][1] + 2.0*vk[0][2] + vk[0][3]) + dt*d(gen);
      //      Ca = Ca + 1.0/6.0*(vk[1][0] + 2.0*vk[1][1] + 2.0*vk[1][2] + vk[1][3]) + dt*d2(gen2);
      Ca = Ca + 1.0/6.0*(vk[1][0] + 2.0*vk[1][1] + 2.0*vk[1][2] + vk[1][3]) + dt*d2(gen2);
      if (Ca < 0){ Ca = 0.01;}
      h2 = h2 + 1.0/6.0*(vk[2][0] + 2.0*vk[2][1] + 2.0*vk[2][2] + vk[2][3]);
      n3 = n3 + 1.0/6.0*(vk[3][0] + 2.0*vk[3][1] + 2.0*vk[3][2] + vk[3][3]);
      h4 = h4 + 1.0/6.0*(vk[4][0] + 2.0*vk[4][1] + 2.0*vk[4][2] + vk[4][3]);
      m5 = m5 + 1.0/6.0*(vk[5][0] + 2.0*vk[5][1] + 2.0*vk[5][2] + vk[5][3]);
      s1 = s1 + 1.0/6.0*(vk[6][0] + 2.0*vk[6][1] + 2.0*vk[6][2] + vk[6][3]);
      x2 = x2 + 1.0/6.0*(vk[7][0] + 2.0*vk[7][1] + 2.0*vk[7][2] + vk[7][3]);
      s2 = s2 + 1.0/6.0*(vk[8][0] + 2.0*vk[8][1] + 2.0*vk[8][2] + vk[8][3]);
      s3 = s3 + 1.0/6.0*(vk[9][0] + 2.0*vk[9][1] + 2.0*vk[9][2] + vk[9][3]);

      if (time > click*20.02){
	Trajfile << Vr << "\t" << Ca << "\t" << time/1000.0 << "\t"<< n3 << "\n";
	click++;
      }

      if (Vr < - 60.0 && sw > 0.5){	
	sw = 0;
	Newfile << time/1000.0 << "\t" << time/1000.0 - Poin1 << "\t" << time/1000.0 - Poin2 << "\t" << Poin2 << "\t" << Poin2 - Poin1 << "\n";
	cout << time << " " << cpoin << endl;
        Poin1 = time/1000.0;         cpoin++;
      }
      if (Vr > 0.0 && sw < 0.5){	
	Poin2 = time/1000.0;
	sw = 1;
      }
     
    }
    cout << cpoin << endl;
  }
} // This ends main
