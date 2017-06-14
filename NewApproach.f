      program NewApproach
      implicit double precision(a-h,o-z)
      parameter (idim=10)
      character*1 relabs
      
      double precision xx(10),w(2500)
      double precision g1,g2,g3,g4,g5,g6,g7,g8,g9,gext1,gext2,gext3
      double precision V1,V2,V3,V4,V5,V6,V7,V8,V9,Vext1,Vext2,Vext3
      double precision I1,I2,I3,I4,I5,I6,I7,I8,I9,Iext1,Iext2,Iext3
      double precision C, A, tau4, K, Ts1, Tx2, Ts2, Ts3, TCa, aCa, Kons
      double precision ccKo, ccKi, ccNao, ccNai
      double precision ccClo, ccCli, ccCao, ccCai,ccMgo
      double precision a1,b1,a2,b2,a3,b3,a4,a5,b5
      double precision m2,m4,m6,m7,m8,h9
      integer iz

      common g1,g2,g3,g4,g5,g6,g7,g8,g9,gext1,gext2,gext3
      common V1,V2,V3,V4,V5,V6,V7,V8,V9,Vext1,Vext2,Vext3
      common I1,I2,I3,I4,I5,I6,I7,I8,I9,Iext1,Iext2,Iext3
      common C, A, tau4, K, Ts1, Tx2, Ts2, Ts3, TCa, aCa, Kons;
      common ccKo, ccKi, ccNao, ccNai, ccClo, ccCli, ccCao, ccCai,ccMgo
      common a1,b1,a2,b2,a3,b3,a4,a5,b5
      common m2,m4,m6,m7,m8,h9
      common k1
      external fcn,d02cjx,d02cjw
      cl = 1;
      tol=10.d0**(- 6.d0)
      relabs='m'
      dt=.1d0
      nsteps=100000


      open(10,file='TransForm.dat')
      do i=0,15

!     Initial values for variables                                                                                            
         xx(1)= -75.0 + 5*15;
         xx(2)=0.51  + 0.000001;
         xx(3)=0.21 + 0.000001;  xx(4)=0.84 + 0.000002;  
         xx(5)=0.31 + 0.000001;  xx(6)=0.11 + 0.000002;
         xx(7)=0.51 + 0.000001;  xx(8)=0.11 + 0.000001;  
         xx(9)=0.21 + 0.000001;
         xx(10) = 5.5 + 0.000001;
         sw = 0
         
!     DO SCRIPT                                                                                                               
         iflag=0
         num=0
         delt=0.d0
         cpoin = 0
         cl = 0
         k1 = i
         do icont=1,nsteps
            ti=dfloat(icont)*dt
            x=ti
            xend=dfloat(icont+1)*dt
            
            call d02cjf(x,xend,idim,xx,fcn,tol,relabs,
     &           d02cjx,d02cjw,w,ifail)
            
            if (xx(1) < - 60.0 .AND. sw > 0.5) then
               sw = 0;
               cpoin=cpoin + 1;
               ti0 = ti
            end if
            
            if (xx(1) > 0.0 .AND. sw < 0.5) then
               sw = 1
            end if
            

            if (icont.gt.cl*10) then
               write(int(100+i),*) ti,xx(1),xx(2),xx(3)
               cl=cl+1
            end if
         end do
         write(6,*) i,cpoin, ti
         write(10,*) cpoin

      end do
      end


      subroutine fcn(t,xx,f)
      implicit double precision(a-h,o-z)
      double precision xx(10),f(10)
      double precision g1,g2,g3,g4,g5,g6,g7,g8,g9,gext1,gext2,gext3
      double precision V1,V2,V3,V4,V5,V6,V7,V8,V9,Vext1,Vext2,Vext3
      double precision I1,I2,I3,I4,I5,I6,I7,I8,I9,Iext1,Iext2,Iext3
      double precision a1,b1,a2,b2,a3,b3,a4,a5,b5,m2,m4,m6,m7,m8,h9
      double precision C, A, tau4, K, Ts1, Tx2, Ts2, Ts3, TCa, aCa, Kons
      double precision ccKo, ccKi, ccNao, ccNai, ccClo
      double precision ccCli, ccCao, ccCai,ccMgo
      
      common g1,g2,g3,g4,g5,g6,g7,g8,g9,gext1,gext2,gext3
      common V1,V2,V3,V4,V5,V6,V7,V8,V9,Vext1,Vext2,Vext3
      common I1,I2,I3,I4,I5,I6,I7,I8,I9,Iext1,Iext2,Iext3
      common C, A, tau4, K, Ts1, Tx2, Ts2, Ts3, TCa, aCa, Kons;
      common ccKo, ccKi, ccNao, ccNai, ccClo, ccCli, ccCao, ccCai,ccMgo;
      common a1,b1,a2,b2,a3,b3,a4,a5,b5
      common m2,m4,m6,m7,m8,h9
      common k1

! Constants                                                                                                               
      C = 1.0;  A = 0.02;  tau4 = 15.0; K = 30.0;
      Ts1 = 2.0; Tx2 = 2.0; Ts2 = 100.0; Ts3 = 10.0; TCa = 121.4;
      aCa = 0.5; Kons = 26.7137;
      
! Gating Parameters                                                                                                       
      g1 = 0.03573; g2 = 2.5; g3 = 2.61868; g4 = 1.79259; g5 = 0.0350135
      g6 = 0.0256867; g7 = 2.5; g8 = 0.0717984; g9 = 0.0166454;
      gext1 = 0.513425; gext2 = 0.00434132; gext3 = 0.00252916;
      
      
!      g7 = 2.5*(1.0 - 0.025*k1);                                                                                             
      g7 = 2.5*(1.0 - 0.027)
      ccKo = (3.5 + 0.1*k1)*1000;
      ccCao = (1.35 - 0.0168*k1)*1000;
      ccMgo = (0.8 - 0.011*k1)*1000;

!      ccKo = (3.5 + 0.1*k1)*1000;
!      ccCao = (1.35 - 0.0167*k1)*1000;
!      ccMgo = (0.8 - 0.0111*k1)*1000;

      ccKi = 140.*1000;
      ccNao= 140.*1000; 
      ccNai= 7.0*1000;
      ccClo= 140.*1000; 
      ccCli= 7.0*1000;
      ccCai = xx(10)

      V1 = Kons*log((ccKo+0.08*ccNao+0.1*ccCli)/
     &(ccKi+0.08*ccNai+0.1*ccClo));
      V2 = Kons*log((ccNao)/(ccNai)); V3=Kons*log((ccKo)/(ccKi));
      V4 = Kons*log((ccKo)/(ccKi)); V5=Kons*log((ccKo)/(ccKi));
      V7 = Kons*log((ccKo)/(ccKi)); V8=Kons*log((ccNao)/(ccNai));
      V9 = Kons*log((ccKo)/(ccKi));
      Vext1=Kons*log((ccKo+ccNao)/(ccKi+ccNai));
      Vext2=Kons*log((ccKo+ccNao+ccCao)/(ccKi+ccNai+ccCai));
      Vext3= - 18.50;


      a1 = 0.1*(xx(1)+33.)/(1-exp(-(xx(1)+33.0)/10.0));
      b1 = 4.0*exp(-(xx(1) + 53.7)/12.0);
      a2 = 0.07*exp(-(xx(1)+50.)/10.0);
      b2 = 1.0/(1.0 + exp(-(xx(1)+20.0)/10.0));
      a3 = 0.01 * (xx(1) + 34.0)/(1 - exp(-(xx(1) + 34.0)/10.0));
      b3 = 0.125*exp(-(xx(1) + 44.0)/25.0);
      a4 = 1.0/( 1.0 + exp( (xx(1) + 80.0)/ 6.0 )  );
      a5 = 1.0/( 1.0 + exp( -( xx(1) + 34.0 )/6.5  )  );
      b5 = 8.0/(exp(-(xx(1)+55.0)/30.0)+exp((xx(1)+55.0)/30.0));
      
      m2 = a1/(a1 + b1);
      m4 = 1.0/( 1.0 + exp(-(xx(1) + 50.0 )/20.0 ) );
      m6 = 1.0/( 1.0 + exp(-( xx(1) + 20.0)/9.0) );
      m7 = 1.0/( 1.0 + (K/xx(10))**(3.5));
      m8 = 1.0/(1.0 + exp(-( xx(1) + 55.7 )/7.7));
      h9 = 1.0/( 1.0 + exp( (xx(1) + 75.0)/4.0));
      fV = 1.0/(1.0 + exp(-(xx(1) - 20.0)/2.0));

      
      V6 = Kons*log( (ccCao) / (xx(10)));
      I1=g1*(xx(1) - V1);
      I2=g2*m2*m2*m2*xx(2)*(xx(1)-V2);
      I3=g3*xx(3)*xx(3)*xx(3)*xx(3)*(xx(1)-V3);
      I4=g4*m4*m4*m4*xx(4)*(xx(1) - V4);
      I5=g5*xx(5)*(xx(1)-V5);
      I6=g6*m6*m6*(xx(1)-V6);
      I7=g7*m7*(xx(1)-V7);
      I8=g8*m8*m8*m8*(xx(1)-V8);
      I9=g9*h9*(xx(1)-V9);
      Iext1=gext1*xx(6)*(xx(1)-Vext1);
      Iext2= 1.1/(1.0 + ccMgo/(8.0*1000))*gext2*xx(8)*(xx(1)-Vext2);
      Iext2= gext2*xx(8)*(xx(1)-Vext2);
      Iext3=gext3*xx(9)*(xx(1)-Vext3);
      
      f(1) = -10.0*A*(I1+I2+I3+I4+I5+I6+I7+I8+I9)/(10.0*A*C)-
     & (Iext1+Iext2+Iext3)/(10.0*A*C);
        f(2) = 4.*a2*(1 - xx(2)) - b2*xx(2);
        f(3) = 4.0*(a3*(1.0 - xx(3)) - b3*xx(3));
        f(4) = (a4 - xx(4))/tau4;
        f(5) = (a5 - xx(5))/b5;
        f(6) =  3.48*fV - xx(6)/Ts1;
        f(7) = 3.48*fV - xx(7)/Tx2;
        f(8) = 0.5*xx(7)*(1.0 - xx(8)) - xx(8)/Ts2;
        f(9) = fV - xx(9)/Ts3;
        f(10) =  -aCa*(10.0*A*I6 + Iext2 ) - xx(10)/TCa;

        return

        end
