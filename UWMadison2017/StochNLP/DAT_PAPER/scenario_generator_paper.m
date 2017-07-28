clc
clear all
close all hidden

% set problem dimensions
   TF=24*3600;   % horizon time - [s] 
   Nt=24;        % number of temporal grid points
 TDEC=10;        % decision time step
   dt=TF/Nt;
    S=50;
 Mnom=41.66;  
        
   Tst=dt*(TDEC+2); % mean start of step [s]
  Tstm=dt*(TDEC-1); 
  Tstp=dt*(TDEC+1);  
   Tdu=5*3600;      % mean duration time [s]
  Tdum=4*3600;
  Tdup=6*3600;
    Ma=Mnom*1.2;   % mean magnitude of demand 
   Mam=Mnom*1.1;
   Map=Mnom*1.3;  
   
 Tstss = Tstm + (Tstp-Tstm).*rand(S,1); 
 Tduss = Tdum + (Tdup-Tdum).*rand(S,1); 
  Mass =  Mam + (Map-Mam).*rand(S,1); 
   tc = 0:dt:TF-dt;
   
 % first segment
 for s = 1:S
 
 ddt=TF/300;    
     
 t1=linspace(0,Tstss(s));
 y1=linspace(Mnom,Mnom);
 
 t2=linspace(Tstss(s)+ddt,Tstss(s)+Tduss(s));
 y2=linspace(Mass(s),Mass(s));
 
 t3=linspace(Tstss(s)+Tduss(s)+ddt,TF);
 y3=linspace(Mnom,Mnom);
 
 tf=[t1 t2 t3];
  y=[y1 y2 y3];
 
 dem(s,:)=interp1(tf,y,tc,'nearest');
 
 end
    
 figure(1) 
 for s=1:S
 stairs(tc,dem(s,:))
 hold on
 end
 
 figure(2)
 plot(tc,dem)
 
 save stoch_scen_paper.dat dem -ascii