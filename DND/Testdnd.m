close all
clear all

% global N
% global alphac

N = 500;
alphac = 0.02;

fixdim = 36;

fname = @f_dnd;
Jname = @J_dnd;
% n=500;
ptsx=50.*(1:N)/(N+1)+(1-(1:N)/(N+1)).*(-23);
sgx=sign(ptsx);
x0=((sgx+1)/2).*exp(-0.8284.*ptsx)+(1-sgx)/2;
clear ptsx;
clear sgx;
IT=[0 0.1];

ttt=fname(0,x0);
ttt=J_dnd(0,x0);
clear ttt;

clear tt;
Ntrials = 1;
optionsS  = odeset('RelTol',1.0e-3,'AbsTol',1.0e-6);  
optionsLL  = llset('RelTol',1.0e-9,'AbsTol',1.0e-12,'dKmax',30,...
    'dKmin',4,'debug',0,'gamma',0.1);
optionsLL2  = llset('RelTol',1.0e-9,'AbsTol',1.0e-12,'dKmax',50,...
    'dKmin',4,'gamma',0.1);
optionsLLe  = odeset('RelTol',1.0e-12,'AbsTol',1.0e-14);
options15s2 = odeset('RelTol',1.0e-9,'AbsTol',1.0e-12,'Jacobian',Jname);
options15s3 = odeset('RelTol',1.0e-9,'AbsTol',1.0e-12);
optionsRK  = odeset('RelTol',1.0e-9,'AbsTol',1.0e-12);
options15s = odeset('RelTol',1.0e-12,'AbsTol',1.0e-14,'Jacobian',Jname);


% tic
% %for i=1:Ntrials
% %  SolLL3Kp = LLDP_Kphi1(fname,Jname,IT,x0,optionsLL);
% SolLL3Kp = LLDP_Kphi1(fname,Jname,IT,x0,optionsLL);
% %end;
% tocLL3Kp=toc;
% % 
% TLL3Kp = SolLL3Kp.x;
% YLL3Kp = real(SolLL3Kp.y)';
% [T,Y] = ode15s(fname,TLL3Kp,x0,options15s);
% % [T,Y] = exactsol(fname,Jname,TLL3Kp, x0);
%   Y = real(Y);
% tocLL3Kp  
% LL3KpRE = RelError(Y,YLL3Kp)
% SolLL3Kp.stats
% clear T
% clear Y
% clear LLDP_Kphi1
% clear SolLL3Kp
% clear phi1LLDP_J_new

% tic
% %for i=1:Ntrials
%  SolLL3Kp = DLLRK45_3_Auto_Kphi1(fname,Jname,IT,x0,optionsLL);
% %end;
% tocLL3Kp=toc;
% 
% TLL3Kp = SolLL3Kp.x;
% YLL3Kp = real(SolLL3Kp.y)';
% [T,Y] = ode15s(fname,TLL3Kp,x0,options15s);
% % [T,Y] = exactsol(fname,Jname,TLL3Kp, x0)
%   Y = real(Y);
% tocLL3Kp  
% LL3KpRE = RelError(Y,YLL3Kp)
% SolLL3Kp.stats
% clear T
% clear Y
% clear DLLRK45_3_Auto_Kphi1
% clear SolLL3Kp

% de 0.04 a 0.1
% [T,Y] = ode15s(fname,IT,x0,options15s);
% x0 = Y(end, :);
% IT = [0.07,0.1];
% 
% tic
% %for i=1:Ntrials
%  SolLL3Kp = DLLRK45_3_Auto_Kphi1(fname,Jname,IT,x0,optionsLL);
% %end;
% tocLL3Kp=toc;
% 
% TLL3Kp = SolLL3Kp.x;
% YLL3Kp = real(SolLL3Kp.y)';
% [T,Y] = ode15s(fname,TLL3Kp,x0,options15s);
% % [T,Y] = exactsol(fname,Jname,TLL3Kp, x0)
%   Y = real(Y);
% tocLL3Kp  
% LL3KpRE = RelError(Y,YLL3Kp)
% SolLL3Kp.stats



tic

 SolLL3Kpj = LLDP_Kphi1_freeJ(fname,IT,x0,optionsLL2);

tocLL3Kpj=toc;

TLL3Kpj = SolLL3Kpj.x;
YLL3Kpj = real(SolLL3Kpj.y)';
[T,Y] = ode15s(fname,TLL3Kpj,x0,options15s);
% [T,Y] = exactsol(fname,Jname,TLL3Kpj,x0,40,optionsLLe);
  Y = real(Y);
tocLL3Kpj  
LL3KpjRE = RelError(Y,YLL3Kpj)
SolLL3Kpj.stats

% 
% tic
% Solode15 = ode15s(fname,IT,x0,options15s2);
% tocode15=toc;
% 
% T = Solode15.x;
% Yode15 = real(Solode15.y)';
% [T,Y] = ode15s(fname,T,x0,options15s);
% % [T,Y] = exactsol(fname,Jname,TLL3Kpj,x0,40,optionsLLe);
%      Y = real(Y);
% tocode15  
% ode15RE = RelError(Y,Yode15)
% Solode15.stats
% 
% tic
% Solode15j2 = ode15sk(fname,IT,x0,options15s2);
% tocode15j2=toc;
% 
% T = Solode15j2.x;
% Yode15 = real(Solode15j2.y)';
% [T,Y] = ode15s(fname,T,x0,options15s);
% % [T,Y] = exactsol(fname,Jname,TLL3Kpj,x0,40,optionsLLe);
%      Y = real(Y);
% tocode15j2
% ode15jRE = RelError(Y,Yode15)
% Solode15j2.stats

% tic
% Solode15j3 = ode15sk_no_pre(fname,IT,x0,options15s2);
% tocode15j3=toc;
% 
% T = Solode15j3.x;
% Yode15 = real(Solode15j3.y)';
% [T,Y] = ode15s(fname,T,x0,options15s);
% % [T,Y] = exactsol(fname,Jname,TLL3Kpj,x0,40,optionsLLe);
%      Y = real(Y);
% tocode15j3
% ode15jRE = RelError(Y,Yode15)
% Solode15j3.stats
