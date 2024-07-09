close all
clear all

N = 500;
alphac = 0.02;

fixdim = 36;

fname = @f_brusselator;
Jname = @J_brusselator;
x0= [1+sin((2*pi/(N+1))*(1:N)),repmat(3,1,N)];
IT=[0 1];


Ntrials = 1;
optionsS  = odeset('RelTol',1.0e-3,'AbsTol',1.0e-6);  
optionsLL  = llset('RelTol',1.0e-6,'AbsTol',1.0e-9,'dKmax',30,...
    'dKmin',4,'debug',1,'gamma',0.05);
optionsLL2  = llset('RelTol',1.0e-9,'AbsTol',1.0e-12,'dKmax',30,...
    'dKmin',5,'debug',1,'gamma',0.05);
optionsLLe  = odeset('RelTol',1.0e-12,'AbsTol',1.0e-14);
options15s2 = odeset('RelTol',1.0e-9,'AbsTol',1.0e-12,'Jacobian',Jname);
options15s3 = odeset('RelTol',1.0e-9,'AbsTol',1.0e-12);
optionsRK  = odeset('RelTol',1.0e-9,'AbsTol',1.0e-12);
options15s = odeset('RelTol',1.0e-12,'AbsTol',1.0e-14,'Jacobian',Jname);


% tic
% % for i=1:Ntrials,
%  SolLL3Kp = LLDP_Kphi1(fname,Jname,IT,x0,optionsLL);
% % end;
% tocLL3Kp=toc;
% 
% TLL3Kp = SolLL3Kp.x;
% YLL3Kp = real(SolLL3Kp.y)';
% [T,Y] = ode15s(fname,TLL3Kp,x0,options15s);
%   Y = real(Y);
% tocLL3Kp  
% LL3KpRE = RelError(Y,YLL3Kp)
% SolLL3Kp.stats

% tic
% % for i=1:Ntrials,
%  SolLL3Kp = LLDP_Kphi1_freeJ(fname,IT,x0,optionsLL2);
% % end;
% tocLL3Kp=toc;
% 
% TLL3Kp = SolLL3Kp.x;
% YLL3Kp = real(SolLL3Kp.y)';
% [T,Y] = ode15s(fname,TLL3Kp,x0,options15s);
% %[T,Y] = DLLRK45_3_Auto_Kphi1(fname,Jname,TLL3Kp,x0,optionsLLe);
%   Y = real(Y);
% tocLL3Kp  
% LL3KpRE = RelError(Y,YLL3Kp)
% SolLL3Kp.stats


tic
% for i=1:Ntrials,
 SolLL3Kp = LLDP_Kphi1_freeJ_old(fname,IT,x0,optionsLL2);
% end;
tocLL3Kp=toc;

TLL3Kp = SolLL3Kp.x;
YLL3Kp = real(SolLL3Kp.y)';
[T,Y] = ode15s(fname,TLL3Kp,x0,options15s);
%[T,Y] = DLLRK45_3_Auto_Kphi1(fname,Jname,TLL3Kp,x0,optionsLLe);
  Y = real(Y);
tocLL3Kp  
LL3KpRE = RelError(Y,YLL3Kp)
SolLL3Kp.stats

% tic
% Solode15 = ode15sk(fname,IT,x0,options15s2);
% tocode15=toc;
% % 
% T = Solode15.x;
% Yode15 = real(Solode15.y)';
% [T,Y] = ode15s(fname,T,x0,options15s);
%      Y = real(Y);
% tocode15  
% ode15RE = RelError(Y,Yode15)
% Solode15.stats

% tic
% Solode15j2 = ode15sk_fj(fname,IT,x0,options15s2);
% tocode15j2=toc;
% % % 
% T = Solode15j2.x;
% Yode15 = real(Solode15j2.y)';
% [T,Y] = ode15s(fname,T,x0,options15s);
%      Y = real(Y);
% tocode15j2  
% ode15jRE = RelError(Y,Yode15)
% Solode15j2.stats


% tic
% Solode15j2 = ode15sk_fj_gmres(fname,IT,x0,options15s2);
% tocode15j2=toc;
% % % 
% T = Solode15j2.x;
% Yode15 = real(Solode15j2.y)';
% [T,Y] = ode15s(fname,T,x0,options15s);
%      Y = real(Y);
% tocode15j2  
% ode15jRE = RelError(Y,Yode15)
% Solode15j2.stats
