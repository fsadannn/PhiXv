close all
clear all

% N = 25000;

fixdim = 36;

fname=@f_burgers;
Jname=@J_burgers;
n=500;
N=n;
x0=((sin((3.0*pi/(N+1)).*(1:N))).^2).*((1.0-1/(N+1).*(1:N)).^(3/2));
IT=[0 0.5];
tt=fname(0,x0(:));
tt=Jname(0,x0(:));
clear tt;
% tic
% ttt=J_bruselator(0,x0);
% cond(ttt)
% toc
% tic
% ttt=J_bruselator(0,x0);
% toc
% clear ttt;

Ntrials = 1;
optionsS  = odeset('RelTol',1.0e-3,'AbsTol',1.0e-6);  
optionsLL  = llset('RelTol',1e-9,'AbsTol',1.0e-12,'dKmax',30,...
    'dKmin',4,'debug',1);
optionsLL2  = llset('RelTol',1.0e-9,'AbsTol',1.0e-12,'dKmax',30,...
    'dKmin',5,'debug',0);
optionsLLe  = odeset('RelTol',1.0e-12,'AbsTol',1.0e-14);
options15s2 = odeset('RelTol',1.0e-9,'AbsTol',1.0e-12,'Jacobian',Jname);
options15s3 = odeset('RelTol',1.0e-9,'AbsTol',1.0e-12);
optionsRK  = odeset('RelTol',1.0e-9,'AbsTol',1.0e-12);
options15s = odeset('RelTol',1.0e-12,'AbsTol',1.0e-14,'Jacobian',Jname);
exp4Options  = odeset('RelTol',1.0e-3,'AbsTol',1.0e-6,'Stats','on');  

% tic
% % for i=1:Ntrials,
%  SolLL3Kp = LLDP_Kphi1_freeJ(fname,IT,x0,optionsLL);
% % end;
% tocLL3Kp=toc;
% 
% TLL3Kp = SolLL3Kp.x;
% YLL3Kp = real(SolLL3Kp.y)';
% [T,Y] = ode15s(fname,TLL3Kp,x0,options15s);
% % [T,Y] = exactsol(fname,Jname,TLL3Kp,x0);
%   Y = real(Y);
% tocLL3Kp  
% LL3KpRE = RelError(Y,YLL3Kp)
% SolLL3Kp.stats

% tic
% % for i=1:Ntrials,
%  SolLL3Kp = DLLRK45_3_Auto_Kphi1_allJ(fname,Jname,IT,x0,optionsLL);
%  %SolLL3Kp = DLLRK45_3_Auto_Kphi1_freeJ(fname,Jname,IT,x0,optionsLL);
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

% tic
% %for i=1:Ntrials,
%  SolLL3Kpj = LLDP_Kphi1_freeJ(fname,IT,x0,optionsLL2);
% %end;
% tocLL3Kpj=toc;
% 
% TLL3Kpj = SolLL3Kpj.x;
% YLL3Kpj = real(SolLL3Kpj.y)';
% [T,Y] = ode15s(fname,TLL3Kpj,x0,options15s);
% % [T,Y] = DLLRK45_3_Auto_Kphi1_exact(fname,Jname,TLL3Kpj,x0,40,optionsLLe);
%   Y = real(Y);
% tocLL3Kpj  
% LL3KpjRE = RelError(Y,YLL3Kpj)
% SolLL3Kpj.stats

% tic
% Solode15 = ode15sk_fj(fname,IT,x0,options15s2);
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
% Solode15j = ode15sk_fj(fname,IT,x0,options15s2);
% tocode15j=toc;
% % % 
% T = Solode15j.x;
% Yode15j = real(Solode15j.y)';
% [T,Y] = exactsol(fname,Jname,T,x0);
%      Y = real(Y);
% tocode15j 
% ode15REj = RelError(Y,Yode15j)
% Solode15j.stats


tic
%for i=1:Ntrials,
[TLL3Kpj,YLL3Kpj,stats]  = exp4_fix_fj(fname,IT,x0,exp4Options);
%end;
tocLL3Kpj=toc;

% TLL3Kpj = SolLL3Kpj.x;
% YLL3Kpj = real(SolLL3Kpj.y)';
[T,Y] = ode15s(fname,TLL3Kpj,x0,options15s);
% [T,Y] = DLLRK45_3_Auto_Kphi1_exact(fname,Jname,TLL3Kpj,x0,40,optionsLLe);
  Y = real(Y);
tocLL3Kpj  
LL3KpjRE = RelError(Y,YLL3Kpj)
stats
