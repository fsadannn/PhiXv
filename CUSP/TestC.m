close all
clear all

global N
global epsc1
global sigmac

N = 500;
epsc1 = 1/(1.0e-4);
sigmac = 1/144;

fixdim = 30;

fname = @f_cusp;
Jname = @J_cusp;
x0= [zeros(1,N),-2*cos((2*pi/(N))*(1:N)),2*sin((2*pi/(N))*(1:N))];
IT=[0 1.0e-4];

tic
tt=J_cusp(0,x0);
%cond(tt)
toc
tic
J_cusp(0,x0);
toc

Ntrials = 1;
optionsS  = odeset('RelTol',1.0e-3,'AbsTol',1.0e-6);  
optionsLL  = llset('RelTol',1.0e-9,'AbsTol',1.0e-12,'dKmax',30,...
    'dKmin',4,'debug',1,'gamma',0.01);
optionsLL2  = llset('RelTol',1.0e-9,'AbsTol',1.0e-12,'dKmax',30,...
    'dKmin',5,'debug',0,'gamma',0.01);
optionsLLe  = odeset('RelTol',1.0e-12,'AbsTol',1.0e-14);
options15s2 = odeset('RelTol',1.0e-9,'AbsTol',1.0e-12,'Jacobian',Jname);
options15s3 = odeset('RelTol',1.0e-9,'AbsTol',1.0e-12);
optionsRK  = odeset('RelTol',1.0e-9,'AbsTol',1.0e-12);
options15s = odeset('RelTol',1.0e-12,'AbsTol',1.0e-14,'Jacobian',Jname);
exp_4 = exp4set('Complex','off',...
    'AbsTol',1.0e-12,'RelTol',1.0e-9,'MatrixFunctions','arnoldi',...
    'CheckDimensions','off',...
    'ClearInternalData','on','OutputFcn',@outf,'KrylovUseExp4StepFun','on');


% tic
% %for i=1:Ntrials,
%  [t, y] = exp4(fname, Jname,IT, x0,exp_4);
% %end;
% tocLL3Kp=toc;
% 
% TLL3Kp = t;
% YLL3Kp = y(:);
% % [T,Y] = ode15s(fname,TLL3Kp,x0,options15s);
% [T,Y] = exactsol(fname,Jname,TLL3Kp,x0,1,30);
% tocLL3Kp  
% LL3KpRE = RelError(Y,YLL3Kp)

% tic
% %for i=1:Ntrials,
% %  SolLL3Kp = LLDP_Kphi1(fname,Jname,IT,x0,optionsLL);
% SolLL3Kp = DLLRK45_3_Auto_Kphi1(fname,Jname,IT,x0,optionsLL);
% %end;
% tocLL3Kp=toc;
% 
% TLL3Kp = SolLL3Kp.x;
% YLL3Kp = real(SolLL3Kp.y)';
% [T,Y] = ode15s(fname,TLL3Kp,x0,options15s);
% % [T,Y] = exactsol(fname,Jname,TLL3Kp,x0,1);
% %[T,Y] = DLLRK45_3_Auto_Kphi1_exact(fname,Jname,TLL3Kp,x0,50,optionsLLe);
%   Y = real(Y);
% tocLL3Kp  
% LL3KpRE = RelError(Y,YLL3Kp)
% SolLL3Kp.stats



% disp('refine');
% tic
% % for i=1:Ntrials,
%  SolLL3Kpj = DLLRK45_3_Auto_Kphi_freeJ(fname,IT,x0,optionsLL);
% % end;
% tocLL3Kpj=toc;
% 
% TLL3Kpj = SolLL3Kpj.x;
% YLL3Kpj = real(SolLL3Kpj.y)';
% [T,Y] = ode15s(fname,TLL3Kpj,x0,options15s);
% %[T,Y] = DLLRK45_3_Auto_Kphi1_exact(fname,Jname,TLL3Kpj,x0,40,optionsLLe);
%      Y = real(Y);
% tocLL3Kpj  
% LL3KpjRE = RelError(Y,YLL3Kpj)
% SolLL3Kpj.stats
% 
% [tt,tt] = size(T');
% t = T(1:tt-1);
% k = SolLL3Kpj.extinf.kdim;
% figure;
% h = scatter(t,k);
% title('phi1LLDPfj Time-Kdim');
% ylabel('Krilov Dimension');
% xlabel('Time');
% print(gcf,'-djpeg','1LLDPfjtkd')


% % disp('mild');
% tic
% % for i=1:Ntrials,
%  SolLL3Kpj = LLDP_Kphi1(fname,Jname,IT,x0,optionsLL);
% % end;
% tocLL3Kpj=toc;
% % 
% TLL3Kpj = SolLL3Kpj.x;
% YLL3Kpj = real(SolLL3Kpj.y)';
% %[T,Y] = ode15s(fname,TLL3Kpj,x0,options15s);
% [T,Y] = exactsol(fname,Jname,TLL3Kpj,x0,1,100);
%      Y = real(Y);
% tocLL3Kpj  
% LL3KpjRE = RelError(Y,YLL3Kpj)
% SolLL3Kpj.stats
% 
% 
% disp('crude');
tic
% for i=1:Ntrials,
 SolLL3Kpj = LLDP_Kphi1_freeJ(fname,IT,x0,optionsLL2);
% end;
tocLL3Kpj=toc;

TLL3Kpj = SolLL3Kpj.x;
YLL3Kpj = real(SolLL3Kpj.y)';
[T,Y] = exactsol(fname,Jname,TLL3Kpj,x0,1,50);
% [T,Y] = ode15s(fname,TLL3Kpj,x0,options15s);
Y = real(Y);
tocLL3Kpj  
LL3KpjRE = RelError(Y,YLL3Kpj)
SolLL3Kpj.stats
% print(gcf,'-djpeg','cp1LLDPfjtkd')
% 
% tic
% Solode15 = ode15sk(fname,IT,x0,options15s2);
% tocode15=toc;
% 
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
