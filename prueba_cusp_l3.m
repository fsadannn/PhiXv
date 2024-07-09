clear all
close all

dim=4;
N = 500;
N2 = 3*N;
nit = 4;
kdmax = ceil(N2/10);
kdmin = 1;

fname = @f_cusp;
Jname = @J_cusp;
x0= [zeros(1,N),-2*cos((2*pi/(N))*(1:N)),2*sin((2*pi/(N))*(1:N))];
x0 = x0(:);
IT=[0 1.0e-4];

% tic
A = Jname(0,x0);
% toc
% tic
% v=ones(N2,1);
v=fname(0,x0);
% toc
V = [2*v,v,v];
r=zeros(N2+3,1);
r(end)=1;
M=[A,V;zeros(3,N2+3)];
M(N2+1,N2+2)=1;
M(N2+2,N2+3)=1;
V2 = [zeros(N2,1),V];
atol = 1e-9;
rtol = 5e-5;
% h=2e-7; % h3
% h=2e-6; % h2
% h=2e-5; % h1
% clear v;
L = [eye(N2,N2),sparse(N2,3)];
hh =  [2e-5,2e-6,2e-7];
rtols = [1e-1,1e-2,1e-3,1e-4,1e-5,1e-6];
atols = [1e-2,1e-3,1e-4,1e-5,1e-6,1e-7];
pbase = "figs/l3/cusp";

for z=1:length(hh)
    h=hh(z);
    % [kaphi,~,kdim,~] = k_phi_adapt(M,M*r,h,r,dim,rtol,atol,kdmax,4);
    % [kaphi2,~,kdim2,~] = k_phi_fj_adapt(fname,V,h,0,x0,dim,rtol,atol,kdmax,4);
    
    % [kaphi,stats] = phipm(h,A,V2,1e-6,false,4 );
    % [kaphi,stats] = phipm(h,M,r,1e-6,false,4 );
    
    exact = L*expm(h.*M)*r;
    % exact = phi1_single(A,v,h,4);
    % exact = expm(A)*v;
    
    % norm(exact-L*(r+kaphi))
    % norm(exact-L*(r+kaphi2))
    
    nhC = h*norm(M,'inf')
    % % scaling calculation
    % [~,e] = log2(nhC);
    % s = max(0,e+1);
    % pd=6;
    % MM = full(h*M);
    %
    % tic
    % M1 = expm64v4(MM,pd,s);
    % toc
    % tic
    % M2 = expm(MM);
    % toc
    %
    % AbsErr(M1,M2);
    
    
    % phimt = zeros(length(rtols),1);
    phifjt = zeros(length(rtols),1);
    phifjt2 = zeros(length(rtols),1);
    % expvt = zeros(length(rtols),1);
    phit = zeros(length(rtols),1);
    
    % phimerr = zeros(length(rtols),1);
    phierr = zeros(length(rtols),1);
    % expverr = zeros(length(rtols),1);
    phifjerr = zeros(length(rtols),1);
    phifjerr2 = zeros(length(rtols),1);
    
    % phimm = zeros(length(rtols),1);
    phifjm = zeros(length(rtols),1);
    phifjm2 = zeros(length(rtols),1);
    % expvt = zeros(length(rtols),1);
    phim = zeros(length(rtols),1);
    
    
    for j=1:length(rtols)
        atol = atols(j);
        rtol = rtols(j);
        init_dim = 1;
        
        exptimes = zeros(nit,1);
        for i=1:nit
            tic
            A = Jname(0,x0);
            M=[A,V;zeros(3,N2+3)];
            M(N2+1,N2+2)=1;
            M(N2+2,N2+3)=1;
            %         M=[A,V;zeros(3,N2+3)];
            [~,~,~,~] = k_phi_adapt(M,M*r,h,r,init_dim,rtol,atol,kdmax,kdmin);
            tt = toc;
            exptimes(i)=tt;
        end
        exptimes = exptimes(2:end);
        exptimes = exptimes(exptimes<mean(exptimes)+2*std(exptimes));
        phit(j) = mean(exptimes);
        [kexp,~,phi_kdim,~] = k_phi_adapt(M,M*r,h,r,init_dim,rtol,atol,kdmax,kdmin);
        phierr(j) = norm(exact-L*kexp);
        phim(j) = phi_kdim;
        clear kexp;
        
        exptimes = zeros(nit,1);
        for i=1:nit
            tic
            [~,~,~,~] = k_phi_fj_adapt(fname,V,h,0,x0,init_dim,rtol,atol,kdmax,kdmin);
            tt = toc;
            exptimes(i)=tt;
        end
        exptimes = exptimes(2:end);
        exptimes = exptimes(exptimes<mean(exptimes)+2*std(exptimes));
        phifjt(j) = mean(exptimes);
        [kexp,~,phi_kdim_fj,~] = k_phi_fj_adapt(fname,V,h,0,x0,init_dim,rtol,atol,kdmax,kdmin);
        phifjerr(j) = norm(exact-L*kexp);
        phifjm(j)=phi_kdim_fj;
        clear kexp;
        
        
        exptimes = zeros(nit,1);
        for i=1:nit
            tic
            [~,~,~,~] = k_phi_fj_adapt2(fname,V,h,0,x0,init_dim,rtol,atol,kdmax,kdmin);
            tt = toc;
            exptimes(i)=tt;
        end
        exptimes = exptimes(2:end);
        exptimes = exptimes(exptimes<mean(exptimes)+2*std(exptimes));
        phifjt2(j) = mean(exptimes);
        [kexp,~,phi_kdim_fj,~] = k_phi_fj_adapt2(fname,V,h,0,x0,init_dim,rtol,atol,kdmax,kdmin);
        phifjerr2(j) = norm(exact-L*kexp);
        phifjm2(j)=phi_kdim_fj;
        clear kexp;
        
        %     phitimes = zeros(nit,1);
        %     for i=1:nit
        %         tic
        %         A = Jname(0,x0);
        %         M=[A,V;zeros(3,N2+3)];
        % M(N2+1,N2+2)=1;
        % M(N2+2,N2+3)=1;
        % %         M=[A,V;zeros(3,N2+3)];
        % %         phi = phipm(h,A,V2,rtol,false,init_dim );
        %         phi = phipm(h,M,r,rtol,false,init_dim );
        %         tt = toc;
        %         phitimes(i)=tt;
        %     end
        %     phitimes = phitimes(2:end);
        %     phitimes = phitimes(phitimes<mean(phitimes)+2*std(phitimes));
        %     phimt(j) = mean(phitimes);
        % %     phi = phipm(h,A,V2,1e-6,false,init_dim );
        %     [phi,stats] = phipm(h,M,r,rtol,false,init_dim );
        %     kdim_phim = stats(3);
        %     phimerr(j) = norm(exact-L*phi);
        %     clear phi;
        %     clear kphi;
        
        %     expv_dim = ceil(mean([kdim_phim,phi_kdim,phi_kdim_fj]));
        %     expvtimes = zeros(nit,1);
        %     for i=1:nit
        %         tic
        %         A = Jname(0,x0);
        %         M=[A,V;zeros(3,N2+3)];
        % M(N2+1,N2+2)=1;
        % M(N2+2,N2+3)=1;
        % %         M=[A,V;zeros(3,N2+3)];
        %         [kexv,~,~] = expv(h,M,r,rtol,expv_dim);
        %         kexv = L*kexv;
        %         tt = toc;
        %         expvtimes(i)=tt;
        %     end
        %     expvtimes = expvtimes(expvtimes<mean(expvtimes)+2*std(expvtimes));
        %     expvt(j) = mean(expvtimes);
        %     [kexv,~,~] = expv(h,M,r,rtol,expv_dim);
        %     kexv = L*kexv;
        %     expverr(j) = norm(exact-kexv);
        %     clear kexv;
    end
    
    % figure
    % hold on
    % plot(log10(phifjt), log10(phifjerr), 'm--o');
    % % plot(log10(expvt), log10(expverr), 'b--+');
    % plot(log10(phit), log10(phierr), 'r--*');
    % % plot(log10(phimt), log10(phimerr), 'g--s');
    % hold off
    % legend('phi-fj','Phi')
    % % legend('phi-fj','Expv','Phi','phim')
    % ayx = xlabel('$\log_{10}(time)$');
    % ayx.Interpreter = 'latex';
    % ayl = ylabel('$\log_{10}(error)$');
    % ayl.Interpreter = 'latex';
    
    figure
    hold on
    plot(log10(rtols), log10(phifjerr), 'm-o');
    plot(log10(rtols), log10(phifjerr2), 'b--s');
    % plot(log10(rtols), log10(expverr), 'b--+');
    plot(log10(rtols), log10(phierr), 'r-*');
    % plot(log10(rtols), log10(phimerr), 'g--s');
    hold off
    legend('JF1-Phi','JF2-Phi','Phi','Location','northwest')
    % legend('phi-fj','Expv','Phi','phim')
    ax1 = xlabel('$\log_{10}(rtol)$');
    ax1.Interpreter = 'latex';
    ayl = ylabel('$\log_{10}(error)$');
    ayl.Interpreter = 'latex';
    ylim(gca,[-17,0])
    box on
    name = sprintf("%s/h%d-e.fig",pbase,z);
    savefig(name)
    close;
    
    % figure
    % hold on
    % plot(log10(rtols), log10(phifjt), 'm--o');
    % % plot(log10(rtols), log10(expvt), 'b--+');
    % plot(log10(rtols),log10( phit), 'r--*');
    % % plot(log10(rtols), log10(phimt), 'g--s');
    % hold off
    % % legend('phi-fj','Expv','Phi','phim')
    % legend('phi-fj','Phi')
    % ax1 = xlabel('$\log_{10}(tols)$');
    % ax1.Interpreter = 'latex';
    % ayl = ylabel('$\log_{10}(time)$');
    % ayl.Interpreter = 'latex';
    
    
    figure
    hold on
    plot(log10(rtols), phifjm, 'm-o');
    plot(log10(rtols), phifjm2, 'b--s');
    % plot(log10(rtols), log10(expvt), 'b--+');
    plot(log10(rtols),phim, 'r-*');
    % plot(log10(rtols), log10(phimt), 'g--s');
    hold off
    % legend('phi-fj','Expv','Phi','phim')
    legend('JF1-Phi','JF2-Phi','Phi')
    ax1 = xlabel('$\log_{10}(rtol)$');
    ax1.Interpreter = 'latex';
    ayl = ylabel('$m$');
    ayl.Interpreter = 'latex';
    box on
    name = sprintf("%s/h%d-m.fig",pbase,z);
    savefig(name)
    close;
end