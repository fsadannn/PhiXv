function [phi, err, dim, nexpo] = k_phi_fj_adapt_adapt(func,A,h,t,y,m,rtol,atol,kdmax,kdmin)
%PHI1LLDP_hJ_f Summary of this function goes here
%   Detailed explanation goes here
kdmax = 100;

fy = func(t,y);
n=sum(size(A));
m=min(m,n);
gap = length(fy);
gap = n-gap;
b = [A(:,end);zeros(gap,1)];
b(end-1)=1;

% the coeficient is in rfacs( m )
% but is expensive call a function
% and coeficient are hardcode
rfacmin=1;
fac=1/log(2);
gamma = 0.001;
% gamma = 0.005;

nexpo=0;
%
% tolr=[1.0e-6 1.0e-9 1.0e-12 1.0e-15];
% Table=[2 3 4 4;3 4 5 6];

btol = 2*eps;
mb = m;
k1 = 3;
ww = atol+rtol.*abs(b);

normy = norm(y);

% hsubspace = h;
delta = sqrt((1+normy)*eps);
idelta = 1/delta;

normb = norm(b);
V = zeros(n,m+1);
H = zeros(m+3,m+3);

rndoff= eps;

tau = h;
tt = 0;

errs = zeros(3,1);
molds = zeros(3,1);
steps = 0;

while tt < h
    
    % build base of A and b
    
    % begin Arnoldi
    V(:,1)=(1/normb).*b;
    % test if norm(b) is 0
    if any(isnan(V(:,1)))
        phi = sparse(n,1);
        err = 0;
        warning('llint:phi1LLDP_hJ_f_na',...
            'First vector of krylov subspace seems to be 0.')
        return;
    end
    
    for j = 1:m
        %     p = A*V(:,j);
        p = [idelta.*(func(t,y+delta.*V(1:end-gap,j))-fy)+A*V(end-gap+1:end,j);V(end-gap+2:end,j);0];
        %     p = [idelta.*(func(t,y+delta.*V(1:end-gap,j))-func(t,y-delta.*V(1:end-gap,j)))+A*V(end-gap+1:end,j);V(end-gap+2:end,j);0];
        s = V(:,1:j);
        H(1:j,j) = s.'*p;
        p = p - s*H(1:j,j);
        s = norm(p);
        if s < btol
            k1 = 0;
            mb = j;
            break
        end
        H(j+1,j) = s;
        V(:,j+1) = (1/s).*p;
    end
    % end Arnoldi
    % using scaling invariance property of Arnoldi
    % rescaling H to have H for the subspace of A and b
    % H=(1/hsubspace).*H;
    % build \hat{H}
    if k1 == 0
        if mb>1
            mb=mb-1;
            warning('llint:phi1LLDP_hJ_f_na',...
                'Breakdown at dimension 1.')
        end
        m=mb;
        hk = H(m+1,m);
        H=[H(1:m,1:m),zeros(m,3);zeros(3,m+3)];
        % V=V(:,1:m+1);
    else
        hk = H(m+1,m);
        H(m+1,m) = 0;
    end
    H(1,m+1) = 1;
    H(m+1,m+2) = 1; H(m+2,m+3) = 1;
    avm1dot = [idelta.*(func(t,y+delta.*V(1:end-gap,m+1))-fy)+A*V(end-gap+1:end,m+1);V(end-gap+2:end,m+1);0];
    
    
    work=1;
    
    while work
        nhC = norm(H,'inf');
        
        if tau*nhc > 600
           tau = nhC/600;
        end
        
        M14 = expm((tau/4).*H);
        M12 = M14*M14;
        M1 = M12*M12;
        nexpo = nexpo + 1;
        
        % calating \hat{E}
        beta = normb;
        %error relative
        % the divsion by h is because Av_{m+1} is in reality
        % hAv{m+1}(avm1dot) and need rescaling
        err_rel=sqrt((1/n)*sum((((hk*M1(m,m+3)*beta/h).*avm1dot)./ww).^2));
        
        if err_rel/gamma>=1 && m<kdmax && k1~=0
            rfacmax=max(1,m/3);
            knew =  log(err_rel/gamma)*fac;
            knew = ceil(m + min(rfacmax,max(knew,rfacmin)));
            
            if m >= 100
                err_rel14=sqrt((1/n)*sum((((hk*M14(m,m+3)*beta/h).*avm1dot)./ww).^2));
                err_rel12=sqrt((1/n)*sum((((hk*M12(m,m+3)*beta/h).*avm1dot)./ww).^2));
            end
            
            
        else
            tt = tt + tau;
           break; 
        end
        
    end
    
end


err = max(err_rel,rndoff);
dim = m;
end

