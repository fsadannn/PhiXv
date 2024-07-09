function dfdt = J_burgers( t, y, theta )
%J_BURGERS Summary of this function goes here
%   Detailed explanation goes here
N = length(y);
nuc = 3e-4;

%persistent K;
persistent ee;

if nargin == 2
    theta = 0;
end

n1 = N+1;
n12 = 0.5*n1;
n2 = nuc*n1*n1;

if length(ee)~= N
    ee = (-2.0*n2).*ones(N,1);
end

[t1,t2] = size(y);
y = n12.*y;
if t1<t2
    y=y(:);
end
%J = spdiags([y sparse(N,1) -y], -1:1, N, N);
%dfdt = J + K;
dfdt = spdiags([y+n2 ee n2-y], -1:1, N, N);

end

