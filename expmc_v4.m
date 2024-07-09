function E = expmc_v4(A,p,s)
% Matrix exponential via Pade approximation (p,p).
% NA = norm(A,'inf');

% Scale A by power of 2 so that its norm is < 1/2 .
% [f,e] = log2(NA);
% s = max(0,e+1);
A = A/2^s;
c1 = 0.5; 
X = A; 
E = eye(size(A)) + c1*A;
D = eye(size(A)) - c1*A;
switch p
    case 2, 
         c2 = 0.083333333333333;
         X = A*X;    E = E + c2 * X;   D = D + c2 * X;
    case 3,
         c2 = 0.100000000000000; c3 = 0.008333333333333;
         X = A*X;    E = E + c2 * X;   D = D + c2 * X;
         X = A*X;    E = E + c3 * X;   D = D - c3 * X;
    case 4,
         c2 = 0.107142857142857; c3 = 0.011904761904762; c4 = 5.952380952380952e-04;
         X = A*X;    E = E + c2 * X;   D = D + c2 * X;
         X = A*X;    E = E + c3 * X;   D = D - c3 * X;
         X = A*X;    E = E + c4 * X;   D = D + c4 * X;    
    case 5,
         c2 = 0.111111111111111; c3 = 0.013888888888889; c4 = 9.920634920634920e-04; c5 = 3.306878306878306e-05;
         X = A*X;    E = E + c2 * X;   D = D + c2 * X;
         X = A*X;    E = E + c3 * X;   D = D - c3 * X;
         X = A*X;    E = E + c4 * X;   D = D + c4 * X;    
         X = A*X;    E = E + c5 * X;   D = D - c5 * X;          
    case 6,
         c2 = 0.113636363636364; c3 = 0.015151515151515; c4 = 0.001262626262626; c5 = 6.313131313131313e-05; c6 = 1.503126503126503e-06;
         X2=A*A;  X4 = X2*X2; X6 = X4*X2;  X3 = X2*A; X5 = X4*A;
         E = c6.*X6;   D = c6.*X6;
         E = E + c5.*X5;  D = D - c5.*X5;
         E = E + c4.*X4;  D = D + c4.*X4;
         E = E + c3.*X3;  D = D - c3.*X3;
         E = E + c2.*X2;  D = D + c2.*X2;
         E = E + c1.*A;  D = D - c1.*A;
         E = E + eye(size(A));  D = D + eye(size(A));
%          X = A*X;    E = E + c2 * X;   D = D + c2 * X;
%          X = A*X;    E = E + c3 * X;   D = D - c3 * X;
%          X = A*X;    E = E + c4 * X;   D = D + c4 * X;    
%          X = A*X;    E = E + c5 * X;   D = D - c5 * X;          
%          X = A*X;    E = E + c6 * X;   D = D + c6 * X;          
end
% Pade approximation for exp(A)
E = D\E;

% Undo scaling by repeated squaring
for k=1:s, E = E*E; end
