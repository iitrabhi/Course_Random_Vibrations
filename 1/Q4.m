clear all;
close all;
clc;

tic

%% Data

fck   = 25;
zhi   = 0.02;                                           % For all floors
wi    = 0;
wf    = 40;
dw    = 0.01;
d_fac = 1/250;
S0    = 1;
Ht    = 3.5;

%% Importing M and K from SAP2000 Model

load('Q4_M');
load('Q4_K');

for i = 1:72
    for j = 1:i
        K(j,i) = K(i,j);
    end
end

for i = 1:72
    for j = 1:i
        M(j,i) = M(i,j);
    end
end

%% Natural Frequencies and Mode Shapes

[Phi,w_s] = eig(K,M);
wn = sqrt(diag(w_s));

% Highest Significant Mode - 12th

a2 = (2*zhi*(wn(1) - wn(12)))/((wn(1).^2) - (wn(12).^2));
a1 = (2*zhi*wn(1)) - (a2*(wn(1).^2));

C = a1*M + a2*K;

%% Input Functiona

w   = wi:dw:wf;
nw  = length(w);

Sgg = S0;

%% Influence Vector

Ilf1(2*3 - 1) = 0.5;
Ilf1(6*3 - 1) = 0.5;
Ilf1(17*3 - 1) = 0.5;
Ilf1(19*3 - 1) = 0.5;
Ilf1(21*3 - 1) = 0.5;
Ilf1(23*3 - 1) = 0.5;
Ilf1(72) = 0;

Ilf = Ilf1.';

%% Frequency Domain Solution

S11 = zeros(1,nw);
S22 = zeros(1,nw);
S12 = zeros(1,nw);
S33 = zeros(1,nw);
S23 = zeros(1,nw);
S13 = zeros(1,nw);

for j = 1:nw,
    Sxx = zeros(3,3);
    Hw  = inv(K-w(j)^2.*M+(i*w(j)).*C); 
    Sxx = Hw*M*Ilf*Sgg*Ilf'*M'*(Hw');
    S11(j) = Sxx(1,1);
    S22(j) = Sxx(2,2);
    S33(j) = Sxx(3,3);
    S12(j) = Sxx(1,2);
    S13(j) = Sxx(1,3);
    S23(j) = Sxx(2,3);
end

figure
semilogy(w,abs(S11),w,abs(S22),'g',w,abs(S33),'r')
legend('S_{11}','S_{22}','S_{33}')
xlabel('\omega (rad/s)');ylabel('S_{XX}(\omega)')
title('Response of the 3-DOF System (Exact Solution):')

Sig_11 = sqrt(trapz(w,abs(S11))).*1000;
Sig_22 = sqrt(trapz(w,abs(S22))).*1000;
Sig_33 = sqrt(trapz(w,abs(S33))).*1000;
Sig_12 = sqrt(trapz(w,abs(S12))).*1000;
Sig_23 = sqrt(trapz(w,abs(S23))).*1000;
Sig_13 = sqrt(trapz(w,abs(S13))).*1000;

d_all  = Ht*d_fac*1000;

%% Inter-Storey Drift

% Xi1 = X1 - X2;
% Xi2 = X2 - X3; 

Sigma1 = sqrt(((Sig_11).^2) + ((Sig_22).^2));                               % Ignoring Covariance
Sigma2 = sqrt(((Sig_22).^2) + ((Sig_33).^2));                               % Ignoring Covariance

Sigma = max(Sigma1, Sigma2);

z = norminv(1 - (0.5*1E-4));
d_obt = z*Sigma;

I0_sqrt = d_all/d_obt;
I0 = I0_sqrt.^2