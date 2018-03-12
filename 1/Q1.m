clear all;
close all;
clc;

tic
%% Data

L     = 10;
B     = 5;
Th    = 0.1;
Ht    = 3.5;
fck   = 25;
bc    = 0.35;
dc    = 0.35;
zhi   = 0.02;                                           % For all floors
S0    = 0.5E-3;
Ilf   = [-1;-1;-1];
wi    = 0;
wf    = 40;
dw    = 0.01;
d_fac = 1/250;
wg    = 20;
zhi_g = 0.2;

% System Properties

E   = 5000*sqrt(fck)*1E6;
I   = bc*dc^3/12;
k  = 4*(12*E*I/Ht^3);
m  = L*B*Th*25*1000/10;

M = [m 0 0; 0 m 0; 0 0 m];
K = [2*k -k 0; -k 2*k -k; 0 -k k];

[Phi,w_s] = eig(K,M);
wn = sqrt(diag(w_s));
a2 = (2*zhi*(wn(1) - wn(3)))/((wn(1).^2) - (wn(3).^2));
a1 = (2*zhi*wn(1)) - (a2*(wn(1).^2));

C = a1*M + a2*K;

%% Input Functiona

w   = wi:dw:wf;
nw  = length(w);

f1  = wg^4+(4*zhi_g^2*wg^2).*w.^2;
f2  = (w.^2-wg^2).^2+(4*zhi_g^2*wg^2).*w.^2;
Sgg = S0.*f1./f2;

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
    Sxx = Hw*M*Ilf*Sgg(j)*Ilf'*M'*(Hw');
    S11(j) = Sxx(1,1);
    S22(j) = Sxx(2,2);
    S33(j) = Sxx(3,3);
    S12(j) = Sxx(1,2);
    S13(j) = Sxx(1,3);
    S23(j) = Sxx(2,3);
end

% figure
% semilogy(w,abs(S11),w,abs(S22),'g',w,abs(S33),'r')
% legend('S_{11}','S_{22}','S_{33}')
% xlabel('\omega (rad/s)');ylabel('S_{XX}(\omega)')
% title('Response of the 3-DOF System (Exact Solution)')

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

Sigma_1 = max(Sigma1, Sigma2);

pf1 = 2*normcdf(-d_all,0,Sigma_1)

%% Base Shear Check

% PSDF of Forces at each level

Sf11 = zeros(1,nw);
Sf22 = zeros(1,nw);
Sf33 = zeros(1,nw);
Sf12 = zeros(1,nw);
Sf23 = zeros(1,nw);
Sf13 = zeros(1,nw);

for j = 1:nw,
    Sff = zeros(3,3);
    Sff = M*Ilf*Sgg(j)*Ilf'*M';
    Sf11(j) = Sxx(1,1);
    Sf22(j) = Sxx(2,2);
    Sf33(j) = Sxx(3,3);
    Sf12(j) = Sxx(1,2);
    Sf13(j) = Sxx(1,3);
    Sf23(j) = Sxx(2,3);
end

Sig_f11 = sqrt(trapz(w,abs(Sf11))).*1E+6;
Sig_f22 = sqrt(trapz(w,abs(Sf22))).*1E+6;
Sig_f33 = sqrt(trapz(w,abs(Sf33))).*1E+6;
Sig_f12 = sqrt(trapz(w,abs(Sf12))).*1E+6;
Sig_f23 = sqrt(trapz(w,abs(Sf23))).*1E+6;
Sig_f13 = sqrt(trapz(w,abs(Sf13))).*1E+6;

Sigma_2 = sqrt(((Sig_f11).^2) + ((Sig_f22).^2) + ((Sig_f33).^2));

F_all = 4*0.29*200*200;

pf2 = 2*normcdf(-F_all,0,Sigma_2)