clear all;
close all;
clc;

tic
%% Data

S0 = 1E-4;
wg = 20;
ng = 0.2;
wi = 0;
wf = 150;
dw = 0.01;
Ilf = [-1,-1,-1,-1,-1].';

%% System Properties

load K
load M
load C
dim = length(K);

%% Input Functiona

w   = wi:dw:wf;
nw  = length(w);

f1  = wg^4+(4*ng^2*wg^2).*w.^2;
f2  = (w.^2-wg^2).^2+(4*ng^2*wg^2).*w.^2;
Sgg = S0.*f1./f2;

figure
semilogy(w,Sgg)
xlabel('\omega (rad/s)')
ylabel('S_{gg}(\omega)')
title('Input PSDF')

%% Exact Solution - Frequency Solution

Sx =  zeros(dim,dim,nw);

for k = 1:nw,
    Sxx = zeros(dim,dim);
    Hw  = inv(K-w(k)^2.*M+(i*w(k)).*C); 
    Sxx = Hw*M*Ilf*Sgg(k)*(Ilf')*(M')*(Hw');
    for l = 1:dim
        for m = 1:dim
            Sx(l,m,k) = Sxx(l,m);
        end
    end
end

figure
semilogy(w,abs(reshape(Sx(1,1,:),[1,nw])), w,abs(reshape(Sx(2,2,:),[1,nw])), w,abs(reshape(Sx(3,3,:),[1,nw])), w,abs(reshape(Sx(4,4,:),[1,nw])), w,abs(reshape(Sx(5,5,:),[1,nw])));
legend('S_{11}','S_{22}','S_{33}','S_{44}','S_{55}')
xlabel('\omega (rad/s)')
ylabel('S_{XX}(\omega)')
title('Response PSDF of the 5-DOF System (Exact Solution)')

for l = 1:dim
    for m = 1:dim
        Sxr = reshape(Sx(l,m,:), [1,nw]);
        Sigma_x_e(l,m) = sqrt(trapz(w,abs(Sxr))).*1E+6;
    end
end

% Top Response RMS

y5_rms_exact = Sigma_x_e(5,5)

% pause

%% Modal Solution

[phi, lam] = eig(K,M);
wn = sqrt(diag(lam));

Mn = phi.'*M*phi;
Kn = phi.'*K*phi;
Cn = phi.'*C*phi;
Fn = phi.'*M*Ilf;

Sz = zeros(dim,dim,nw);

for k = 1:nw,
    Szz = zeros(dim,dim);
    Hw  = inv(Kn-w(k)^2.*Mn+(i*w(k)).*Cn); 
    Szz = Hw*M*Ilf*Sgg(k)*(Ilf')*(M')*(Hw');
    for l = 1:dim
        for m = 1:dim
            Sz(l,m,k) = Szz(l,m);
        end
    end
end

figure
semilogy(w,abs(reshape(Sz(1,1,:),[1,nw])), w,abs(reshape(Sz(2,2,:),[1,nw])), w,abs(reshape(Sz(3,3,:),[1,nw])), w,abs(reshape(Sz(4,4,:),[1,nw])), w,abs(reshape(Sz(5,5,:),[1,nw])));
legend('Sz_{11}','Sz_{22}','Sz_{33}','Sz_{44}','Sz_{55}')
xlabel('\omega (rad/s)')
ylabel('S_{ZZ}(\omega)')
title('Modal Response PSDF of the 5-DOF System')

for l = 1:dim
    for m = 1:dim
        Szr = reshape(Sz(l,m,:), [1,nw]);
        Sigma_z(l,m) = sqrt(trapz(w,abs(Szr))).*1E+6;
    end
end

for l = 1:dim
    for m = 1:dim
        Cz(l,m) = Sigma_z(l,m).^2;
    end
end

Var_x_m = phi*Cz*phi.'';

for l = 1:dim
    for m = 1:dim
        Sigma_x_m(l,m) = sqrt(Var_x_m(l,m));
    end
end

y5_rms_no_combo = Sigma_x_m(5,5)                         % No Modal Combinations

%% Global Responses per mode

Syy = zeros(dim,dim,dim);

for i = 1:dim
    Syy(i,:,:) = phi(:,i).*Cz(5,5)*phi(:,1).';
end

Syy1 = reshape(Syy(:,5,5), [1,dim]);

% SRSS

y5_srss = 0;

for i = 1:dim
    y5_srss = y5_srss + (Syy1(i).^2);
end

y5_rms_srss = sqrt(y5_srss)                          % SRSS

% CQC

y5_cqc = 0;

tau = 0.05;

for i = 1:dim
    for j = 1:dim
        b(i,j) = wn(j)/wn(i);
        corr(i,j) = ((8*(tau.^2)*(1+b(i,j))*(b(i,j).^(1.5)))/(((1-(b(i,j).^2)).^2)+(4*(tau.^2)*b(i,j)*((1+b(i,j)).^2))));
        
        y5_cqc = y5_cqc + (Syy1(i)*corr(i,j)*Syy1(j));
    end
end

y5_rms_cqc = sqrt(y5_cqc)