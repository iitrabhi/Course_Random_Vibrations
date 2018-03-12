clear;
close all;
clc;

tic
%% Beam Stiffness and Mass

syms x xa xb y1 y2 th1 th2 c1 c2 c3 c4 le E I md;
xb = le + xa;

A = [1 xa (xa^2) (xa^3); 0 1 (2*xa) (3*(xa^2)); 1 xb (xb^2) (xb^3); 0 1 (2*xb) (3*(xb^2))];
B = inv(A);
ph1(x) = B(1,1) + B(2,1)*x + B(3,1)*(x^2) + B(4,1)*(x^3);
ph2(x) = B(1,2) + B(2,2)*x + B(3,2)*(x^2) + B(4,2)*(x^3);
ph3(x) = B(1,3) + B(2,3)*x + B(3,3)*(x^2) + B(4,3)*(x^3);
ph4(x) = B(1,4) + B(2,4)*x + B(3,4)*(x^2) + B(4,4)*(x^3);

ph1(x) = simplify(ph1(x));
ph2(x) = simplify(ph2(x));
ph3(x) = simplify(ph3(x));
ph4(x) = simplify(ph4(x));

y(x) = ph1(x)*y1 + ph2(x)*th1 + ph3(x)*y2 + ph4(x)*th2;
yd(x) = diff(y(x));
ydd(x) = diff(yd(x));

w1(x) = ph1(x);
wd1(x) = diff(w1(x));
wdd1(x) = diff(wd1(x));

k1(x) = int((E*I)*ydd(x)*wdd1(x), x, xa, xb);
c1 = coeffs(k1(x), y1);
Ke(1,1) = c1(1,2);
c2 = coeffs(k1(x), th1);
Ke(1,2) = c2(1,2);
c3 = coeffs(k1(x), y2);
Ke(1,3) = c3(1,2);
c4 = coeffs(k1(x), th2);
Ke(1,4) = c4(1,2);

w2(x) = ph2(x);
wd2(x) = diff(w2(x));
wdd2(x) = diff(wd2(x));

k2(x) = int((E*I)*ydd(x)*wdd2(x), x, xa, xb);
c1 = coeffs(k2(x), y1);
Ke(2,1) = c1(1,2);
c2 = coeffs(k2(x), th1);
Ke(2,2) = c2(1,2);
c3 = coeffs(k2(x), y2);
Ke(2,3) = c3(1,2);
c4 = coeffs(k2(x), th2);
Ke(2,4) = c4(1,2);

w3(x) = ph3(x);
wd3(x) = diff(w3(x));
wdd3(x) = diff(wd3(x));

k3(x) = int((E*I)*ydd(x)*wdd3(x), x, xa, xb);
c1 = coeffs(k3(x), y1);
Ke(3,1) = c1(1,2);
c2 = coeffs(k3(x), th1);
Ke(3,2) = c2(1,2);
c3 = coeffs(k3(x), y2);
Ke(3,3) = c3(1,2);
c4 = coeffs(k3(x), th2);
Ke(3,4) = c4(1,2);

w4(x) = ph4(x);
wd4(x) = diff(w4(x));
wdd4(x) = diff(wd4(x));

k4(x) = int((E*I)*ydd(x)*wdd4(x), x, xa, xb);
c1 = coeffs(k4(x), y1);
Ke(4,1) = c1(1,2);
c2 = coeffs(k4(x), th1);
Ke(4,2) = c2(1,2);
c3 = coeffs(k4(x), y2);
Ke(4,3) = c3(1,2);
c4 = coeffs(k4(x), th2);
Ke(4,4) = c4(1,2);

Me(1,1) = int(md*w1(x)*w1(x), x, xa, xb);
Me(1,2) = int(md*w1(x)*w2(x), x, xa, xb);
Me(1,3) = int(md*w1(x)*w3(x), x, xa, xb);
Me(1,4) = int(md*w1(x)*w4(x), x, xa, xb);
Me(2,1) = int(md*w2(x)*w1(x), x, xa, xb);
Me(2,2) = int(md*w2(x)*w2(x), x, xa, xb);
Me(2,3) = int(md*w2(x)*w3(x), x, xa, xb);
Me(2,4) = int(md*w2(x)*w4(x), x, xa, xb);
Me(3,1) = int(md*w3(x)*w1(x), x, xa, xb);
Me(3,2) = int(md*w3(x)*w2(x), x, xa, xb);
Me(3,3) = int(md*w3(x)*w3(x), x, xa, xb);
Me(3,4) = int(md*w3(x)*w4(x), x, xa, xb);
Me(4,1) = int(md*w4(x)*w1(x), x, xa, xb);
Me(4,2) = int(md*w4(x)*w2(x), x, xa, xb);
Me(4,3) = int(md*w4(x)*w3(x), x, xa, xb);
Me(4,4) = int(md*w4(x)*w4(x), x, xa, xb);

n = 2;                                              % Number of DoF for 1 point

%% Solution - Finite Element Part

% Properties of Section

E1 = 69E+9;                                         % Modulus of Rigidity - Al
ro = 2700;

bf = 0.05;
tf = 0.01;
bw = 0.01;
dw = 0.08;

Ar = (bf*tf*2) + (bw*dw);
cg = 0.5*((2*tf)+dw);
I1 = (2*((bf*(tf.^(3))*(1/12))+(bf*tf*((cg-(0.5*tf)).^(2))))) + (bw*(dw.^(3))*(1/12));  % MI

syms F1 F2;                                         
L = 1;                                              % Span of Beam
md1 = Ar*ro;                                         % Distributed Mass

% Elements

N = 3;                                              % Number of Elements
RN = [1 ((N+1)*2) - 1];                             % Restrained DoF
n_RN = length(RN);

for i = 1:N
    EL(i,:) = [i,i+1];                              % Elements defined by nodes
end

for i = 1:N
    for k = 1:n;
        Dof(i,2*k) = 2*(EL(i,k));
        Dof(i,2*k - 1) = 2*(EL(i,k)) - 1;
    end
end

le1 = L/N;

% Global Stiffness and Mass Matrix

Ke1 = double(subs(Ke, [E,I,le], [E1,I1,le1]));

for i = 1:N
    for j = 1:2*n
        for k = 1:2*n
            KeG(Dof(i,j),Dof(i,k)) = Ke1(j,k);
        end
    end
end

Me1 = double(subs(Me, [md,le], [md1,le1]));

for i = 1:N
    for j = 1:2*n
        for k = 1:2*n
            MeG(Dof(i,j),Dof(i,k)) = Me1(j,k);
        end
    end
end

% Natural Frequencies

A1 = inv(MeG)*KeG;
[Phi,w_s] = eig(A1);
wn = sqrt(diag(w_s));

% Force Influence Vector

i1 = 1 + (L/3)/le1;
i2 = 1 + ((2*L)/3)/le1;

Ilf(i1*2 - 1,1) = 1;
Ilf(i2*2 - 1,1) = 0.75;
Ilf(Dof(N,2*n),1) = 0;

w   = 0:0.01:40;
nw  = length(w);

Sgg = 1;

for j = 1:nw,
    Sxx = zeros(2,2);
    Hw  = inv(KeG-w(j)^2.*MeG); 
    Sxx = Hw*Ilf*Sgg*Ilf'*(Hw');
    S11(j) = Sxx(1,1);
    S22(j) = Sxx(2,2);
    S33(j) = Sxx(3,3);
    S44(j) = Sxx(4,4);
    S55(j) = Sxx(5,5);
    S66(j) = Sxx(6,6);
    S77(j) = Sxx(7,7);
    S88(j) = Sxx(8,8);
end

Sig_11 = sqrt(trapz(w,abs(S11))).*1000;
Sig_22 = sqrt(trapz(w,abs(S22))).*1000;
Sig_33 = sqrt(trapz(w,abs(S33))).*1000;
Sig_44 = sqrt(trapz(w,abs(S44))).*1000;
Sig_55 = sqrt(trapz(w,abs(S55))).*1000;
Sig_66 = sqrt(trapz(w,abs(S66))).*1000;
Sig_77 = sqrt(trapz(w,abs(S77))).*1000;
Sig_88 = sqrt(trapz(w,abs(S88))).*1000;

Sig = [Sig_11; Sig_33; Sig_55; Sig_77];

Sigma = max(Sig);

d_all = L/250;

z = norminv(1 - (0.5*(1E-4)));
d_obt = z*Sigma;

I0_sqrt = d_all/d_obt;
I0 = I0_sqrt.^2

%% Bending And Shear Stress

% Cant Do It