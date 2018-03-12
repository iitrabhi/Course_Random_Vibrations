clear all;
close all;
clc;

tic % Takes 452 secs
%% Evaluation of Autocorrelation 

tmax = 50;
dt = 0.01;
t = 0:dt:tmax;                                                      % Time Interval
tau = -tmax:dt:tmax;
n = length(t);

y1 = cos(5*t);                                                      % First function

for i = 1:n
    y2(i) = (exp((-0.5)*t(i)))*(sin(10*t(i)) + cos(15*t(i)));       % Second function
end

y3 = wgn(n,1,0);                                                    % Third function

y4 = exp((-0.5)*t);

% Auto-Correlation - Matlab Function

Ryy1 = xcorr(y1);

Ryy2 = xcorr(y2);

Ryy3 = xcorr(y3);

Ryy4 = xcorr(y4);

% Auto-Correlation - Own Function

for i = 1:n
    jmax = 1 + ((tmax - abs(tau(i)))/dt);
    dj = single(abs(tau(i))/dt);
    
    for j = 1:jmax
        R1(i,j) = y1(n-(j-1)).*y1((n-(j-1))-dj);
        R2(i,j) = y2(n-(j-1)).*y2((n-(j-1))-dj);
        R3(i,j) = y3(n-(j-1)).*y3((n-(j-1))-dj);
        R4(i,j) = y4(n-(j-1)).*y4((n-(j-1))-dj);
    end
end

for k = 1:n
    r1 = R1(k,:);
    Ry1(k) = sum(r1);
    Ry1(2*n - k) = Ry1(k);
    
    r2 = R2(k,:);
    Ry2(k) = sum(r2);
    Ry2(2*n - k) = Ry2(k);
    
    r3 = R3(k,:);
    Ry3(k) = sum(r3);
    Ry3(2*n - k) = Ry3(k);
    
    r4 = R4(k,:);
    Ry4(k) = sum(r4);
    Ry4(2*n - k) = Ry4(k);
end


% Plots - Matlab Function

figure(1)
subplot(2,2,1)
plot(tau, Ryy1)
xlabel('\tau(s)')
ylabel('Auto-Correlation for y1')

subplot(2,2,2)
plot(tau, Ryy2)
xlabel('\tau(s)')
ylabel('Auto-Correlation for y2')

subplot(2,2,3)
plot(tau, Ryy3)
xlabel('\tau(s)')
ylabel('Auto-Correlation for y3')

subplot(2,2,4)
plot(tau, Ryy4)
xlabel('\tau(s)')
ylabel('Auto-Correlation for y4')


figure(2)
subplot(2,2,1)
plot(tau, Ry1)
xlabel('\tau(s)')
ylabel('Auto-Correlation for y1')

subplot(2,2,2)
plot(tau, Ry2)
xlabel('\tau(s)')
ylabel('Auto-Correlation for y2')

subplot(2,2,3)
plot(tau, Ry3)
xlabel('\tau(s)')
ylabel('Auto-Correlation for y3')

subplot(2,2,4)
plot(tau, Ry4)
xlabel('\tau(s)')
ylabel('Auto-Correlation for y4')