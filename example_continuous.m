clc
clear
close all
reset(symengine)

tic

%Uniform
a = 2;
b = 4;
fCDF = @(x) (x-a)/(b-a) * heaviside(x-a)*heaviside(b-x) + heaviside(x-b);
X = 0:6;

CDF = zeros(1,numel(X));
syms t
CF = (exp(1i*t*b) - exp(1i*t*a))/(1i*t*(b-a));

parfor j=1:numel(X)
    CDF(j) = Davies_inversion_continuous( CF, X(j) );
end

figure;
hold on
plot(X,arrayfun(fCDF,X),'DisplayName','Analytical');
plot(X,CDF,'x:','DisplayName','Numerical');
title('Uniform');
legend('show');
hold off

%Gamma
k = 9;
theta = 0.5;
fCDF = @(x) gammainc(x/theta,k,'lower');
X = 0:0.1:10;

CDF = zeros(1,numel(X));
syms t
CF = (1-theta*1i*t)^-k;
parfor j=1:numel(X)
    CDF(j) = Davies_inversion_continuous( CF, X(j) );
end

figure;
hold on
plot(X,arrayfun(fCDF,X),'DisplayName','Analytical');
plot(X,CDF,'x:','DisplayName','Numerical');
legend('show');
title('Gamma');
hold off

elapsedTime = toc;
s = seconds(elapsedTime);
s.Format = 'hh:mm:ss.SSS';
disp('Elapsed time (hms): ');
disp(s);