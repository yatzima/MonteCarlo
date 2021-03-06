%%Power production of a wind turbine
clc
clear
close all
load('powercurve_V112.mat');

lambda = [10.6 9.7 9.2 8.0 7.8 8.1 7.8 8.1 9.1 9.9 10.6 10.6];
k = [2.0 2.0 2.0 1.9 1.9 1.9 1.9 1.9 2.0 1.9 2.0 2.0];
lamdaG = [5.8 6.5 6.5 6.5 6.5 6.5 6.5 6.5 6.5 6.5 6.5 6.5];
kappaG = [3 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5];
month = 1;
N = 1000;

%Defining the stochastic wind speed V for different months
f = @(v, month) wblpdf(v, lambda(month), k(month)); 
%Defining the g function
g = @(v, month) gampdf(v, lamdaG(month), kappaG(month));

%Defining the random-generator for different months
Frand = @(month) wblrnd(lambda(month), k(month), 1, N);
Grand = @(month) gamrnd(lamdaG(month), kappaG(month), 1, N);

%Defining the inverse
a = 25;
b = 3;
FV = @(u, month) u*(wblcdf(a, lambda(month), k(month)) - wblcdf(b, lambda(month), k(month))) + wblcdf(b, lambda(month), k(month));
FU = @(U, month) wblinv(U, lambda(month), k(month));

%Crude monte carlo
draw1 = Frand(1);
y1 = P(draw1);
tau1 = mean(P(draw1));
std1 = std(P(draw1));

%Truncated
u = rand(1, N);
we = FV(u, month);
draw2 = FU(we, month);
comp = wblcdf(25, lambda(month), k(month)) - wblcdf(3, lambda(month), k(month));
tau2 = mean(P(draw2))*comp;
std2 = std(P(draw2)*comp) ; %Ska man multiplicera?

%Importance Sampling
draw3 = Grand(month);
y3 = g(draw3, month);
phiomega = P(draw3)*(f(draw3, month)/g(draw3, month));
tau3 = mean(phiomega);
std3 = std(phiomega);

%Antithetic Sampling
u = rand(1, N);
draw41 = FU(u, month);
draw42 = FU(1-u, month);
V1 = P(draw41);
V2 = P(draw42);
W = (V1+V2)/2;
tau4 = mean(W) * comp;
std4 = std(W * comp);

%Calculating the effiency and the mean probability
rho = 1.225;
d = 112;
Ptot = @(v) (1/2) .* rho .* pi .* ((d.^2)/4) .* v.^3;
E1 = zeros(N,12);
E2 = zeros(N,12);
prob = zeros(1, 12);
counter = 1;
for month = 1:12
    X = Frand(month);
    E1(:, counter) = P(X);
    E2(:, counter) = Ptot(X);
    counter = counter + 1;
    prob(1, month) = length(find(X>=3 & X<=25))/N;
end
E1 = mean(P(X));
E2 = mean(Ptot(X));
prob = mean(prob);
effi = E1/E2;

%Plot the quote and the dist for f and g
figure(1)
hold on
draw1 = zeros(N, 12);
y1 = zeros(N, 12);
counter = 1;

for month = 1:12
    draw1(:, counter) = Frand(month)';
    y1(:, counter) = P(draw1(:, counter)).*f(draw1(:, counter), month);
    plot(draw1(:,counter), y1(:,counter), '*')
    counter = counter + 1;
end
plot(draw3, 3e6*y3, 'o')
legend('Function f','Function g')

%Plot the quote and the dist for f and g
% figure(1)
% hold on
% draw = zeros(N, 12);
% y = zeros(N, 12);
% counter = 1;
% for month = 1:12
%     draw(:, counter) = Grand(month)';
%     y(:, counter) = P(draw(:, counter)).*(f(draw(:, counter), month)./g(draw(:, counter), month));
%     plot(draw1(:,counter), y(:,counter), '*')
%     counter = counter + 1;
% end
% title('Quote')

%Defining the random-generator for different months
Frand = @(month, N) wblrnd(lambda(month), k(month), 1, N);
Grand = @(month, N) gamrnd(lamdaG(month), kappaG(month), 1, N);

N = 1:50:2000;
counter = 1;
tau1 = zeros(size(N));
tau2 = zeros(size(N));
tau3 = zeros(size(N));
tau4 = zeros(size(N));
lambda95 = norminv(0.975);

ci1 = zeros(length(N), 2);
ci2 = zeros(length(N), 2);
ci3 = zeros(length(N), 2);
ci4 = zeros(length(N), 2);

month = 1;
comp = wblcdf(25, lambda(month), k(month)) - wblcdf(3, lambda(month), k(month));
for N = 1:50:2000;
    %Crude Monte Carlo
    draw1 = Frand(month, N);
    y1 = P(draw1);
    tau1(counter) = mean(y1);
    std1 = std(y1);
    
    %Truncated
    u = rand(1, N);
    we = FV(u, month);
    draw2 = FU(we, month);
    tau2(counter) = mean(P(draw2)) * comp;
    std2 = std(P(draw2) * comp);
    
    %Importance Sampling
    draw3 = Grand(month, N);
    phiomega = P(draw3)*(f(draw3, month)/g(draw3, month));
    tau3(counter) = mean(phiomega);
    std3 = std(phiomega);
    
    %Antithetic Sampling
    uniform1 = rand(1, N);
    uniform2 = 1 - rand(1, N);
    draw41 = FU(uniform1, month);
    draw42 = FU(uniform2, month);
    V1 = P(draw41);
    V2 = P(draw42);
    W = (V1+V2)/2;
    tau4(counter) = mean(W);
    std4 = std(W);
    
    %CI1 - Crude MC
    ci1(counter, :) = tau1(counter) + [1 -1] * lambda95*(std1/sqrt(N));

    %CI2 - Truncated
    ci2(counter, :) = tau2(counter) + [1 -1] * lambda95*(std2/sqrt(N));

    %CI3 - Importance sampling
    ci3(counter, :) = tau3(counter) + [1 -1] * lambda95*(std3/sqrt(N));

    %CI4 - Antithetic sampling
    ci4(counter, :) = tau4(counter) + [1 -1] * lambda95*(std4/sqrt(N));
    
    counter = counter + 1;
end

N = 1:50:2000;

figure(3)
plot(N,tau1, 'r');
hold on
plot(N,tau2, 'g');
plot(N,tau3, 'b');
plot(N,tau4, 'k');
title('Expected value based on difference MC-methods for March')
legend('Crude MC','Truncated Crude MC', 'Importance Sampling', 'Antithetic Sampling')
xlabel('N')
ylabel('Power output')

figure(4)
hold on
plot(N, ci1(:,1), 'r');
plot(N, ci1(:,2), 'r');
plot(N, ci2(:,1), 'g');
plot(N, ci2(:,2), 'g');
plot(N, ci3(:,1), 'b');
plot(N, ci3(:,2), 'b');
plot(N, ci4(:,1), 'k');
plot(N, ci4(:,2), 'k');
title('Confidence interval based on difference MC-methods for January')
legend('Crude MC 1', 'Crude MC 2', 'Truncated Crude MC 1', 'Truncated Crude MC 2', 'Importance Sampling 1', 'Importance Sampling 2', 'Antithetic Sampling 1', 'Antithetic Sampling 2')
xlabel('N')
ylabel('Something')

%% Combined power production of two wind turbines
clear
close all
load('powercurve_V112.mat');

%Defining sample size.
N = 1e4;

%Defining the given constants.
alpha = 0.638; pe = 3; q = 1.5; lambda = 9.13; k = 1.96;

%Defining shape and scale parameters for the gamma distribution used in IS
%for one dimensional case.
lamdaG = 6.7;
kappaG = 2.5;

%Defining mu and sigma for the normal distribution used as instrumental
%function in IS for bivariate case.
mu = 11;
sigma = 3;

%Defining functions that are used to formulate the joint PDF
F = @(v) wblcdf(v, lambda, k);
f = @(v) wblpdf(v, lambda, k);
g = @(v) gampdf(v, lamdaG, kappaG);

%Fjoint = @(v1, v2) F(v1) .* F(v2) .* (1 + alpha.*(1 - F(v1).^pe).^q .* (1 - F(v2).^pe).^q);
fjoint = @(v1, v2) f(v1) .* f(v2) .* (1 + alpha .* (1 - F(v1).^pe).^(q-1) .* (1 - F(v2).^pe).^(q-1) .* (F(v1).^pe .* (1+pe*q)-1) .* (F(v2).^pe .* (1+pe*q) - 1));
pjoint = @(v1, v2) P(v1) .* P(v2);

%Defining 
Grand = @(N) gamrnd(lamdaG, kappaG, N, 1); 
%gjoint = @(v1, v2) gampdf(v1, lamdaG, kappaG) .* gampdf(v2, lamdaG, kappaG) ;
%gjoint = @(v1, v2) mvnpdf(v1, mu, sigma) .* mvnpdf(v2, mu, sigma);
gjoint = @(v1, v2) normpdf(v1, mu, sigma) .* normpdf(v2, mu, sigma);

x = linspace(0, 25, 100)';
y = linspace(0, 25, 100)';
p1 = pjoint(x, y');
z1 = fjoint(x, y');
%z2 = mvnpdf(x, mu, sigma) * mvnpdf(y, mu, sigma)';
%z2 = gjoint(x, y');
z2 = normpdf(x, mu, sigma) * normpdf(y, mu, sigma)';
z3 = z1./z2;

figure(1)
surf(x, y, z1.*p1);
title('The objective function multiplied with the target function')
xlabel('x')
ylabel('y')
zlabel('z')

figure(2)
surf(x, y, z2);
title('The instrumental function chosen as a normal distribution')
xlabel('x')
ylabel('y')
zlabel('z')

figure(3)
surf(x, y, z3);
title('The phi-omega quota')
xlabel('x')
ylabel('y')
zlabel('z')

%Calculating the expected value of the combined 
draw = Grand(N);
omega = f(draw)./g(draw);
phiomega1 = P(draw).* omega;
tau = mean(phiomega1);
tau2 = 2 * tau;

%Calculating the covariance
%drawV1 = Grand(N);
%drawV2 = Grand(N);
drawV1 = normrnd(mu, sigma, N, 1);
drawV2 = normrnd(mu, sigma, N, 1);
%drawV1 = mvnrnd(mu, sigma, N);
%drawV2 = mvnrnd(mu, sigma, N);

phiobjective = pjoint(drawV1, drawV2);
ftarget = fjoint(drawV1, drawV2);
ginstrumental = gjoint(drawV1, drawV2);

phiomega2 = phiobjective .* (ftarget./ginstrumental);
cov = mean(phiomega2) - tau.^2;

%Calculating the variance and standard deviation
vari = 2*var(phiomega1) + 2*cov;
stde = sqrt(vari);

%Calculating the probability
Pind1 = @(v1, v2) P(v1) + P(v2) > 3.075e6;
%Pind1 = @(v1, v2) heaviside(P(v1) + P(v2) - 3.075e6);
Pind2 = @(v1, v2) P(v1) + P(v2) < 3.075e6;
N = 1e6;
%drawV1 = Grand(N);
%drawV2 = Grand(N);
%drawV1 = mvnrnd(mu, sigma, N);
%drawV2 = mvnrnd(mu, sigma, N);
drawV1 = normrnd(mu, sigma, N, 1);
drawV2 = normrnd(mu, sigma, N, 1);

omega = fjoint(drawV1, drawV2)./gjoint(drawV1,drawV2);
phi1 = Pind1(drawV1, drawV2);
phi2 = Pind2(drawV1, drawV2);
sum1 = phi1 .* omega;
sum2 = phi2 .* omega;
prob1 = mean(sum1);
prob2 = mean(sum2);
lambda95 = norminv(0.975);
std1 = std(sum1);
ci1 = prob1 + [-1 1] * lambda95*(std1/sqrt(N));
std2 = std(sum2);
ci2 = prob2 + [-1 1] * lambda95*(std2/sqrt(N));
probtot = prob1 + prob2;