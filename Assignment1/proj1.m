%%Power production of a wind turbine
clc
clear
close all
load('powercurve_V112.mat');

lambda = [10.6 9.7 9.2 8.0 7.8 8.1 7.8 8.1 9.1 9.9 10.6 10.6];
k = [2.0 2.0 2.0 1.9 1.9 1.9 1.9 1.9 2.0 1.9 2.0 2.0];
const1 = [5.8 6.5 6.5 6.5 6.5 6.5 6.5 6.5 6.5 6.5 6.5 6.5];
const2 = [3 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5];

N = 1000;

%Defining the stochastic wind speed V for different months
f = @(v, month) wblpdf(v, lambda(month), k(month)); 
%Defining the g function
g = @(v, month) gampdf(v, const1(month), const2(month));

%Defining the random-generator for different months
Frand = @(month) wblrnd(lambda(month), k(month), 1, N);
Grand = @(month) gamrnd(const1(month), const2(month), 1, N);

%Defining the inverse
a = 25;
b = 3;
FV = @(u, month) u*(wblcdf(a, lambda(month), k(month)) - wblcdf(b, lambda(month), k(month))) + wblcdf(b, lambda(month), k(month));
FU = @(U, month) wblinv(U, lambda(month), k(month));
%% Basics
figure(5) 
lin = linspace(0,30);
plot(lin,P(lin))
xlabel('Wind Speed (m/s)')
ylabel('Power Output (W)')
title('Power curve of a V112 turbine')
figure(6)
hold on
for month = 1:12
    plot(lin,f(lin,month))
end
title('Distribution of wind speeds each month')
xlabel('Wind Speeds (m/s)')
ylabel('Distribution of wind speed')
legend('Jan','Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
%% Plotting for 2a)-2c)
close all
ci1 = zeros(12,2);
ci2 = zeros(12,2);
ci3 = zeros(12,2);
ci4 = zeros(12,2);
width1 = zeros(12,1);
width2 = zeros(12,1);
width3 = zeros(12,1);
width4 = zeros(12,1);
lambda95 = 1.96;

for month = 1:12
%Crude monte carlo
draw1 = wblrnd(lambda(month), k(month), 1, N);
%y1 = P(draw1);
tau1 = mean(P(draw1));
std1 = std(P(draw1));

%Truncated
u = rand(1, N);
we = FV(u, month);
draw2 = FU(we, month);
comp = wblcdf(25, lambda(month), k(month)) - wblcdf(3, lambda(month), k(month));
tau2 = mean(P(draw2))*comp;
std2 = std(P(draw2)*comp);

%Importance Sampling
draw3 = Grand(month);
y3 = g(draw3, month);
phiomega = P(draw3)*(f(draw3, month)/g(draw3, month));
tau3 = mean(phiomega);
std3 = std(phiomega);

%Antithetic Sampling
uniform1 = rand(1, N/2);
uniform2 = 1 - rand(1, N/2);
draw41 = FU(uniform1, month);
draw42 = FU(uniform2, month);
V1 = P(draw41);
V2 = P(draw42);
W = (V1+V2)/2;
tau4 = mean(W);
std4 = std(W);

ci1(month,:) = tau1 + [-1,1]*lambda95*std1/sqrt(N);
ci2(month,:) = tau2 + [-1,1]*lambda95*std2/sqrt(N);
ci3(month,:) = tau3 + [-1,1]*lambda95*std3/sqrt(N);
ci4(month,:) = tau4 + [-1,1]*lambda95*std4/sqrt(N);

width1(month) = 2*lambda95*std1/sqrt(N);
width2(month) = 2*lambda95*std2/sqrt(N);
width3(month) = 2*lambda95*std3/sqrt(N);
width4(month) = 2*lambda95*std4/sqrt(N);
end

avCi1 = [mean(ci1(:,1)), mean(ci1(:,2))];
avCi2 = [mean(ci2(:,1)), mean(ci2(:,2))];
avCi3 = [mean(ci3(:,1)), mean(ci3(:,2))];
avCi4 = [mean(ci4(:,1)), mean(ci4(:,2))];

avWidth1 = mean(width1);
avWidth2 = mean(width2);
avWidth3 = mean(width3);
avWidth4 = mean(width4);

figure(9)
hold on
p1 = plot(linspace(1,12,12), ci1(:,1),'r');
plot(linspace(1,12,12), ci1(:,2),'r')
p2 = plot(linspace(1,12,12), ci2(:,1),'b');
plot(linspace(1,12,12), ci2(:,2),'b')
p3 = plot(linspace(1,12,12), ci3(:,1), 'y');
plot(linspace(1,12,12), ci3(:,2), 'y')
p4 = plot(linspace(1,12,12), ci4(:,1),'g');
plot(linspace(1,12,12), ci4(:,2),'g')
legend([p1,p2,p3,p4],'Crude MC', 'Truncated', 'Importance Sampling', 'Antithetic')
xlim([1,12])

%% 2d)
%Calculating the effiency and the mean probability
rho = 1.225;
d = 112;
Ptot = @(v) (1/2) .* rho .* pi .* ((d.^2)/4) .* v.^3;
E1 = zeros(N,12);
E2 = zeros(N,12);
prob = zeros(1, 12);
for month = 1:12
    X = Frand(month);
    E1(:, month) = P(X);
    E2(:, month) = Ptot(X);
    prob(1, month) = length(find(X>=3 & X<=25))/N;
end
avE1 = mean(P(X));
avE2 = mean(Ptot(X));
avProb = mean(prob);
effi = E1/E2;
%% 2e)
pvec = zeros(1,12);
EP = zeros(3,12);
EPtot = zeros(1,12);
rho = 1.225;
lambda95 = 1.96;
d = 112;
ci = zeros(2,12);
tau = zeros(1,12);
EPtotFunc = @(month) 1/8*rho*pi*d^2*lambda(month)^3*gamma(3/k(month) + 1);
for month = 1:12
    %Importance Sampling
    draw = Grand(month);
    phiomega = P(draw)*(f(draw, month)/g(draw, month));
    EP(2,month) = mean(phiomega);
    diff = lambda95*(std(phiomega)/sqrt(N));
    EP(1,month) = EP(2,month) + diff;
    EP(3,month) = EP(2,month) - diff;
    EPtot(month) = EPtotFunc(month);
end
eff = EP(2,:)./EPtot;
stde = EP(1,:)-EP(3,:);
stdPost = stde./EPtot;
mean(stdPost);
mean(eff);
%% 2f)
N = 1000;
Pmeans = zeros(1,12);
for month = 1:12
    draw = Frand(month);
    Ptemp = P(draw);
    Pmeans(month) = mean(Ptemp);
end
capfac = sum(Pmeans)/(12*3.075*10^6);
avfac = mean(prob);
%% Plot the quota of Pf/g
close all
 figure(2)
 N = 1000;
 c1 = 6.5;
 c2 = 2.3;
 hold on
 draw = zeros(N, 12);
 y = zeros(N, 12);
 for month = 1:12
    draw(:, month) = gamrnd(lambda(month), k(month),1,N);
    y(:, month) = P(draw(:, month)).*(f(draw(:, month), month)./gampdf(draw(:,month), c1,c2));
    plot(draw(:,month), y(:,month), '.')
    xlim([3,25])
 end
title('Quota of P*f/g')
xlabel('Windspeed (m/s)')
ylabel('P*f/g')
legend('Jan','Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
%% Defining the random-generator for different months
Frand = @(month, N) wblrnd(lambda(month), k(month), 1, N);
Grand = @(month, N) gamrnd(const1(month), const2(month), 1, N);

N = 1:50:4000;
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

for N = 1:50:4000
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
    %phiomega = P(draw3)*(f(draw3, month)/g(draw3, month));
    phiomega = P(draw3)*(f(draw3, month)/gampdf(draw3,c1,c2));
    tau3(counter) = mean(phiomega);
    std3 = std(phiomega);
    
    %Antithetic Sampling
    uniform1 = rand(1, round(N/2));
    uniform2 = 1 - rand(1, round(N/2));
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

N = 1:50:4000;

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
p1 = plot(N, ci1(:,1), 'r');
plot(N, ci1(:,2), 'r');
p2 = plot(N, ci2(:,1), 'g');
plot(N, ci2(:,2), 'g');
p3 = plot(N, ci3(:,1), 'b');
plot(N, ci3(:,2), 'b');
p4 = plot(N, ci4(:,1), 'k');
plot(N, ci4(:,2), 'k');
title('Confidence interval based on difference MC-methods for January')
%legend('Crude MC 1', 'Crude MC 2', 'Truncated Crude MC 1', 'Truncated Crude MC 2', 'Importance Sampling 1', 'Importance Sampling 2', 'Antithetic Sampling 1', 'Antithetic Sampling 2')
legend([p1 p2 p3 p4], 'Crude MC', 'Truncated Crude MC', 'Importance Sampling', 'Antithetic');
xlabel('N')
ylabel('Power Output')

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
gjoint = @(v1, v2) normpdf(v1, mu, sigma) .* normpdf(v2, mu, sigma);

x = linspace(0, 25, 100)';
y = linspace(0, 25, 100)';
p1 = pjoint(x, y');
z1 = fjoint(x, y');
z2 = normpdf(x, mu, sigma) * normpdf(y, mu, sigma)';
z3 = z1./z2;

figure(1)
surf(x, y, z1.*p1);
set(gca,'fontsize', 18)
title('The objective function multiplied with the target function')
xlabel('x [m/s]','fontsize', 18)
ylabel('y [m/s)','fontsize', 18)
zlabel('z [W]','fontsize', 18)

figure(2)
surf(x, y, z2);
set(gca,'fontsize', 18)
title('The instrumental function chosen as a normal distribution')
xlabel('x [m/s]','fontsize', 18)
ylabel('y [m/s)','fontsize', 18)
zlabel('z')

figure(3)
surf(x, y, z3);
title('The phi-omega quota')
xlabel('x [m/s]')
ylabel('y [m/s)')
zlabel('z [W]')

%Calculating the expected value of the combined 
draw = Grand(N);
omega = f(draw)./g(draw);
phiomega1 = P(draw).* omega;
tau = mean(phiomega1);
tau2 = 2 * tau;

%Calculating the covariance
drawV1 = normrnd(mu, sigma, N, 1);
drawV2 = normrnd(mu, sigma, N, 1);

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
Pind2 = @(v1, v2) P(v1) + P(v2) < 3.075e6;
N = 1e6;

drawV1 = normrnd(mu, sigma, N, 1);
drawV2 = normrnd(mu, sigma, N, 1);

omega = fjoint(drawV1, drawV2)./gjoint(drawV1,drawV2);
phi1 = Pind1(drawV1, drawV2);
phi2 = Pind2(drawV1, drawV2);
sum1 = phi1 .* omega;
sum2 = phi2 .* omega;
prob1 = mean(sum1);
prob2 = mean(sum2);

%Creating the confidence interval for the two probabilities
lambda95 = norminv(0.975);
std1 = std(sum1);
ci1 = prob1 + [-1 1] * lambda95*(std1/sqrt(N));
std2 = std(sum2);
ci2 = prob2 + [-1 1] * lambda95*(std2/sqrt(N));
probtot = prob1 + prob2;