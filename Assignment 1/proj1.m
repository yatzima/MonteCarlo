%%Power production of a wind turbine
clc
clear
close all
load('powercurve_V112.mat');

lambda = [10.6 9.7 9.2 8.0 7.8 8.1 7.8 8.1 9.1 9.9 10.6 10.6];
k = [2.0 2.0 2.0 1.9 1.9 1.9 1.9 1.9 2.0 1.9 2.0 2.0];

month = 1;
N = 1000;
lin = linspace(3,25);
phi = P;
const = 6;

%Defining the stochastic wind speed V for different months
f = @(v, month) wblpdf(v, lambda(month), k(month)); 
%Defining the g function
g = @(v, month) gampdf(v, lambda(month) - const, k(month));

%Defining the random-generator for different months
Frand = @(month) wblrnd(lambda(month), k(month), 1, N);
Grand = @(month) gamrnd(lambda(month) - const, k(month), 1, N);

%Defining the inverse
FV = @(v, month) 1 - exp(-(v/lambda(month))^k(month)); 
FU = @(U, month) lambda(month)*(-log(U*(FV(25,month)-FV(3,month)) + FV(3,month))).^(1/k(month));  % inverse truncated

%Crude monte carlo
draw1 = Frand(1);
y1 = P(draw1);
tau1 = mean(P(draw1));
std1 = std(P(draw1));

%Truncated
u = rand(1, N);
draw2 = FU(u, month);
tau2 = mean(P(draw2)) * (FV(25,month) - FV(3, month));
std2 = std(P(draw2)) * (FV(25,month) - FV(3, month)); %Ska man multiplicera?

%Importance Sampling
draw3 = Grand(month);
y3 = g(draw3, month);
phiomega = phi(draw3)*(f(draw3, month)/g(draw3, month));
tau3 = mean(phiomega);
std3 = std(phiomega);

%Antithetic Sampling
u = rand(N, 1);
draw41 = FU(u, month);
draw42 = FU(1-u, month);
V1 = P(draw41);
V2 = P(draw42);
W = (V1+V2)/2;
tau4 = mean(W) * (FV(25,month) - FV(3, month));
std4 = std(W) * (FV(25,month) - FV(3, month));

%CI1 - Crude
lambda95 = norminv(0.975);
ci1 = [tau1 - lambda95*(std1/sqrt(N)), tau1 + lambda95*(std1/sqrt(N))];

%CI2 - Truncated
ci2 = [tau2 - lambda95*(std2/sqrt(N)), tau2 + lambda95*(std2/sqrt(N))];

%CI3 - Importance sampling
ci3 = [tau3 - lambda95*(std3/sqrt(N)), tau3 + lambda95*(std3/sqrt(N))];

%CI4 - Antithetic sampling
ci4 = [tau4 - lambda95*(std4/sqrt(N)), tau4 + lambda95*(std4/sqrt(N))];

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

%%
%Plot the quote and the dist for f and g
figure(1)
hold on
draw1 = zeros(N, 12);
y1 = zeros(N, 12);
counter = 1;

for month = 1:12
    draw1(:, counter) = Frand(month)';
    y1(:, counter) = f(draw1(:, counter), month);
    plot(draw1(:,counter), y1(:,counter), '*')
    counter = counter + 1;
end
plot(draw3, y3, 'o')
legend('Function f','Function g')

%Defining the random-generator for different months
Frand = @(month, N) wblrnd(lambda(month), k(month), 1, N);
Grand = @(month, N) gamrnd(lambda(month) - 5, k(month), 1, N);

N = 1:50:10000;
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

for N = 1:50:10000;
    %Crude monte carlo
    draw1 = Frand(1, N);
    y1 = P(draw1);
    tau1(counter) = mean(P(draw1));
    std1 = std(P(draw1));
    
    %Truncated
    u = rand(1, N);
    draw2 = FU(u, month);
    tau2(counter) = mean(P(draw2)) * (FV(25,month) - FV(3, month));
    std2 = std(P(draw2)) * (FV(25,month) - FV(3, month));
    
    %Importance Sampling
    draw3 = Grand(month, N);
    y3 = g(draw3, month);
    phiomega = phi(draw3)*(f(draw3, month)/g(draw3, month));
    tau3(counter) = mean(phiomega);
    std3 = std(phiomega);
    
    %Antithetic Sampling
    u = rand(N, 1);
    draw41 = FU(u, month);
    draw42 = FU(1-u, month);
    V1 = P(draw41);
    V2 = P(draw42);
    W = (V1+V2)/2;
    tau4(counter) = mean(W) * (FV(25,month) - FV(3, month));
    std4 = std(W) * (FV(25,month) - FV(3, month));
    
    %CI1 - Crude MC
    ci1(counter, :) = [tau1(counter) - lambda95*(std1/sqrt(N)), tau1(counter) + lambda95*(std1/sqrt(N))];

    %CI2 - Truncated
    ci2(counter, :) = [tau2(counter) - lambda95*(std2/sqrt(N)), tau2(counter) + lambda95*(std2/sqrt(N))];

    %CI3 - Importance sampling
    ci3(counter, :) = [tau3(counter) - lambda95*(std3/sqrt(N)), tau3(counter) + lambda95*(std3/sqrt(N))];

    %CI4 - Antithetic sampling
    ci4(counter, :) = [tau4(counter) - lambda95*(std4/sqrt(N)), tau4(counter) + lambda95*(std4/sqrt(N))];

    
    counter = counter + 1;
end

N = 1:50:10000;

figure(3)
plot(N,tau1, 'r');
hold on
plot(N,tau2, 'g');
plot(N,tau3, 'b');
plot(N,tau4, 'k');
title('Expected value based on difference MC-methods for January')
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
