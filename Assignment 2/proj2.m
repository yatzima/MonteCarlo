%% 3 A Naive approach - Random Walk
%Defining 
close all
d = 2;
N = 1e3;
walk = zeros(N, 2);
n = [-1,0; 0,1; 1,0; 0,-1];

for i = 3 : N + 1
    newDir = datasample(n,1);
    next = newDir + walk(i-2, :);
        
    if ismember(next, walk, 'rows') == ones(1, 2)
        walk(i-1,:) = next;
        break
    end 
    
    walk(i - 1, :) = next;
end
plot(walk(1:i-2, 1), walk(1:i-2, 2),'-*')
hold on
plot(walk(i-1,1),walk(i-1,2),'-*r')
walk(1:i-1,:)
%%
%Defining the instrumental distribution gn as a standard random walk in Z2

%Plot the SAW
%uniform distr
%plot vector on 
%% 4 Improving - actual self avoiding walk
close all
d = 2;
N = 1e3;
walk = zeros(N, 2);
n = [-1,0; 0,1; 1,0; 0,-1];
%for k  = 1:100
for i = 3 : N + 1
    % n is the neighbour-vector
    neigh(1,:) = walk(i-2,:) + [-1,0];
    neigh(2,:) = walk(i-2,:) + [0,1];
    neigh(3,:) = walk(i-2,:) + [1,0];
    neigh(4,:) = walk(i-2,:) + [0,-1];
    
    if ismember(neigh,walk,'rows') == ones(4,1)
        break
    end
    
    newDir = datasample(n,1);
    
    while ismember(newDir + walk(i-2,:),walk,'rows') == ones(1,2)
        newDir = datasample(n,1);
    end
    
    next = newDir + walk(i-2, :);
    walk(i - 1, :) = next;
end
%if i == N+1
%    counter = counter + 1;
%end
%end 
plot(walk(1:i-2, 1), walk(1:i-2, 2),'-*')
hold on
plot(walk(i-2, 1), walk(i-2,2), '*r')
%% 4 SIS Sample
N = 1e3;
X_0 = zeros(N,2);

for i = 1:N
    X_0(i) = 
end