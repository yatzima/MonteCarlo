%% 3 A Naive approach - Random Walk
%Defining 
close all
d = 2;
N = 1e2;
noWalk = zeros(1, 2);
walk = zeros(N, 2);
     
for i = 3 : N + 1
    newDir = unidrnd(3, 1, 2) - 2;
     
    while newDir == noWalk
       newDir = unidrnd(3, 1, 2) - 2;
    end
     
    next = newDir + walk(i-2, :);
     
    if ismember(next,walk,'rows') == ones(1, 2)
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
n = zeros(8,2);
counter = 0;
%for k  = 1:100
for i = 3 : N + 1
    % n is the neighbour-vector
    n(1,:) = walk(i-2,:) + [-1,0];
    n(2,:) = walk(i-2,:) + [-1,1];
    n(3,:) = walk(i-2,:) + [0,1];
    n(4,:) = walk(i-2,:) + [1,1];
    n(5,:) = walk(i-2,:) + [1,0];
    n(6,:) = walk(i-2,:) + [1,-1];
    n(7,:) = walk(i-2,:) + [0,-1];
    n(8,:) = walk(i-2,:) + [-1,-1];
     
    if ismember(n,walk,'rows') == ones(8,1)
        break
    end
     
    newDir = unidrnd(3, 1, 2) - 2;
    while ismember(newDir + walk(i-2,:),walk,'rows') == ones(1,2) 
       newDir = unidrnd(3, 1, 2) - 2;
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