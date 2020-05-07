clc
clear all
close all

N = 10;                 % number of agents
m = 2;                  % m = 2 for 2-D plane
% initial positions uniformly
q = sqrt(11)*rand(N,m); %initial positions uniformly in [0,4]^2 box
p_min = 0;
p_max = 4;
p = p_max+(p_min-p_max)*rand(N,m);
cell_pos = [(p_max-p_min)/2 (p_max-p_min)/2];

plot(q(:,1),q(:,2),'b*')
for i = 1:length(q)
   text(q(i,1), q(i,2), sprintf('node %d', i)); 
end
hold on

nei = {};
dt = 0.01;
endt = 2;
t = 0:dt:endt;
F1 = 50;
c_v = 0.01;
R = 1.6;
rs = R*ones(N,1);
q_bar = (1/N)*sum(q);

for i = 1:N
    V1(i) = ((norm(q(i,:)-q_bar).^2)+c_v)/(rs(i)^2);
    n1(i) = sqrt(V1(i))*(randn/10);                     % Guasssian noise N(0,V1(i)) 
    nei{1,i} = N_i(i,q,rs(i));
    x(i,1) = F1 + n1(i);                                % Initial measurement
end

% Plot neighbor network.
for i = 1:N
    node_x = q(i,1);
    node_y = q(i,2);
    sprintf('node %d:', i);

    for j = 1:length(nei{1,i})
        nei_x = q(nei{1,i}(j),1);
        nei_y = q(nei{1,i}(j),2);
        line([node_x, nei_x], [node_y, nei_y]);
        
        hold on;
    end
end
hold off

% ML is the average of the measurements.
theta_ML = 1/N * ones(N,1) * sum (x(:,1));
err_thresh = 0.0001;

% For Max Degree
err = 100000;
t=1;
x1(:,t) = x(:,t);

while(err >= err_thresh) 
    w = weight_design_max_degree(N, nei);       % Using equation (8)
    x1(:,t+1) = update_x(x1(:,t), w);           % Using equation (4)
    E1(:,t+1) = x1(:,t+1) - theta_ML;
    err = max(E1(:,t+1));
    t = t+1;
    disp(t);
end

% Plot values
figure;
plot(1:10,x1(:,t),'*',1:10,x1(:,1),'r');
dim = [.6 .6 .3 .3];
str = sprintf('Number of Iteration: %d', t);
annotation('textbox',dim,'String',str,'FitBoxToText','on');

% For Metropolis
err = 100000;
t=1;
x2(:,t) = x(:,t);

while(err >= err_thresh) 
    w = weight_design_metropolis(N, nei);       % Using equation (9)
    x2(:,t+1) = update_x(x2(:,t), w);           % Using equation (4)
    E2(:,t+1) = x2(:,t+1) - theta_ML;
    err = max(E2(:,t+1));
    t = t+1;
    disp(t);
end

figure;
plot(1:10,x2(:,t),'*',1:10,x2(:,1),'r');
dim = [.6 .6 .3 .3];
str = sprintf('Number of Iteration: %d', t);
annotation('textbox',dim,'String',str,'FitBoxToText','on');