clc
clear all
close all

% From provided dataset
dxy = 0.5;
field_min = -6;
field_max = 6;
x = (field_min:dxy:field_max)';
y = (field_min:dxy:field_max)';
F = 12*12;
load('Scalar_field_data_open_by_MatLab_only.mat');

% Plot the Provided Scalar field
surf(x,y,-F);
xlabel('x1'); 
ylabel('x2'); 
zlabel('Probability Density for given sample');

% Initialization
N = 30;     % Number of agents
m = 2;
K = 25*25;  % Number of cells

% Field, F setup
for i = 1:length(x)
    for j = 1:length(y)
        q_c((i-1)*25+j,:)=[x(i) y(j)];
    end
end

% Initial nodes positions
q = field_min + (field_max-field_min) * rand(N,m);
figure;
plot(q(:,1),q(:,2),'k*')
hold on

% From provided value
dt = 0.01;
endt = 2;
t = 0:dt:endt;
c_v = 0.01;
R = 5;
rs = R*ones(N,1);

plot_sensing_range(q(:,1), q(:,2), rs);

% Setup neighbor network
some_really_small_noise = 0.000001;
for k = 1:K
    for i = 1:N
        dist = norm(q(i,:) - q_c(k,:));
        if dist <= rs(i)
            O(i, k) = 1;
            V(i, k) = ((norm(q(i,:)-q_c(k,:)).^2)+c_v)/(rs(i)^2);
        else
            O(i, k) = 0;                                    % Not observable
            V(i, k) = 1000000;                              % For handling wrong calculation.
        end
        n1(i, k) = sqrt(V(i, k)) * (rand/4) ;               % Guasssian noise N(0,V1(i)) 
        m1(i, k) = O(i, k) * (F(k) + n1(i, k)) + some_really_small_noise;
        x1{1}(i, k) = m1(i, k);
        nei{i} = N_i(i, q, rs(i));
        c1ww(i) = (2*c_v) / ((rs(i)^2) * (length(nei{i})));
    end
end

% Plot neighbor network.
figure;
plot(q(:,1),q(:,2),'b*')
hold on
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

err_thresh = 0.0001;
err = 10000;
l=1;
cw1 = min(c1ww);
while(err >= err_thresh && l < 40000)
    for k = 1:K
        w = weight_design1( cw1, V(:,k) );
        x1{l+1}(:,k) = update_x(x1{l}(:,k), w);
        E(k,l+1) = update_E(m1,V(:,k),w); 
        EE(k,l+1) = max(x1{l+1}(:,k)-E(k,l+1));
        EE_error_calc(l+1) = max(EE(k,l+1));
    end
    err = min(wkeep(EE_error_calc,length(EE_error_calc)-1,'r'));
    if mod(l,2000) == 0
        disp(err)
    end
    l = l+1;
end

disp(l);

for i = 1:length(x)
    for j = 1:length(y)
        F_1((i-1)*25+j)=F(i,j);
    end
end

F_hat = (1/N) * sum(x1{1,l-1}(:,1:K)); 

figure;
plot(1:625, F_1-F_hat);

EE = F_1 - F_hat;
for i = 1:length(x)
    for j = 1:length(y)
        ff(j,i) = F_hat((i-1)*25 + j);
        ee(j,i) = EE((i-1)*25 + j);
    end
end

figure;
surf(x,y,-ff);
xlabel('x1'); ylabel('x2'); zlabel('Probability Density'); 

figure;
surf(x,y,-F);
xlabel('x1'); ylabel('x2'); zlabel('Probability Density');

figure;
surf(x,y,-ee);
xlabel('x1'); ylabel('x2'); zlabel('Probability Density');