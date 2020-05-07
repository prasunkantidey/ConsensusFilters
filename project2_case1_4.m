clc
clear all
close all
weight_design = 1;

N = 10;                 % number of agents
m = 2;                  % m = 2 for 2-D plane
q = sqrt(11)*rand(N,m); %initial positions uniformly in [0,4]^2 box

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
    n1(i) = sqrt(V1(i))*(randn/10);                         % Guasssian noise N(0,V1(i)) 
    m1(i) = F1+n1(i);                                       % Initial measurement
    x1(i,1) = m1(i);
    nei{1,i} = N_i(i,q,rs(i));
    if length(nei{1,i}) == 0
        c1ww(i) = 0;
    else
        c1ww(i) = (2*c_v)/((rs(i)^2)*(length(nei{1,i})));   % c1ww
    end
    c2ww(i) =c_v/(rs(i)^2);                                 % c2ww
end

% Identify Maximum and Minimum neighbored node.
minimum_nei = 1;
maximum_nei = 1;
for i = 1:N
    if length(nei{1,i}) < length(nei{1,minimum_nei})
        minimum_nei = i;
    end
    if length(nei{1,i}) > length(nei{1,maximum_nei})
        maximum_nei = i;
    end
end

priv_minimum_nei_count = length(nei{1,minimum_nei});
priv_maximum_nei_count = length(nei{1,maximum_nei});

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

err_thresh = 0.0001;
err = 100000;
cw1 = min(c1ww);
cw2 = min(c2ww);
if weight_design == 1
    cw = cw1;
else
    cw = cw2;
end

l=1;

while(err>=err_thresh) 
    if weight_design == 1
        w = weight_design1(cw, V1);                 % Using equation (16)
    else
        w = weight_design2(cw, V1, nei);            % Using equation (22)
    end
    
    x1(:,l+1) = update_x(x1(:,l), w);               % Using the equation (11)
    E(l+1) = update_E(m1,V1,w);                     % Estimate value
    EE(l+1) = max(x1(:,l+1))-E(l+1);                % Store err to check the max error.
    err = min(wkeep(EE, length(EE)-1, 'r'));
    if mod(l,100) == 0
        disp(err)
    end
    l = l+1;
    
    if l == 300                                     % Update location
        q = sqrt(8)*rand(N,m);                     
        
        figure;
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
            n1(i) = sqrt(V1(i))*(randn/10);                         % Guasssian noise N(0,V1(i))
            m1(i) = F1+n1(i);                                       % Initial measurement
            x1(i,l) = m1(i);
            nei{1,i} = N_i(i,q,rs(i));
            if length(nei{1,i}) == 0
                c1ww(i) = 0;
            else
                c1ww(i) = (2*c_v)/((rs(i)^2)*(length(nei{1,i})));   % c1ww
            end
            c2ww(i) = c_v/(rs(i)^2);                                % c2ww
        end

        % Plot new neighbor network
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

        cw1 = min(c1ww);
        cw2 = min(c2ww);
        if weight_design == 1
            cw = cw1;
        else
            cw = cw2;
        end
        
    end
    disp(l);
end

% Plot values
figure; 
plot(1:10,E(:,l),'*',1:10,x1(:,1),'r');

consensus = x1 - E(l);

% Plot consensus
figure;
for i = 1:N
    plot(consensus(i,:));
    hold on;
end

% Plot smallest and highest neighbor nodes convergence.
figure;
plot(consensus(minimum_nei,:));
hold on;
plot(consensus(maximum_nei,:));
legend(sprintf('Smallest Neighbor (Node %d), Priv %d Neighbors, Current %d Neighbors',minimum_nei,priv_minimum_nei_count,length(nei{1,minimum_nei})),...
    sprintf('Highest Neighbor (Node %d), Priv %d Neighbors, Current %d Neighbors',maximum_nei,priv_maximum_nei_count,length(nei{1,maximum_nei})));