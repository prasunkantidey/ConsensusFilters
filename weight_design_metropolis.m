function [ w ] = weight_design_metropolis( N, nei )
%WEIGHT_DESIGN_METROPOLICE Summary of this function goes here
%   Detailed explanation goes here

    w = zeros(N, N);
    
    for i = 1:N
        for k = 1:length(nei{1,i})
            for j = 1:N
                if i ~= j
                    if j == nei{1,i}(k)
                        w(i, j) = 1/(1 + max(length(nei{1,i}), length(nei{1,j})));
                    end
                end
            end
        end
        w(i,i) = 1 - sum(w(i,:));
    end
    
end

