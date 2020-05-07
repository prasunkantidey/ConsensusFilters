function [ w ] = weight_design_max_degree( N, nei )
%WEIGHT_DESIGN_MAX_DEGREE Summary of this function goes here
%   Detailed explanation goes here

    w = zeros(N, N);
    for i = 1:N
        for k = 1:length(nei{1,i})
            for j = 1:N
                if i == j
                    w(i,i) = 1 - length(nei{1,i})/N;
                else
                    if nei{1,i}(k) == j
                        w(i,j) = 1/N;
                    end
                end
            end
        end
    end
end

