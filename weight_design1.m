function [ w ] = weight_design1( cw1, V )
%WEIGHT_DESIGN1 Summary of this function goes here
%   Detailed explanation goes here

    w = zeros(length(V), length(V));
    for i = 1:length(V)
        for j = 1:length(V)
            if i~=j
                w(i,j) = cw1/(V(i) + V(j));
            end
        end
        w(i,i) = 1 - sum(w(i,:));
    end

end
