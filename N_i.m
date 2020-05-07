function [ nei ] = N_i(i,q,r)
%UPDATE_E Summary of this function goes here
%   Detailed explanation goes here

    nei = [];
    for j=1:length(q)
        temp = [q(i,1),q(i,2); q(j,1),q(j,2)];
        dist = pdist(temp,'euclidean');
        
        if (abs(dist) < r && dist ~= 0)
            nei = [nei j];
        end
    end
end
