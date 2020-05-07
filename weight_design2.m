function [ w ] = weight_design2( cw2, V, nei )
%WEIGHT_DESIGN1 Summary of this function goes here
%   Detailed explanation goes here

    w = zeros(length(V), length(V));
    for i = 1:length(V)
       w(i,i) = cw2/V(i);
       for k = 1:length(nei{1,i})
           for j = 1:length(V)
               if i ~= j
                   if (nei{1,i}(k)==j)
                       w(i,j) = (1-w(i,i))/length(nei{1,i});
                   end
               end
           end
       end
    end
    
end
