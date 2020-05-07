function [ E ] = update_E( m, V, w)
%UPDATE_E Summary of this function goes here
%   Detailed explanation goes here

    sum_w_ii = 0;
    numerator = 0;
    
    for i = 1:length(V)
        numerator = numerator + w(i,i)*m(i);
        sum_w_ii = sum_w_ii + w(i,i);        
    end

    E = numerator/sum_w_ii;    
end

