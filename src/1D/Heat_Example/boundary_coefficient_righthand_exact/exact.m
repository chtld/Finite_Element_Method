function [ result ] = exact( x,t )
%EXACT 此处显示有关此函数的摘要
%   此处显示详细说明
result = exp(-pi^2*t)*cos(pi*x)+(1-cos(t));


end

