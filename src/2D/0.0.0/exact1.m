function [ result ] = exact1( x,y )
%EXACT1 此处显示有关此函数的摘要
%   此处显示详细说明
result = (y-y^2-x*y+x*y^2)*exp(x+y)+(x*y-x*y^2-0.5*x^2*y+0.5*x^2*y^2)*exp(x+y);

end

