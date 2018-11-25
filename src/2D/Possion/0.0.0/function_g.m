function result=function_g(x,y)
if x==-1
    result = -1.5*y*(1-y)*exp(-1+y);
elseif x==1
    result = 0.5*y*(1-y)*exp(1+y);
elseif y==-1
    result = -2*x*(1-x/2)*exp(x-1);
elseif y==1
    result = 0;
else
    error('the given x is not a boundary nodes');
end
end