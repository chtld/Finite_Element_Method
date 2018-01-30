function result=function_g(x)
if x==0
    result = 0;
elseif x==1
    result = cos(1);
else
    warning='the given x is not a boundary nodes';
end
end