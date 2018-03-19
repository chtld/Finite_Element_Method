function result = function_f(x,y)
result = -exp(x+y)*(y*(1-y)*(1-x-x^2/2)+x*(1-x/2)*(-3*y-y^2));
end