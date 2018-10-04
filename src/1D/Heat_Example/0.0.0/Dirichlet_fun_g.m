function result=Dirichlet_fun_g(x,t)
if x==0
    result = 1 + t;
elseif x==1
    result = 2 + t;
else
    error('the given x is not a boundary nodes');
end
end