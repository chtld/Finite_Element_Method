function [ A, b ] = treat_Robin_boundary( A, b, boundarynodes,coef_fun,Robin_fun_pb,Robin_fun_qb,Pb)
nbn=size(boundarynodes,2);
for k=1:nbn
   if boundarynodes(1,k) ==-3
      i = boundarynodes(2,k);
      b(i) = b(i)+feval(Robin_fun_pb,Pb(i))*feval(coef_fun,Pb(i));
      A(i,i)=A(i,i)+feval(Robin_fun_qb,Pb(i))*feval(coef_fun,Pb(i));
   end
end

end

function result = Robin_fun_pb(x)
result = 1;
end
function result = Robin_fun_qb(x)
result = 1;
end