function result=FE_function_1D(x, uh_local, vertices, basis_type, der)
  
number_of_local_basis=size(uh_local,1);
result=0;

for k=1:number_of_local_basis
    
    result=result+uh_local(k)*FE_local_basis_1D(x, vertices, basis_type, k, der);
    
end


end