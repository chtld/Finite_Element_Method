function result=FE_function_2D(x,y, uh_local, vertices, basis_type, der_x,der_y)
  
number_of_local_basis=size(uh_local,1);
result=0;

for k=1:number_of_local_basis
    
    result=result+uh_local(k)*FE_local_basis_2D( x,y,vertices,basis_type,k,der_x,der_y );
    
end


end