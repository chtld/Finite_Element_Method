function int_value=Guass_int_for_error_1D (exact_solution_fun,Guass_weight, Guass_point, uh_local, vertices, basis_type, der)
   %p63
int_value=0;
Gpn=size(Guass_point,2);

for i=1:Gpn
    
    int_value=int_value+Guass_weight(i)*(feval(exact_solution_fun,Guass_point(i))-FE_function_1D(Guass_point(i), uh_local, vertices, basis_type, der))^2;
    
end
