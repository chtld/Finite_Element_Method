function int_value=Guass_int_for_error_2D (exact_solution_fun,Gauss_weight, Gauss_point, uh_local, vertices, basis_type, der_x,der_y, time)
   %p63
int_value=0;
Gpn=size(Gauss_point,1);

for i=1:Gpn
    
    int_value=int_value+Gauss_weight(i)*(feval(exact_solution_fun,Gauss_point(i,1),Gauss_point(i,2), time)-FE_function_2D(Gauss_point(i,1),Gauss_point(i,2), uh_local, vertices, basis_type, der_x,der_y))^2;
    
end
