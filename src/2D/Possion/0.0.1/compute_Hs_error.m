function error=compute_Hs_error(exact_solution_fun, solution,P,T, Tb, basis_type, der_x,der_y,Gauss_type)

error=0;
number_of_elments = size(Tb,2);
for n=1:number_of_elments
    vertices = P(:,T(:,n));
    uh_local=solution(Tb(:,n));
    [Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(Gauss_type);
    [Gauss_weight,Gauss_point]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices);
    
    local_error=Guass_int_for_error_2D(exact_solution_fun,Gauss_weight, Gauss_point, uh_local, vertices, basis_type, der_x,der_y);
    error=error+local_error;
    
end
error=sqrt(error);
end
