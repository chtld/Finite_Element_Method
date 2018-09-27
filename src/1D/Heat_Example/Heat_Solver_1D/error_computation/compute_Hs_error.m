function error=compute_Hs_error(exact_solution_fun, solution,P,T, Tb, basis_type, der,Gauss_type,t)

error=0;
number_of_elments = size(Tb,2);
for n=1:number_of_elments
    vertices = P(:,T(:,n));
    uh_local=solution(Tb(:,n));
    [Gauss_weight,Gauss_point]=generate_Gauss_formula(vertices, Gauss_type);
    
    local_error=Guass_int_for_error_1D(exact_solution_fun,Gauss_weight, Gauss_point, uh_local, vertices, basis_type, der, t);
    error=error+local_error;
    
end
error=sqrt(error);
