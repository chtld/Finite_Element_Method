function error=compute_inf_error(exact_solution_fun, solution,P,T,Tb, basis_type, der_x, der_y, time)
number_of_elements = size(T,2);
for n=1:number_of_elements    
    vertices=P(:,T(:,n));
    uh_local=solution(Tb(:,n));
    evaluation_points=select_evalution_points(vertices);
    number_of_evaluation_points=size(evaluation_points,2);
    
    for i=1:number_of_evaluation_points
       error(i)=abs(feval(exact_solution_fun, evaluation_points(1,i), evaluation_points(2,i), time)-FE_function_2D(evaluation_points(1,i),evaluation_points(2,i), uh_local, vertices, basis_type, der_x,der_y)); 
    end
    
end
end

function points=select_evalution_points(vertices)
    [Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(4);
    [Gauss_weight,Gauss_point]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices);
    points=[vertices Gauss_point' ];
end
