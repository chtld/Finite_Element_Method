function error=compute_inf_error(exact_solution_fun, solution, FE, der, t)
%º∆À„Œﬁ«Ó∑∂ ˝ŒÛ≤Ó
T = FE.T;
P = FE.P;
Tb = FE.Tb;
basis_type_trial = FE.basis_type_trial;
number_of_elements = size(T,2);
k=1;
for n=1:number_of_elements
    
    vertices=P(:,T(:,n));
    uh_local=solution(Tb(:,n));
    evaluation_points=select_evalution_points(vertices);
    number_of_evaluation_points=size(evaluation_points,2);
    for i=1:number_of_evaluation_points
       error(k)=abs(feval(exact_solution_fun,evaluation_points(i),t)...
                    -FE_function_1D(evaluation_points(i), uh_local,...
                    vertices, basis_type_trial, der)); 
       k = k+1;
    end
    
end