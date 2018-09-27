function [solution] = compute_exact_solution(exact_solution_fun, Pb,Tb)

k = 1;
number_of_points = size(Pb,2);
for n=1:number_of_points    
    solution(n, 1) = feval(exact_solution_fun, Pb(1, n), Pb(2, n));
end

end