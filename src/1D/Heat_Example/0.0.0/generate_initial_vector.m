function x = generate_initial_vector(initial_fun, Pb, initial)
number_of_nodes = size(Pb, 2);
x = zeros(number_of_nodes, 1);
for i=1: number_of_nodes
    x(i) = feval(initial_fun, Pb(i), initial);
end
end