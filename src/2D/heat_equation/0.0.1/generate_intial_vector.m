function result = generate_intial_vector(intial, FE)
result = feval(intial, FE.Pb(1, :), FE.Pb(1, :))';
end