function plot_figure(Mesh, FE, result)
if FE.basis_type_test == 202
    Mesh.hx = Mesh.hx / 2;
    Mesh.hy = Mesh.hy / 2;
end
x = Mesh.left: Mesh.hx: Mesh.right;
y = Mesh.bottom: Mesh.hy: Mesh.top;
[X, Y] = meshgrid(x, y);
exact = reshape(result.exact_solution, size(X));
numerical = reshape(result.solution, size(X));
subplot(121);
mesh(X, Y, exact);
title('exact solution');
subplot(122);
mesh(X, Y, numerical);
title('numerical solution');
end