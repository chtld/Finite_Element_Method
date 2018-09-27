clc
clear
close all
format long
%%
%网格信息
mesh.type = 'triangular';
mesh.left = -1;
mesh.right = 1;
mesh.bottom = -1;
mesh.top = 1;
mesh.hx = 0.25;
mesh.hy = 0.25;

%%
%方程信息
pde.exact_sol.u = @(x, y) x .* y .* (1 - x ./ 2) .* (1 - y) .* exp(x + y);
pde.exact_sol.u_x = @(x, y) (y - y.^2 - x .* y + x .* y.^2) .* exp(x + y)+(x .* y-x .* y.^2-0.5 .* x.^2 .* y+0.5 .* x.^2 .* y.^2) .* exp(x + y);
pde.exact_sol.u_y = @(x, y) (x - 2 .* x .* y - 0.5 .* x.^2 + x.^2 .* y) .* exp(x + y) + (x .* y - x .* y.^2 - 0.5 .* x.^2 .* y + 0.5 .* x.^2 .* y.^2) .* exp(x + y);
pde.rhs = @(x, y) -exp(x + y) .* (y .* (1 - y) .* (1 - x - x.^2 ./ 2) + x .* (1 - x ./ 2) .* (-3 .* y - y.^2));
pde.coef = @(x, y) 1;
pde.bdry.num = 1;
s1 = @(x, y) -1.5 .* y .* (1 - y) .* exp(-1 + y);
s2 = @(x, y) 0.5 .* y .* (1-y) .* exp(1 + y);
s3 = @(x, y) -2 .* x .* (1-x./2) .* exp(x - 1);
s4 = @(x, y) 0;
b1 = @(x, y) ((x == mesh.left && y ~= mesh.bottom && y ~= mesh.top) .* s1(x, y));
b2 = @(x, y) ((x == mesh.right && y ~= mesh.bottom && y ~= mesh.top) .* s2(x,y));
b3 = @(x, y) ((y == mesh.bottom) .* s3(x, y));
b4 = @(x, y) ((y == mesh.top) .* s4(x, y));
pde.bdry.u_dirichlet = @(x, y) (b1(x, y) + b2(x, y) + b3(x, y) + b4(x, y));

%%
FE.basis_type_test = 201;
FE.basis_type_trial = 201;
FE.Gauss_type = 4;

%%
ns = [4, 8, 16, 32, 64];
 fprintf('h        L_inf_err         L2_err        H1_err\n');
for i = 1: length(ns)
    mesh.hx = 1.0 / ns(i);
    mesh.hy = 1.0 / ns(i);
    [result] = Poisson_solver_2D(pde, mesh, FE);
    fprintf('1/%d  \t%e   \t%e   \t%e\n', ns(i), result.error.L_inf, result.errorl.L2, result.error.H1);
end
