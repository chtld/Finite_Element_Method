clc
clear
close all
format long
%%
%网格信息
mesh.type = 'triangular';
mesh.left = 0;
mesh.right = 1;
mesh.bottom = 0;
mesh.top = 1;
mesh.hx = 0.25;
mesh.hy = 0.25;

%%
%方程信息
pde.exact_sol.u = @(x, y, t) exp(x + y + t);
pde.exact_sol.u_x = @(x, y, t) exp(x + y + t);
pde.exact_sol.u_y = @(x, y, t) exp(x + y + t);
pde.intial = @(x, y, t) exp(x + y);
pde.rhs = @(x, y, t) -3 * exp(x + y + t);
pde.coef1 = @(x, y, t) 1;
pde.coef2 = @(x, y, t) 2;
pde.bdry.num = 1;
pde.initialtime = 0;
pde.endtime = 1;
pde.dt = 0.00001;
pde.theta = 0;
s1 = @(x, y, t) exp(y + t);
s2 = @(x, y, t) exp(1 + y + t);
s3 = @(x, y, t) exp(x + t);
s4 = @(x, y, t) exp(1 + x + t);
b1 = @(x, y, t) ((x == mesh.left && y ~= mesh.bottom && y ~= mesh.top) .* s1(x, y, t));
b2 = @(x, y, t) ((x == mesh.right && y ~= mesh.bottom && y ~= mesh.top) .* s2(x, y, t));
b3 = @(x, y, t) ((y == mesh.bottom) .* s3(x, y, t));
b4 = @(x, y, t) ((y == mesh.top) .* s4(x, y, t));
pde.bdry.u_dirichlet = @(x, y, t) (b1(x, y, t) + b2(x, y, t) + b3(x, y, t) + b4(x, y, t));

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
%     pde.dt = 1.0 / ns(i)^2;
    [result] = heat_solver_2D(pde, mesh, FE);
    fprintf('1/%d  \t%e   \t%e   \t%e\n', ns(i), result.error.L_inf, result.errorl.L2, result.error.H1);
end
