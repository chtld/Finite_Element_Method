function  error = compute_error( pde, FE, final )
%COMPUTE_ERROR 此函数为误差计算模块
%   此处显示详细说明
%   误差计算
inf_error1 = compute_inf_error(pde.exact_solution, pde.solution, FE, 0, final);      %计算无穷范数误差
inf_error = max(inf_error1);
L2_error = compute_Hs_error(pde.exact_solution, pde.solution, FE, 0, final);         %计算L2范数误差
H1_semi_error = compute_Hs_error(pde.gradient, pde.solution, FE, 1, final);          %计算H1半范误差
error = [inf_error, L2_error, H1_semi_error];
fprintf('%7.4e  %7.4e  %7.4e\n',inf_error,L2_error,H1_semi_error);
end

