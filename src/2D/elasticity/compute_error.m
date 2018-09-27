function [result] = compute_error(exact, exact1, exact2, mesh, FE, solution)
P = mesh.P;
T = mesh.T;
Tb = FE.Tb;
basis_type = FE.basis_type_trial;
error=compute_inf_error(exact, solution,P,T,Tb, basis_type, 0,0);
inf_error = max(error);
L2_error=compute_Hs_error(exact, solution,P,T, Tb, basis_type, 0,0,FE.Gauss_type);
H1_error1=compute_Hs_error(exact1, solution,P,T, Tb, basis_type, 1,0,FE.Gauss_type);
H1_error2=compute_Hs_error(exact2, solution,P,T, Tb, basis_type, 0,1,FE.Gauss_type);
H1_error=sqrt(H1_error1^2+H1_error2^2);
result.error.L_inf = inf_error;
result.errorl.L2 = L2_error;
result.error.H1 = H1_error;
end