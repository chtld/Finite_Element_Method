function int_value = Gauss_vol_int_trial_test(coe_fun,Gauss_weight,Gauss_point,vertices,basis_type_trial,basis_index_trial,der_trial,...
                                                                                        basis_type_test,basis_index_test,der_test,t)
   int_value = 0;
   Gpn = size(Gauss_point,2);  %得到高斯点的个数
   for i = 1: Gpn              %循环所有的高斯点
      int_value = int_value +Gauss_weight(i)*feval(coe_fun,Gauss_point(i),t)...
          *FE_local_basis_1D(Gauss_point(i),vertices,basis_type_trial,basis_index_trial,der_trial)...
          *FE_local_basis_1D(Gauss_point(i),vertices,basis_type_test,basis_index_test,der_test); 
    
   end
end