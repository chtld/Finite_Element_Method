function int_value = Gauss_int_vector_test(rightside,Gauss_weight,Gauss_point,vertices,basis_type_test,basis_index_test,der_test,t)
%GAUSS_INT_VECTOR_TEST 此处显示有关此函数的摘要
%   此处显示详细说明
   int_value = 0;
   Gpn = size(Gauss_point,2);
   for i = 1: Gpn
      int_value = int_value +Gauss_weight(i)*feval(rightside,Gauss_point(i),t)...
          *FE_local_basis_1D(Gauss_point(i),vertices,basis_type_test,basis_index_test,der_test); 
   end

end

