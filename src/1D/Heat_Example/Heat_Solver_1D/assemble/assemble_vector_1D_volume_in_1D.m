function b = assemble_vector_1D_volume_in_1D(rightside,P,T,Pb,Tb,vector_size,basis_type_test,Gauss_type, t)
%ASSEMBLE_VECTOR_1D_VOLUME_IN_1D 此处显示有关此函数的摘要
%   此处显示详细说明
    b = zeros(vector_size,1);
    number_of_elements=size(T,2);
    number_of_local_basis=get_number_of_local_basis(basis_type_test);
    for n=1: number_of_elements
        vertices=P(:,T(:,n));
        [Gauss_weight,Gauss_point]=generate_Gauss_formula(vertices,Gauss_type);
        for beta = 1:number_of_local_basis
            int_value=Gauss_int_vector_test(rightside,Gauss_weight,Gauss_point,vertices,basis_type_test,beta,0,t);
             b(Tb(beta,n))=b(Tb(beta,n))+int_value;
        end
    end


end

