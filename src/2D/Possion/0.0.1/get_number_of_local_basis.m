function [ result ] = get_number_of_local_basis( basis_type )
%GET_NUMBER_OF_LOCAL_BASIS 此处显示有关此函数的摘要
%   此处显示详细说明
if basis_type==201
    result=3;
elseif basis_type==202
    result=6;
else
    warning='no such basis type'
end
end

