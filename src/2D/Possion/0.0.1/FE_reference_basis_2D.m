function [ result ] = FE_reference_basis_2D( x,y,basis_type,basis_index,der_x,der_y )
%FE_REFERENCE_BASIS_2D 此处显示有关此函数的摘要
%   此处显示详细说明
if basis_type==201
    if basis_index==1
       if der_x==0&&der_y==0
           result=-x-y+1;
       elseif der_x==1&&der_y==0
           result=-1;
       elseif der_y==1&&der_x==0
           result=-1;
       elseif der_x+der_y>=2
           result =0;
       else
           error('no such der');
       end
    elseif basis_index==2
        if der_x==0&&der_y==0
            result=x;
        elseif der_x==1&&der_y==0
            result=1;
        elseif der_x==0&&der_y==1
            result=0;
        elseif der_x+der_y>=2
            result=0;
        else
            error('you are wrong');
        end
    elseif basis_index==3
        if der_x==0&&der_y==0
            result=y;
        elseif der_x==1&&der_y==0
            result=0;
        elseif der_x==0&&der_y==1
            result=1;
        elseif der_x+der_y>=2
            result=0;
        else
            error('you are wrong');
        end
            
    else
        error('you are wrong');
    end
elseif basis_type==202
     if basis_index==1
       if der_x==0&&der_y==0
           result=2*x^2+2*y^2+4*x*y-3*y-3*x+1;    
       elseif der_x==1&&der_y==0
           result=4*x+4*y-3;
       elseif der_y==1&&der_x==0
           result=4*y+4*x-3;
       elseif der_x==1&&der_y==1
           result=4;
       elseif der_x==2&&der_y==0
           result=4;
       elseif der_y==2&&der_x==0
           result=4;
       elseif der_x+der_y>=3
           result =0;
       else
           error('no such basis fun');
       end
     elseif basis_index==2
       if der_x==0&&der_y==0
           result=2*x^2-x;    
       elseif der_x==1&&der_y==0
           result=4*x-1;
       elseif der_x==2&&der_y==0
           result=4;
       elseif der_y>=1||der_x>=3
           result=0;
       else
           error('no such basis fun');
       end
     elseif basis_index==3
       if der_x==0&&der_y==0
           result=2*y^2-y;    
       elseif der_y==1&&der_x==0
           result=4*y-1;
       elseif der_y==2&&der_x==0
           result=4;
       elseif der_x>=1||der_y>=3
           result=0;
       else
           error('no such basis fun');
       end
     elseif basis_index==4
       if der_x==0&&der_y==0
           result=-4*x^2-4*x*y+4*x;    
       elseif der_x==1&&der_y==0
           result=-8*x-4*y+4;
       elseif der_y==1&&der_x==0
           result=-4*x;
       elseif der_x==1&&der_y==1
           result=-4;
       elseif der_x==2&&der_y==0
           result=-8;
       elseif der_y==2&&der_x==0
           result=0;
       elseif der_x+der_y>=3
           result =0;
       else
           error('no such basis fun');
       end
     elseif basis_index==5
       if der_x==0&&der_y==0
           result=4*x*y;    
       elseif der_x==1&&der_y==0
           result=4*y;
       elseif der_y==1&&der_x==0
           result=4*x;
       elseif der_x==1&&der_y==1
           result=4;
       elseif der_x>=2||der_y>=2
           result=0;
       else
           error('no such basis fun');
       end
       elseif basis_index==6
       if der_x==0&&der_y==0
           result=-4*x^2-4*x*y+4*x;    
       elseif der_x==1&&der_y==0
           result=-8*x-4*y+4;
       elseif der_y==1&&der_x==0
           result=-4*x;
       elseif der_x==1&&der_y==1
           result=-4;
       elseif der_x==2&&der_y==0
           result=-8;
       elseif der_y==2&&der_x==0
           result=0;
       elseif der_x+der_y>=3
           result =0;
       else
           error('no such basis fun');
       end
     else
       error('no such basis function');
     end
    
else

end

