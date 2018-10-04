% L形区域
% [-1, 1] x [-1, 1] 的正方形区域除去右下角[-1, 0] x [0, 1]的正方形区域
node = [1, 0; 1, 1; 0, 1; -1, 1; -1, 0; -1, -1; 0, -1; 0, 0]; %网格点坐标矩阵，P矩阵
elem = [1, 2, 8; 3, 8, 2; 8, 3, 5; 4, 5, 3; 7, 8, 6; 5, 6, 8]; %网格，连接关系，T矩阵
showmesh(node, elem); axis on;
findelem(node, elem); %绘制三角形编号
findnode(node); %绘制节点编号

% 对网格进行3次一致的细分
for i = 1: 3
    [node, elem] = uniformrefine(node, elem);
end
showmesh(node, elem);