% ����������
% [-1, 1] x [-1, 1] x [-1, 1]��������������������ʷ�
node  = [-1, -1, -1; 1, -1, -1; 1, 1, -1; -1, 1, -1;
         -1, -1 , 1; 1, -1, 1; 1, 1, 1; -1, 1, 1]; %�ڵ�����
elem = [1, 2, 3, 7; 1, 6, 2, 7; 1, 5 , 6, 7; 1, 8, 5, 7; 1, 4, 8, 7; 1, 3, 4, 7];
clf;
showmesh3(node, elem, [], 'FaceAlpha', 0.25);
view([-53, 8]); axis on; % ѡ�񿴵ĽǶ�
findelem3(node, elem); % ������������
findnode3(node); % ���ƽڵ���