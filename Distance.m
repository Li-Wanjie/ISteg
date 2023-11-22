function dist = Distance(cover,stego)
% ŷ�Ͼ�����Ϊ�׺ͶȺ���
% Input��cover,stego
% Output:dist��ŷ�Ͼ���

  H_kv = [
    -1, +2, -1;
    +2, -4, +2;
    -1, +2, -1
  ]; % KB filter

[A,B] = size(cover);
coverH = imfilter(cover, H_kv, 'symmetric', 'conv', 'same');
stegoH = imfilter(stego, H_kv, 'symmetric', 'conv', 'same');

% coverF = spam686(cover);     % ��ȡcover��SPAM����
% stegoF = spam686(stego);     % ��ȡstego��SPAM����
% x = (coverF)';              % ������ת��Ϊ������
% y = (stegoF)';

x = reshape(coverH,[1,A*B]);
y = reshape(stegoH,[1,A*B]);

dist = (sum((x - y).^2) )^(1/2);
% dist = [x;y];
% dist = pdist(dist);         % �����������֮���ŷ�Ͼ���
end

