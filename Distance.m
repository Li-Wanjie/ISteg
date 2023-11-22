function dist = Distance(cover,stego)
% 欧氏距离作为亲和度函数
% Input：cover,stego
% Output:dist：欧氏距离

  H_kv = [
    -1, +2, -1;
    +2, -4, +2;
    -1, +2, -1
  ]; % KB filter

[A,B] = size(cover);
coverH = imfilter(cover, H_kv, 'symmetric', 'conv', 'same');
stegoH = imfilter(stego, H_kv, 'symmetric', 'conv', 'same');

% coverF = spam686(cover);     % 提取cover的SPAM特征
% stegoF = spam686(stego);     % 提取stego的SPAM特征
% x = (coverF)';              % 列向量转化为行向量
% y = (stegoF)';

x = reshape(coverH,[1,A*B]);
y = reshape(stegoH,[1,A*B]);

dist = (sum((x - y).^2) )^(1/2);
% dist = [x;y];
% dist = pdist(dist);         % 计算各行向量之间的欧氏距离
end

