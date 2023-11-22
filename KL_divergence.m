function score_KL = KL_divergence(cover,stego)
% KL_divergence 作为亲和度函数
% Input：cover,stego
% Output:score_KL：KL-divergence

% Make sure cover and stego sum to 1
if any(cover(:))
    cover = cover / sum(cover(:));
end

if any(stego(:))
    stego = stego / sum(stego(:));
end

% compute Kullback-Leibler Divergence
score_KL = sum(sum(cover.* log(eps + cover./(stego + eps))));

if cover == stego
    score_KL = 0;
end
end
