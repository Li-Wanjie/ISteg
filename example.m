clear all;
clc;
% %cover = [0.2 0.4 0.4];
% % stego = [0.3 0.2 0.5];
cover = double(imread('1.pgm'));
stego = double(imread('s1.pgm'));
% 
% enhanced_stego = stego_post_process(cover, stego);
% % imshow(enhanced_stego,[]);
% % figure;
% % imshow((double(cover)-double(enhanced_stego)+1)/2);
% % 
% % figure;
% % imshow((double(stego)-double(enhanced_stego)+1)/2);
% 
% score_KL1 = KL_divergence(cover, stego);
% score_KL2 = KL_divergence(cover, enhanced_stego);

