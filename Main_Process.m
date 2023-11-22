clear all;
clc;
cover_path = 'G:\ExperimentCodes\BOSSbase_1.01(256@256)\';
stego_path = 'G:\stego\bossbase0.1\';
save_stego = 'G:\Results\IA_process\IPDstego_BOSSBase(256@256)\POP_30_NCL_40_0.1bpp\stego_IA_Post_Fast_0.1bpp\';
for i = 1:10000
    cover_image = imread([cover_path,num2str(i),'.pgm']);
    stego_image = imread([stego_path,num2str(i),'.pgm']);
    post_image = uint8(IA_Post(cover_image,stego_image));
    imwrite(post_image,[save_stego,num2str(i),'.pgm']);
    fprintf('IA_Post后处理图片(256@256)序号:%f\n',i);
end


