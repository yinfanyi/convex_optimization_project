b=8; % block size
R=4; % redundancy factor  冗余因素
K=R*b^2; % number of atoms in the dictionary

sigma = 30; 

load McM01_noise.mat u_n %加载噪声图像

%以下等价于sample code中的imgsc函数
%u_n = uint8(u_n);
%imagesc(u_n);

Image = double(imread('McM01.tif'));%读取原图的正确姿势！

 %[Image,pp]=imread('McM02.tif');%读取原图信息
%Image=im2double(Image);%将原图信息中uint8转换成double

R = Image(:, :, 1);%原图RGB中R通道，下同
G = Image(:, :, 2);
B = Image(:, :, 3);

R_noise = u_n(:, :, 1);%噪声图RGB中R通道，下同
G_noise = u_n(:, :, 2);
B_noise = u_n(:, :, 3);

%以下计算噪声图相比原图的PSNR值，采用直接公式计算以及sample code中psnr函数计算两种方式

PSNR_before_Denoising = psnr(Image, u_n);

%PSNR_before_Denoising_R = 20*log10(255/sqrt(mean((R_noise(:)-R(:)).^2)));
%PSNR_before_Denoising_G = 20*log10(255/sqrt(mean((G_noise(:)-G(:)).^2)));
%PSNR_before_Denoising_B = 20*log10(255/sqrt(mean((B_noise(:)-B(:)).^2)));

%PSNR_before_Denoising_R = psnr(R, R_noise);
%PSNR_before_Denoising_G = psnr(G, G_noise);
%PSNR_before_Denoising_B = psnr(B, B_noise);


tic%计时
%KSVD降噪，采用噪点图片信息生成字典
[Image_Out_R,~] = Denoise_Image_KSVD(R_noise, K);
[Image_Out_G,~] = Denoise_Image_KSVD(G_noise, K);
[Image_Out_B,output] = Denoise_Image_KSVD(B_noise, K);

Image_Out_RGB = cat(3, Image_Out_R, Image_Out_G, Image_Out_B);

%以下计算字典学习后图像相比原图的PSNR值，采用直接公式计算以及sample code中psnr函数计算两种方式

PSNR_after_Denoising = psnr(Image, Image_Out_RGB);

%PSNR_after_Denoising_R = 10*log10(255*255/(mean((Image_Out_R(:)-R(:)).^2)));
%PSNR_after_Denoising_G = 20*log10(255/sqrt(mean((Image_Out_G(:)-G(:)).^2)));
%PSNR_after_Denoising_B = 20*log10(255/sqrt(mean((Image_Out_B(:)-B(:)).^2)));

%PSNR_after_Denoising_R = psnr(R, Image_Out_R);
%PSNR_after_Denoising_G = psnr(G, Image_Out_G);
%PSNR_after_Denoising_B = psnr(B, Image_Out_B);

%需要补充PSNR的具体数值
%可以用subplot将原图、噪点图和字典学习后获得的图像绘制在一起
imgsc(Image_Out_RGB);
title('RGB PSNR = ');
toc