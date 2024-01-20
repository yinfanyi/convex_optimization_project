b=8; % block size
R=4; % redundancy factor  ��������
K=R*b^2; % number of atoms in the dictionary

sigma = 30; 

load McM01_noise.mat u_n %��������ͼ��

%���µȼ���sample code�е�imgsc����
%u_n = uint8(u_n);
%imagesc(u_n);

Image = double(imread('McM01.tif'));%��ȡԭͼ����ȷ���ƣ�

 %[Image,pp]=imread('McM02.tif');%��ȡԭͼ��Ϣ
%Image=im2double(Image);%��ԭͼ��Ϣ��uint8ת����double

R = Image(:, :, 1);%ԭͼRGB��Rͨ������ͬ
G = Image(:, :, 2);
B = Image(:, :, 3);

R_noise = u_n(:, :, 1);%����ͼRGB��Rͨ������ͬ
G_noise = u_n(:, :, 2);
B_noise = u_n(:, :, 3);

%���¼�������ͼ���ԭͼ��PSNRֵ������ֱ�ӹ�ʽ�����Լ�sample code��psnr�����������ַ�ʽ

PSNR_before_Denoising = psnr(Image, u_n);

%PSNR_before_Denoising_R = 20*log10(255/sqrt(mean((R_noise(:)-R(:)).^2)));
%PSNR_before_Denoising_G = 20*log10(255/sqrt(mean((G_noise(:)-G(:)).^2)));
%PSNR_before_Denoising_B = 20*log10(255/sqrt(mean((B_noise(:)-B(:)).^2)));

%PSNR_before_Denoising_R = psnr(R, R_noise);
%PSNR_before_Denoising_G = psnr(G, G_noise);
%PSNR_before_Denoising_B = psnr(B, B_noise);


tic%��ʱ
%KSVD���룬�������ͼƬ��Ϣ�����ֵ�
[Image_Out_R,~] = Denoise_Image_KSVD(R_noise, K);
[Image_Out_G,~] = Denoise_Image_KSVD(G_noise, K);
[Image_Out_B,output] = Denoise_Image_KSVD(B_noise, K);

Image_Out_RGB = cat(3, Image_Out_R, Image_Out_G, Image_Out_B);

%���¼����ֵ�ѧϰ��ͼ�����ԭͼ��PSNRֵ������ֱ�ӹ�ʽ�����Լ�sample code��psnr�����������ַ�ʽ

PSNR_after_Denoising = psnr(Image, Image_Out_RGB);

%PSNR_after_Denoising_R = 10*log10(255*255/(mean((Image_Out_R(:)-R(:)).^2)));
%PSNR_after_Denoising_G = 20*log10(255/sqrt(mean((Image_Out_G(:)-G(:)).^2)));
%PSNR_after_Denoising_B = 20*log10(255/sqrt(mean((Image_Out_B(:)-B(:)).^2)));

%PSNR_after_Denoising_R = psnr(R, Image_Out_R);
%PSNR_after_Denoising_G = psnr(G, Image_Out_G);
%PSNR_after_Denoising_B = psnr(B, Image_Out_B);

%��Ҫ����PSNR�ľ�����ֵ
%������subplot��ԭͼ�����ͼ���ֵ�ѧϰ���õ�ͼ�������һ��
imgsc(Image_Out_RGB);
title('RGB PSNR = ');
toc