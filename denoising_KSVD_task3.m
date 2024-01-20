b=8; % block size
R=4; % redundancy factor  冗余因素
K=R*b^2; % number of atoms in the dictionary

sigma = 30;
%%
%for number = 1:18
%    if number < 10
%        number_str = ['0', num2str(number)];
%    else
%        number_str = num2str(number);
%    end

    load McM06_noise.mat u_n %加载噪声图像
    %tmp_dir = ['McM images/McM',number_str,'_noise.mat'];
    %load(tmp_dir, 'u_n');

    %以下等价于sample code中的imgsc函数
    %u_n = uint8(u_n);
    %imagesc(u_n);

    Image = double(imread('McM06.tif'));%读取原图的正确姿势！
    %filename = ['McM images/McM',number_str,'.tif']; % 需要进入sample code工作目录
    %Image = double(imread(filename));

    %[Image,pp]=imread('McM02.tif');%读取原图信息
    %Image=im2double(Image);%将原图信息中uint8转换成double

    %R = Image(:, :, 1);%原图RGB中R通道，下同
    %G = Image(:, :, 2);
    %B = Image(:, :, 3);

    %R_noise = u_n(:, :, 1);%噪声图RGB中R通道，下同
    %G_noise = u_n(:, :, 2);
    %B_noise = u_n(:, :, 3);

    %以下计算噪声图相比原图的PSNR值，采用直接公式计算以及sample code中psnr函数计算两种方式

    PSNR_before_Denoising = psnr(Image, u_n);

    %PSNR_before_Denoising_R = 20*log10(255/sqrt(mean((R_noise(:)-R(:)).^2)));
    %PSNR_before_Denoising_G = 20*log10(255/sqrt(mean((G_noise(:)-G(:)).^2)));
    %PSNR_before_Denoising_B = 20*log10(255/sqrt(mean((B_noise(:)-B(:)).^2)));

    %PSNR_before_Denoising_R = psnr(R, R_noise);
    %PSNR_before_Denoising_G = psnr(G, G_noise);
    %PSNR_before_Denoising_B = psnr(B, B_noise);


    tic%计时
    %KSVD降噪，采用原图片生成字典
    [Image_Out_R,output_R] = Denoise_Image_KSVD_task3(R_noise, K, R);
    [Image_Out_G,output_G] = Denoise_Image_KSVD_task3(G_noise, K, G);
    [Image_Out_B,output_B] = Denoise_Image_KSVD_task3(B_noise, K, B);
    toc
    Image_Out_RGB = cat(3, Image_Out_R, Image_Out_G, Image_Out_B);
    D_R = output_R.D;
    D_G = output_G.D;
    D_B = output_B.D;
    %%
    %以下计算字典学习后图像相比原图的PSNR值，采用直接公式计算以及sample code中psnr函数计算两种方式
    %figure(1);clf;
    %subplot(1,3,1);
    %Dictionary2Image(D_R);
    %title('red');
    %subplot(1,3,2);
    %Dictionary2Image(D_G);
    %title('green');
    %subplot(1,3,3);
    %Dictionary2Image(D_B);
    %title('blue');

    PSNR_after_Denoising = psnr(Image, Image_Out_RGB);

    q = 0.1 .* PSNR_after_Denoising;
    Q = 10 .^ q;
    mse = 255 * 255 ./ (Q);
    MSE = mean(mse);
    PSNR_after_Denoising_Average = 10*log10(255*255/MSE);

    DR = kron(D_R, [1; 1; 1; 1]);
    DG = kron(D_G, [1; 1; 1; 1]);
    DB = kron(D_B, [1; 1; 1; 1]);
    D_RGB = cat(3, DR, DG, DB);
    plot_color_dictionary(D_RGB);
    tmp_dir = ['C:\Users\nonarne\Desktop\code\data\task2\task2_2\','dictionary', number_str, '.png'];
    saveas(gcf, tmp_dir);
    %PSNR_after_Denoising_R = 10*log10(255*255/(mean((Image_Out_R(:)-R(:)).^2)));
    %PSNR_after_Denoising_G = 20*log10(255/sqrt(mean((Image_Out_G(:)-G(:)).^2)));
    %PSNR_after_Denoising_B = 20*log10(255/sqrt(mean((Image_Out_B(:)-B(:)).^2)));

    %PSNR_after_Denoising_R = psnr(R, Image_Out_R);
    %PSNR_after_Denoising_G = psnr(G, Image_Out_G);
    %PSNR_after_Denoising_B = psnr(B, Image_Out_B);

    %需要补充PSNR的具体数值
    %可以用subplot将原图、噪点图和字典学习后获得的图像绘制在一起
    %figure(2);clf;
    %imgsc(Image_Out_RGB);
    %title('RGB PSNR = ',sprintf('%.2f ', PSNR_after_Denoising, PSNR_after_Denoising_Average));
    %tmp_dir = ['C:\Users\nonarne\Desktop\code\data\task3_ksvd\','recover_', number_str, '.png'];
    %saveas(gcf, tmp_dir);
    
end