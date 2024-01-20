function p = psnr(g, f)
% g: ground truth image
% f: noisy/restored image

g = double(g); % in case of data format is unit8,12,16
f = double(f);

if ndims(f)<ndims(g)
    error ('Dimesion of two images don''t match!');
end

max_ = 255;
min_ = 0;

[m,n, ~] = size(g);
N = m*n;

s = size(g,3);
if s==1 % single channel

    mse = norm(f-g, 'fro')^2/N;
    p = 10.*log10((max_-min_)^2/mse);

else % multi-channel

    p = zeros(s,1);
    for i = 1:s
        mse = norm(f(:,:,i)-g(:,:,i), 'fro')^2/N;
        p(i) = 10.*log10((max_-min_)^2/mse);
    end

end
