function [D_all, ek_all] = color_image_dictionary(varargin)
Parser = inputParser;
addOptional(Parser, 'filename', 'lena_512.png');
addOptional(Parser, 'size_of_atom', [1, 1] * 9);
addOptional(Parser, 'num_of_atom', 81);
addOptional(Parser, 'overlapping_size', 4);
addOptional(Parser, 'regularization_parameter', 50);
addOptional(Parser, 'tolerance', 1e-4);
addOptional(Parser, 'max_iterations', 1e4);
parse(Parser, varargin{:});

filename = Parser.Results.filename;
p = Parser.Results.size_of_atom;
nA = Parser.Results.num_of_atom;
gap = Parser.Results.overlapping_size;
lambda = Parser.Results.regularization_parameter;
tol = Parser.Results.tolerance;
maxits = Parser.Results.max_iterations;

np = prod(p); % number of pixels in each patch
D_all = zeros(np, nA, 3);
ek_all = zeros(maxits, 3); %% record residual error

% 读取图像
u = double(imread(filename));

% 将图像的rgb分别分为三个图像，储存在临时文件夹里
r = u(:,:,1); 
g = u(:,:,2); 
b = u(:,:,3); 
filename_temp = cell(3,1);
filename_temp{1} = 'temp/temp_r.png';
filename_temp{2} = 'temp/temp_g.png';
filename_temp{3} = 'temp/temp_b.png';
imwrite(uint8(r), filename_temp{1}); 
imwrite(uint8(g), filename_temp{2}); 
imwrite(uint8(b), filename_temp{3});

for i=1:3
    [D, ek] = gray_image_dictionary('filename', filename_temp{i}, ...
        'size_of_atom', p, 'num_of_atom', nA, 'overlapping_size', gap, ...
        'regularization_parameter', lambda, 'tolerance', tol, ...
        'max_iterations', maxits);
    D_all(:, :, i) = D;
    ek_all(1:size(ek, 1), i) = ek;
    delete(filename_temp{i}) % 删除临时文件夹的文件
end


end

