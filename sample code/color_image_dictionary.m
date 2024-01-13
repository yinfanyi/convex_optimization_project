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
   
% 将图像的rgb分别分为三个图像，储存在临时文件夹里

% 读取临时文件夹中的文件
for i=1:3
    [D, ek] = gray_image_dictionary(varargin);
    D_all(:, :, i) = D;
    ek_all(:, i) = ek;
end

end

