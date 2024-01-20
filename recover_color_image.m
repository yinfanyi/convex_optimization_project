function [u_denoise_RGB] = recover_color_image(varargin)
Parser = inputParser;
addParameter(Parser, 'noise_image_RGB', zeros(500, 500, 3));
addOptional(Parser, 'size_of_atom', [1, 1] * 9);
addOptional(Parser, 'num_of_atom', 81);
addOptional(Parser, 'overlapping_size', 4);
addOptional(Parser, 'tolerance', 1e-4);
addOptional(Parser, 'max_iterations', 1e4);
addParameter(Parser, 'dictionary_RGB', zeros(81,81,3));
addOptional(Parser, 'regularization_parameter', 100)
addParameter(Parser, 'sparse_coefficient_RGB', zeros(81,15376, 3))
addParameter(Parser, 'mP_noise_RGB', zeros(15376, 3))
parse(Parser, varargin{:});

noise_image_RGB = Parser.Results.noise_image_RGB;
p = Parser.Results.size_of_atom;
gap = Parser.Results.overlapping_size;
nA = Parser.Results.num_of_atom;
tol = Parser.Results.tolerance;
maxits = Parser.Results.max_iterations;
D_RGB = Parser.Results.dictionary_RGB;
mu = Parser.Results.regularization_parameter;
C_RGB = Parser.Results.sparse_coefficient_RGB;
mP_noise_RGB = Parser.Results.mP_noise_RGB;


u_denoise_R = recover_gray_image('noise_image', noise_image_RGB(:,:,1),...
    'size_of_atom',p,'num_of_atom',nA,'overlapping_size',gap, ...
    'tolerance',tol, 'max_iterations',maxits, ...
    'dictionary',D_RGB(:,:,1),'regularization_parameter', mu,...
    'sparse_coefficient',C_RGB(:,:,1),'mP_noise', mP_noise_RGB(:,1));

u_denoise_G = recover_gray_image('noise_image', noise_image_RGB(:,:,2),...
    'size_of_atom',p,'num_of_atom',nA,'overlapping_size',gap, ...
    'tolerance',tol, 'max_iterations',maxits, ...
    'dictionary',D_RGB(:,:,2),'regularization_parameter', mu,...
    'sparse_coefficient',C_RGB(:,:,2),'mP_noise', mP_noise_RGB(:,2));

u_denoise_B = recover_gray_image('noise_image', noise_image_RGB(:,:,3),...
    'size_of_atom',p,'num_of_atom',nA,'overlapping_size',gap, ...
    'tolerance',tol, 'max_iterations',maxits, ...
    'dictionary',D_RGB(:,:,3),'regularization_parameter', mu,...
    'sparse_coefficient',C_RGB(:,:,3),'mP_noise', mP_noise_RGB(:,3));

u_denoise_RGB = cat(3, u_denoise_R, u_denoise_G, u_denoise_B);
end

