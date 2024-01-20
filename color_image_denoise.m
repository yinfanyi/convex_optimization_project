function [C_RGB, ek_lasso_RGB, mP_noise_RGB] = color_image_denoise(varargin)
Parser = inputParser;
addParameter(Parser, 'noise_image', zeros(500, 500, 3));
addOptional(Parser, 'size_of_atom', [1, 1] * 9);
addOptional(Parser, 'num_of_atom', 81);
addOptional(Parser, 'overlapping_size', 4);
addOptional(Parser, 'tolerance', 1e-4);
addOptional(Parser, 'max_iterations', 1e4);
addParameter(Parser, 'dictionary', zeros(81,81,3));
addOptional(Parser, 'regularization_parameter', 100)
parse(Parser, varargin{:});

noise_image = Parser.Results.noise_image;
p = Parser.Results.size_of_atom;
gap = Parser.Results.overlapping_size;
nA = Parser.Results.num_of_atom;
tol = Parser.Results.tolerance;
maxits = Parser.Results.max_iterations;
D = Parser.Results.dictionary;
mu = Parser.Results.regularization_parameter;


[C_R, ek_lasso_R, mP_noise_R] = gray_image_denoise('noise_image', noise_image(:, :, 1), ...
    'size_of_atom', p, 'num_of_atom', nA, 'overlapping_size', gap, ...
    'tolerance', tol, 'max_iterations', maxits, 'dictionary', D(:, :, 1), ...
    'regularization_parameter', mu);
[C_G, ek_lasso_G, mP_noise_G] = gray_image_denoise('noise_image', noise_image(:, :, 2), ...
    'size_of_atom', p, 'num_of_atom', nA, 'overlapping_size', gap, ...
    'tolerance', tol, 'max_iterations', maxits, 'dictionary', D(:, :, 2), ...
    'regularization_parameter', mu);
[C_B, ek_lasso_B, mP_noise_B] = gray_image_denoise('noise_image', noise_image(:, :, 3), ...
    'size_of_atom', p, 'num_of_atom', nA, 'overlapping_size', gap, ...
    'tolerance', tol, 'max_iterations', maxits, 'dictionary', D(:, :, 3), ...
    'regularization_parameter', mu);

C_RGB = zeros(size(C_R,1),size(C_R,2),3);
C_RGB(:,:,1) = C_R;
C_RGB(:,:,2) = C_G;
C_RGB(:,:,3) = C_B;

ek_lasso_RGB = zeros(maxits,3);
ek_lasso_RGB(1:size(ek_lasso_R,1), 1) = ek_lasso_R;
ek_lasso_RGB(1:size(ek_lasso_G, 1),2) = ek_lasso_G;
ek_lasso_RGB(1:size(ek_lasso_B, 1),3) = ek_lasso_B;

mP_noise_RGB = zeros(size(mP_noise_R, 1), 3);
mP_noise_RGB(:,1) = mP_noise_R;
mP_noise_RGB(:,2) = mP_noise_G;
mP_noise_RGB(:,3) = mP_noise_B;

end

