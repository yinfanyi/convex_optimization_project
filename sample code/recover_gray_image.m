function [u_denoise] = recover_gray_image(varargin)
Parser = inputParser;
addOptional(Parser, 'noise_image', zeros(500, 500, 3));
addOptional(Parser, 'size_of_atom', [1, 1] * 9);
addOptional(Parser, 'num_of_atom', 81);
addOptional(Parser, 'overlapping_size', 4);
addOptional(Parser, 'tolerance', 1e-4);
addOptional(Parser, 'max_iterations', 1e4);
addOptional(Parser, 'dictionary', zeros(81,81,3));
addOptional(Parser, 'regularization_parameter', 100)
parse(Parser, varargin{:});

noise_image = Parser.Results.noise_image;
p = Parser.Results.size_of_atom;
gap = Parser.Results.overlapping_size;
nA = Parser.Results.num_of_atom;
tol = Parser.Results.tolerance;
maxits = Parser.Results.max_iterations;
D_rgb = Parser.Results.dictionary;
mu = Parser.Results.regularization_parameter;


np = prod(p); % number of pixels in each patch
n = size(noise_image, 1:2);
I = [1:gap:n(1)-p(1), n(1)-p(1)+1];
J = [1:gap:n(2)-p(2), n(2)-p(2)+1];
nP = numel(I)* numel(J); % total number of patches
lP = zeros(nP, 5); % location of each patch
mP = zeros(nP, 1);
X = zeros(np, nP);
D0 = eye(np, nA) /1;
A0 = repmat(D0(:, 1:np)\ X, nA/np, 1);
A = A0;

y = D_McM12(:,:,1)* C;
Y = 0* McM12_noise(:,:,1);
W = 0* McM12_noise(:,:,1);

for i=1:nP
    p_location = lP(i, 2:5);
    patch = y(:, i) + mP_noise(i); % - mean(patch(:));

    range_i = p_location(1):p_location(2);
    range_j = p_location(3):p_location(4);

    Y(range_i, range_j) = Y(range_i, range_j) + reshape(patch, p);
    W(range_i, range_j) = W(range_i, range_j) + 1;
end

u_denoise = Y./W;

end
