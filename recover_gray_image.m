function [u_denoise] = recover_gray_image(varargin)
Parser = inputParser;
addParameter(Parser, 'noise_image', zeros(500, 500));
addOptional(Parser, 'size_of_atom', [1, 1] * 9);
addOptional(Parser, 'num_of_atom', 81);
addOptional(Parser, 'overlapping_size', 4);
addOptional(Parser, 'tolerance', 1e-4);
addOptional(Parser, 'max_iterations', 1e4);
addOptional(Parser, 'dictionary', zeros(81,81));
addOptional(Parser, 'regularization_parameter', 100)
addParameter(Parser, 'sparse_coefficient', zeros(81,15376))
addParameter(Parser, 'mP_noise', zeros(15376, 1))
parse(Parser, varargin{:});

noise_image = Parser.Results.noise_image;
p = Parser.Results.size_of_atom;
gap = Parser.Results.overlapping_size;
nA = Parser.Results.num_of_atom;
tol = Parser.Results.tolerance;
maxits = Parser.Results.max_iterations;
D = Parser.Results.dictionary;
mu = Parser.Results.regularization_parameter;
C = Parser.Results.sparse_coefficient;
mP_noise = Parser.Results.mP_noise;



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

c = 1;
for j = 1:numel(J)
    range_j = J(j): J(j)+p(2) -1;
    for i=1:numel(I)
        range_i = I(i): I(i)+p(1)-1;

        lP(c, :) = [c, range_i(1), range_i(end), range_j(1), range_j(end)];

        patch = noise_image(range_i, range_j);
        patch_ = patch(:) - 1* mean(patch(:)); % remove the DC part

        mP(c) = mean(patch(:)); % mean value of each patch is recorded for reconstruction
        X(:, c) = patch_; 

        c = c + 1;
    end
end

y = D * C;
Y = 0 * noise_image;
W = 0 * noise_image;

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

