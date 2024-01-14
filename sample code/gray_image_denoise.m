function [C, ek_lasso, mP_noise] = gray_image_denoise(varargin)
Parser = inputParser;
addParameter(Parser, 'noise_image', zeros(500, 500));
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

% 将noise_image分为RBG通道
% 三个通道的各自生成patch
% 现在只有R通道
u_noise = noise_image;
X_noise = zeros(np, nP);
mP_noise = zeros(nP, 1);

c = 1;
for j = 1:numel(J)
    range_j = J(j): J(j)+p(2) -1;
    for i=1:numel(I)
        range_i = I(i): I(i)+p(1)-1;

        lP(c, :) = [c, range_i(1), range_i(end), range_j(1), range_j(end)];

        patch = u_noise(range_i, range_j);
        patch_ = patch(:) - 1* mean(patch(:)); % remove the DC part

        mP_noise(c) = mean(patch(:));
        % mean value of each patch is recorded for reconstruction
        X_noise(:, c) = patch_; % - mean(patch(:));

        c = c + 1;
    end
end

% the different in meanvalue is actually not large.
fprintf(['Max difference between the values of patch mean value: ' ...
    '%.2f...\n'], max(abs(mP-mP_noise)));

% 使用LASSO方法降噪
L = norm(D)^2; % Lipschitz constant of the gradient of 0.5*||X-DC||_F^2
gamma = 1/ L;
C = A; % initialize C with A above

ek_lasso = zeros(maxits, 1);

its = 1;
while its<=maxits
    C_old = C;

    % gradient descent on the quadratic part
    grad = D'*(D*C - X_noise);
    w = C - gamma* grad;

    % this is for dealing the l0-norm, called hard-thresholding shrinkage
    % which is similar to OMP (orthogonal matching pursuit)
    C = wthresh(w, 'h', gamma*mu);

    res = C - C_old;
    ek_lasso(its) = norm(res, 'fro');

    %%% stopping criterion is not used here...
    if ek_lasso(its) < tol || ek_lasso(its) > 1e10
        if its>1e1-1
            ek_lasso = ek_lasso(1:its-1);
            break;
        end
    end

    if mod(its,5e2)==0
        fprintf('%06d, rel_error: %05.2e, obj_value: %07.4e, l1_norm: %05.2e...\n',...
            its, ek_lasso(its), norm(X-D*C, 'fro')^2/2+mu*norm(C(:), 1), mu*norm(C(:), 1));
    end
    its = its + 1;
end
end

