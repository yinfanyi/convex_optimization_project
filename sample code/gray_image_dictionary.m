function [D, ek] = gray_image_dictionary(varargin)
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

% 读取图像
u = double(imread(filename));
u0 = u;
n = size(u, 1:2);
figure(001); clf;
imgsc(u0);

% 字典学习步骤的基本设置和准备
np = prod(p); % number of pixels in each patch
I = [1:gap:n(1)-p(1), n(1)-p(1)+1];
J = [1:gap:n(2)-p(2), n(2)-p(2)+1];
nP = numel(I)* numel(J); % total number of patches
lP = zeros(nP, 5); % location of each patch

% 提取图像块
X = zeros(np, nP);
mP = zeros(nP, 1);

c = 1;
for j = 1:numel(J)
    range_j = J(j): J(j)+p(2) -1;
    for i=1:numel(I)
        range_i = I(i): I(i)+p(1)-1;

        lP(c, :) = [c, range_i(1), range_i(end), range_j(1), range_j(end)];

        patch = u(range_i, range_j);
        patch_ = patch(:) - 1* mean(patch(:)); % remove the DC part

        mP(c) = mean(patch(:)); % mean value of each patch is recorded for reconstruction
        X(:, c) = patch_; 

        c = c + 1;
    end
end

% TODO: 适应彩色图片
D0 = eye(np, nA) /1;
A0 = repmat(D0(:, 1:np)\ X, nA/np, 1);



D = D0;
A = A0;

ek = zeros(maxits, 1); %% record residual error

tic;
its = 1;
while its<=maxits-1
    D_old = D;
    % A_old = A;

    % gradient descent on D
    L_D = norm(A*A') + .1; % estimated Lipschitz constant
    % gamma_D = 1.9/L_D; % the step-size

    grad_D = (-X + D*A)* (A'); % gradient respect to D
    D = D - (1.9/L_D)* grad_D; 

    % enforcing unit length constraint on to the atoms
    for j=1:nA
        D(:, j) = D(:, j) / norm(D(:,j));
    end

    % gradient descent on A
    L_A = norm(D*D') + .1; % estimated Lipschitz constant
    gamma_A = 1.9/L_A; % the step-size

    grad_A = (D')* (-X + D*A); % gradient respect to A
    w = A - gamma_A* grad_A;

    % this is for dealing the l1-norm, called soft-thresholding shrinkage
    A = sign(w) .* max(abs(w)-lambda*gamma_A, 0);

    ek(its) = norm(D(:)-D_old(:)) /norm(D(:));

    if mod(its, 5e2)==0; lambda = max(lambda/1.5, 1/2); end

    if ek(its) < tol || ek(its) > 1e10 
        ek = ek(1:its-1);
        break; 
    end

    if mod(its,5e2)==0
        fprintf('%06d, time: %05.2fs, rel_error: %05.2e, obj_value: %05.2e...\n',...
            its, toc, ek(its), norm(X-D*A, 'fro')^2);

        imgsc(D);
        pause(.01);
    end

    its = its + 1;
end

end

