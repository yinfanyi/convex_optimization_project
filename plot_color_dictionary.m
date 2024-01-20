
function plot_color_dictionary(D_RGB, varargin)
%PLOT_COLOR_DICTIONARY 绘制彩色图片的字典
%   此处显示详细说明

Parser = inputParser;
addOptional(Parser, 'size_of_atom', [1, 1] * 16);
addOptional(Parser, 'is_save_image_data', false);
addOptional(Parser, 'save_dir', [])
parse(Parser, varargin{:});

p = Parser.Results.size_of_atom;
is_save_image_data = Parser.Results.is_save_image_data;
tmp_dir = Parser.Results.save_dir;

i_size = prod(p) + p(1) - 1;
imD = zeros(i_size, i_size);
title_name = {'red', 'green', 'blue'};
figure; clf;
hold on
for k=1:3
    c = 0;
    for i=1: p(1)
        range_i = (i-1)*(p(1)+1) + 1 : i*(p(1)+1)-1;
        for j = 1: p(2)
            c = c + 1;

            range_j = (j-1)*(p(2)+1) + 1 : j*(p(2)+1)-1;

            atom = D_RGB(:, c, k);
            imD(range_i, range_j) = reshape(atom, p) - min(atom(:));
        end
    end
    subplot(3,1,k);
    imagesc(imD);
    colormap(jet);
    axis off;
    colorbar;
    title(title_name{1,k});
end
hold off
if is_save_image_data
    saveas(gcf, tmp_dir);
end
end

