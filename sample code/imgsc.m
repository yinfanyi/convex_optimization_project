function imgsc(u)
% 处理的u都是double形式，但显示的u需要uint8形式
if size(u,3)>1
    u = uint8(u);
end

if ~isreal(u)
    u = abs(u);
end

% figure,
imagesc(u);
colormap(gray);
% colormap(jet);
axis('image');
axis off;

colorbar;
