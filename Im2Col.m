function [Blocks,Idx] = Im2Col(I,Block_Size,Sliding_Distance);

if (Sliding_Distance==1)
    Blocks = im2col(I,Block_Size,'sliding');
    Idx = [1:size(Blocks,2)];
    return
end

Idx_Mat = zeros(size(I)-Block_Size+1);
Idx_Mat([[1:Sliding_Distance:end-1],end],[[1:Sliding_Distance:end-1],end]) = 1;
Idx = find(Idx_Mat);
[Rows,Cols] = ind2sub(size(Idx_Mat),Idx);
Blocks = zeros(prod(Block_Size),length(Idx));
for i = 1:length(Idx)
    Current_Block = I(Rows(i):Rows(i)+Block_Size(1)-1,Cols(i):Cols(i)+Block_Size(2)-1);
    Blocks(:,i) = Current_Block(:);
end
