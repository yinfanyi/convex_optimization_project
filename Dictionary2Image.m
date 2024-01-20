function I = Dictionary2Image(D, Num_Rows, Num_Cols,X,Y,SVF)

Border_Size = 1;
Column_Scan_Flag = 1;
Strech_Each_Vector_Flag = 1;
Show_Image_Flag = 1;

if (isempty(who('X')))
    X = 8;
    Y = 8;
end
if (isempty(who('sortVarFlag')))
    a=1
    SVF = 1;
end

Num_Elements = size(D,2);
if (isempty(who('numRows')))
    Num_Rows = floor(sqrt(Num_Elements));
    Num_Cols = Num_Rows;
end
if (isempty(who('strechEachVecFlag'))) 
    Strech_Each_Vector_Flag = 0;
end
if (isempty(who('showImFlag'))) 
    Show_Image_Flag = 1;
end

Size_for_Each_Image = sqrt(size(D,1))+Border_Size;
I = zeros(Size_for_Each_Image*Num_Rows+Border_Size,Size_for_Each_Image*Num_Cols+Border_Size,3);

I(:,:,1) = 0;
I(:,:,2) = 0;
I(:,:,3) = 1;

if (Strech_Each_Vector_Flag)
    for Counter = 1:size(D,2)
        D(:,Counter) = D(:,Counter)-min(D(:,Counter));
        if (max(D(:,Counter)))
            D(:,Counter) = D(:,Counter)./max(D(:,Counter));
        end
    end
end

if (SVF)
    vars = var(D);
    [~,Indices] = sort(vars');
    Indices = fliplr(Indices);
    D = [D(:,1:SVF-1),D(:,Indices+SVF-1)];
    Signs = sign(D(1,:));
    Signs(find(Signs==0)) = 1;
    D = D.*repmat(Signs,size(D,1),1);
    D = D(:,1:Num_Rows*Num_Cols);
end

Counter=1;
for j = 1:Num_Rows
    for i = 1:Num_Cols
            I(Border_Size+(j-1)*Size_for_Each_Image+1:j*Size_for_Each_Image,Border_Size+(i-1)*Size_for_Each_Image+1:i*Size_for_Each_Image,1)=reshape(D(:,Counter),X,Y);
            I(Border_Size+(j-1)*Size_for_Each_Image+1:j*Size_for_Each_Image,Border_Size+(i-1)*Size_for_Each_Image+1:i*Size_for_Each_Image,2)=reshape(D(:,Counter),X,Y);
            I(Border_Size+(j-1)*Size_for_Each_Image+1:j*Size_for_Each_Image,Border_Size+(i-1)*Size_for_Each_Image+1:i*Size_for_Each_Image,3)=reshape(D(:,Counter),X,Y);
        Counter = Counter+1;
    end
end

if (Show_Image_Flag) 
    I = I-min(min(min(I)));
    I = I./max(max(max(I)));
    imshow(I,[]);
end
