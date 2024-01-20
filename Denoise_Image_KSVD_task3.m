function [ImageOut,output] = Denoise_Image_KSVD_task3(Image,K, image_original)

% 输入：Image：噪声图（灰度或单通道）
%       K：字典中原子个数              
%其他参数：Block_Size：算法采用的blocks大小，所有blocks都是方阵，因此只需给出矩阵的行或列数即可，默认值为8
%         Error_Factor：允许误差，默认值为1.15
%         Max_Blocks_Considered：允许被处理的blocks最大值，取决于设备存储情况和性能，若可用blocks数量大于此参数，则blocks之间的距离将增加
%         Sliding_Factor：处理的blocks之间的距离，默认值为1，若图像很大，则此数值将自动增加，该参数越大，计算速度也越快
%         KSVD_Iterations：噪声图的KSVD迭代次数，若此参数小于图像中blocks的个数，则将从所有可用blocks中随机选取，默认值为10
%         Max_Blocks_to_Train_On：训练的最大blocks数量，默认值为65000，若图片过大，该值可能略小
%         Display_Flag：若此开关状态为开，则每次迭代后会给出系数平均值，以表征算法进度，默认值为1（开）
%         WaitBar_On：可设置为0或1，若设置为1，则运行脚本时将出现一个进度条，显示算法进度
% 输出：ImageOut：一个二维数组，其大小与输入的图像一致，包含复原后的图像
%       output.D：训练得到的字典

% 首先从原图中训练字典

Reduce_DC = 1;
[N1,N2] = size(image_original);
WaitBar_On = 1;
KSVD_Iterations = 10;
Error_Factor = 1.15;
Max_Blocks_Considered = 260000;
Sliding_Factor = 1;
Block_Size = 8;
Max_Blocks_to_Train_On = 65000;
Display_Flag = 1;

if(prod([N1,N2]-Block_Size+1)> Max_Blocks_to_Train_On)
    Rand_Permutation =  randperm(prod([N1,N2]-Block_Size+1));
    Selected_Blocks = Rand_Permutation(1:Max_Blocks_to_Train_On);

    Block_Matrix = zeros(Block_Size^2,Max_Blocks_to_Train_On);
    for i = 1:Max_Blocks_to_Train_On
        [Row,Col] = ind2sub(size(image_original)-Block_Size+1,Selected_Blocks(i));
        Current_Block = image_original(Row:Row+Block_Size-1,Col:Col+Block_Size-1);
        Block_Matrix(:,i) = Current_Block(:);
    end
else
    Block_Matrix = im2col(image_original,[Block_Size,Block_Size],'sliding');
end

param.K = K;
param.Num_Iteration = KSVD_Iterations ;
param.Error_Flag = 1;
param.Error_Goal = 20*Error_Factor;
param.Preserve_DCAtom = 0;

Alpha=ceil(sqrt(K));
Beta=zeros(Block_Size,Alpha);
for k=0:1:Alpha-1
    Theta=cos([0:1:Block_Size-1]'*k*pi/Alpha);
    if k>0
        Theta=Theta-mean(Theta); 
    end
    Beta(:,k+1)=Theta/norm(Theta);
end
Beta=kron(Beta,Beta);

param.Initial_Dictionary = Beta(:,1:param.K );
param.Initialization_Method =  'GivenMatrix';

if (Reduce_DC)
    Vector_of_Means = mean(Block_Matrix);
    Block_Matrix = Block_Matrix-ones(size(Block_Matrix,1),1)*Vector_of_Means;
end

if (WaitBar_On)
    %param.Num_Iteration = param.Num_Iteration + 1;
    Counter_for_WaitBar = param.Num_Iteration+1;
    Title = waitbar(0,'Denoising In Process ...');
    param.WaitBar_Handle = Title;
    param.Counter_for_WaitBar = Counter_for_WaitBar;
end

param.Display_Progress = Display_Flag;
[Dictionary,output] = KSVD(Block_Matrix,param);
output.D = Dictionary;

if (Display_Flag)
    disp('finished Trainning dictionary');
end

% 使用得到的字典对图片进行降噪
Error_T = 20*Error_Factor;

while (prod(floor((size(Image)-Block_Size)/Sliding_Factor)+1)>Max_Blocks_Considered)
    Sliding_Factor = Sliding_Factor+1;
end

[Blocks,Idx] = Im2Col(Image,[Block_Size,Block_Size],Sliding_Factor);

if (WaitBar_On)
    NewCounter_for_WaitBar = (param.Num_Iteration+1)*size(Blocks,2);
end

for p = 1:30000:size(Blocks,2)
    if (WaitBar_On)
        waitbar(((param.Num_Iteration*size(Blocks,2))+p)/NewCounter_for_WaitBar);
    end
    Jump_Size = min(p+30000-1,size(Blocks,2));
    if (Reduce_DC)
        Vector_of_Means = mean(Blocks(:,p:Jump_Size));
        Blocks(:,p:Jump_Size) = Blocks(:,p:Jump_Size) - repmat(Vector_of_Means,size(Blocks,1),1);
    end
    
    Coeffiecients = OMP_Error(Dictionary,Blocks(:,p:Jump_Size),Error_T);
    if (Reduce_DC)
        Blocks(:,p:Jump_Size)= Dictionary*Coeffiecients + ones(size(Blocks,1),1) * Vector_of_Means;
    else
        Blocks(:,p:Jump_Size)= Dictionary*Coeffiecients;
    end
end

Count = 1;
Weight = zeros(N1,N2);
Image_Out = zeros(N1,N2);
[Rows,Cols] = ind2sub(size(Image)-Block_Size+1,Idx);
for i  = 1:length(Cols)
    Col = Cols(i); Row = Rows(i);        
    Block =reshape(Blocks(:,Count),[Block_Size,Block_Size]);
    Image_Out(Row:Row+Block_Size-1,Col:Col+Block_Size-1)=Image_Out(Row:Row+Block_Size-1,Col:Col+Block_Size-1)+Block;
    Weight(Row:Row+Block_Size-1,Col:Col+Block_Size-1)=Weight(Row:Row+Block_Size-1,Col:Col+Block_Size-1)+ones(Block_Size);
    Count = Count+1;
end

if (WaitBar_On)
    close(Title);
end
ImageOut = (Image+0.034*20*Image_Out)./(1+0.034*20*Weight);

