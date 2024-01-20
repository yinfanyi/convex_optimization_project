function [ImageOut,output] = Denoise_Image_KSVD(Image,K,varargin)

% ���룺Image������ͼ���ҶȻ�ͨ����
%       K���ֵ���ԭ�Ӹ���              
%����������Block_Size���㷨���õ�blocks��С������blocks���Ƿ������ֻ�����������л��������ɣ�Ĭ��ֵΪ8
%         Error_Factor��������Ĭ��ֵΪ1.15
%         Max_Blocks_Considered�����������blocks���ֵ��ȡ�����豸�洢��������ܣ�������blocks�������ڴ˲�������blocks֮��ľ��뽫����
%         Sliding_Factor�������blocks֮��ľ��룬Ĭ��ֵΪ1����ͼ��ܴ������ֵ���Զ����ӣ��ò���Խ�󣬼����ٶ�ҲԽ��
%         KSVD_Iterations������ͼ��KSVD�������������˲���С��ͼ����blocks�ĸ������򽫴����п���blocks�����ѡȡ��Ĭ��ֵΪ10
%         Max_Blocks_to_Train_On��ѵ�������blocks������Ĭ��ֵΪ65000����ͼƬ���󣬸�ֵ������С
%         Display_Flag�����˿���״̬Ϊ������ÿ�ε���������ϵ��ƽ��ֵ���Ա����㷨���ȣ�Ĭ��ֵΪ1������
%         WaitBar_On��������Ϊ0��1��������Ϊ1�������нű�ʱ������һ������������ʾ�㷨����
% �����ImageOut��һ����ά���飬���С�������ͼ��һ�£�������ԭ���ͼ��
%       output.D��ѵ���õ����ֵ�

% ���ȴ�����ͼ��ѵ���ֵ�

Reduce_DC = 1;
[N1,N2] = size(Image);
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
        [Row,Col] = ind2sub(size(Image)-Block_Size+1,Selected_Blocks(i));
        Current_Block = Image(Row:Row+Block_Size-1,Col:Col+Block_Size-1);
        Block_Matrix(:,i) = Current_Block(:);
    end
else
    Block_Matrix = im2col(Image,[Block_Size,Block_Size],'sliding');
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

% ʹ�õõ����ֵ��ͼƬ���н���
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
end;

if (WaitBar_On)
    close(Title);
end
ImageOut = (Image+0.034*20*Image_Out)./(1+0.034*20*Weight);

