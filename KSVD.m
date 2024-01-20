function [Dictionary,output] = KSVD(...
    Data,...
    param)
%K-SVD
%输入：
%Data：一个n*N矩阵，包含N个维度为n的向量 
%param：包含算法所需全部参数的一个struct
%       包含：K, ...：需要训练的字典元素个数
%             Num_Iteration,...：迭代次数
%             Error_Flag...：若为0，则使用固定个数的系数来表示每个信号，若为1，则使用任意数量的系数，直到达到给定表示误差。如果是这样，必须指定允许误差为param.Error_Goal
%             Preserve_DCAtom...：若为1，则字典中的第一个原子将被设置为恒常量。
%             Initialization_Method,...：初始化字典的方式
%             Display_Progress, ...：若为1，则显示进度，若param.Error_Flag为0，则显示RMSE，如果 param.Error_Flag==1 则显示表示每个信号所需的系数平均数。
%输出：
%Dictionary：大小为n*(param.K)
%output：包含实时迭代信息的struct
%        包含：Coefficient_Matrix：最终的系数矩阵
%              ratio：若定义了真正的字典，这一参数保存一个长度为param.Num_Iteration 的向量，其中包括每次迭代的检测率
%              Total_Error：每次迭代后的累计误差
%              Num_Coefficient：长度为param.Num_Iteration的向量，表示系数平均数
% =========================================================================

if (~isfield(param,'displayProgress'))
    param.displayProgress = 1;
end

Total_Error(1) = 99999;%累积误差

if (isfield(param,'errorFlag')==0)
    param.Error_Flag = 1;
end

if (isfield(param,'TrueDictionary'))
    Display_Error_with_True_Dictionary = 1;
    Error_Between_Dictionaries = zeros(param.Num_Iteration+1,1);
    ratio = zeros(param.Num_Iteration+1,1);
else
    Display_Error_with_True_Dictionary = 0;
	ratio = 0;
end

if (param.Preserve_DCAtom>0)
    Fixed_Dictionary_Element(1:size(Data,1),1) = 1/sqrt(size(Data,1));
else
    Fixed_Dictionary_Element = [];
end

if (size(Data,2) < param.K)
    disp('Size of data is smaller than the dictionary size. Trivial solution...');
    Dictionary = Data(:,1:size(Data,2));
    return;
elseif (strcmp(param.Initialization_Method,'DataElements'))
    Dictionary(:,1:param.K-param.Preserve_DCAtom) = Data(:,1:param.K-param.Preserve_DCAtom);
elseif (strcmp(param.Initialization_Method,'GivenMatrix'))
    Dictionary(:,1:param.K-param.Preserve_DCAtom) = param.Initial_Dictionary(:,1:param.K-param.Preserve_DCAtom);
end

if (param.Preserve_DCAtom)
    tmp_Mat = Fixed_Dictionary_Element \ Dictionary;
    Dictionary = Dictionary - Fixed_Dictionary_Element*tmp_Mat;
end

%字典归一化
Dictionary = Dictionary*diag(1./sqrt(sum(Dictionary.*Dictionary)));
%Dictionary = Dictionary.*repmat(sign(Dictionary(1,:)),size(Dictionary,1),1); 
Total_Error = zeros(1,param.Num_Iteration);

for Iter_Num = 1:param.Num_Iteration
    if (param.Error_Flag==0)   
        Coefficient_Matrix = OMP([Fixed_Dictionary_Element,Dictionary],Data, param.L);
    else  
        Coefficient_Matrix = OMP_Error([Fixed_Dictionary_Element,Dictionary],Data, param.Error_Goal);
        param.L = 1;
    end
    
    Replaced_Vector_Counter = 0;
	Random_Perm = randperm(size(Dictionary,2));
    for j = Random_Perm
        [Better_Dictionary_Element,Coefficient_Matrix,Added_New_Vector] = I_findBetterDictionaryElement(Data,...    
            [Fixed_Dictionary_Element,Dictionary],j+size(Fixed_Dictionary_Element,2),...
            Coefficient_Matrix,param.L);
        Dictionary(:,j) = Better_Dictionary_Element;
        if (param.Preserve_DCAtom)
            Tmp_Coefficient = Fixed_Dictionary_Element\Better_Dictionary_Element;
            Dictionary(:,j) = Better_Dictionary_Element - Fixed_Dictionary_Element*Tmp_Coefficient;
            Dictionary(:,j) = Dictionary(:,j)./sqrt(Dictionary(:,j)'*Dictionary(:,j));
        end
        Replaced_Vector_Counter = Replaced_Vector_Counter+Added_New_Vector;
    end
    
    if (Iter_Num>1 && param.Display_Progress)
        if (param.Error_Flag==0)
            output.Total_Error(Iter_Num-1) = sqrt(sum(sum((Data-[Fixed_Dictionary_Element,Dictionary]*Coefficient_Matrix).^2))/numel(Data));
            disp(['Iteration   ',num2str(Iter_Num),'   Total error is: ',num2str(output.Total_Error(Iter_Num-1))]);
        else
            output.Num_Coefficient(Iter_Num-1) = length(find(Coefficient_Matrix))/size(Data,2);
            disp(['Iteration   ',num2str(Iter_Num),'   Average number of coefficients: ',num2str(output.Num_Coefficient(Iter_Num-1))]);
        end
    end
    if (Display_Error_with_True_Dictionary )
        [ratio(Iter_Num+1),Error_Between_Dictionaries(Iter_Num+1)] = I_findDistanseBetweenDictionaries(param.True_Dictionary,Dictionary);
        disp(strcat(['Iteration  ', num2str(Iter_Num),' ratio of restored elements: ',num2str(ratio(Iter_Num+1))]));
        output.ratio = ratio;
    end
    
   Dictionary = I_clearDictionary(Dictionary,Coefficient_Matrix(size(Fixed_Dictionary_Element,2)+1:end,:),Data);
       
    if (isfield(param,'waitBarHandle'))
        waitbar(Iter_Num/param.Counter_for_Wait_Bar);
    end
end

output.Coefficient_Matrix = Coefficient_Matrix;
Dictionary = [Fixed_Dictionary_Element,Dictionary];

function [Better_Dictionary_Element,Coefficient_Matrix,New_Vector_Added] = I_findBetterDictionaryElement(Data,Dictionary,j,Coefficient_Matrix,~)
if (isempty(who('numCoefUsed')))
    Num_Coefficient_Used = 1;
end
Relevant_Data_Indices = find(Coefficient_Matrix(j,:)); 
if (length(Relevant_Data_Indices)<1) 
    Error_Mat = Data-Dictionary*Coefficient_Matrix;
    Error_Norm_Vector = sum(Error_Mat.^2);
    [~,i] = max(Error_Norm_Vector);
    Better_Dictionary_Element = Data(:,i);
    Better_Dictionary_Element = Better_Dictionary_Element./sqrt(Better_Dictionary_Element'*Better_Dictionary_Element);
    Better_Dictionary_Element = Better_Dictionary_Element.*sign(Better_Dictionary_Element(1));
    Coefficient_Matrix(j,:) = 0;
    New_Vector_Added = 1;
    return;
end

New_Vector_Added = 0;
Tmp_Coefficient_Matrix = Coefficient_Matrix(:,Relevant_Data_Indices); 
Tmp_Coefficient_Matrix(j,:) = 0;
Errors =(Data(:,Relevant_Data_Indices) - Dictionary*Tmp_Coefficient_Matrix); 
[Better_Dictionary_Element,Singular_Value,Beta_Vector] = svds(Errors,1);
Coefficient_Matrix(j,Relevant_Data_Indices) = Singular_Value*Beta_Vector';

function [ratio,Total_Distances] = I_findDistanseBetweenDictionaries(Original,New)
CatchCounter = 0;
Total_Distances = 0;
for i = 1:size(New,2)
    New(:,i) = sign(New(1,i))*New(:,i);
end
for i = 1:size(Original,2)
    d = sign(Original(1,i))*Original(:,i);
    distances =sum ( (New-repmat(d,1,size(New,2))).^2);
    [~,index] = min(distances);
    Error_of_Element = 1-abs(New(:,index)'*d);
    Total_Distances = Total_Distances+Error_of_Element;
    Catch_Counter = Catch_Counter+(Error_of_Element<0.01);
end
ratio = 100*Catch_Counter/size(Original,2);

function Dictionary = I_clearDictionary(Dictionary,Coeffiecient_Matrix,Data)
T2 = 0.99;
T1 = 3;
K=size(Dictionary,2); 
Error=sum((Data-Dictionary*Coeffiecient_Matrix).^2,1); 
G=Dictionary'*Dictionary; 
G = G-diag(diag(G));
for j=1:1:K
    if max(G(j,:))>T2 || length(find(abs(Coeffiecient_Matrix(j,:))>1e-7))<=T1 
        [~,pos]=max(Error);
        Clear_Dictionary=1;
        Error(pos(1))=0;
        Dictionary(:,j)=Data(:,pos(1))/norm(Data(:,pos(1)));
        G=Dictionary'*Dictionary;
        G = G-diag(diag(G));
    end
end