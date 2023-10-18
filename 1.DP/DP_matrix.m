clear
clc

% 矩阵设置 可选自设或随机
% DP_MA=[5  13 13 2  6  10;
%        8  4  30 1  12 7 ;
%        12 3  7  14 7  3 ;
%        6  15 6  17 13 2 ;
%        3  27 7  12 20 7 ];
DP_MA=randi([1,50],9,9);
   
% 需要用到的矩阵初始化
rows=size(DP_MA,1);
cols=size(DP_MA,2);

DPmid=zeros(rows,1);
DPres=zeros(rows,1);
path =zeros(rows,cols);
path2=zeros(rows,cols);
fin_path=zeros(cols,2);
fin_path2=zeros(cols,2);
fin_num=zeros(1,cols);
fin_num2=zeros(1,cols);
pathmid=zeros(3,1);
resmid=zeros(3,1);

% 相加得到最小值的数相对于当前数的位置，-1代表上一行
dips=[-1;0;1];
% 第二个最优解是否存在判断值，0为无，1为有
multi_set=0;

% 将第一列的数据写入中继矩阵
for i=1:rows
    DPmid(i,1)=DP_MA(i,1);
end

% 从第二列开始对每一个元素进行遍历
for j=2:cols
    for i= 1:rows
        m=1;
        % 将当前遍历到的数与前一列数据进行相加计算
        for dip=-1:1
            k=i+dip;
            l=j-1;
            % 若第一行或是最后一行进行运算，自动更改索引
            if k==0
                k=rows;
            end
            if k>rows
                k=1;
            end
            % 将计算结果暂存到resmid矩阵
            resmid(m,1)=DPmid(k,1)+DP_MA(i,j);
            m=m+1;
        end
        % 将resmid矩阵中的最小值写入到最终的代价矩阵
        % 并算出得到最小值是与上一列哪一个数相加，结果有-1、0、+1，存在path矩阵
        DPres(i,1)=min(resmid);
        [p,~]=find(resmid==min(resmid));
        x = p(1,1);
        
        path(i,j-1)=dips(x,1);
        path2(i,j-1)=dips(x,1);
        resmid=zeros(3,1);
        % 查证是否有第二个最优解
        if size(p,1)>1
            multi_set=1;
            path2(i,j-1)=dips(p(2,1),1);
        end
    end
    DPmid=DPres;
end

% 查询最优解所在相对路径
DPres;
[p,~]=find(DPres==min(DPres));
x = p(1,1);
fin_path(cols,1)=x;
fin_path2(cols,1)=x;
n=1:cols-1;
k=x;
for m=n(:,end:-1:1)
    % 若第一行或是最后一行进行运算，自动更改索引
    if k<1 
        k=cols;
    end
    if k>cols 
        k=1;
    end
    
    %将该列索引值加上前一列索引值得到每一列的索引
    path_idx=fin_path(m+1,1)+ path(k,m);
    if path_idx<1 
        path_idx=cols;
    end
    if path_idx>cols 
        path_idx=1;
    end
    fin_path(m,1)=path_idx;
    
    k=k+ path(k,m);
end
 % 每个索引值后加入列数
for i=1:cols
    fin_path(i,2)=i;
    fin_path2(i,2)=i;
end
% 根据索引查询过程中的每一个数字
for j=1:cols
    j;
    we=fin_path(j,1);
    us=fin_path(j,2);
    fin_num(1,j)=DP_MA(fin_path(j,1),fin_path(j,2));
end

% 第二个最优解求解 同上
if multi_set==1
    n=1:cols-1;
    k=x;
    for m=n(:,end:-1:1)
        if k<1 
            k=cols;
        end
        if k>cols 
            k=1;
        end
        
        path_idx=fin_path2(m+1,1)+ path2(k,m);
        if path_idx<1 
            path_idx=cols;
        end
        if path_idx>cols 
            path_idx=1;
        end
            fin_path2(m,1)=path_idx;
            k=k+ path2(k,m);
    end
    for j=1:cols
    fin_num2(1,j)=DP_MA(fin_path2(j,1),fin_path2(j,2));
    end
end

% 结果可视化输出
disp('代价矩阵为：')
disp(DP_MA)
disp('最小代价为：')
disp(min(DPres))
disp('最终最优化后数字顺序为：')
disp(fin_num)
disp('数字位置为：')
disp(fin_path)

%第二个最优解可视化输出
if multi_set==1
    disp('该问题有第二个最优解！')
    disp('第二个最优解最小代价为：')
    disp(min(DPres))
    disp('第二个最优解最终最优化后数字顺序为：')
    disp(fin_num2)
    disp('第二个最优解数字位置为：')
    disp(fin_path2)
end

