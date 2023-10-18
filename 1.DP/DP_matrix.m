clear
clc

% �������� ��ѡ��������
% DP_MA=[5  13 13 2  6  10;
%        8  4  30 1  12 7 ;
%        12 3  7  14 7  3 ;
%        6  15 6  17 13 2 ;
%        3  27 7  12 20 7 ];
DP_MA=randi([1,50],9,9);
   
% ��Ҫ�õ��ľ����ʼ��
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

% ��ӵõ���Сֵ��������ڵ�ǰ����λ�ã�-1������һ��
dips=[-1;0;1];
% �ڶ������Ž��Ƿ�����ж�ֵ��0Ϊ�ޣ�1Ϊ��
multi_set=0;

% ����һ�е�����д���м̾���
for i=1:rows
    DPmid(i,1)=DP_MA(i,1);
end

% �ӵڶ��п�ʼ��ÿһ��Ԫ�ؽ��б���
for j=2:cols
    for i= 1:rows
        m=1;
        % ����ǰ������������ǰһ�����ݽ�����Ӽ���
        for dip=-1:1
            k=i+dip;
            l=j-1;
            % ����һ�л������һ�н������㣬�Զ���������
            if k==0
                k=rows;
            end
            if k>rows
                k=1;
            end
            % ���������ݴ浽resmid����
            resmid(m,1)=DPmid(k,1)+DP_MA(i,j);
            m=m+1;
        end
        % ��resmid�����е���Сֵд�뵽���յĴ��۾���
        % ������õ���Сֵ������һ����һ������ӣ������-1��0��+1������path����
        DPres(i,1)=min(resmid);
        [p,~]=find(resmid==min(resmid));
        x = p(1,1);
        
        path(i,j-1)=dips(x,1);
        path2(i,j-1)=dips(x,1);
        resmid=zeros(3,1);
        % ��֤�Ƿ��еڶ������Ž�
        if size(p,1)>1
            multi_set=1;
            path2(i,j-1)=dips(p(2,1),1);
        end
    end
    DPmid=DPres;
end

% ��ѯ���Ž��������·��
DPres;
[p,~]=find(DPres==min(DPres));
x = p(1,1);
fin_path(cols,1)=x;
fin_path2(cols,1)=x;
n=1:cols-1;
k=x;
for m=n(:,end:-1:1)
    % ����һ�л������һ�н������㣬�Զ���������
    if k<1 
        k=cols;
    end
    if k>cols 
        k=1;
    end
    
    %����������ֵ����ǰһ������ֵ�õ�ÿһ�е�����
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
 % ÿ������ֵ���������
for i=1:cols
    fin_path(i,2)=i;
    fin_path2(i,2)=i;
end
% ����������ѯ�����е�ÿһ������
for j=1:cols
    j;
    we=fin_path(j,1);
    us=fin_path(j,2);
    fin_num(1,j)=DP_MA(fin_path(j,1),fin_path(j,2));
end

% �ڶ������Ž���� ͬ��
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

% ������ӻ����
disp('���۾���Ϊ��')
disp(DP_MA)
disp('��С����Ϊ��')
disp(min(DPres))
disp('�������Ż�������˳��Ϊ��')
disp(fin_num)
disp('����λ��Ϊ��')
disp(fin_path)

%�ڶ������Ž���ӻ����
if multi_set==1
    disp('�������еڶ������Ž⣡')
    disp('�ڶ������Ž���С����Ϊ��')
    disp(min(DPres))
    disp('�ڶ������Ž��������Ż�������˳��Ϊ��')
    disp(fin_num2)
    disp('�ڶ������Ž�����λ��Ϊ��')
    disp(fin_path2)
end

