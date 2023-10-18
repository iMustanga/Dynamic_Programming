clear
clc

% x_ran=3;
% y_ran=7;
x_ran=randi([2,9]);
y_ran=randi([2,9]);
ground=zeros(x_ran,y_ran);
rows=size(ground,1);
cols=size(ground,2);
ground(rows,cols)=Inf;

n=1:cols;
m=1:rows;

for i=1:rows
   for j= 1:cols
       if i==1 || j==1
           ground(i,j)=1;
       else
           ground(i,j)=ground(i-1,j)+ground(i,j-1);
       end
       
   end
end

if ground(rows,cols)~=Inf
    disp('·�����������ɣ�')
    disp(['������Ϊһ�� ',num2str(x_ran),' x ',num2str(y_ran),' �ľ���'])
    disp('������·������Ϊ��')
    disp(ground(rows,cols))
    disp('���������·����Ϊ��')
    disp(ground)
end