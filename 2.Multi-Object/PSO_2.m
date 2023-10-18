clear all;
close all;
clc;
t1=tic;
N=100;                         %Ⱥ�����Ӹ���
D=2;                           %����ά��
T=200;                         %����������
c1=1.5;                        %ѧϰ����1
c2=1.5;                        %ѧϰ����2
Wmax=0.8;                      %����Ȩ�����ֵ
Wmin=0.4;                      %����Ȩ����Сֵ
Xmax=4;                        %λ�����ֵ
Xmin=-4;                       %λ����Сֵ
Vmax=1;                        %�ٶ����ֵ
Vmin=-1;                       %�ٶ���Сֵ

x=-4:0.02:4;
y=-4:0.02:4;
N=size(x,2);
for i=1:N
    for j=1:N
        z(i,j)=3*cos(x(i)*y(j))+x(i)+y(j)*y(j);
    end
end
mesh(x,y,z)
xlabel('x')
ylabel('y')

%��ʼ������
x=rand(N,D)*(Xmax-Xmin)+Xmin;
v=rand(N,D)*(Vmax-Vmin)+Vmin;
%��ʼ����������λ�ú�����ֵ
p=x;
pbest=ones(1,D);
for i=1:N
    pbest(i)=func2(x(i,:));
end
%��ʼ��ȫ������λ�ú�����ֵ
g=ones(1,D);
gbest=inf;
for i=1:N
    if (pbest(i)<gbest)
        g=p(i,:);
        gbest=pbest(i);
    end
end
gb=ones(1,T);
%���չ�ʽ���ε���ֱ�����㾫�Ȼ��ߵ�������
for i=1:T
    for j=1:N
        %���¸�������λ�ú�����ֵ
        if (func2(x(j,:))<pbest(j))
            p(j,:)=x(j,:);
            pbest(j)=func2(x(j,:));
        end
        %����ȫ������λ�ú�����ֵ
        if (pbest(j)<gbest)
            g=p(j,:);
            gbest=pbest(j);
        end
        %���㶯̬����Ȩ��ֵ
        w=Wmax-(Wmax-Wmin)*i/T;
        %����λ�ú��ٶ�
        v(j,:)=w*v(j,:)+c1*rand*(p(j,:)-x(j,:))+c2*rand*(g-x(j,:));
        x(j,:)=x(j,:)+v(j,:);
        %�߽���������
        for ii=1:D
            if (v(j,ii)<Vmin)||(v(j,ii)>Vmax)
               v(j,ii)=rand*(Vmax-Vmin)+Vmin;
            end
            if (x(j,ii)<Xmin)||(x(j,ii)>Xmax)
                x(j,ii)=rand*(Xmax-Xmin)+Xmin;
            end
        end
    end
    gb(i)=gbest;
end
g;%���Ÿ���
gb(end);%����ֵ
figure
plot(gb)
xlabel('��������')
ylabel('��Ӧ��ֵ')
title('��Ӧ�Ƚ�������')
disp(['����õ�����ֵΪ��', num2str(gb(end))])
disp(['����ֵʱx(1)ֵΪ:',num2str(g(1)),'�� x(2)ֵΪ:',num2str(g(2))])
toc(t1)
%��Ӧ�Ⱥ���
function value=func2(x)
value=3*cos(x(1)*x(2))+x(1)+x(2)^2;
end

