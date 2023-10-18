% ����Ⱥ�㷨
%pop������Ⱥ����
%dim��������ά��
%ub���������Ͻ磬[1,dim]����
%lb���������½磬[1,dim]����
%fobj������Ӧ�Ⱥ�����ָ�룩
%MaxIter��������������
%Best_Pos����x�����ֵ 
%Best_Score����������Ӧ��
  
clc;
clear all;
close all;
t1=tic;
pop=50;  %��Ⱥ����
dim=2;   %����ά��

% ub=[10,10]; % �����Ͻ磬[1,dim]����
% lb=[-10,-10];% �����½磬[1,dim]����

ub=[10,10]; % �����Ͻ磬[1,dim]����
lb=[3,3];% �����½磬[1,dim]����

% vmax=[2,2];
% vmin=[-2,-2];

vmax=[2,2];
vmin=[-2,-2];

maxIter=100; % ����������
fobj=@(X)fun(X); % ��Ӧ�Ⱥ�����ָ�룩
%Best_Pos����x�����ֵ Best_Score����������Ӧ��
[Best_Pos,Best_fitness,IterCurve]=pso(pop,dim,ub,lb,fobj,vmax,vmin,maxIter);
figure
plot(IterCurve,'b','linewidth',1.5);
xlabel('��������')
ylabel('Ŀ�꺯��ֵ')
grid on;
toc(t1)
disp(['���õ���x1��x2��:',num2str(Best_Pos(1)),' ',num2str(Best_Pos(2))]);
disp(['���Ž��Ӧ�ĺ���:',num2str(Best_fitness)]);
 
function [X]=initialization(pop,ub,lb,dim)
    for i=1:pop
        for j=1:dim
            X(i,j)=(ub(j)-lb(j))*rand()+lb(j); 
        end 
    end
end
 
function fitness=fun(x)
    fitness=sum(x.^2);
end
 
function [X]=BoundaryCheck(X,ub,lb,dim)
    for i=1:dim
        if X(i)>ub(i)
            X(i)=ub(i);
        end
        if X(i)<lb(i)
            X(i)=lb(i);
        end
    end
end
 
function [Best_Pos,Best_fitness,IterCurve]=pso(pop,dim,ub,lb,fobj,vmax,vmin,maxIter)
c1=2.0;
c2=2.0;
V=initialization(pop,vmax,vmin,dim);
X=initialization(pop,ub,lb,dim);
fitness=zeros(1,pop);
for i=1:pop
    fitness(i)=fobj(X(i,:));
end
pBest=X;
pBestFitness=fitness;
[~,index]=min(fitness);
gBestFitness=fitness(index);
gBest=X(index,:);
Xnew=X;
fitnessNew=fitness;
for t=1:maxIter
    for i=1:pop
        r1=rand(1,dim);
        r2=rand(1:dim);
        V(i,:)=V(i,:)+c1.*r1.*(pBest(i,:)-X(i,:))+c2.*r2.*(gBest-X(i,:));
        V(i,:)=BoundaryCheck(V(i,:),vmax,vmin,dim);
        Xnew(i,:)=X(i,:)+V(i,:);
        fitnessNew(i)=fobj(Xnew(1,:));
        if fitnessNew(i)<pBestFitness(i)
            pBest(i,:)=Xnew(i,:);
            pBestFitness(i)=fitnessNew(i);
        end
        if fitnessNew(i)<gBestFitness
            gBestFitness=fitnessNew(i);
            gBest=Xnew(i,:);
        end
    end
    X=Xnew;
    fitness=fitnessNew;
    Best_Pos=gBest;
    Best_fitness=gBestFitness;
    IterCurve(t)=gBestFitness;
end
end