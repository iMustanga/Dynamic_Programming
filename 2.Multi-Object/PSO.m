% 粒子群算法
%pop——种群数量
%dim——问题维度
%ub——变量上界，[1,dim]矩阵
%lb——变量下界，[1,dim]矩阵
%fobj——适应度函数（指针）
%MaxIter——最大迭代次数
%Best_Pos——x的最佳值 
%Best_Score——最优适应度
  
clc;
clear all;
close all;
t1=tic;
pop=50;  %种群数量
dim=2;   %问题维度

% ub=[10,10]; % 变量上界，[1,dim]矩阵
% lb=[-10,-10];% 变量下界，[1,dim]矩阵

ub=[10,10]; % 变量上界，[1,dim]矩阵
lb=[3,3];% 变量下界，[1,dim]矩阵

% vmax=[2,2];
% vmin=[-2,-2];

vmax=[2,2];
vmin=[-2,-2];

maxIter=100; % 最大迭代次数
fobj=@(X)fun(X); % 适应度函数（指针）
%Best_Pos——x的最佳值 Best_Score——最优适应度
[Best_Pos,Best_fitness,IterCurve]=pso(pop,dim,ub,lb,fobj,vmax,vmin,maxIter);
figure
plot(IterCurve,'b','linewidth',1.5);
xlabel('迭代次数')
ylabel('目标函数值')
grid on;
toc(t1)
disp(['求解得到的x1，x2是:',num2str(Best_Pos(1)),' ',num2str(Best_Pos(2))]);
disp(['最优解对应的函数:',num2str(Best_fitness)]);
 
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