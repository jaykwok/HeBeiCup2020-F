clear all;
clc;

global H;
global P;
global sig2;
global K;
global pos

%���û�������
H=150;%�������˻����и߶�H,(m)
l0=10;%�������˻������ٶ�l0,(m/s)
P=10;%���û�վ���͹���,(W)
sig2=2*1E-5;%�����������ʦ�2 = 2��10-5,(w)
rho=3.1623;%�������ަ�1����2
max=100000;%���������������ޣ�ʵ���Ͽ���������ڴﵽ��������֮ǰ�ҵ�����ѵ㣩
t=0.5;%���嵥������ʱ��
left=[-l0*t,0];
right=[l0*t,0];
up=[0,l0*t];
down=[0,-l0*t];

K=3;%���蹲��K����վ
gamma=0.5;%�����ۿ����Ӧ�
a0=100;%�������˻���ʼ������
b0=100;%�������˻���ʼ������

%���û�������
x=zeros(max,2);%���˻���i������֮�������,���е�һ��Ϊ�����꣬�ڶ���Ϊ������

pos=zeros(K,2);%��վ������,���е�һ��Ϊ�����꣬�ڶ���Ϊ������
pos(1,:)=[1350, 1950];
pos(2,:)=[1600, 1790];
pos(3,:)=[1590, 2090];

%���ɿ���������뾶
r_safe=sqrt(P/(rho*sig2)-H^2);%�Ի�վΪԲ�ĵĿ�����������İ뾶
theta = linspace(0,2*pi,10000);%�������������
A(:,1) = pos(1,1)+r_safe*cos(theta);
A(:,2) = pos(1,2)+r_safe*sin(theta);
B(:,1) = pos(2,1)+r_safe*cos(theta);
B(:,2) = pos(2,2)+r_safe*sin(theta);
C(:,1) = pos(3,1)+r_safe*cos(theta);
C(:,2) = pos(3,2)+r_safe*sin(theta);
D = union(A,B,'rows');%�ϲ����л�վ�Ŀ���������
D = union(C,D,'rows');

R=zeros(max,1);%��i��ʱ�̵ĺ���Ϣ�������ۺ���

%�����һ�������㣺

ins=left;
best=getR([a0,b0]+left);
    
if (best<getR([a0,b0]+right))
    ins=right;
    best=getR([a0,b0]+right);
end
    
if (best<getR([a0,b0]+up))
    ins=up;
    best=getR([a0,b0]+up);
end
    
if (best<getR([a0,b0]+down))
    ins=down;
    best=getR([a0,b0]+down);
end
R(1)=best;
x(1,:)=[a0,b0]+ins;

t_1st=0;%��ʼ����һ�ν�������������ʱ��

for i=1:1:max
    
    if(t_1st==0)%�жϵ�һ�ν�������������ʱ��
        if(((x(i,1)-pos(1,1))^2+(x(i,2)-pos(1,2))^2<=r_safe^2)||((x(i,1)-pos(2,1))^2+(x(i,2)-pos(2,2))^2<=r_safe^2)||((x(i,1)-pos(3,1))^2+(x(i,2)-pos(3,2))^2<=r_safe^2))
            t_1st=i*t;
        end
    end
    
    if (best<getR(x(i,:)+left))
        ins=left;
        best=getR(x(i,:)+left);
    end
    
    if (best<getR(x(i,:)+right))
        ins=right;
        best=getR(x(i,:)+right);
    end
    
    if (best<getR(x(i,:)+up))
        ins=up;
        best=getR(x(i,:)+up);
    end
    
    if (best<getR(x(i,:)+down))
        ins=down;
        best=getR(x(i,:)+down);
    end
    
    if(best==R(i))
        index=i;
        break;
    end
    

    R(i+1)=best;
    x(i+1,:)=x(i,:)+ins;
end

E=0;
for i=1:1:max
    E=E+power(gamma,i-1)*R(i);
end

plot(D(:,1),D(:,2));hold on;%���ƿɱ���������

plot(x(1:index,1),x(1:index,2));%���ƺ���ͼ
s=scatter(pos(:,1),pos(:,2),"*");%���ƻ�վλ��
s.LineWidth = 5;
s.MarkerEdgeColor = 'black';

disp(['��һ�ν�������������ʱ��Ϊ:',num2str(t_1st)]);

disp(['��������������ʱ��Ϊ:',num2str(t*index)]);


axis square;
title(['�˺�������󻯳����ۿ�Ԥ��Ϊ',num2str(E)]);
xlabel('x������') 
ylabel('y������') 
legend({'��վ�Ŀ���������','���˻����й켣','��վλ��'},'Location','southeast')
hold off;

function Rtemp=getR(x)%�õ�����Ϣ��������
global H;
global P;
global sig2;
global K;
global pos

d2=zeros(1,K);%���˻������վ�ľ���
h=zeros(1,K);%��վ�����˻�֮����ŵ���������
Rk=zeros(1,K);%���˻���ͨ�Ż�վ����Ϣ�������ۺ�������
sum=0;

for j=1:1:K
    d2(j)=H^2+(x(1)-pos(j,1))^2+(x(2)-pos(j,2))^2;
    h(j)=1/d2(j);
    Rk(j)=log2(1+P*h(j)/sig2);
    sum=sum+Rk(j);
end

Rtemp=sum;
end
