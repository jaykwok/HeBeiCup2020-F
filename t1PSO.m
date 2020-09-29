clear all;
clc;

%{
����Ⱥ�㷨���Ĳ��֣�
v(k)=wv(k-1)+c1r1(pbest-x(k-1))+c2r2(gbest-x(k-1))
x(k+1)=x(k)+v(k)
%}

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
max=1000;%���������������ޣ�ʵ���Ͽ���������ڴﵽ��������֮ǰ�ҵ�����ѵ㣩
t=1/2;%���嵥������ʱ��
left=[-l0*t,0];
right=[l0*t,0];
up=[0,l0*t];
down=[0,-l0*t];
ins=[up;down;left;right];%�ֱ�����������ָ��


K=3;%���蹲��K����վ
gamma=0.5;%�����ۿ����Ӧ�
a0=100;%�������˻���ʼ������
b0=100;%�������˻���ʼ������

%���û�������
drone=30;%����һ�ηų������˻�����
x=zeros(max,2,drone);%���˻���i������֮�������,���е�һ��Ϊ�����꣬�ڶ���Ϊ������
v=zeros(max,2,drone);%�����ٶȸ�������
pbest=zeros(drone,2);%�������������ʷλ��

%��������Ⱥ�㷨���Ĳ���
w=0.8;
c1=1.49445;
c2=1.49445;

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


R=zeros(max,drone);%��i��ʱ�̵ĺ���Ϣ�������ۺ���

%�����һ�������㣺
for j=1:1:drone
    pbest(j,:)=[a0,b0];
end

gbest=[a0,b0];

for j=1:1:drone
    x(1,:,j)=[a0,b0]+ins(1+floor(4*rand),:);
    v(1,:,j)=ins(1+floor(4*rand),:);
    if (getR(pbest(j,:))<getR(x(1,:,j)))
        pbest(j,:)=x(1,:,j);
    end
    if (pbest(j,:)>gbest)
        gbest=pbest(j,:);
    end
    R(1,j)=getR(x(1,:,j));
end



t_1st=0;%��ʼ����һ�ν�������������ʱ��

index=zeros(drone,1);

for i=2:1:max
    for j=1:1:drone
        
        if(t_1st==0)%�жϵ�һ�ν�������������ʱ��
            if(norm(x(i,:,1)-pos(1,:))<=r_safe||norm((x(i,:,1)-pos(2,:)))<=r_safe||norm((x(i,:,1)-pos(3,:)))<=r_safe)
                t_1st=i*t;
            end
        end
        

        v(i,:,j)=w*v(i-1,:,j)+c1*rand*(pbest(j,:)-x(i-1,:,j))+c2*rand*(gbest-x(i-1,:,j));%�����ٶ�v(k)
        if(abs(v(i,1,j))>abs(v(i,2,j))&&v(i,1,j)>0)%�ж��ٶ�V���Ӷ�������������
            v(i,:,j)=right;
        elseif (abs(v(i,1,j))>abs(v(i,2,j))&&v(i,1,j)<0)
            v(i,:,j)=left;
        elseif (abs(v(i,1,j))<abs(v(i,2,j))&&v(i,2,j)>0)
            v(i,:,j)=up;
        elseif (abs(v(i,1,j))<abs(v(i,2,j))&&v(i,2,j)<0)
            v(i,:,j)=down;
        else
            v(i,:,j)=0;
        end
        x(i,:,j)=x(i-1,:,j)+v(i,:,j);
        
        if (getR(pbest(j,:))<getR(x(i,:,j)))
            pbest(j,:)=x(i,:,j);
        end
        if (getR(gbest)<getR(pbest(j,:)))
            gbest=pbest(j,:);
        end
        
        R(i,j)=getR(x(i,:,j));
        
    end
end
    
    


for i=1:1:max
    for j=1:1:drone
        if (~index(j))
            if(x(i,:,j)==[1505,1955])
                index(j)=i;
            end
        end
    end
    
end

E=0;

for i=1:1:index(drone)
    E=E+power(gamma,i-1)*R(i,drone);
end

plot(D(:,1),D(:,2));hold on;%���ƿɱ���������
s=scatter(pos(:,1),pos(:,2),"*");%���ƻ�վλ��
s.LineWidth = 5;
s.MarkerEdgeColor = 'black';

for j=1:1:drone

plot(x(1:index(j),1,j),x(1:index(j),2,j));%���ƺ���ͼ

axis square;
title(['�˺�������󻯳����ۿ�Ԥ��Ϊ',num2str(E)]);
xlabel('x������') 
ylabel('y������') 
end
legend({'��վ�Ŀ���������','��վλ��','���˻����й켣'},'Location','southeast')
hold off;

disp(['��һ�ν�������������ʱ��Ϊ:',num2str(t_1st)]);

disp(['��������������ʱ��Ϊ:',num2str(t*index(drone))]);

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
