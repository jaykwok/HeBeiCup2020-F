clc;

%����һЩ��Ҫ����Ϊȫ�ֱ���
global raw;
global freq;
global num;
global T;
global T0;
global nack_count;
global nack_freq;
global nack_index;

freq_num=6;%������Ż��ɵ�Ƶ���
f=1:1:freq_num;%������Ż��ĵ�Ƶ����
raw=[5,1,4,2,3,6];%����һ�����������˻��ķɿ��źŵ�Ƶ��,����Ƶ��ͼδ֪����˴˴�ͨ�����Ƶ��������
phase=rand;%ͨ���������ʵ�����˻��ɿ��ź���λ��
freq=217;%���˻�ÿ������ƵƵ��
num=6;%����ɿ��ź�һ����������Ƶ����
samp_freq=1000;%���ø��Ż�����/��ƵƵ��
samp_t=100;%���ø��Ż�������ʱ��
max=samp_freq*samp_t;%����������

f_simu=zeros(max,1);%���ø��Ż��ź�Ƶ�ʾ���
f_real=zeros(max,1);%�������˻��ɿ��ź�Ƶ�ʾ���

j=1;%��ʼ��������
T=floor(num*samp_freq/freq);%һ���ɿ��ź����ڵĲ�������
T0=floor(samp_freq/freq);%һ����Ƶ����Ĳ�������

nack_count=0;%��ʼ������NACK�źŴ���
nack_index=1;%��ʼ��NACK�źŵ�λ��
nack_freq=0;%��ʼ������NACK�ź�ʱ���Ƶ��

f_real(1)=f_raw(phase+1/samp_freq);%����ʵ�ʷɿ��ź�Ƶ�ʵĵ�һ��������

if (ifmatch(f(1+mod(j,freq_num)),f_real(1)))  %���õ�һ��
    f_simu(1)=f(1+mod(j,freq_num));
    nack_count=nack_count+1;
    nack_index=1;
    nack_freq=f_simu(1);
    if(1+T<=max)
        f_simu(1+T)=f_simu(1);
    end
end
j=j+1;


for i=2:1:max
    
    f_real(i)=f_raw(phase+i/samp_freq);%����ʵ�ʷɿ��źŸ��ݲ���Ƶ�ʲ���
    
    if (f_simu(i)==0)%�����Ż�ģ��Ƶ��Ϊ�յ�ʱ��
        if (ifmatch(f_simu(i-1),f_real(i)))  %�����һ���Ƶ������ƥ�䵱ǰʵ�ʷɿ��ź�Ƶ�ʣ�����������ϵ��ź�Ƶ��
            
            f_simu(i)=f_simu(i-1);
            
            if (ifsos(i,f_simu(i)))%���SOS�ź�
                break;
            end
            
            continue;
        end
        if (ifmatch(f(1+mod(j,freq_num)),f_real(i)))  %�����Ƶ��ȷ
            f_simu(i)=f(1+mod(j,freq_num));
            
            if (ifsos(i,f_simu(i)))%���SOS�ź�
                break;
            end
            
            if(i+T<=max)%����һ�����ڵĸõ㶼Ϊ��Ƶ��
                f_simu(i+T)=f_simu(i);
            end
            continue;
        else
            j=j+1;%���ʧ�ܣ���������Ż�Ƶ��
        end
    else
        if (ifmatch(f_simu(i),f_real(i)))  %��֤������ź�Ƶ���Ƿ���ʵ�����
            
            if (ifsos(i,f_simu(i)))%���SOS�ź�
                break;
            end
            if(i+T<=max)%�������ϣ����¸����ڵĸõ�Ƶ�ʵ�Ϊ��ǰƵ��
                f_simu(i+T)=f_simu(i);
            end
            continue;
        else
            f_simu(i)=f_simu(i-1);
        end
    end
end

%{
fileID = fopen('data_t21.txt','a');
fprintf(fileID,'%d\r\n',nack_index/samp_freq);
fclose(fileID);
%}

plot(f_simu(nack_index-T:nack_index));hold on;
plot(f_real(nack_index-T:nack_index));hold off;
title(['������֪������£����Ż���ƵƵ��Ϊ',num2str(samp_freq),'ʱ�ķ���Ƶ�ʺͷɿ��źŶԱ�ͼ']);
xlabel('����SOS��ʱ�����ڵĲ�������');
ylabel('Ƶ��');
axis square;
legend({'���Ż�ģ���ź�Ƶ��','�ɿ��ź�ʵ��Ƶ��'},'Location','southeast')

disp(['���Ƴɹ������ѵ�ʱ��Ϊ:',num2str(nack_index/samp_freq)]);

function f=f_raw(t)%���˻���tʱ�̵ķɿ��ź�Ƶ��
global raw;
global freq;
global num;

t=mod(t,(num/freq));
index=floor(t*freq);
f=raw(index+1);
end

function NACK=ifmatch(f,f_raw)%���췢��NACK�źŵĺ���
if (f==f_raw)
    NACK=1;
else
    NACK=0;
end
end

function sos=ifsos(i,f_simu)
global nack_count;
global nack_freq;
global nack_index;
global num;
global T0;
sos=0;

if (nack_freq==0)
    nack_freq=f_simu;
end

if (f_simu~=nack_freq)
    if(i-nack_index>=T0)
        nack_count=0;
        nack_index=i;
    else
        nack_count=nack_count+1;
        nack_index=i;
        nack_freq=f_simu;
    end
else
    nack_index=i;
end
if (nack_count==num)
    sos=1;
end
end