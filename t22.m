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

raw=[5,1,4,2,3,6];%����һ�����������˻��ķɿ��źŵ�Ƶ��,����Ƶ��ͼδ֪����˴˴�ͨ�����Ƶ��������
phase=rand;%ͨ���������ʵ�����˻��ɿ��ź���λ��
freq=217;%���˻�ÿ������ƵƵ��
num=6;%����ɿ��ź�һ����������Ƶ����
freq_num=2*num;%������Ż��ɵ�Ƶ���
f=1:1:freq_num;%������Ż��ĵ�Ƶ����
samp_freq=1000;%���ø��Ż�����/��ƵƵ��
samp_t=100;%���ø��Ż�������ʱ��
max=samp_freq*samp_t;%����������

f_simu=zeros(max,1);%���ø��Ż��ź�Ƶ�ʾ���
f_real=zeros(max,1);%�������˻��ɿ��ź�Ƶ�ʾ���

j=1;%��ʼ��������
T=floor(num*samp_freq/freq);%һ���ɿ��ź����ڵĲ�������
T0=floor(samp_freq/freq);%һ����Ƶ����Ĳ�������

T_simu=0;%��ʼ�����Ż�̽��ķɿ��ź�����Ƶ��
T_f=0;%��ʼ�����Ż�̽������ʱ��Ĳο�Ƶ��
T_index1=0;%��ʼ�����Ż�̽��ɿ��ź����ڵĿ�ʼ���
T_index2=0;%��ʼ�����Ż�̽��ɿ��ź����ڵĽ������

nack_count=0;%��ʼ������NACK�źŴ���
nack_index=1;%��ʼ��NACK�źŵ�λ��
nack_freq=0;%��ʼ������NACK�ź�ʱ���Ƶ��

f_real(1)=f_raw(phase+1/samp_freq);
if (ifmatch(f(1+mod(j,freq_num)),f_real(1)))  %���õ�һ����������
    T_f=f(1+mod(j,freq_num));
    f_simu(1)=T_f;
    nack_count=nack_count+1;
    nack_index=1;
    nack_freq=f_simu(1);
end
if(~T_f)%�����Ż�̽������ʱ��Ĳο�Ƶ��Ϊ0��ʱ��ŶԵ�����+1
    j=j+1;
end

for i=2:1:max
    f_real(i)=f_raw(phase+i/samp_freq);%����ʵ�ʷɿ��źŸ��ݲ���Ƶ�ʲ���
    
    if (T_simu<=0)%��δ����ɿ��ź����ڵ�����£�
        
        f_simu(i)=T_f;%��ģ���ź�Ƶ�ʵ��ڸ��Ż�̽������ʱ��Ĳο�Ƶ��
        if (ifmatch(f(1+mod(j,freq_num)),f_real(i)))  %�����Ƶ��ȷ
            if(~T_f)%�����Ż�̽������ʱ��Ĳο�Ƶ��Ϊ0��ʱ��
                T_f=f(1+mod(j,freq_num));
                f_simu(i)=T_f;
            end
            
            continue;
            
        elseif (ifmatch(f_simu(i-1),f_real(i-1)))%�����������һ���Ƶ��ƥ�䵫��һ�㲻ƥ���ʱ��
            
            if (T_index1==0)%������Ż�̽��ɿ��ź����ڵĿ�ʼ��ǻ�δ���ã��������ÿ�ʼ���
                T_index1=i-1;
                
            else
                T_index2=i-1;%���ø��Ż�̽��ɿ��ź����ڵĽ������
            end
            
        end
        
        if(~T_f)%�����Ż�̽������ʱ��Ĳο�Ƶ��Ϊ0��ʱ��
            j=j+1;
        end
        
        T_simu=T_index2-T_index1-1;%���ݸ��Ż�̽��ɿ��ź����ڵĿ�ʼ������������ƶϷɿ��ź�һ�����ڵĲ�������
        
    else%������˷ɿ��ź�һ�����ڲ�������֮��
        
        if (f_simu(i)==0)%�����Ż�ģ��Ƶ��Ϊ�յ�ʱ��
            if (ifmatch(f_simu(i-1),f_real(i)))  %�����һ���Ƶ������ƥ�䵱ǰʵ�ʷɿ��ź�Ƶ�ʣ�����������ϵ��ź�Ƶ��
            
                f_simu(i)=f_simu(i-1);
                
                if(i+T_simu<=max)%����һ�����ڵĸõ㶼Ϊ��Ƶ��
                    f_simu(i+T_simu)=f_simu(i);
                end
                
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
            
                if(i+T_simu<=max)%����һ�����ڵĸõ㶼Ϊ��Ƶ��
                    f_simu(i+T_simu)=f_simu(i);
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
                if(i+T_simu<=max)%�������ϣ����¸����ڵĸõ�Ƶ�ʵ�Ϊ��ǰƵ��
                    f_simu(i+T_simu)=f_simu(i);
                end
                continue;
            else
                f_simu(i)=f_simu(i-1);
            end
        end
    end
end


fileID = fopen('data_t22.txt','a');
fprintf(fileID,'%d\r\n',nack_index/samp_freq);
fclose(fileID);

%{
plot(f_simu(nack_index-T:nack_index));hold on;
plot(f_real(nack_index-T:nack_index));hold off;
title(['����δ֪������£����Ż���ƵƵ��Ϊ',num2str(samp_freq),'ʱ�ķ���Ƶ�ʺͷɿ��źŶԱ�ͼ']);
xlabel('����SOS��ʱ�����ڵĲ�������');
ylabel('Ƶ��');
axis square;
legend({'���Ż�ģ���ź�Ƶ��','�ɿ��ź�ʵ��Ƶ��'},'Location','southeast')

disp(['���Ƴɹ������ѵ�ʱ��Ϊ:',num2str(nack_index/samp_freq)]);
%}
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

function sos=ifsos(i,f_simu)%���췢��SOS�źź���
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