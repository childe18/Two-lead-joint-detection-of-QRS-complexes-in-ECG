%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  ���ܣ��Ľ�PT�㷨��ECG�źŴ����Ľ�3��
%%%%%%  �����ˣ���1�����ͨ��1������peak���ȹ�С�����ͨ��2�������������ж�
%%%%%%          ��2��ֱ�ӴӰ����ź����жϣ������Ӵ�ͨ�˲��������ж�
%%%%%%          (3)��imprvoe2�Ļ����ϣ���������������
%%%%%%          (4)��ͨ�˲������޸�Ϊ5-20Hz
%%%%%%          (5)�޸��˼������޲���
%%%%%%          ��6���޸���ͨ��2��peak�ж�ʱ������ޣ�ԭ����4������ˮƽ����Ϊ3��
%%%%%%           (7)3.1�������116�����������ܽϺã�3.2����������������peak��⣬���ۺ��������½���3.3�ۺ���������������,Ч���Ϻ�
%%%%%%  edited by xws��2020.7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%6%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear ; close all;
load('D:\WorkFile\Matlab\ECGdata\232.mat');

ECGdataM1=M(:,1);       %% ȡ��һ·����
ECGdataM2=M(:,2);
Fs=360;                 %% ������ 

SAMPLES2READ=4000;      %% ÿ�δ����������
SampleStart=1+round(360*1353.2)-200;   %% ��ȡ������ʼλ��,ע�������Ŵ�1��ʼ
needPlot=1;             %% �Ƿ���Ҫ��ͼ

% SAMPLES2READ=650000;      %% ÿ�δ����������
% SampleStart=1;   %% ��ȡ������ʼλ��,ע�������Ŵ�1��ʼ
% needPlot=0;


data_Org=ECGdataM1(SampleStart:SampleStart+SAMPLES2READ-1);

%% ��ͼ��ʱ��̶ȼ���
TIME=(SampleStart-1)/Fs+(0:(SAMPLES2READ-1))/Fs;    %%%��ͼʱ������ʱ��̶���Ҫ���ϵ���
[ANNOT_QRS,ATRTIME_QRS]=fun_collect_QRS_notes_from_matFile(ANNOT,ATRTIME);
ind_QRS= find(ATRTIME_QRS <= TIME(end) & ATRTIME_QRS>=TIME(1)); %%%��Ǽ�ʱ��
ATRTIMED_QRS= ATRTIME_QRS(ind_QRS);
ANNOT_QRS=round(ANNOT_QRS);
ANNOTD_QRS= ANNOT_QRS(ind_QRS);

%% ��ͼ��ԭʼ�ź�(�������ЧR����) 
if needPlot==1
    figure(1);subplot(211);
    plot(TIME,data_Org,'b');
    
    for k=1:length(ATRTIMED_QRS)
        text(ATRTIMED_QRS(k),0,num2str(ANNOTD_QRS(k)));
    end;
    title('ԭʼECG�źţ�I��');
    
    figure(8);subplot(221);
    plot(TIME,data_Org,'b');
    title('116 record��lead I');
end


%% ����CZT�任����RR��ʼֵ
f1=0.8;f2=2.2;m=1024; % ��RR����
w = exp(-j*2*pi*(f2-f1)/(m*Fs));
a = exp(j*2*pi*f1/Fs);
z = abs(czt(data_Org(1:10*Fs),m,w,a));
fz = ((0:length(z)-1)'*(f2-f1)/length(z)) + f1;
[max_val,max_index]=max(z);
max_index_f=max_index*(f2-f1)/length(z)+f1; 
max_index_f_half=max_index_f/2; %��߷��ȴ��п��ܶ�Ӧ����Ƶ��
if max_index_f_half>f1
    max_index_f_half_index=ceil((max_index_f_half-f1)*length(z)/(f2-f1));  %��Ƶ�ʷ��ư�Ƶ��index
    if max_index_f_half_index-30<1
        half_f_begin=1;
    else
        half_f_begin=max_index_f_half_index-30;
    end
    half_f_end=max_index_f_half_index+30;
    [z_half_f_max,half_f_max_pos]=max(z(half_f_begin:half_f_end));
    if z_half_f_max>=0.4*max_val
        max_index_f=(half_f_begin+half_f_max_pos-1)*(f2-f1)/length(z)+f1; 
        max_val=z_half_f_max;
    end
end
RR_initial=1/max_index_f;

if needPlot==1 
    figure(2);
    plot(fz,z);
    hold on; 
    plot(max_index_f,max_val,'o','color','r');
    title('CZT�任��RR��ʼֵ');
end


%% ��ͨ�˲� (Filter in between 5-15 Hz)
f1=5; %cuttoff low frequency to get rid of baseline wander
f2=20; %cuttoff frequency to discard high frequency noise
Wn=[f1 f2]*2/Fs; % cutt off based on fs
N = 3; % order of 3 less processing
[a,b] = butter(N,Wn); %bandpass filtering
data_temp1 = filtfilt(a,b,data_Org);
normalization_Ref=max( abs(data_temp1));
data_BPF = data_temp1/ normalization_Ref;
clear data_temp1;
if needPlot==1
    figure(3);subplot(211);
    plot(TIME,data_BPF,'b');
    for k=1:length(ATRTIMED_QRS)
        hold on;
        plot([ATRTIMED_QRS(k),ATRTIMED_QRS(k)],[min(data_BPF),max(data_BPF)],'linestyle',':','color','g');
    end;
    title('��ͨ�˲�ECG�źţ�I��');
end


%% ��ǰ��һ�׵���
data_Diff=diff(data_BPF);


%% ���źŰ���
data_Env=abs(hilbert(data_Diff));
if needPlot==1
   figure(4);subplot(211);
    plot(TIME(1:SAMPLES2READ-1),data_Env,'b');
    for k=1:length(ATRTIMED_QRS)
        hold on;
        plot([ATRTIMED_QRS(k),ATRTIMED_QRS(k)],[min(data_Env),max(data_Env)],'linestyle',':','color','r');
    end;
    title('ϣ�����ذ��磨I��');
end


%% ��ʼ�������ź���ֵ��������peak
THR_sig_Env=max(data_Env(1:10*Fs))*0.25;
THR_noise_Env=mean(data_Env(1:10*Fs));
[pks,locs] = findpeaks(data_Env,'minpeakheight',THR_sig_Env,'MINPEAKDISTANCE',round(0.25*Fs));

pkstime=(locs+SampleStart-2)/Fs;
if needPlot
    hold on;
    plot(pkstime,pks,'o','color','r'); 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ͨ��II����
dataM2_Org=ECGdataM2(SampleStart:SampleStart+SAMPLES2READ-1);
data_temp2=filtfilt(a,b,dataM2_Org);
dataM2_BPF=data_temp2/ normalization_Ref;
clear data_temp2;
dataM2_Diff=diff(dataM2_BPF);
dataM2_Env=abs(hilbert(dataM2_Diff));

if needPlot==1
    figure(1);subplot(212);
    plot(TIME,dataM2_Org,'b');
    
    for k=1:length(ATRTIMED_QRS)
        text(ATRTIMED_QRS(k),0,num2str(ANNOTD_QRS(k)));
    end;
    title('ԭʼECG�źţ�II��');
    
    figure(8);subplot(222);
    plot(TIME,dataM2_Org,'b');
    title('116 record��lead II');
end

if needPlot==1 
    figure(3);    subplot(212);
    plot(TIME,dataM2_BPF,'b');
    for k=1:length(ATRTIMED_QRS)
        hold on;
        plot([ATRTIMED_QRS(k),ATRTIMED_QRS(k)],[min(data_BPF),max(data_BPF)],'linestyle',':','color','g');
    end;
    title('��ͨ�˲�ECG�źţ�II��');
end

if needPlot==1 
   figure(4);subplot(212);
    plot(TIME(1:SAMPLES2READ-1),dataM2_Env,'b');
    for k=1:length(ATRTIMED_QRS)
        hold on;
        plot([ATRTIMED_QRS(k),ATRTIMED_QRS(k)],[min(data_Env),max(data_Env)],'linestyle',':','color','r');
    end;
    title('ϣ�����ذ��磨II��');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ���peak
% ��ʼ��RR���ڣ� �жϼ�ֵ��
THR_RR_Env=RR_initial*Fs;
update_THR_Amp_cof=0.33;               % peak�����ȸ���ϵ��
update_THR_RR_cof=0.1;               % peak���RR�������ϵ��
QRS_Amp_array=[];           %����QRS��ֵ����
QRS_Index_array=[];         %����QRS��ֵλ������
for i = 1 : length(pks)
    peak_locs=locs(i)+SampleStart-1;
    if length(QRS_Amp_array) >= 2
        if locs(i)-5*Fs<1
            bgn=1;
        else
            bgn=locs(i)-5*Fs;
        end
        THR_noise_Env=mean(data_Env(bgn:locs(i)));
        latestRR = peak_locs-QRS_Index_array(end);
        if latestRR <=round(0.28*Fs) %%%%%%%%%%%%%%%%%%%%%%%%%200ms��Ӧ���ڣ�ȡ���Ƚϸ���
            if pks(i) > QRS_Amp_array(end) & pks(i)>2*THR_noise_Env
                QRS_Amp_array(end) =pks(i);
                QRS_Index_array(end)=peak_locs;
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
            end  
        elseif latestRR < round(0.6*THR_RR_Env)      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 200ms��Ӧ��------0.4*THR_RR_Env
            [peak_amp,peak_index_Env, peak_exist ] = fun_check_M2dataEnv_peak_exist_at_peaklocs_improve3( peak_locs-SampleStart,dataM2_Env,0.15*Fs,Fs );
            if  pks(i)>=0.3*THR_sig_Env  & peak_exist==1 
                % ���peak���Ⱥܸߣ�����Ϊ���µ�R����
                QRS_Amp_array =[QRS_Amp_array   pks(i)];
                QRS_Index_array=[QRS_Index_array    peak_locs];
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
            end
            
        elseif latestRR < round(0.8*THR_RR_Env)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 0.4*THR_RR_Env---0.8*THR_RR_Env
            % ���peak���Ⱥܸߣ�����Ϊ���µ�R����
            if pks(i)>=0.75*THR_sig_Env  
                % ���peak���Ⱥܸߣ�����Ϊ���µ�R����
                QRS_Amp_array =[QRS_Amp_array   pks(i)];
                QRS_Index_array=[QRS_Index_array    peak_locs];
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
            elseif pks(i)>=0.25*THR_sig_Env  & pks(i)>1.5*THR_noise_Env    %ͨ��1���в��ߵ�peak����ͨ��2���нϸߵ�peak
                 [peak_amp,peak_index_Env, peak_exist ] = fun_check_M2dataEnv_peak_exist_at_peaklocs_improve3( peak_locs-SampleStart,dataM2_Env,0.15*Fs,Fs );
                if peak_exist==1
                    QRS_Amp_array =[QRS_Amp_array   pks(i)];
                    QRS_Index_array=[QRS_Index_array    peak_locs];
                    [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
                end   
            end
        elseif latestRR < round(1.2*THR_RR_Env)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 0.8*THR_RR_Env--1.2*THR_RR_Env    �������,ֱ�Ӽ�¼       
            QRS_Amp_array =[QRS_Amp_array   pks(i)];
            QRS_Index_array=[QRS_Index_array    peak_locs];
            [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
        elseif   latestRR < round(1.6*THR_RR_Env)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1.2*THR_RR_Env--1.6*THR_RR_Env    �п���©��,���¼������Ϊ����һ����Чpeak+0.4*THR_RR_Env----����peak-0.4*THR_RR_Env
            begin_pos_Env_temp=QRS_Index_array(end)+round(0.4*Fs)-SampleStart+1;
            end_pos_Env_temp=peak_locs-round(0.4*Fs)-SampleStart+1;
            [max_sig_temp,max_index_temp]=max(data_Env(begin_pos_Env_temp:end_pos_Env_temp));
            if max_sig_temp>0.6*THR_sig_Env  & pks(i)>2*THR_noise_Env   %�����޼��
                % �ȱ��潵�����¼��������Ǹ�peak
                qrs_index_temp=max_index_temp;
                QRS_Amp_array =[QRS_Amp_array   max_sig_temp];
                QRS_Index_array=[QRS_Index_array    SampleStart+begin_pos_Env_temp+max_index_temp-2];
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
                %�ٱ���ԭ��RR��϶�ϴ���Ǹ�peak
                QRS_Amp_array =[QRS_Amp_array   pks(i)];
                QRS_Index_array=[QRS_Index_array    peak_locs];
                %��������
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
            elseif max_sig_temp>0.15*THR_sig_Env   & pks(i)>2*THR_noise_Env  %�����޼�⣬ͬʱ���ͨ��2peak
                centre_pos=round(0.5*(peak_locs+QRS_Index_array(end)));
                win=peak_locs-QRS_Index_array(end)-round(0.45*Fs)*2;
                [peak_amp,peak_index_Env, peak_exist ] = fun_check_M2dataEnv_peak_exist_at_peaklocs_improve3(centre_pos-SampleStart,dataM2_Env,win,Fs);
                if peak_exist==1
                    % ������ͨ��2���������Ǹ�peak
                    QRS_Amp_array =[QRS_Amp_array   max_sig_temp];
                    QRS_Index_array=[QRS_Index_array    peak_index_Env+SampleStart];
                    [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
                end
                % �ٱ���ԭ��RR��϶�ϴ���Ǹ�peak
                QRS_Amp_array =[QRS_Amp_array   pks(i)];
                QRS_Index_array=[QRS_Index_array    peak_locs];
                %��������
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
            else
                 %û���¼��peak��ֻ����ԭ��RR��϶�ϴ���Ǹ�peak
                QRS_Amp_array =[QRS_Amp_array   pks(i)];
                QRS_Index_array=[QRS_Index_array    peak_locs];
                %��������
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
            end
        elseif    latestRR < round(2.2*THR_RR_Env) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RR>1.6*THR_RR_Env    �п���©����peak,���¼������Ϊ����һ����Чpeak+0.4*THR_RR_Env----����peak-0.4*THR_RR_Env
            begin_pos_Env_temp=QRS_Index_array(end)+round(0.4*Fs)-SampleStart+1;
            end_pos_Env_temp=peak_locs-round(0.4*Fs)-SampleStart+1;
            [max_sig_temp,max_index_temp]=max(data_Env(begin_pos_Env_temp:end_pos_Env_temp));
            if max_sig_temp>0.4*THR_sig_Env  & pks(i)>2*THR_noise_Env   %�����޼��
                % �ȱ��潵�����¼��������Ǹ�peak
                qrs_index_temp=max_index_temp;
                QRS_Amp_array =[QRS_Amp_array   max_sig_temp];
                QRS_Index_array=[QRS_Index_array    SampleStart+begin_pos_Env_temp+max_index_temp-2];
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
                %�ٱ���ԭ��RR��϶�ϴ���Ǹ�peak
                QRS_Amp_array =[QRS_Amp_array   pks(i)];
                QRS_Index_array=[QRS_Index_array    peak_locs];
                %��������
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
            elseif max_sig_temp>0.1*THR_sig_Env   & pks(i)>2*THR_noise_Env  %�����޼�⣬ͬʱ���ͨ��2peak
                centre_pos=round(0.5*(peak_locs+QRS_Index_array(end)));
                win=peak_locs-QRS_Index_array(end)-round(0.45*Fs)*2;
                [peak_amp,peak_index_Env, peak_exist ] = fun_check_M2dataEnv_peak_exist_at_peaklocs_improve3(centre_pos-SampleStart,dataM2_Env,win,Fs);
                if peak_exist==1
                    % ������ͨ��2���������Ǹ�peak
                    QRS_Amp_array =[QRS_Amp_array   max_sig_temp];
                    QRS_Index_array=[QRS_Index_array    peak_index_Env+SampleStart];
                    [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
                end
                % �ٱ���ԭ��RR��϶�ϴ���Ǹ�peak
                QRS_Amp_array =[QRS_Amp_array   pks(i)];
                QRS_Index_array=[QRS_Index_array    peak_locs];
                %��������
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
            else
                 %û���¼��peak��ֻ����ԭ��RR��϶�ϴ���Ǹ�peak
                QRS_Amp_array =[QRS_Amp_array   pks(i)];
                QRS_Index_array=[QRS_Index_array    peak_locs];
                %��������
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
            end
        elseif  latestRR < round(5*THR_RR_Env)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1.6-2.2 *THR_RR_Env    ����ʴ���©��,���¼������Ϊ����һ����Чpeak+0.4*THR_RR_Env----����peak-0.4*THR_RR_Env
            centre_pos=round(0.5*(peak_locs+QRS_Index_array(end)))-SampleStart;
            win=peak_locs-QRS_Index_array(end)-round(0.4*Fs)*2;
            [peak_ampA,peak_indexA_Env, peak_existA ] = fun_check_M2dataEnv_peak_exist_at_peaklocs_improve3(centre_pos,dataM2_Env,win,Fs);
            if peak_existA==1
                % ������ͨ��2���������Ǹ�peakA
                QRS_Amp_array =[QRS_Amp_array   data_Env(peak_indexA_Env)];
                QRS_Index_array=[QRS_Index_array    peak_indexA_Env+SampleStart];
                
                % peakA��ǰ�����һ��peak֮����ܻ���©���peakB
                if QRS_Index_array(end)-QRS_Index_array(end-1)>1.6*THR_RR_Env
                    centre_pos=round(0.5*(QRS_Index_array(end-1)+QRS_Index_array(end)))-SampleStart;
                    win=QRS_Index_array(end)-QRS_Index_array(end-1)-round(0.4*Fs)*2;
                    [peak_ampB,peak_indexB_Env, peak_existB ] = fun_check_M2dataEnv_peak_exist_at_peaklocs_improve3(centre_pos,dataM2_Env,win,Fs);
                    if peak_existB==1
                        % ������ͨ��2���������Ǹ�peakB,ע���ǲ���
                        QRS_Amp_array=[QRS_Amp_array(1:end-1),data_Env(peak_indexB_Env),QRS_Amp_array(end)];
                        QRS_Index_array=[QRS_Index_array(1:end-1),peak_indexB_Env+SampleStart,QRS_Index_array(end)];
                    end
                end
                
                % peakA�����ոռ�������һ��peak֮����ܻ���©���peakC
                if  peak_locs-QRS_Index_array(end)> 1.6*THR_RR_Env
                    centre_pos=round(0.5*(peak_locs+QRS_Index_array(end)))-SampleStart;
                    win=peak_locs-QRS_Index_array(end)-round(0.4*Fs)*2;
                    [peak_ampC,peak_indexC_Env, peak_existC ] = fun_check_M2dataEnv_peak_exist_at_peaklocs_improve3(centre_pos,dataM2_Env,win,Fs);
                    if peak_existC==1
                        % ������ͨ��2���������Ǹ�peakC
                        QRS_Amp_array =[QRS_Amp_array   data_Env(peak_indexC_Env)];
                        QRS_Index_array=[QRS_Index_array    peak_indexC_Env+SampleStart];
                    end
                end
                %��������
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
            end
            %û���¼��peak��ֻ����ԭ��RR��϶�ϴ���Ǹ�peak
            QRS_Amp_array =[QRS_Amp_array   pks(i)];
            QRS_Index_array=[QRS_Index_array    peak_locs];
            %��������
            [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
        else    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% >3.2 *THR_RR_Env    ���ܴ���©����peak
            begin_pos=QRS_Index_array(end)+round(0.4*Fs)-SampleStart+1;
            end_pos=peak_locs-round(0.4*Fs)-SampleStart+1;
            [peak_amp, peak_index,peak_num] = fun_check_M2dataEnv_peaks_num_in_area_improve3( begin_pos,end_pos,dataM2_Env,Fs);
            if peak_num>0
                for kk=1:peak_num
                    QRS_Amp_array =[QRS_Amp_array   data_Env(peak_index(kk))];
                    QRS_Index_array=[QRS_Index_array    peak_index(kk)+SampleStart];
                end
                %��������
                 [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
            end
            % �ٱ���ԭ��RR��϶�ϴ���Ǹ�peak
                QRS_Amp_array =[QRS_Amp_array   pks(i)];
                QRS_Index_array=[QRS_Index_array    peak_locs];
                %��������
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
        end     

   
    else  %%    QRS_Amp_array<5,��Ϊ��ʼ���׶�
        if length(QRS_Index_array)==1     %��2��peak��ȡ�ϴ���
            latestRR = peak_locs-QRS_Index_array(end);
            if latestRR <=0.4*THR_RR_Env %�����RR<0.4*THR_RR_Env,z����Ϊ��ͬһ��peak��ȡ���Ƚϸ���
                if pks(i) > QRS_Amp_array(end) 
                    QRS_Amp_array(end) =pks(i);
                    QRS_Index_array(end)=peak_locs;
                end 
            else
                QRS_Amp_array =[QRS_Amp_array   pks(i)];
                QRS_Index_array=[QRS_Index_array    peak_locs];
                THR_sig_Env=THR_sig_Env*(1-update_THR_Amp_cof)+update_THR_Amp_cof*QRS_Amp_array(end);
                THR_RR_Env=THR_RR_Env*(1-update_THR_RR_cof)+update_THR_RR_cof*(QRS_Index_array(end)-QRS_Index_array(end-1));
            end
        else     %��1��peak,ֱ�Ӽ�¼
            QRS_Amp_array =[QRS_Amp_array   pks(i)];
            QRS_Index_array=[QRS_Index_array    peak_locs];
        end
    end
end

%% 
if needPlot
    figure(6);
    plot(TIME(1:SAMPLES2READ-1),data_Env,'b');
    hold on;
    plot(QRS_Index_array/Fs,QRS_Amp_array,'o','color','r'); 
    hold on;
    for k=1:length(ATRTIMED_QRS)
        plot([ATRTIMED_QRS(k),ATRTIMED_QRS(k)],[0,1],'linestyle',':','color','g');
    end;
    title('����ʶ����');
end

[correct_judge_QRS,missing_judge_QRS,error_judge_QRS]=fun_compare_QRS_detect_statistics(ATRTIME_QRS,QRS_Index_array/Fs);
if length(missing_judge_QRS)==1 && missing_judge_QRS(1)==0
    miss_num=0;
else
    miss_num=length(missing_judge_QRS);
end
if length(error_judge_QRS)==1 && error_judge_QRS(1)==0
    error_num=0;
else
    error_num=length(error_judge_QRS);
end

info=['missing_num=',num2str(miss_num),'     error_num=',num2str(error_num)];
disp(info);