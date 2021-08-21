%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  功能：改进PT算法对ECG信号处理（改进3）
%%%%%%  增加了：（1）如果通道1的数据peak幅度过小，则从通道2的数据来辅助判断
%%%%%%          （2）直接从包络信号上判断，而不从带通滤波数据上判断
%%%%%%          (3)在imprvoe2的基础上，增加了噪声门限
%%%%%%          (4)带通滤波区间修改为5-20Hz
%%%%%%          (5)修改了间期门限参数
%%%%%%          （6）修改了通道2中peak判断时候的门限，原大于4倍噪声水平，改为3倍
%%%%%%           (7)3.1程序除了116号数据外性能较好，3.2程序增加了区域多个peak检测，但综合性能有下降。3.3综合了上述两个程序,效果较好
%%%%%%  edited by xws，2020.7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%6%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear ; close all;
load('D:\WorkFile\Matlab\ECGdata\232.mat');

ECGdataM1=M(:,1);       %% 取第一路数据
ECGdataM2=M(:,2);
Fs=360;                 %% 抽样率 

SAMPLES2READ=4000;      %% 每次处理的样点数
SampleStart=1+round(360*1353.2)-200;   %% 所取样点起始位置,注意样点编号从1开始
needPlot=1;             %% 是否需要画图

% SAMPLES2READ=650000;      %% 每次处理的样点数
% SampleStart=1;   %% 所取样点起始位置,注意样点编号从1开始
% needPlot=0;


data_Org=ECGdataM1(SampleStart:SampleStart+SAMPLES2READ-1);

%% 画图的时间刻度计算
TIME=(SampleStart-1)/Fs+(0:(SAMPLES2READ-1))/Fs;    %%%画图时，横轴时间刻度需要加上底数
[ANNOT_QRS,ATRTIME_QRS]=fun_collect_QRS_notes_from_matFile(ANNOT,ATRTIME);
ind_QRS= find(ATRTIME_QRS <= TIME(end) & ATRTIME_QRS>=TIME(1)); %%%标记及时间
ATRTIMED_QRS= ATRTIME_QRS(ind_QRS);
ANNOT_QRS=round(ANNOT_QRS);
ANNOTD_QRS= ANNOT_QRS(ind_QRS);

%% 画图，原始信号(仅标记有效R波峰) 
if needPlot==1
    figure(1);subplot(211);
    plot(TIME,data_Org,'b');
    
    for k=1:length(ATRTIMED_QRS)
        text(ATRTIMED_QRS(k),0,num2str(ANNOTD_QRS(k)));
    end;
    title('原始ECG信号（I）');
    
    figure(8);subplot(221);
    plot(TIME,data_Org,'b');
    title('116 record，lead I');
end


%% 根据CZT变换设置RR初始值
f1=0.8;f2=2.2;m=1024; % 看RR间期
w = exp(-j*2*pi*(f2-f1)/(m*Fs));
a = exp(j*2*pi*f1/Fs);
z = abs(czt(data_Org(1:10*Fs),m,w,a));
fz = ((0:length(z)-1)'*(f2-f1)/length(z)) + f1;
[max_val,max_index]=max(z);
max_index_f=max_index*(f2-f1)/length(z)+f1; 
max_index_f_half=max_index_f/2; %最高幅度处有可能对应二倍频处
if max_index_f_half>f1
    max_index_f_half_index=ceil((max_index_f_half-f1)*length(z)/(f2-f1));  %由频率反推半频的index
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
    title('CZT变换求RR初始值');
end


%% 带通滤波 (Filter in between 5-15 Hz)
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
    title('带通滤波ECG信号（I）');
end


%% 求前向一阶导数
data_Diff=diff(data_BPF);


%% 求信号包络
data_Env=abs(hilbert(data_Diff));
if needPlot==1
   figure(4);subplot(211);
    plot(TIME(1:SAMPLES2READ-1),data_Env,'b');
    for k=1:length(ATRTIMED_QRS)
        hold on;
        plot([ATRTIMED_QRS(k),ATRTIMED_QRS(k)],[min(data_Env),max(data_Env)],'linestyle',':','color','r');
    end;
    title('希尔伯特包络（I）');
end


%% 初始化包络信号阈值，检测包络peak
THR_sig_Env=max(data_Env(1:10*Fs))*0.25;
THR_noise_Env=mean(data_Env(1:10*Fs));
[pks,locs] = findpeaks(data_Env,'minpeakheight',THR_sig_Env,'MINPEAKDISTANCE',round(0.25*Fs));

pkstime=(locs+SampleStart-2)/Fs;
if needPlot
    hold on;
    plot(pkstime,pks,'o','color','r'); 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 通道II处理
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
    title('原始ECG信号（II）');
    
    figure(8);subplot(222);
    plot(TIME,dataM2_Org,'b');
    title('116 record，lead II');
end

if needPlot==1 
    figure(3);    subplot(212);
    plot(TIME,dataM2_BPF,'b');
    for k=1:length(ATRTIMED_QRS)
        hold on;
        plot([ATRTIMED_QRS(k),ATRTIMED_QRS(k)],[min(data_BPF),max(data_BPF)],'linestyle',':','color','g');
    end;
    title('带通滤波ECG信号（II）');
end

if needPlot==1 
   figure(4);subplot(212);
    plot(TIME(1:SAMPLES2READ-1),dataM2_Env,'b');
    for k=1:length(ATRTIMED_QRS)
        hold on;
        plot([ATRTIMED_QRS(k),ATRTIMED_QRS(k)],[min(data_Env),max(data_Env)],'linestyle',':','color','r');
    end;
    title('希尔伯特包络（II）');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 检测peak
% 初始化RR间期， 判断极值点
THR_RR_Env=RR_initial*Fs;
update_THR_Amp_cof=0.33;               % peak检测幅度更新系数
update_THR_RR_cof=0.1;               % peak检测RR间隔更新系数
QRS_Amp_array=[];           %保存QRS极值幅度
QRS_Index_array=[];         %保存QRS极值位置索引
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
        if latestRR <=round(0.28*Fs) %%%%%%%%%%%%%%%%%%%%%%%%%200ms不应期内，取幅度较高者
            if pks(i) > QRS_Amp_array(end) & pks(i)>2*THR_noise_Env
                QRS_Amp_array(end) =pks(i);
                QRS_Index_array(end)=peak_locs;
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
            end  
        elseif latestRR < round(0.6*THR_RR_Env)      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 200ms不应期------0.4*THR_RR_Env
            [peak_amp,peak_index_Env, peak_exist ] = fun_check_M2dataEnv_peak_exist_at_peaklocs_improve3( peak_locs-SampleStart,dataM2_Env,0.15*Fs,Fs );
            if  pks(i)>=0.3*THR_sig_Env  & peak_exist==1 
                % 如果peak幅度很高，则认为是新的R波峰
                QRS_Amp_array =[QRS_Amp_array   pks(i)];
                QRS_Index_array=[QRS_Index_array    peak_locs];
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
            end
            
        elseif latestRR < round(0.8*THR_RR_Env)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 0.4*THR_RR_Env---0.8*THR_RR_Env
            % 如果peak幅度很高，则认为是新的R波峰
            if pks(i)>=0.75*THR_sig_Env  
                % 如果peak幅度很高，则认为是新的R波峰
                QRS_Amp_array =[QRS_Amp_array   pks(i)];
                QRS_Index_array=[QRS_Index_array    peak_locs];
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
            elseif pks(i)>=0.25*THR_sig_Env  & pks(i)>1.5*THR_noise_Env    %通道1上有不高的peak。但通道2上有较高的peak
                 [peak_amp,peak_index_Env, peak_exist ] = fun_check_M2dataEnv_peak_exist_at_peaklocs_improve3( peak_locs-SampleStart,dataM2_Env,0.15*Fs,Fs );
                if peak_exist==1
                    QRS_Amp_array =[QRS_Amp_array   pks(i)];
                    QRS_Index_array=[QRS_Index_array    peak_locs];
                    [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
                end   
            end
        elseif latestRR < round(1.2*THR_RR_Env)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 0.8*THR_RR_Env--1.2*THR_RR_Env    正常情况,直接记录       
            QRS_Amp_array =[QRS_Amp_array   pks(i)];
            QRS_Index_array=[QRS_Index_array    peak_locs];
            [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
        elseif   latestRR < round(1.6*THR_RR_Env)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1.2*THR_RR_Env--1.6*THR_RR_Env    有可能漏检,重新检测区间为：上一次有效peak+0.4*THR_RR_Env----本次peak-0.4*THR_RR_Env
            begin_pos_Env_temp=QRS_Index_array(end)+round(0.4*Fs)-SampleStart+1;
            end_pos_Env_temp=peak_locs-round(0.4*Fs)-SampleStart+1;
            [max_sig_temp,max_index_temp]=max(data_Env(begin_pos_Env_temp:end_pos_Env_temp));
            if max_sig_temp>0.6*THR_sig_Env  & pks(i)>2*THR_noise_Env   %降门限检测
                % 先保存降门限新检测出来的那个peak
                qrs_index_temp=max_index_temp;
                QRS_Amp_array =[QRS_Amp_array   max_sig_temp];
                QRS_Index_array=[QRS_Index_array    SampleStart+begin_pos_Env_temp+max_index_temp-2];
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
                %再保存原来RR间隙较大的那个peak
                QRS_Amp_array =[QRS_Amp_array   pks(i)];
                QRS_Index_array=[QRS_Index_array    peak_locs];
                %更新门限
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
            elseif max_sig_temp>0.15*THR_sig_Env   & pks(i)>2*THR_noise_Env  %降门限检测，同时检测通道2peak
                centre_pos=round(0.5*(peak_locs+QRS_Index_array(end)));
                win=peak_locs-QRS_Index_array(end)-round(0.45*Fs)*2;
                [peak_amp,peak_index_Env, peak_exist ] = fun_check_M2dataEnv_peak_exist_at_peaklocs_improve3(centre_pos-SampleStart,dataM2_Env,win,Fs);
                if peak_exist==1
                    % 保存由通道2检测出来的那个peak
                    QRS_Amp_array =[QRS_Amp_array   max_sig_temp];
                    QRS_Index_array=[QRS_Index_array    peak_index_Env+SampleStart];
                    [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
                end
                % 再保存原来RR间隙较大的那个peak
                QRS_Amp_array =[QRS_Amp_array   pks(i)];
                QRS_Index_array=[QRS_Index_array    peak_locs];
                %更新门限
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
            else
                 %没有新检出peak，只保存原来RR间隙较大的那个peak
                QRS_Amp_array =[QRS_Amp_array   pks(i)];
                QRS_Index_array=[QRS_Index_array    peak_locs];
                %更新门限
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
            end
        elseif    latestRR < round(2.2*THR_RR_Env) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RR>1.6*THR_RR_Env    有可能漏检多个peak,重新检测区间为：上一次有效peak+0.4*THR_RR_Env----本次peak-0.4*THR_RR_Env
            begin_pos_Env_temp=QRS_Index_array(end)+round(0.4*Fs)-SampleStart+1;
            end_pos_Env_temp=peak_locs-round(0.4*Fs)-SampleStart+1;
            [max_sig_temp,max_index_temp]=max(data_Env(begin_pos_Env_temp:end_pos_Env_temp));
            if max_sig_temp>0.4*THR_sig_Env  & pks(i)>2*THR_noise_Env   %降门限检测
                % 先保存降门限新检测出来的那个peak
                qrs_index_temp=max_index_temp;
                QRS_Amp_array =[QRS_Amp_array   max_sig_temp];
                QRS_Index_array=[QRS_Index_array    SampleStart+begin_pos_Env_temp+max_index_temp-2];
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
                %再保存原来RR间隙较大的那个peak
                QRS_Amp_array =[QRS_Amp_array   pks(i)];
                QRS_Index_array=[QRS_Index_array    peak_locs];
                %更新门限
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
            elseif max_sig_temp>0.1*THR_sig_Env   & pks(i)>2*THR_noise_Env  %降门限检测，同时检测通道2peak
                centre_pos=round(0.5*(peak_locs+QRS_Index_array(end)));
                win=peak_locs-QRS_Index_array(end)-round(0.45*Fs)*2;
                [peak_amp,peak_index_Env, peak_exist ] = fun_check_M2dataEnv_peak_exist_at_peaklocs_improve3(centre_pos-SampleStart,dataM2_Env,win,Fs);
                if peak_exist==1
                    % 保存由通道2检测出来的那个peak
                    QRS_Amp_array =[QRS_Amp_array   max_sig_temp];
                    QRS_Index_array=[QRS_Index_array    peak_index_Env+SampleStart];
                    [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
                end
                % 再保存原来RR间隙较大的那个peak
                QRS_Amp_array =[QRS_Amp_array   pks(i)];
                QRS_Index_array=[QRS_Index_array    peak_locs];
                %更新门限
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
            else
                 %没有新检出peak，只保存原来RR间隙较大的那个peak
                QRS_Amp_array =[QRS_Amp_array   pks(i)];
                QRS_Index_array=[QRS_Index_array    peak_locs];
                %更新门限
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
            end
        elseif  latestRR < round(5*THR_RR_Env)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1.6-2.2 *THR_RR_Env    大概率存在漏检,重新检测区间为：上一次有效peak+0.4*THR_RR_Env----本次peak-0.4*THR_RR_Env
            centre_pos=round(0.5*(peak_locs+QRS_Index_array(end)))-SampleStart;
            win=peak_locs-QRS_Index_array(end)-round(0.4*Fs)*2;
            [peak_ampA,peak_indexA_Env, peak_existA ] = fun_check_M2dataEnv_peak_exist_at_peaklocs_improve3(centre_pos,dataM2_Env,win,Fs);
            if peak_existA==1
                % 保存由通道2检测出来的那个peakA
                QRS_Amp_array =[QRS_Amp_array   data_Env(peak_indexA_Env)];
                QRS_Index_array=[QRS_Index_array    peak_indexA_Env+SampleStart];
                
                % peakA与前面最后一个peak之间可能还有漏检的peakB
                if QRS_Index_array(end)-QRS_Index_array(end-1)>1.6*THR_RR_Env
                    centre_pos=round(0.5*(QRS_Index_array(end-1)+QRS_Index_array(end)))-SampleStart;
                    win=QRS_Index_array(end)-QRS_Index_array(end-1)-round(0.4*Fs)*2;
                    [peak_ampB,peak_indexB_Env, peak_existB ] = fun_check_M2dataEnv_peak_exist_at_peaklocs_improve3(centre_pos,dataM2_Env,win,Fs);
                    if peak_existB==1
                        % 保存由通道2检测出来的那个peakB,注意是插入
                        QRS_Amp_array=[QRS_Amp_array(1:end-1),data_Env(peak_indexB_Env),QRS_Amp_array(end)];
                        QRS_Index_array=[QRS_Index_array(1:end-1),peak_indexB_Env+SampleStart,QRS_Index_array(end)];
                    end
                end
                
                % peakA与后面刚刚检测出来的一个peak之间可能还有漏检的peakC
                if  peak_locs-QRS_Index_array(end)> 1.6*THR_RR_Env
                    centre_pos=round(0.5*(peak_locs+QRS_Index_array(end)))-SampleStart;
                    win=peak_locs-QRS_Index_array(end)-round(0.4*Fs)*2;
                    [peak_ampC,peak_indexC_Env, peak_existC ] = fun_check_M2dataEnv_peak_exist_at_peaklocs_improve3(centre_pos,dataM2_Env,win,Fs);
                    if peak_existC==1
                        % 保存由通道2检测出来的那个peakC
                        QRS_Amp_array =[QRS_Amp_array   data_Env(peak_indexC_Env)];
                        QRS_Index_array=[QRS_Index_array    peak_indexC_Env+SampleStart];
                    end
                end
                %更新门限
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
            end
            %没有新检出peak，只保存原来RR间隙较大的那个peak
            QRS_Amp_array =[QRS_Amp_array   pks(i)];
            QRS_Index_array=[QRS_Index_array    peak_locs];
            %更新门限
            [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
        else    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% >3.2 *THR_RR_Env    可能存在漏检多个peak
            begin_pos=QRS_Index_array(end)+round(0.4*Fs)-SampleStart+1;
            end_pos=peak_locs-round(0.4*Fs)-SampleStart+1;
            [peak_amp, peak_index,peak_num] = fun_check_M2dataEnv_peaks_num_in_area_improve3( begin_pos,end_pos,dataM2_Env,Fs);
            if peak_num>0
                for kk=1:peak_num
                    QRS_Amp_array =[QRS_Amp_array   data_Env(peak_index(kk))];
                    QRS_Index_array=[QRS_Index_array    peak_index(kk)+SampleStart];
                end
                %更新门限
                 [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
            end
            % 再保存原来RR间隙较大的那个peak
                QRS_Amp_array =[QRS_Amp_array   pks(i)];
                QRS_Index_array=[QRS_Index_array    peak_locs];
                %更新门限
                [ THR_sig_Env,THR_RR_Env ] = fun_update_THR(THR_sig_Env,QRS_Amp_array,update_THR_Amp_cof,THR_RR_Env,QRS_Index_array,update_THR_RR_cof );
        end     

   
    else  %%    QRS_Amp_array<5,即为初始化阶段
        if length(QRS_Index_array)==1     %第2个peak，取较大者
            latestRR = peak_locs-QRS_Index_array(end);
            if latestRR <=0.4*THR_RR_Env %若间隔RR<0.4*THR_RR_Env,z则认为是同一个peak，取幅度较高者
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
        else     %第1个peak,直接记录
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
    title('最终识别结果');
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