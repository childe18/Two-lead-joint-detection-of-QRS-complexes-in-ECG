function [peak_amp, peak_index,peak_exist ] = fun_check_M2dataEnv_peak_exist_at_peaklocs_improve3( peaklocs,M2data,search_win,Fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% �������ã�������M2data�е�peaklocsλ�ã�Ѱ���Ƿ����peak���������ڿ��Ϊsearch_win
%%%% peaklocs������peak������λ��
%%%% search_win���������ڿ��
%%%% peak_exist���������peak��peak_exist=1����������peak��peak_exist=0
%%%% peak_amp��������peak��peak_amp��Ϊ��peak�ķ���
%%%% peak_index��������peak��peak_index��Ϊ��peak������,ע���index��peaklocsͬ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%��Ѱ��4�������ֵ�����÷����ж�����
begin_pos=peaklocs-round(2*Fs);
end_pos=peaklocs+round(2*Fs);
if begin_pos<1
    begin_pos=1;
    end_pos=peaklocs+round(4*Fs);
end
if end_pos>length(M2data)
    begin_pos=peaklocs-round(4*Fs);
    end_pos=length(M2data);
end

% ���������ϵ����ֵ����ֵ,��Ϊpeak�ο�ֵ
max_val=max(M2data(begin_pos:end_pos));
mean_val=mean(M2data(begin_pos:end_pos));
THR_peak1=max(max_val*0.4,mean_val*3);

%%%��Ѱ��peaklocs�������ڿ��Ϊsearch_win��Χ������peak
begin_pos=peaklocs-round(0.5*search_win);
end_pos=peaklocs+round(0.5*search_win);
if begin_pos<1
    begin_pos=1;
end
if end_pos>length(M2data)
    end_pos=length(M2data);
end
% Ѱ�����������ϵ����ֵ
tempdata=M2data(begin_pos:end_pos);
[peak_amp,peak_index]=max(tempdata);
peak_exist=0;
if peak_amp <= 2*mean_val % ��Ϊû��peak
    peak_exist=0;
    return;
elseif peak_amp>=THR_peak1
    peak_exist=1;
    peak_index=peaklocs-round(0.5*search_win)+peak_index;
    return;
else 
    %�Ƚ�����������peak��ƽ
    for i=peak_index:length(tempdata)
        if tempdata(i)> 2*mean_val
            tempdata(i)=mean_val;
        else
            break;
        end
    end
    
    if peak_index>1
        for i=peak_index-1:-1:1
            if tempdata(i)> 1.5*mean_val
                tempdata(i)=mean_val;
            else
                break;
            end 
        end
    end
    % ��������ʣ���������ֵ��������ֵ��2.5��������ΪҲ���ڽϵ͵�peak
    max_val2=max(tempdata);
    mean_val2=mean(tempdata);
    if max_val2 < 2.5*mean_val2
        peak_exist=1;
        peak_index=peaklocs-round(0.5*search_win)+peak_index;
        return;
    end
end


end
