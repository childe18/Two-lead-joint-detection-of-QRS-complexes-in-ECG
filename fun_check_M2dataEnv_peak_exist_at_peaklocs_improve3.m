function [peak_amp, peak_index,peak_exist ] = fun_check_M2dataEnv_peak_exist_at_peaklocs_improve3( peaklocs,M2data,search_win,Fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 函数作用：在数据M2data中的peaklocs位置，寻找是否存在peak，搜索窗口宽度为search_win
%%%% peaklocs：搜索peak的中心位置
%%%% search_win：搜索窗口宽度
%%%% peak_exist：如果存在peak，peak_exist=1，若不存在peak，peak_exist=0
%%%% peak_amp：若存在peak，peak_amp即为该peak的幅度
%%%% peak_index：若存在peak，peak_index即为该peak的索引,注意该index与peaklocs同起点
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%先寻找4秒内最大值，设置幅度判定门限
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

% 扩大区间上的最大值、均值,作为peak参考值
max_val=max(M2data(begin_pos:end_pos));
mean_val=mean(M2data(begin_pos:end_pos));
THR_peak1=max(max_val*0.4,mean_val*3);

%%%再寻找peaklocs附近窗口宽度为search_win范围内有无peak
begin_pos=peaklocs-round(0.5*search_win);
end_pos=peaklocs+round(0.5*search_win);
if begin_pos<1
    begin_pos=1;
end
if end_pos>length(M2data)
    end_pos=length(M2data);
end
% 寻找搜索区间上的最大值
tempdata=M2data(begin_pos:end_pos);
[peak_amp,peak_index]=max(tempdata);
peak_exist=0;
if peak_amp <= 2*mean_val % 认为没有peak
    peak_exist=0;
    return;
elseif peak_amp>=THR_peak1
    peak_exist=1;
    peak_index=peaklocs-round(0.5*search_win)+peak_index;
    return;
else 
    %先将搜索区间内peak削平
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
    % 如果削峰后剩余数据最大值不超过均值的2.5倍，则认为也存在较低的peak
    max_val2=max(tempdata);
    mean_val2=mean(tempdata);
    if max_val2 < 2.5*mean_val2
        peak_exist=1;
        peak_index=peaklocs-round(0.5*search_win)+peak_index;
        return;
    end
end


end
