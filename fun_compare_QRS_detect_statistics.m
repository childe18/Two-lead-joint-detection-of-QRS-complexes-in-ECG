%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% 函数：fun_compare_QRS_detect_error_rate
%%%%%%  功能：识别结束后,比较识别正确率
%%%%%%  函数使用：
%%%%%%          correct_judge_QRS:正确识别索引（ANNOT_QRS中的索引）
%%%%%%          missing_judge_QRS：遗漏识别索引（ANNOT_QRS中的索引）
%%%%%%          error_judge_QRS：误判索引（pks中的索引）
%%%%%%          ATRTIME_QRS：MIT-BIH给出的有效识别时间点(已剔除部分)
%%%%%%          pkstime：本算法识别的时间点
%%%%%%  edited by xws，2020.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [correct_judge_QRS,missing_judge_QRS,error_judge_QRS]=fun_compare_QRS_detect_statistics(ATRTIME_QRS,pkstime)
%%%%%%%首先计算正确识别的个数，150毫秒内算作有效识别
correct_index=1;
missing_index=1;
for i=1:length(ATRTIME_QRS)
    for j=1:length(pkstime)
        matching=0; %初始化为不匹配
        if abs(ATRTIME_QRS(i)-pkstime(j))<=0.15
            correct_judge_QRS(correct_index)=i;
            correct_index=correct_index+1;
            matching=1;
            break;
        end
    end
    if matching==0
        missing_judge_QRS(missing_index)=i;
        missing_index=missing_index+1;
    end
end

if missing_index==1 %没有遗漏的情况下才执行
    missing_judge_QRS(missing_index)=0;
end



error_judge_index=1;
for i=1:length(pkstime)
    matching=0;
    for j=1:length(ATRTIME_QRS)
        if abs(pkstime(i)-ATRTIME_QRS(j))<=0.15
            matching=1;
            break;
        end
    end
    if matching==0
        error_judge_QRS(error_judge_index)=i;
        error_judge_index=error_judge_index+1;
    end
end

if error_judge_index==1 %没有误判的情况下才执行
    error_judge_QRS(error_judge_index)=0;
end




