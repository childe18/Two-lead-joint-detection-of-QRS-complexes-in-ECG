%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% ������fun_compare_QRS_detect_error_rate
%%%%%%  ���ܣ�ʶ�������,�Ƚ�ʶ����ȷ��
%%%%%%  ����ʹ�ã�
%%%%%%          correct_judge_QRS:��ȷʶ��������ANNOT_QRS�е�������
%%%%%%          missing_judge_QRS����©ʶ��������ANNOT_QRS�е�������
%%%%%%          error_judge_QRS������������pks�е�������
%%%%%%          ATRTIME_QRS��MIT-BIH��������Чʶ��ʱ���(���޳�����)
%%%%%%          pkstime�����㷨ʶ���ʱ���
%%%%%%  edited by xws��2020.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [correct_judge_QRS,missing_judge_QRS,error_judge_QRS]=fun_compare_QRS_detect_statistics(ATRTIME_QRS,pkstime)
%%%%%%%���ȼ�����ȷʶ��ĸ�����150������������Чʶ��
correct_index=1;
missing_index=1;
for i=1:length(ATRTIME_QRS)
    for j=1:length(pkstime)
        matching=0; %��ʼ��Ϊ��ƥ��
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

if missing_index==1 %û����©������²�ִ��
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

if error_judge_index==1 %û�����е�����²�ִ��
    error_judge_QRS(error_judge_index)=0;
end




