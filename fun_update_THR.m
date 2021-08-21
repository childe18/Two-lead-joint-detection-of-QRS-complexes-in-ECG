function [ THR_Amp,THR_RR ] = fun_update_THR(old_THR_Amp,Amp_array,Amp_cof,old_THR_RR,Index_array,RR_cof )
%% 1
THR_Amp=old_THR_Amp*(1-Amp_cof)+Amp_array(end)*Amp_cof;
THR_RR=old_THR_RR*(1-RR_cof)+(Index_array(end)-Index_array(end-1))*RR_cof;

%% 2
% THR_Amp=old_THR_Amp*(1-Amp_cof)+(Amp_array(end)*0.7+Amp_array(end-1)*0.3)*Amp_cof;
% THR_RR=old_THR_RR*(1-RR_cof)+(Index_array(end)-Index_array(end-1))*RR_cof;

end

