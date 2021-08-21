%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   fun_collect_QRS_notes_from_matFile
%%%%作用： 将原始记录数据提中的QRS复合波部分提取出来
%%%%       注意：只提取QRS复合波的部分标记文件
%%%%使用：      ANNOT：原始数据记录中的标记
%%%%            ATRTIME：原始数据记录中的时间标记
%%%%            ANNOT_QRS：提取QRS的标记
%%%%            ATRTIME_QRS：提取QRS的时间标记
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% 只需要统计数据库中1，2，3，4，5，6，7，8，9，10，11，12，13，34，38，符合
%%%% 《动态心电自动分析中QRS复合波检测算法研究（浙大博）》


% clear all; close all;
% load('D:\software\MatlabWorkFile\ECGdata\102.mat');
% QRS_index=find(ANNOT==1 | ANNOT==2 | ANNOT==3 | ANNOT==4 | ANNOT==5 | ANNOT==6 | ANNOT==7 | ANNOT==8 | ANNOT==9 | ANNOT==10 | ANNOT==11 | ANNOT==12 | ANNOT==13 | ANNOT==34 | ANNOT==38 );
% for i=1:length(QRS_index)
%     ANNOT_QRS(i,1)=ANNOT(QRS_index(i));
%     ATRTIME_QRS(i,1)=ATRTIME(QRS_index(i));
% end


function [ANNOT_QRS,ATRTIME_QRS]=fun_collect_QRS_notes_from_matFile(ANNOT,ATRTIME)
QRS_index=find(ANNOT==1 | ANNOT==2 | ANNOT==3 | ANNOT==4 | ANNOT==5 | ANNOT==6 | ANNOT==7 | ANNOT==8 | ANNOT==9 | ANNOT==10 | ANNOT==11 | ANNOT==12 | ANNOT==13 | ANNOT==34 | ANNOT==38 );
for i=1:length(QRS_index)
    ANNOT_QRS(i,1)=ANNOT(QRS_index(i));
    ATRTIME_QRS(i,1)=ATRTIME(QRS_index(i));
end

end


% clear all; close all;
% dataPath='D:\software\MatlabWorkFile\ECGdata\';
% %dataNum='100';
% dataType='.mat';
% %dataFile=strcat(dataPath,dataNum,dataType);
% %load(datafile);
% %load('D:\software\MatlabWorkFile\ECGdata\118.mat');
% 
% dataNumList={   '100','101','102','103','104','105','106','107','108','109',...
%                 '111','112','113','114','115','116','117','118','119','121', ...
%                 '122','123','124','200','201','202','203','205','207','208',...
%                 '209','210','212','213','214','215','217','219','220','221',...
%                 '222','223','228','230','231','232','233','234'};
% %%%% 只需要统计数据库中1，2，3，4，5，6，7，8，9，10，11，12，13，34，38，符合
% %%%% 《动态心电自动分析中QRS复合波检测算法研究（浙大博）》
% for i=1:length(dataNumList)
%     dataNum=dataNumList(i);
%     dataFile=strcat(dataPath,dataNum,dataType);
%     dataFile=dataFile{1};
%     load(dataFile);

