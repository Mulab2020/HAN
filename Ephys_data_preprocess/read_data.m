function data = read_data(FileName,ch_num)
%READ_DATA reads data from a binary file in float precision
    disp('Loading data...');
    read_length = [];
    h = fopen(FileName,'r');
    X = fread(h,Inf,'float');
    fclose(h);




% for DET_VRorLS_shy_v4
    if ch_num == 16
        data.ch1 = X(1:ch_num:end);
        data.ch2 = X(2:ch_num:end);        
        data.camera =  X(3:ch_num:end);
        data.trial_num = X(5:ch_num:end);
        data.CurrObsDur = X(6:ch_num:end);
        data.fishx = X(7:ch_num:end);
        data.fishy = X(8:ch_num:end);
        data.stimParam1 = X(9:ch_num:end);
        data.stimParam2 = X(10:ch_num:end);
        data.fgain = X(11:ch_num:end);
        data.lgain = X(12:ch_num:end);
        data.rgain = X(13:ch_num:end);
        data.piezo =  X(14:ch_num:end);
        data.optoOngoing = X(15:ch_num:end);
        data.LRD =  X(16:ch_num:end);
    elseif ch_num==13
        data.ch1 = X(1:ch_num:end);
        data.ch2 = X(2:ch_num:end);
        data.camera = X(3:ch_num:end);
        data.mode = X(4:ch_num:end);
        data.trial_num = X(5:ch_num:end);
        data.CurrObsDur = X(6:ch_num:end);
        data.fishx = X(7:ch_num:end);
        data.fishy = X(8:ch_num:end);
        data.fgain = X(11:ch_num:end);
        data.lgain = X(12:ch_num:end);
        data.rgain = X(13:ch_num:end);
    elseif ch_num == 17
        data.ch1 = X(1:ch_num:end);
        data.ch2 = X(2:ch_num:end);
        data.camera =  X(3:ch_num:end);
        data.trial_num = X(5:ch_num:end);
        data.CurrObsDur = X(6:ch_num:end);
        data.fishx = X(7:ch_num:end);
        data.fishy = X(8:ch_num:end);
        data.stimParam1 = X(9:ch_num:end);
        data.stimParam2 = X(10:ch_num:end);
        data.fgain = X(11:ch_num:end);
        data.lgain = X(12:ch_num:end);
        data.rgain = X(13:ch_num:end);
        data.galvo_1 =  X(14:ch_num:end);
        data.galvo_2 = X(15:ch_num:end);
        data.galvo_3 =  X(16:ch_num:end);
        data.galvo_4 = X(17:ch_num:end);
%     elseif ch_num==13 % for 20230418_1_1_7d_VITL
%         data.ch1 = X(1:ch_num:420625*780);
%         data.ch2 = X(2:ch_num:420625*780);
%         data.camera = X(3:ch_num:420625*780);
%         data.mode = X(4:ch_num:420625*780);
%         data.trial_num = X(5:ch_num:420625*780);
%         data.CurrObsDur = X(6:ch_num:420625*780);
%         data.fishx = X(7:ch_num:420625*780);
%         data.fishy = X(8:ch_num:420625*780);
%         data.fgain = X(11:ch_num:420625*780);
%         data.lgain = X(12:ch_num:420625*780);
%         data.rgain = X(13:ch_num:420625*780);
    end
    clear X
    fieldNames = fieldnames(data);
    for i = 1:length(fieldNames)
        data_i = getfield(data,fieldNames{i});
        read_length = [read_length length(data_i)];
    end
    if length(unique(read_length))==1
        disp('Data loaded');
    else
        disp('Please check the ch_num');
    end
end