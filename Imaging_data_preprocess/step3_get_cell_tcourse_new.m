clear all;close all;
tic
% img_dir='Z:\YuMu\SPIM\20140722\20140722_2_2_cy74_6d_givingup_replay_0vel_20140722_203337\Registered';
% input_dirs=struct;
% input_dirs(1).dir='X:\YuMu\20140916_fish2_12\20140916_2_1_cy74_6d_GU_stationary_OMR_20140916_175006\Registered';
folder_num = 1;
input_dir = 'E:\shy\20250507_3_1_6dpf_VODL2_20250507_133458\registered';
motion_name = [input_dir, '\motion_tcourse.tif'];
% finf = imfinfo(motion_name);
% motion_fig = figure('position',[100 100 1700 1200]);
% for i = 1:(min(numel(finf),32))
%     tmp = imread(motion_name,i);
%     imagesc(tmp);
%      set(gca,'visible','off')
%      waitforbuttonpress
% end
ending_frame =input('input the last frame:');
frame_rate =input('input frame rate:');

%%
%sample_len=190;   %frequency for cauculating the baseline
sample_len = round(frame_rate * 60* 3); % cauculating baseline within 3 mintutes window, considering using 5 min if slow responses appears
corr_thre=0.7;
move_thre=1;
poolnumber=6;
adapting_frame = 1;
offset = adapting_frame - 1;
save(fullfile(input_dir,'offset'),'offset');

frame_name = [input_dir, '\frame_info.txt'];
fileID = fopen(frame_name,'w');
fprintf(fileID, '%6.2f  %6.2f',adapting_frame,ending_frame);
fclose(fileID);

%%
%for dd=1:length(folder_num)
    
%     input_dir=input_dirs(dd).dir;

    load(fullfile(input_dir,'cell_info.mat'));
%     load(fullfile(input_dir,'motion_param.mat'));
%     ave_stack=readtiff(fullfile(input_dir,'ave.tif'));
    info = imfinfo([input_dir, '\ave.tif']);
    numberOfPages = length(info);
    for k = 1 : numberOfPages
          ave_stack(:,:,k)  = imread([input_dir, '\ave.tif'],k);
    end
%     dim=read_LSstack_size(fullfile(img_dir,'Stack dimensions.log'));
%     ave_stack=read_LSstack_fast_float(fullfile(img_dir,'ave.stackf'),dim);
    background = imread(fullfile(input_dir,'Background_1.tif'));  

    numcell=length(cell_info);
    dim=size(ave_stack);
    if length(dim)==2
        dim = [dim, 1];
    end

    cell_zlist=[cell_info.slice];

    zplane_list=min(cell_zlist):max(cell_zlist);


    outputName = [input_dir, '/Plane' num2str(1, '%.2d') '.stack'];
    dim2=read_LSstack_info(outputName,dim(1:2));
    totlen=dim2(3);

    zplane_inds=struct;
    for i=zplane_list
        zplane_inds(i).cellinds=find(cell_zlist==i);
    end

    disp('Getting raw time course and applying low cut filter...');
    cell_resp=zeros(numcell,totlen,'single');

    cell_resp1_cells=cell(1,length(zplane_list));

   parpool(poolnumber);

    parfor zplane=zplane_list
        cellinds=zplane_inds(zplane).cellinds;  
        outputName = [input_dir, '\Plane' num2str(zplane, '%.2d') '.stack'];
        stack=read_LSstack_fast1(outputName,dim);
        disp(['Getting time course of plane ',num2str(zplane)]);
        cellresp1=zeros(length(cellinds),totlen);
        cinfo=cell_info(cellinds);
        sliceinds=int64(0:(totlen-1))*int64(dim(1)*dim(2));

        for i=1:length(cellinds)    
            inds=cinfo(i).inds;
            tcourse=single(get_cell_tcourse_mex64(stack,sliceinds,int64(inds)));
            cellresp1(i,:)=tcourse;
        end
        %clear(stack);
        cell_resp1_cell{zplane}=cellresp1;

    end

    delete(gcp('nocreate'));

    for i=zplane_list
        cellinds=zplane_inds(i).cellinds;  
        cell_resp(cellinds,:)=cell_resp1_cell{i};
    end
    
    if ending_frame == 0
    cell_resp=cell_resp(:,adapting_frame:end);
    else
    cell_resp=cell_resp(:,adapting_frame:ending_frame);
    end
    cell_resp_dim=size(cell_resp);
    write_LSstack_fast_float(fullfile(input_dir,'cell_resp.stackf'),cell_resp);
    save(fullfile(input_dir,'cell_resp_dim'),'cell_resp_dim');
    
    %% low-cut timecourse of cells exponential fit
    disp('low-cut timecourse of cells exponential fit....');
    tmp=sort(double(background(:)),'ascend');
    back_value=mean(tmp(1:round(length(tmp)/20)));
    
   
    
    totcell=size(cell_resp,1);
    totaltime=size(cell_resp,2);
    nrep=floor(totaltime/sample_len);
    baselines=zeros(totcell,nrep);


    bottomlen=round(sample_len/5); %calculating F0, bottom 20%

    for j=1:nrep
        tmp=cell_resp(:,(j-1)*sample_len+(1:sample_len));    
        tmp=sort(tmp,2);  
        baselines(:,j)=mean(tmp(:,1:bottomlen), 2)-back_value;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%
    baselines(baselines<0) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%

    fitbaseline=log(baselines./repmat(mean(baselines(:,1:2),2),[1 nrep]));
    xs=round(sample_len/2)+(0:sample_len:((nrep-1)*sample_len));

    baseline2=zeros(totcell,totaltime);
    baseline3 = zeros(totcell,totaltime);

    cs=zeros(totcell,1);
    % maybe parfor
    % keep this chunked style, change into sliding window later
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!@!
    for i=1:totcell
        cf=polyfit(xs,fitbaseline(i,:),1);
        cs(i)=cf(1);
        baseline2(i,:)=exp(cf(1)*(1:totaltime))*mean(baselines(i,1:2));
        
%        baseline3(i,:) = interp1((1:66)/66*36000,baselines(i,:),(1:36000));
    end

    cell_resp2=(cell_resp-back_value)./baseline2;
%     cell_resp2=(cell_resp-back_value)./baseline3;
%     save(fullfile(input_dir,'baseline_fit.mat'),'xs','baselines','baseline2','back_value');
    

    cell_resp_dim=size(cell_resp2);
    write_LSstack_fast_float(fullfile(input_dir,'cell_resp_lowcut.stackf'),cell_resp2);
    save(fullfile(input_dir,'cell_resp_dim_lowcut'),'cell_resp_dim');
%   cell_resp2_ori=cell_resp2;
    %% detrending or loading
    disp('detrending....');
    tic
    cell_resp_old= cell_resp2;
    CR = zeros(size(cell_resp2));
    tmax=size(cell_resp2,2);
    ncells=size(cell_resp2,1);
    poolnumber = 40;
   parpool(poolnumber);
   
    parfor i=1:ncells
        CR(i,:) = cell_resp(i,:) ./ rolling_percentile_filter(cell_resp(i,:),600,15);
%         CR(i,:) = cell_resp2(i,:) ./ rolling_percentile_filter(cell_resp2(i,:),300,15);
%         cr = cell_resp2(i,:);
%         crd = 0*cr;
%         for j=1:100:tmax
%             if j<=150
%                 tlim1 = 1;
%                 tlim2 = 300;
%             elseif j>tmax-150
%                 tlim1 = tmax-300;
%                 tlim2 = tmax;
%             else
%                 tlim1 = j-150;
%                 tlim2 = j+150;
%             end
%             crr = cr(tlim1:tlim2);
%             crd(max(1,j-50):min(tmax,j+50)) = prctile(crr,15);
%         end
%         crds = crd;  %smooth(crd,120); % using smooth is bad, because of boundary effects
%         if mod(i,5000)==0   % for manual control! to make sure it works well
%             plot(cr)
%             hold on
%             plot(crd,'r')
%             plot(crds,'g')
%             plot(cr-crds,'k')
%             hold off
%             title(num2str(i))
%             %         waitforbuttonpress
%             drawnow
%         end
%         if mod(i,1000)==0,i,end
%         CR(i,:) = cr-crds+1;  % take the z-score - optional. Can also do CR(i,:) = cr-crds;
%         %     CR(i,:) = zscore(cr-crds);  % take the z-score - optional. Can also do CR(i,:) = cr-crds;
    end
 %   delete(gcp('nocreate'));
    cell_resp2 = single(CR);
%     write_LSstack_fast_float(fullfile(input_dir,'CR_detrend.stackf'),cell_resp2);
    toc
    %% removing double_counted cells
    
%     if dim(3)>1
%         totrepeats=floor(cell_resp_dim(2)/ sample_len);
%         cell_resp2=cell_resp2(:,1:totrepeats*sample_len);
%         %cell_resp3=cell_resp2(:,1:20*sample_len);
% 
% 
%         disp('Removing multiple-counted cells...');
%         %cell_resp2_ncorr=cell_resp3-repmat(squeeze(mean(reshape(cell_resp3,[cell_resp_dim(1) sample_len 20]),3)),[1 20]);
%         cell_resp2_ncorr=cell_resp2-repmat(squeeze(mean(reshape(cell_resp2,[cell_resp_dim(1) sample_len totrepeats]),3)),[1 totrepeats]);
% 
%         [r,c]=find(ones(5,5));
%         r=r-3;c=c-3;
%         moveinds=c*dim(1)+r;
% 
%         cell_remove=zeros(1,length(cell_info));
%         for z=1:dim(3)-1
% 
%             zplane1=zeros(dim(1),dim(2));
%             zplane2=zeros(dim(1),dim(2));
% 
%             cinds1=find([cell_info.slice]==z  & cell_remove==0);
%             cinds2=find([cell_info.slice]==z+1 & cell_remove==0);
% 
%             for i=1:length(cinds1)
%                 cellcenter=dim(1)*(cell_info(cinds1(i)).center(2)-1)+cell_info(cinds1(i)).center(1);
%                 mask=find(cellcenter+moveinds>0 & cellcenter+moveinds<dim(1)*dim(2));
%                 zplane1(cellcenter+moveinds(mask))=cinds1(i);
%             end
% 
%             for i=1:length(cinds2)
%                 cellcenter=dim(1)*(cell_info(cinds2(i)).center(2)-1)+cell_info(cinds2(i)).center(1);
%                 mask=find(cellcenter+moveinds>0 & cellcenter+moveinds<dim(1)*dim(2));
%                 zplane2(cellcenter+moveinds(mask))=cinds2(i);
%             end 
% 
%             zplane3=(zplane1>0).*(zplane2>0);
% 
%             CC=bwconncomp(zplane3);
%             rlist=zeros(1,CC.NumObjects);
%             clist=zeros(1,CC.NumObjects);
% 
%             for i=1:CC.NumObjects
%                 tmp=CC.PixelIdxList{i};
%                 cnum1=max(zplane1(tmp));
%                 cnum2=max(zplane2(tmp));
%                 clist(i)=cnum2;
% 
%                 r=corrcoef_pair_mex(double(cell_resp2_ncorr(cnum1,:)),double(cell_resp2_ncorr(cnum2,:)));
%                 rlist(i)=r;
% 
%             end
% 
%             removelist=find(rlist>corr_thre);
%             cell_remove(clist(removelist))=1;
%         end
%         keep_inds=find(cell_remove==0);
%         remove_inds=find(cell_remove>0);
%         disp([num2str(length(remove_inds)), ' cells removed by multi-counting'])
% 
%         cell_info=cell_info(keep_inds);
%         cell_resp2=cell_resp2_ori(keep_inds,:);
% 
%     end
    
    %%  calculating motion confidence

    disp('Removing deformed cells....');
%     dists=zeros(1,length(cell_info));
% 
%     for c=1:length(cell_info)
% 
%         z=cell_info(c).slice;
%         yx=cell_info(c).center;
% 
%         grid_y=mod(motion_param(z).indslist,dim(1));grid_y(grid_y==0)=dim(1);
%         grid_x=ceil(motion_param(z).indslist/dim(1));
% 
%         distance=sqrt((grid_y-yx(1)).^2 + (grid_x-yx(2)).^2);
%         [V,inds]=sort(distance,'ascend');
% 
%         cell_info(c).motion=mean(motion_param(z).tilt_med(inds(1:min([5,length(inds)])),:),1);  
% 
%     end
% 
%     motion_matrix=[cell_info.motion];
%     motion_y=motion_matrix(1:3:end);
%     motion_x=motion_matrix(2:3:end);
%     motion_z=motion_matrix(3:3:end);
%     remove_inds=find(abs(motion_y)>move_thre | abs(motion_x)>move_thre | abs(motion_z)>move_thre);
%     disp([num2str(length(remove_inds)), ' cells removed by motion_assessment'])
% 
%     keep_inds=setdiff(1:length(cell_info),remove_inds);
% 
%     cell_info=cell_info(keep_inds);
%     cell_resp2=cell_resp2(keep_inds,:);
% 
% 
    cell_resp_dim=size(cell_resp2);
    
%    imwrite(ave_stack,fullfile(input_dir,'ave.tif'));
    
    save(fullfile(input_dir,'cell_resp_dim_processed'),'cell_resp_dim');
    write_LSstack_fast_float(fullfile(input_dir,'cell_resp_processed.stackf'),cell_resp2);
    save(fullfile(input_dir,'cell_info_processed.mat'),'cell_info');
    
    delete(gcp('nocreate'));
%end
toc