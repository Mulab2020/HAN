clear all;close all;
tic
img_dir='E:\shy\20250507_3_1_6dpf_VODL2_20250507_133458\registered';

dim=read_LSstack_size(fullfile(img_dir,'\Stack dimensions.log'));
    
%%stack_av=read_LSstack_fast_float(fullfile(img_dir,'\ave.stackf'),dim);

info = imfinfo([img_dir, '\ave.tif']);
numberOfPages = length(info);
for k = 1 : numberOfPages
      stack_av(:,:,k)  = imread([img_dir, '\ave.tif'],k);
end



%%
button = 2;

while button == 2
    
    br_threshold=input('input the brightness threshold:');
    
    dim=size(stack_av);
    if length(dim) == 2
        dim = [dim, 1];
    end
    cell_info=struct;
    cell_color=zeros([dim(1) dim(2)*2+1 dim(3) 3],'uint8');
    cellnum=0;
    
    cont_threshold=4;
    cell_rad=5;%6 for cytosolic; 5 for nuclear; 12 for single plane
    
    
    zlist=1:size(stack_av,3);
    %%  setset masking filter;
    
    ave_rad=round(cell_rad/2)+1;
    [avedisk,  ave_se, r1, c1, maskinds]=make_recog_disk(round(cell_rad/2)+1,dim);
    [maxdisk,  max_se, r2, c2, maxinds]=make_recog_disk(cell_rad+2,dim);
    onedisk=makeDisk2(3,7);
    one_se=strel(onedisk);
    
    %% for rank calculation template
    
    r=cell_rad*2;
    dimp=[dim(1)+r*2 dim(2)+r*2];
    oop=zeros(dimp);
    oop(r+1:end-r,r+1:end-r)=1;
    one_inds=find(oop);
    
    [mdisk,  ~, ~, ~, rankinds]=make_recog_disk(r,dimp);
    rank_ones=double(maskones2D_mex(int32([dim(1) dim(2)]),int32(mdisk),int32(size(mdisk))))';
    
    %%  recognize cells
    
    allmask=zeros(dim(1),dim(2));
    imlen=dim(1)*dim(2);
    
    for z=zlist
        
        im=stack_av(:,:,z);
        allmask(:)=0;
        contimage = local_contrast_mex(single(im),int32(32),single(cont_threshold));
        contimage = imdilate(imerode(contimage,one_se),one_se);
        contimage = contimage.*uint8((im>br_threshold));
        
        candidates=find(contimage);
        
        if ~isempty(candidates)
            %% recognizing cells in the first round
            
            imrank = calc_rank_simple2(im,rank_ones,cell_rad*2,oop,one_inds,rankinds,candidates);
            aveimg = double(local_average_mex(single(imrank),int32(c1),int32(r1),int32(candidates)));
            maximg = double(local_max_mex(single(aveimg),int32(c2),int32(r2),int32(candidates)));
            
            inds=find(maximg(candidates)>0 & aveimg(candidates) >0.4);
            mask2=zeros(dim(1),dim(2));
            
            for i=1:length(inds)
                cinds=candidates(inds(i))+maskinds;
                cinds(cinds > imlen | cinds<1)=[];
                mask2(cinds)=1;
            end
            
            allmask=mask2;
            
            %% recognizing cells in the second round
            
            mask3=ones(size(im),'uint8')-imdilate(uint8(allmask),max_se);
            mask3 = imdilate(imerode(mask3,one_se),one_se);
            candidates2=candidates(mask3(candidates)>0);
            
            imrank2 = calc_rank_simple2(im,rank_ones,cell_rad*2,oop,one_inds,rankinds,candidates2);
            aveimg2 = double(local_average_mex(single(imrank2),int32(c1),int32(r1),int32(candidates2)));
            maximg2 = double(local_max_mex(single(aveimg2),int32(c1),int32(r1),int32(candidates2)));
            
            inds=find(maximg2(candidates2)>0 & aveimg2(candidates2) >0.4);
            mask2=zeros(dim(1),dim(2));
            
            for i=1:length(inds)
                cinds=candidates2(inds(i))+maskinds;
                cinds(cinds > imlen | cinds<1)=[];
                mask2(cinds)=1;
            end
            
            allmask=allmask+mask2;
            
            %% create each cell ROIs
            
            [celllabel, totcell]=bwlabel(allmask,8);
            if totcell>0
                cell_info=create_cell_info_fish(cell_info,celllabel, totcell,z);
            end
        else
            
            totcell=0;
            celllabel=zeros(size(im));
        end
        
        cell_color(:,:,z,:)=reshape(imMask2D_fish(im,celllabel,candidates),[dim(1) dim(2)*2+1 1 3]);
        
        cellnum=cellnum+totcell;
        disp(num2str(z));
        
    end
    
    figure('position',[300 300 1000 800]);
    for z = zlist
        imagesc(squeeze(cell_color(:,:,z,:)));
        [~,~,button] = ginput(1);
        if button ~= 1
            break;
        end
    end
end


%%
cellmask_name = ['cellmask_',num2str(br_threshold),'_',num2str(cont_threshold),'_',num2str(cell_rad),'.tif'];
write_colortiff_mex(fullfile(img_dir,cellmask_name),cell_color, int32(size(cell_color)));
save(fullfile(img_dir,'cell_info.mat'),'cell_info');
toc
disp(['total cellnum: ',num2str(cellnum)]);
