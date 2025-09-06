clear all;close all;
tic

input_dirs=struct;

input_dirs(1).fname = 'E:\shy\20250507_3_1_6dpf_VODL2_20250507_133458\registered';
ending_frame = 0;

for fnum=1:length(input_dirs)
    
    input_dir=input_dirs(fnum).fname;
    disp(input_dir);
    sdim=double(read_LSstack_size(fullfile(input_dir,'Stack dimensions.log')));
%%%%%%%%%%%
    %sdim(3) = 1;
%%%%%%%%%%%
    radius=30;
%         radius=50;
    interval=30;
    numpool=4;
    sample_rate =input('input the sample_rate:');
    zcycle=20;
%     zcycle=120*sample_rate;
    xypixeldist=0.406;
    zpixeldist=5;
%     brightness_thre=104;
    brightness_thre =input('input the brightness_thre:');

    rwidth=radius*2+1;
    zlist=1:sdim(3);
    ref_stacks=struct;
    refstack=zeros(sdim(1),sdim(2),sdim(3));
    
%    pool_num = min(floor(sdim(3)/4)+3,12);

%%%%%%%
    pool_num = sdim(3);
%%%%%%%%
    parpool(pool_num);
    for i=1:sdim(3)
        ref_stacks(i).sdim=sdim;
    end
    
    %% reading files from PlaneXX.stack files

    parfor zz=1:length(zlist)
    %for zz=1:length(zlist)
        fname=['Plane',num2str(zlist(zz),'%.2d'),'.stack'];
        d=ref_stacks(zz).sdim;
        fstack=read_LSstack_fast2(fullfile(input_dir,fname),[d(1) d(2)],[1 zcycle]);
        refstack(:,:,zz) = create_zslice_ave_mex(uint16(fstack), int32([d(1) d(2)]),int32((1:zcycle)));  
    end
    

    for i=1:sdim(3)    
        k=(-2:2)+i;
        inds=find(k>0 & k<=sdim(3));
        ref_stacks(i).refstack=refstack(:,:,k(inds));
    end

    %%


    output=struct;
    motion_param=struct;
    move_tcourse=struct;
    
    parfor zz=zlist
    %for zz=zlist

        fname=['Plane',num2str(zlist(zz),'%.2d'),'.stack'];
        disp(['Calculating deformation of plane',num2str(zz)]);
        d=ref_stacks(zz).sdim;
        d2=double(read_LSstack_info(fullfile(input_dir,fname),[d(1) d(2)]));
        
        if ending_frame ~= 0
            totcycle=floor(ending_frame/zcycle);
        else
        totcycle=floor(d2(3)/zcycle);
        end
        stack=zeros(d(1),d(2),totcycle);
        move_tcourse(zz).tcourse=((1:totcycle)-1)*zcycle+round(zcycle/2);

        for z=1:totcycle
            fstack=read_LSstack_fast2(fullfile(input_dir,fname),[d(1) d(2)],[((z-1)*zcycle+1) z*zcycle]);
            stack(:,:,z)=create_zslice_ave_mex(uint16(fstack),int32([d(1) d(2)]),int32((1:zcycle)));
        end


        sampleimg=double(mean(stack(:,:,1),3));
        tcourse=1:totcycle;

        regimg=repmat(imNormalize99(sampleimg),[1 1 3]);
        regimg2=regimg;
        pgrid=radius;
        ygrid_num=floor((d(1)-rwidth*2)/pgrid);
        xgrid_num=floor((d(2)-rwidth*2)/pgrid);
        npoints=ygrid_num*xgrid_num;
        indslist=zeros(npoints,2);

        [r, c]=find(ones(rwidth)>0);
        r=r-radius-1;c=c-radius-1;
        inds=(c-1)*d(1)+r;

        cc=1;
        for x=1:xgrid_num
            for y=1:ygrid_num
                rinds=(rwidth+(x-1)*pgrid)*d(1)+rwidth+(y-1)*pgrid+1;
                indslist(cc,1)=rinds;
                indslist(cc,2)=mean(sampleimg(rinds+inds));
                cc=cc+1;
            end
        end
        
        %%
        
        thre_inds=find(indslist(:,2)>brightness_thre);
        indslist2=indslist(thre_inds,1);

        for i=1:length(indslist2)
            rinds=indslist2(i);
            regimg(rinds+inds)=1;
            regimg(rinds+inds+d(1)*d(2))=0;
            regimg(rinds+inds+d(1)*d(2)*2)=0;             
        end


        simg=double(stack(:,:,1));
        rs=zeros(length(indslist2),3,totcycle);  
        tilt=zeros(length(indslist2),3);    
        tilt_med=zeros(length(indslist2),3); 
        zmove=-2:2;

        for i=1:length(indslist2)
            rinds=indslist2(i);
            narray=zeros(length(r),1);
            source=fft2(reshape(simg(rinds+inds),[rwidth rwidth]));

            for j=1:totcycle
                target=reshape(stack(rinds+inds+(j-1)*d(1)*d(2)),[rwidth rwidth]);

                buff=fftshift(ifft2(source.*conj(fft2(target))));
                [~,maxinds] = max(abs(buff(:)));

                a=-(mod(maxinds,rwidth)-radius-1);
                if a==-radius-1; a=radius;end

                rs(i,1,j)=a; 
                rs(i,2,j)=-(ceil(maxinds/rwidth)-radius-1);   

                moveinds=-rs(i,2,j)*d(1)-rs(i,1,j);

                tt=zeros(1,5);
                zinds=find(zmove+zz > 0 & zmove+zz<=d(3) );
                for k=1:length(zinds)                  
                         ztargetlist=ref_stacks(zz).refstack(rinds+moveinds+inds+(k-1)*d(1)*d(2));                     
                         tt(zinds(k))=corrcoef_pair_mex(target(:),ztargetlist(:));                   
                end



                [~,zmaxinds] = max(tt(zinds));
                rs(i,3,j)=zmove(zinds(zmaxinds));

            end

            py=polyfit(tcourse',squeeze(rs(i,1,:)),1);
            px=polyfit(tcourse',squeeze(rs(i,2,:)),1);
            pz=polyfit(tcourse',squeeze(rs(i,3,:)),1);

            tilt(i,1)=py(1)*(length(tcourse)-1);
            tilt(i,2)=px(1)*(length(tcourse)-1);
            tilt(i,3)=pz(1)*(length(tcourse)-1);

        end

        %%
        move_tcourse(zz).rs_ave_xy = mean(sqrt(squeeze(rs(:,1,:)).^2 + squeeze(rs(:,2,:)).^2));
        move_tcourse(zz).rs_std_xy =  std(sqrt(squeeze(rs(:,1,:)).^2 + squeeze(rs(:,2,:)).^2));
        
        move_tcourse(zz).rs_ave_z  = mean(squeeze(rs(:,3,:)));
        move_tcourse(zz).rs_std_z  = std(squeeze(rs(:,3,:)));
        
        
        
        
        [r, c]=find(ones(3,3));
        median_inds=(c-2)*ygrid_num+(r-2);

        for jj=1:3

            move_matrix=zeros(ygrid_num,xgrid_num);
            ones_matrix=zeros(ygrid_num,xgrid_num);

            for ii=1:length(indslist2)
                move_matrix(thre_inds(ii))=tilt(ii,jj);
                ones_matrix(thre_inds(ii))=1;
            end

            for ii=1:length(indslist2)
                tmp = thre_inds(ii)+median_inds;
                inds1=thre_inds(ii)+median_inds(tmp>0 & tmp<=xgrid_num*ygrid_num);
                inds2=inds1(ones_matrix(inds1)>0);
                tilt_med(ii,jj)=median(move_matrix(inds2));
            end
        end
       
        

        output(zz).masks=regimg;
        output(zz).tilt=tilt;
        output(zz).tilt_med=tilt_med;
        output(zz).indslist2=indslist2;
        output(zz).regimg2=regimg2;

        motion_param(zz).tilt=tilt;
        motion_param(zz).tilt_med=tilt_med;
        motion_param(zz).indslist=indslist2;
        motion_param(zz).xymove_av=mean(sqrt(tilt_med(:,1).^2 + tilt_med(:,2).^2))*xypixeldist;
        motion_param(zz).xymove_sd=std(sqrt(tilt_med(:,1).^2 + tilt_med(:,2).^2)*xypixeldist,[],1);
        motion_param(zz).zmove_av=mean(abs(tilt_med(:,3)))*zpixeldist;
        motion_param(zz).zmove_sd=std(abs(tilt_med(:,3))*zpixeldist);


   

    end
    delete(gcp('nocreate'));


    %%
    
    h1=figure(1);set(h1,'Position',[300 400 500 500]);
  
    movegraph=struct;
    xtcourse=move_tcourse(1).tcourse;

    for zz=1:length(zlist)
        if length(xtcourse)==size(move_tcourse(zz).rs_ave_xy,2)
            clf(h1);
            errorbar(xtcourse,move_tcourse(zz).rs_ave_xy*xypixeldist, move_tcourse(zz).rs_std_xy*xypixeldist,'mo-','linewidth',2);hold on;
            errorbar(xtcourse,move_tcourse(zz).rs_ave_z*zpixeldist , move_tcourse(zz).rs_std_z*zpixeldist ,'co-','linewidth',2);hold off;
            ylim([-10 10]);xlim([0 max(xtcourse)]);
            title({['Plane',num2str(zz),': motion timecourse'],'magenta=xy, cyan=z'});
            CC=getframe(h1);    
            movegraph(zz).graph=CC.cdata;
        end
    end

    odim=size(movegraph(1).graph);
%     graphstack=zeros(odim(1),odim(2),length(zlist),3);

%     for i=1:length(zlist)
%         if length(xtcourse)==size(move_tcourse(i).rs_ave_xy,2);
%         graphstack(:,:,i,:)=reshape(movegraph(i).graph,[odim(1) odim(2) 1 3]);
%         end
%     end
% 
%     writetiff(uint8(graphstack),fullfile(input_dir,'motion_tcourse.tif'));
    for i = 1:length(zlist)
        imwrite(squeeze(reshape(movegraph(i).graph,[odim(1) odim(2) 1 3])),fullfile(input_dir,'motion_tcourse.tif'),'WriteMode','append');
    end
    
    %%
    h2=figure(2);set(h2,'Position',[150 150 400 800]);
    dim=size(output(1).masks);
    colorlist=[0 0 1;
               0 1 1;
               0 1 0;
               1 1 0;
               1 0 0;];

    zmovelist=[-2 -1 0 1 2];


    for zz=1:length(zlist)
        clf(h2);
        
        indslist2=output(zz).indslist2;
        tilt=output(zz).tilt_med;
        regimg2=output(zz).regimg2;
        ha=axes; set(ha,'Position',[0 0 1 1]);
        
        image(regimg2);hold on;    
        hq=quiver(ceil(indslist2/dim(1)),mod(indslist2,dim(1)),tilt(:,2),tilt(:,1),0,'linewidth',2,'Color',[1 0 0]);

        hold off;

        hU = get(hq,'UData') ;
        hV = get(hq,'VData') ;
        set(hq,'UData',10*hU,'VData',10*hV)
        axis off;

        CC=getframe(h2);    
        output(zz).arrows=CC.cdata;
        cla;
        
        set(ha,'Position',[0 0 1 1]);
        image(regimg2);hold on;

        for i=1:5
            inds=find(round(tilt(:,3))==zmovelist(i));
            if ~isempty(inds)
                for j=1:length(inds)
                    rectangle('Position',[ceil(indslist2(inds(j))/dim(1)),mod(indslist2(inds(j)),dim(1)), 20, 20],'FaceColor',colorlist(i,:));
                end
            end
        end

        hold off;
        axis off;

        CC=getframe(h2);    
        output(zz).zmotion=CC.cdata;

    end

    odim=size(output(1).arrows);
%     outstack=zeros(odim(1),odim(2)*2+1,length(zlist),3);

%     for i=1:length(zlist)
%         outstack(:,1:odim(2),i,:)=reshape(output(i).arrows,[odim(1) odim(2) 1 3]);
%         outstack(:,(1:odim(2))+odim(2)+1,i,:)=reshape(output(i).zmotion,[odim(1) odim(2) 1 3]);
%     end
% 
%     writetiff(uint8(outstack),fullfile(input_dir,'motion.tif'));
    for i = 1:length(zlist)
        imwrite(squeeze(reshape(output(i).zmotion,[odim(1) odim(2) 1 3])),fullfile(input_dir,'motion.tif'),'WriteMode','append');
    end


    save(fullfile(input_dir,'motion_param.mat'),'motion_param');
    %%

    h3=figure(3);
    set(h3,'Position',[100 100 1500 750]);
    subplot(1,2,1);
    for zz=1:length(zlist)
        indslist=output(zz).indslist2;
        tilt=output(zz).tilt_med;
        plot(ones(1,length(output(zz).indslist2))*zz,sqrt(tilt(:,1).^2+tilt(:,2).^2)*xypixeldist,'.');hold on;
    end
    errorbar([motion_param(zlist).xymove_av],[motion_param(zlist).xymove_sd],'r','LineWidth',2,'LineStyle','none');
    scatter(zlist,[motion_param.xymove_av],'ro','fill');hold off;
    xlim([min(zlist)-1 max(zlist)+1]);ylim([-1 5]);title('XY motion');

    subplot(1,2,2);
    for zz=1:length(zlist)
        indslist=output(zz).indslist2;
        tilt=output(zz).tilt_med;
        plot(ones(1,length(output(zz).indslist2))*zz,abs(tilt(:,3))*zpixeldist,'.');hold on;
    end
    errorbar([motion_param.zmove_av],[motion_param.zmove_sd],'r','LineWidth',2,'LineStyle', 'none');
    scatter(zlist,[motion_param.zmove_av],'ro','fill');hold off;
    xlim([min(zlist)-1 max(zlist)+1]);ylim([-1 10]);title('Z motion');


    set(h3,'PaperPositionMode','auto'); 
    saveas(h3,fullfile(input_dir,['motion_graph.tif']),'tif');
    saveas(h3,fullfile(input_dir,['motion_graph.eps']),'eps');
    
end
toc










