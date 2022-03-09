function [EvaLum] = GPconstancy_GI(input_im,param)
%GPCONSTANCY_GI Summary of this function goes here
    mask=param.mask;
    numGPs=param.numGPs;
    delta_threshold=param.delta_threshold;
    real_rgb=param.real_rgb;
    i=param.runtime.i;
    name_img=param.runtime.name_img;
    

    img_column=reshape(input_im,[],3);
    r=input_im(:,:,1); g=input_im(:,:,2); b=input_im(:,:,3);
    
    %denoise
    hh = fspecial('average',[7 7]);
    r = imfilter(r,hh,'circular');g = imfilter(g,hh,'circular');b = imfilter(b,hh,'circular');
    
    %mask 0 elements
    mask=mask | (r==0) | (g==0) | (b==0);
    r(r==0)=eps;  g(g==0)=eps;  b(b==0)=eps;norm1=r+g+b;

    %mask low contrast pixels     
    delta_r=DerivGauss(r,.5); delta_g=DerivGauss(g,.5); delta_b=DerivGauss(b,.5);
    mask=mask | delta_r<=delta_threshold & delta_g<=delta_threshold & delta_b<=delta_threshold;

   
    log_r=log(r)-log(norm1); log_b=log(b)-log(norm1);
    
    delta_log_r=DerivGauss(log_r,.5); 
    delta_log_b=DerivGauss(log_b,.5);
    mask=mask | (delta_log_r==Inf) | (delta_log_b==Inf);
    

    data=[(delta_log_r(:)),(delta_log_b(:))];%,Mr(:)-norm1_M(:),Mg(:)-norm1_M(:),Mb(:)-norm1_M(:)];
    mink_norm=2;
    norm2_data=power(sum(power(data,mink_norm),2),1/mink_norm);
    map_uniquelight=reshape(norm2_data,size(delta_log_r));
    
    %(1)mask out some pixels
    map_uniquelight(mask==1)=max(map_uniquelight(:));
    
    %(2)denoise 
    hh = fspecial('average',[7 7]);
    map_uniquelight = imfilter(map_uniquelight,hh,'circular');
    

    %filter by map_uniquelight
    Greyidx_unique=map_uniquelight;
    
    sort_unique=sort(Greyidx_unique(:));
    Gidx_unique = zeros(size(Greyidx_unique));
    
    
    Gidx_unique(Greyidx_unique<=sort_unique(floor(numGPs))) = 1;
    choosen_pixels=img_column(Gidx_unique==1,:);
    EvaLum=normr(mean((choosen_pixels),1));
    i=1;
    if param.prior.use
        uv_chosen=RgbToUv(choosen_pixels')';
        u_ind=floor((uv_chosen(:,1)-(-1))/.04);
        v_ind=floor((uv_chosen(:,2)-(-1))/.04);
        u_ind=int16(u_ind); v_ind=int16(v_ind);
        
        ind_camera=param.prior.ind_camera(param.runtime.i);
        prior=param.prior.prior_map{ind_camera};
        sum_posterior=zeros(1,3);
        for ind=1:1:length(u_ind)
           
           sum_posterior=sum_posterior+(choosen_pixels(ind,:))...
               .*prior(max(min(u_ind(ind),size(prior,1)),1),...
               max(min(v_ind(ind),size(prior,2)),1)); 
        end
        EvaLum=normr(sum_posterior);

    end
    

    if param.visualization.greypixel_comparison
        input_im=input_im/max(input_im(:));
        img_column=reshape(input_im,[],3);
        im_gt=input_im;
        im_gt(:,:,1)=im_gt(:,:,1)/(real_rgb(i,1)*sqrt(3));
        im_gt(:,:,2)=im_gt(:,:,2)/(real_rgb(i,2)*sqrt(3));
        im_gt(:,:,3)=im_gt(:,:,3)/(real_rgb(i,3)*sqrt(3));
        im_gt=im_gt/max(im_gt(:));
        I_input = uint8(255 * rgb2srgb(max(0, min(1, input_im)),false));
        I_true = uint8(255 * rgb2srgb(max(0, min(1, im_gt)),false));
        
        Gidx_dgp=Gidx_unique;
        Gidx_dgp=repmat(Gidx_dgp,[1,1,3]);
        choosen_dgp=uint8(double(I_true).*Gidx_dgp);
        
        %old gray pixel method
        map_old_greyness = GetGreyidx(input_im,'GPedge',0.5);
        sort_old=sort(map_old_greyness(:));
        Gidx_gp = zeros(size(map_old_greyness));
        Gidx_gp(map_old_greyness<=sort_old(floor(numGPs))) = 1;
        
        Gidx_gp=repmat(Gidx_gp,[1,1,3]);
        choosen_gp=uint8(double(I_true).*Gidx_gp);
        
        comp=cat(2,I_input,255*ones(size(I_input,1),10,3),I_true,255*ones(size(I_input,1),10,3),...
            choosen_dgp,255*ones(size(I_input,1),10,3),choosen_gp);
        name_comp=sprintf(['comp_' name_img]);
        imwrite(comp,fullfile(param.visualization.comp_dir,name_comp));
    end
    
    if param.visualization.bigsequenceimage
        %ColormapUsa;
        Map       = jet(255);
        sequence_1=rgb2srgb(input_im,false);
        
%         sequence_2=(delta_log_r+0*max(delta_log_r(:)))/(1*max(delta_log_r(:)));
%         sequence_2=uint8(sequence_2*255.0);
%         sequence_2=ind2rgb(sequence_2,Map);
%         
%         sequence_3=(delta_log_b+0*max(delta_log_b(:)))/(1*max(delta_log_b(:)));
%         sequence_3=uint8(sequence_3*255.0);
%         sequence_3=ind2rgb(sequence_3,Map);
        
        norm2_data=power(sum(power(data,mink_norm),2),1/mink_norm);
        sequence_4=reshape(norm2_data,size(delta_log_r));
        sequence_4=(sequence_4+0*max(sequence_4(:)))/(1*max(sequence_4(:)));
        sequence_4=uint8(sequence_4*255.0);
        sequence_4=ind2rgb(sequence_4,Map);
        
        norm2_data(mask==1)=max(norm2_data(:));
        sequence_5=reshape(norm2_data,size(delta_log_r));
        sequence_5=(sequence_5+0*max(sequence_5(:)))/(1*max(sequence_5(:)));
        sequence_5=uint8(sequence_5*255.0);
        sequence_5=ind2rgb(sequence_5,Map);
        
        sequence_6=(map_uniquelight)/(max(map_uniquelight(:)));
        sequence_6=uint8(sequence_6*255.0);
        sequence_6= ind2rgb(sequence_6, Map);
        
        input_im=input_im/max(input_im(:));
        img_column=reshape(input_im,[],3);
        im_gt=input_im;
        im_gt(:,:,1)=im_gt(:,:,1)/(real_rgb(i,1)*sqrt(3));
        im_gt(:,:,2)=im_gt(:,:,2)/(real_rgb(i,2)*sqrt(3));
        im_gt(:,:,3)=im_gt(:,:,3)/(real_rgb(i,3)*sqrt(3));
        im_gt=im_gt/max(im_gt(:));
        I_input = uint8(255 * rgb2srgb(max(0, min(1, input_im)),false));
        I_true = uint8(255 * rgb2srgb(max(0, min(1, im_gt)),false));
        
        Gidx_unique(Greyidx_unique<=sort_unique(floor(10*numGPs))) = 1;
        Gidx_dgp=Gidx_unique;
        Gidx_dgp=repmat(Gidx_dgp,[1,1,3]);
        choosen_dgp=uint8(double(I_input).*Gidx_dgp);
        
        pred_bar=ones(size(sequence_1,1),200,3);
        pred_bar(:,:,1)=EvaLum(1); pred_bar(:,:,2)=EvaLum(2);pred_bar(:,:,3)=EvaLum(3);
        
        gt_bar=ones(size(sequence_1,1),200,3);
        gt_bar(:,:,1)=real_rgb(i,1); gt_bar(:,:,2)=real_rgb(i,2);gt_bar(:,:,3)=real_rgb(i,3);

        bigsequence=cat(2,sequence_1,ones(size(sequence_1,1),10,3),...
            sequence_4,ones(size(sequence_1,1),10,3),sequence_5,ones(size(sequence_1,1),10,3),...
            sequence_6,ones(size(sequence_1,1),10,3),...
            double(choosen_dgp)/255.0...
             ,ones(size(sequence_1,1),20,3), pred_bar, ones(size(sequence_1,1),20,3),gt_bar...
            );

        name_sequence=sprintf(['sequence_' name_img])
        imwrite(bigsequence, name_sequence);
    end

end

