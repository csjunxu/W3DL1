clear;
Original_image_dir  =    'C:\Users\csjunxu\Desktop\Projects\WODL\20images\';
% Original_image_dir  =    'C:\Users\csjunxu\Desktop\Projects\WODL\20newimages\';
Sdir = regexp(Original_image_dir, '\', 'split');
fpath = fullfile(Original_image_dir, '*.png');
im_dir  = dir(fpath);
im_num = length(im_dir);


for nSig = 40 % [20 40 60 80 100]
    if nSig<=20
        Par.ps = 8;
        Par.win = 30;
        Par.step = 3;
        Par.outerIter = 8;
        Par.nlspini = 70;
        Par.nlspgap = 10;
        Par.delta = 0.05;
        Par.lambda2 = 0.78;
    elseif      nSig<=40
        Par.ps = 8;
        Par.win = 35;
        Par.step = 3;
        Par.outerIter = 12;
        Par.nlspini = 90;
        Par.nlspgap = 10;
        Par.delta = 0.06;
        Par.lambda2 = 0.8;
    elseif      nSig<=60
        Par.ps = 8;
        Par.win = 35;
        Par.step = 3;
        Par.outerIter = 12;
        Par.nlspini = 120;
        Par.nlspgap = 15;
        Par.delta = 0.05;
        Par.lambda2 = 0.7;
    elseif      nSig<=80
        Par.ps = 9;
        Par.win = 40;
        Par.step = 4;
        Par.outerIter = 14;
        Par.nlspini = 140;
        Par.nlspgap = 15;
        Par.delta = 0.05;
        Par.lambda2 = 0.7;
    else
        Par.ps = 9;
        Par.win = 40;
        Par.step = 4;
        Par.outerIter = 14;
        Par.nlspini = 140;
        Par.nlspgap = 15;
        Par.delta = 0.05;
        Par.lambda2 = 0.7;
    end
    Par.lambda1 = 0.1;
    Par.maxIter = 100;
    Par.innerIter = 2;
    Par.model    = 1;
    % record all the results in each iteration
    Par.PSNR = zeros(Par.outerIter, im_num, 'double');
    Par.SSIM = zeros(Par.outerIter, im_num, 'double');
    T512 = [];
    T256 = [];
    for i = 4% 1:im_num
        Par.nlsp = Par.nlspini;  % number of non-local patches
        Par.image = i;
        Par.nSig = nSig/255;
        Par.I =  im2double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
        S = regexp(im_dir(i).name, '\.', 'split');
        randn('seed',0);
        Par.nim =   Par.I + Par.nSig*randn(size(Par.I));
        %
        fprintf('%s :\n',im_dir(i).name);
        PSNR =   csnr( Par.nim*255, Par.I*255, 0, 0 );
        SSIM      =  cal_ssim( Par.nim*255, Par.I*255, 0, 0 );
        fprintf('The initial value of PSNR = %2.4f, SSIM = %2.4f \n', PSNR,SSIM);
        %
        time0 = clock;
        [im_out, Par]  =  WLSWSC_Sigma_WA(Par);
        if size(Par.I,1) == 512
            T512 = [T512 etime(clock,time0)];
            fprintf('Total elapsed time = %f s\n', (etime(clock,time0)) );
        elseif size(Par.I,1) ==256
            T256 = [T256 etime(clock,time0)];
            fprintf('Total elapsed time = %f s\n', (etime(clock,time0)) );
        end
        im_out(im_out>1)=1;
        im_out(im_out<0)=0;
        % calculate the PSNR
        Par.PSNR(Par.outerIter, Par.image)  =   csnr( im_out*255, Par.I*255, 0, 0 );
        Par.SSIM(Par.outerIter, Par.image)      =  cal_ssim( im_out*255, Par.I*255, 0, 0 );
        imname = sprintf('C:/Users/csjunxu/Desktop/NIPS2017/W3Results/TWSC/TWSC_Sigma_WAG_nSig%d_%s', nSig, im_dir(i).name);
        imwrite(im_out,imname);
        fprintf('%s : PSNR = %2.4f, SSIM = %2.4f \n',im_dir(i).name, Par.PSNR(Par.outerIter, Par.image),Par.SSIM(Par.outerIter, Par.image)     );
    end
    mPSNR=mean(Par.PSNR,2);
    [~, idx] = max(mPSNR);
    PSNR =Par.PSNR(idx,:);
    SSIM = Par.SSIM(idx,:);
    mSSIM=mean(SSIM,2);
    mT512 = mean(T512);
    sT512 = std(T512);
    mT256 = mean(T256);
    sT256 = std(T256);
    fprintf('The best PSNR result is at %d iteration. \n',idx);
    fprintf('The average PSNR = %2.4f, SSIM = %2.4f. \n', mPSNR(idx),mSSIM);
    name = sprintf(['C:/Users/csjunxu/Desktop/NIPS2017/W3Results/TWSC_Sigma_WAG_' Sdir{end-1} '_nSig' num2str(nSig) '_step' num2str(Par.step) '_nlspini' num2str(Par.nlspini) '_nlspgap' num2str(Par.nlspgap) '_delta' num2str(Par.delta) '_l1' num2str(Par.lambda1) '_l2' num2str(Par.lambda2) '.mat']);
    save(name,'nSig','PSNR','SSIM','mPSNR','mSSIM','mT512','sT512','mT256','sT256');
end