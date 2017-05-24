clear;
Original_image_dir  =    'C:\Users\csjunxu\Desktop\Projects\WODL\20images\';
% Original_image_dir  =    'C:\Users\csjunxu\Desktop\Projects\WODL\20newimages\';
Sdir = regexp(Original_image_dir, '\', 'split');
fpath = fullfile(Original_image_dir, '*.png');
im_dir  = dir(fpath);
im_num = length(im_dir);


nSig = 40;
Par.ps = 8;
Par.win = 30;
Par.step = 3;

Par.outerIter = 12;
Par.innerIter = 2;

for nlspini = [70]
    Par.nlspini = nlspini;
    for lambda1 = [0]
        Par.lambda1 = lambda1;
        for nlspgap = [10]
            Par.nlspgap = nlspgap;
            for delta = [0.05]
                Par.delta = delta;
                for lambda2 = [0.8 0.7 0.6 0.5]
                    Par.lambda2 = lambda2;
                    % record all the results in each iteration
                    Par.PSNR = zeros(Par.outerIter, im_num, 'double');
                    Par.SSIM = zeros(Par.outerIter, im_num, 'double');
                    T512 = [];
                    T256 = [];
                    for i = 1:im_num
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
                        [im_out, Par]  =  W3DL1_Sigma_WA(Par);
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
                        %             imname = sprintf('nSig%d_clsnum%d_delta%2.2f_lambda%2.2f_%s', nSig, cls_num, delta, lambda, im_dir(i).name);
                        %             imwrite(im_out,imname);
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
                    name = sprintf(['W3DL1_Sigma_WAG_' Sdir{end-1} '_nSig' num2str(nSig) '_ps' num2str(Par.ps) '_step' num2str(Par.step) '_nlspini' num2str(Par.nlspini) '_nlspgap' num2str(Par.nlspgap) '_delta' num2str(delta) '_l1' num2str(lambda1) '_l2' num2str(lambda2) '.mat']);
                    save(name,'nSig','PSNR','SSIM','mPSNR','mSSIM','mT512','sT512','mT256','sT256');
                end
            end
        end
    end
end