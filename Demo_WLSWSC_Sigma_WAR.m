clear;
% GT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2017\cc_Results\Real_MeanImage\';
% GT_fpath = fullfile(GT_Original_image_dir, '*.png');
% TT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2017\cc_Results\Real_NoisyImage\';
% TT_fpath = fullfile(TT_Original_image_dir, '*.png');
% GT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2017\cc_Results\Real_ccnoise_denoised_Part\';
% GT_fpath = fullfile(GT_Original_image_dir, '*mean.png');
% TT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2017\cc_Results\Real_ccnoise_denoised_Part\';
% TT_fpath = fullfile(TT_Original_image_dir, '*real.png');
GT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2017\1_Results\Real_NoisyImage\';
GT_fpath = fullfile(GT_Original_image_dir, '*.png');
TT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2017\1_Results\Real_NoisyImage\';
TT_fpath = fullfile(TT_Original_image_dir, '*.png');
GT_im_dir  = dir(GT_fpath);
TT_im_dir  = dir(TT_fpath);
im_num = length(TT_im_dir);

% Parameters
Par.ps = 6;        % patch size
Par.step = 3;       % the step of two neighbor patches
Par.win = 20;   % size of window around the patch

Par.Outerloop = 8;
Par.Innerloop = 2;
Par.maxIter = 10;
Par.maxrho = 100;
Par.nlspini = 70;

Par.model = 1;
Par.display = 1;
Par.delta = 0;
Par.nlspgap = 10;
for mu = [1.1]
    Par.mu = mu;
    for rho = [0.5]
        Par.rho = rho;
        for lambda1 = [0]
            Par.lambda1 = lambda1;
            for lambda2 = [2.4:0.2:3]
                Par.lambda2 = lambda2;
                % set Parameters
                % record all the results in each iteration
                Par.PSNR = zeros(Par.Outerloop, im_num, 'double');
                Par.SSIM = zeros(Par.Outerloop, im_num, 'double');
                CCPSNR = [];
                CCSSIM = [];
                for i = [4] %1 : im_num
                    Par.nlsp = Par.nlspini;  % number of non-local patches
                    Par.image = i;
                    Par.nim = double(imread(fullfile(TT_Original_image_dir, TT_im_dir(i).name) ));
                    Par.I = double(imread(fullfile(GT_Original_image_dir, GT_im_dir(i).name)));
                    S = regexp(TT_im_dir(i).name, '\.', 'split');
                    IMname = S{1};
                    [h,w,ch] = size(Par.nim);
                    % noise estimation
                    for c = 1:ch
                        %                     Par.nSig0(c) = NoiseLevel(Par.nim(:, :, c));
                        Par.nSig(c) = NoiseEstimation(Par.nim(:, :, c), Par.ps);
                    end
                    % initial PSNR and SSIM
                    fprintf('%s: \n', TT_im_dir(i).name);
                    CCPSNR = [CCPSNR csnr( Par.nim, Par.I, 0, 0 )];
                    CCSSIM = [CCSSIM cal_ssim( Par.nim, Par.I, 0, 0 )];
                    fprintf('The initial PSNR = %2.4f, SSIM = %2.4f. \n', CCPSNR(end), CCSSIM(end));
                    % denoising
                    t1=clock;
                    [IMout, Par]  =  WLSWSC_Sigma_WAR(Par);
                    t2=clock;
                    etime(t2,t1)
                    alltime(Par.image)  = etime(t2, t1);
                    % calculate the PSNR
                    Par.PSNR(Par.Outerloop, Par.image)  =   csnr( IMout, Par.I, 0, 0 );
                    Par.SSIM(Par.Outerloop, Par.image)      =  cal_ssim( IMout, Par.I, 0, 0 );
                    %% output
                    fprintf('The final PSNR = %2.4f, SSIM = %2.4f. \n', Par.PSNR(Par.Outerloop, Par.image), Par.SSIM(Par.Outerloop, Par.image));
                    %% output
                    imname = sprintf(['C:/Users/csjunxu/Desktop/NIPS2017/W3Results/TWSC/TWSC_1AR_' num2str(lambda2) '_' TT_im_dir(i).name]);
                    imwrite(IMout/255, imname);
                end
%                 mPSNR=mean(Par.PSNR,2);
%                 [~, idx] = max(mPSNR);
%                 PSNR =Par.PSNR(idx,:);
%                 SSIM = Par.SSIM(idx,:);
%                 mSSIM=mean(SSIM,2);
%                 mtime  = mean(alltime);
%                 mCCPSNR = mean(CCPSNR);
%                 mCCSSIM = mean(CCSSIM);
%                 save(['C:/Users/csjunxu/Desktop/NIPS2017/W3Results/WLSWSC_WAR_delta' num2str(Par.delta) '_ps' num2str(Par.ps) '_step' num2str(Par.step) '_nlspini' num2str(Par.nlspini) '_nlspgap' num2str(Par.nlspgap) '_rho' num2str(Par.rho) '_mu' num2str(Par.mu) '_l1' num2str(Par.lambda1) '_l2' num2str(Par.lambda2) '.mat'],'alltime','mtime','PSNR','mPSNR','SSIM','mSSIM','CCPSNR','mCCPSNR','CCSSIM','mCCSSIM');
            end
        end
    end
end