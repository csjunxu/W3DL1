clear;
Original_image_dir  =    'C:\Users\csjunxu\Desktop\Projects\WODL\20images\';
% Original_image_dir  =    'C:\Users\csjunxu\Desktop\Projects\WODL\20newimages\';
Sdir = regexp(Original_image_dir, '\', 'split');
fpath = fullfile(Original_image_dir, '*.png');
im_dir  = dir(fpath);
im_num = length(im_dir);


nSig                      =   20;         % The standard variance of the additive Gaussian noise;
for sp                        =  [0.3]        % salt and pepper
    Type = 0;
    
    if sp < 0.3
        Par.ps = 7; % patch size
        Par.step = 3;
        Par.Outerloop = 10;
        Par.nlspini = 50;
    else
        Par.ps = 8;
        Par.step = 3;
        Par.Outerloop = 12;
        Par.nlspini = 60;
    end
    
    Par.win = 30;
    Par.model = 1;
    Par.display = 1;
    Par.Innerloop = 2;
    Par.maxIter = 200;
    Par.maxrho = 100;
    
    for delta = [0]
        Par.delta = delta;
        for nlspgap = [0 5]
            Par.nlspgap = nlspgap;
            for rho = [0.5]
                Par.rho = rho;
                for mu = [1.1]
                    Par.mu = mu;
                    for lambda1 = [0.001 0.0005]
                        Par.lambda1 = lambda1;
                        for lambda2 = [1.1 1.2]
                            Par.lambda2 = lambda2;
                            % record all the results in each iteration
                            Par.PSNR = zeros(Par.Outerloop, im_num, 'double');
                            Par.SSIM = zeros(Par.Outerloop, im_num, 'double');
                            T512 = [];
                            T256 = [];
                            for i = 1:im_num
                                Par.nlsp = Par.nlspini;  % number of non-local patches
                                Par.image = i;
                                Par.nSig = nSig;
                                Par.I =  double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
                                S = regexp(im_dir(i).name, '\.', 'split');
                                %                                 randn('seed',0);
                                %                                 Par.nim =   Par.I + Par.nSig*randn(size(Par.I));
                                %                                 randn('seed',0);
                                %                                 [Par.nim,Narr]          =   impulsenoise(Par.nim,sp,Type);
                                if Type == 0
                                    Par.nim = double( imread(['images/G' num2str(nSig) '_SPIN' num2str(sp) '_' im_dir(i).name]));
                                elseif Type == 1
                                    Par.nim = double( imread(['images/G' num2str(nSig) '_RVIN' num2str(sp) '_' im_dir(i).name]));
                                else
                                    break;
                                end
                                [Par.pim,ind]           =   adpmedft(Par.nim,19);
                                
                                fprintf('%s :\n',im_dir(i).name);
                                PSNR =   csnr( Par.nim, Par.I, 0, 0 );
                                SSIM      =  cal_ssim( Par.nim, Par.I, 0, 0 );
                                fprintf('The initial value of PSNR = %2.4f, SSIM = %2.4f \n', PSNR,SSIM);
                                
                                PSNR =   csnr( Par.pim, Par.I, 0, 0 );
                                SSIM      =  cal_ssim( Par.pim, Par.I, 0, 0 );
                                fprintf('The initial value of PSNR = %2.4f, SSIM = %2.4f \n', PSNR,SSIM);
                                Par.nSig = NoiseEstimation(Par.pim, Par.ps)
                                time0 = clock;
                                [im_out, Par]  =  WLSWSC_Sigma_W1W2GINcw(Par);
                                if size(Par.I,1) == 512
                                    T512 = [T512 etime(clock,time0)];
                                    fprintf('Total elapsed time = %f s\n', (etime(clock,time0)) );
                                elseif size(Par.I,1) ==256
                                    T256 = [T256 etime(clock,time0)];
                                    fprintf('Total elapsed time = %f s\n', (etime(clock,time0)) );
                                end
                                im_out(im_out>255)=255;
                                im_out(im_out<0)=0;
                                % calculate the PSNR
                                Par.PSNR(Par.Outerloop, Par.image)  =   csnr( im_out, Par.I, 0, 0 );
                                Par.SSIM(Par.Outerloop, Par.image)      =  cal_ssim( im_out, Par.I, 0, 0 );
                                %                                 imname = sprintf(['GSPIN/W1W2GIN_' Sdir{end-1} '_nSig' num2str(nSig) '_sp' num2str(sp) '_ps' num2str(Par.ps) '_step' num2str(Par.step) '_nlspini' num2str(Par.nlspini) '_nlspgap' num2str(Par.nlspgap) '_delta' num2str(delta) '_l1' num2str(lambda1) '_l2' num2str(lambda2) '_rho' num2str(rho) '_mu' num2str(mu) im_dir(i).name]);
                                %                                 imwrite(im_out/255,imname);
                                fprintf('%s : PSNR = %2.4f, SSIM = %2.4f \n',im_dir(i).name, Par.PSNR(Par.Outerloop, Par.image),Par.SSIM(Par.Outerloop, Par.image)     );
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
                            name = sprintf(['GSPIN/W1W2GIN_p_' Sdir{end-1} '_nSig' num2str(nSig) '_sp' num2str(sp) '_ps' num2str(Par.ps) '_step' num2str(Par.step) '_nlspini' num2str(Par.nlspini) '_nlspgap' num2str(Par.nlspgap) '_delta' num2str(delta) '_l1' num2str(lambda1) '_l2' num2str(lambda2) '_rho' num2str(rho) '_mu' num2str(mu) '.mat']);
                            save(name,'nSig','sp','PSNR','SSIM','mPSNR','mSSIM','mT512','sT512','mT256','sT256');
                        end
                    end
                end
            end
        end
    end
end