clear;
Original_image_dir  =    'C:\Users\csjunxu\Desktop\Projects\WODL\20images\';
% Original_image_dir  =    'C:\Users\csjunxu\Desktop\Projects\WODL\20newimages\';
Sdir = regexp(Original_image_dir, '\', 'split');
fpath = fullfile(Original_image_dir, '*.png');
im_dir  = dir(fpath);
im_num = length(im_dir);

nSig0 = 40;
par.ps = 8; % patch size
par.step = 3; % the step of two neighbor patches
par.win = 30;

par.outerIter = 12;
par.innerIter = 2;

for step = [3 2]
    par.step = step;
    for nlspini = [90 100]
        par.nlspini = nlspini;
        for nlspgap = [10 5 0]
            par.nlspgap = nlspgap;
            for delta = 0.06:-0.02:0
                par.delta = delta;
                for lambdals = [0.85 0.84 0.83 0.82 0.81 0.8]
                    par.lambdals = lambdals;
                    % record all the results in each iteration
                    par.PSNR = zeros(par.outerIter, im_num, 'single');
                    par.SSIM = zeros(par.outerIter, im_num, 'single');
                    T512 = [];
                    T256 = [];
                    for i = 1:im_num
                        par.nlsp = par.nlspini;  % number of non-local patches
                        par.image = i;
                        par.nSig0 = nSig0/255;
                        par.I =  im2double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
                        S = regexp(im_dir(i).name, '\.', 'split');
                        randn('seed',0);
                        par.nim =   par.I + par.nSig0*randn(size(par.I));
                        %
                        fprintf('%s :\n',im_dir(i).name);
                        PSNR =   csnr( par.nim*255, par.I*255, 0, 0 );
                        SSIM      =  cal_ssim( par.nim*255, par.I*255, 0, 0 );
                        fprintf('The initial value of PSNR = %2.4f, SSIM = %2.4f \n', PSNR,SSIM);
                        %
                        time0 = clock;
                        [im_out, par]  =  WLSWSC_SigmaN_WA(par);
                        if size(par.I,1) == 512
                            T512 = [T512 etime(clock,time0)];
                            fprintf('Total elapsed time = %f s\n', (etime(clock,time0)) );
                        elseif size(par.I,1) ==256
                            T256 = [T256 etime(clock,time0)];
                            fprintf('Total elapsed time = %f s\n', (etime(clock,time0)) );
                        end
                        im_out(im_out>1)=1;
                        im_out(im_out<0)=0;
                        % calculate the PSNR
                        par.PSNR(par.outerIter, par.image)  =   csnr( im_out*255, par.I*255, 0, 0 );
                        par.SSIM(par.outerIter, par.image)      =  cal_ssim( im_out*255, par.I*255, 0, 0 );
                        %             imname = sprintf('nSig%d_clsnum%d_delta%2.2f_lambda%2.2f_%s', nSig, cls_num, delta, lambda, im_dir(i).name);
                        %             imwrite(im_out,imname);
                        fprintf('%s : PSNR = %2.4f, SSIM = %2.4f \n',im_dir(i).name, par.PSNR(par.outerIter, par.image),par.SSIM(par.outerIter, par.image)     );
                    end
                    mPSNR=mean(par.PSNR,2);
                    [~, idx] = max(mPSNR);
                    PSNR =par.PSNR(idx,:);
                    SSIM = par.SSIM(idx,:);
                    mSSIM=mean(SSIM,2);
                    mT512 = mean(T512);
                    sT512 = std(T512);
                    mT256 = mean(T256);
                    sT256 = std(T256);
                    fprintf('The best PSNR result is at %d iteration. \n',idx);
                    fprintf('The average PSNR = %2.4f, SSIM = %2.4f. \n', mPSNR(idx),mSSIM);
                    name = sprintf(['WLSWSC_SigmaN_WAG_' Sdir{end-1} '_nSig' num2str(nSig0) '_step' num2str(par.step) '_nlspini' num2str(par.nlspini) '_nlspgap' num2str(par.nlspgap) '_delta' num2str(delta) '_lls' num2str(lambdals) '.mat']);
                    save(name,'nSig0','PSNR','SSIM','mPSNR','mSSIM','mT512','sT512','mT256','sT256');
                end
            end
        end
    end
end