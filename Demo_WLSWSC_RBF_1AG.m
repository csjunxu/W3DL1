clear;
Original_image_dir  =    'C:\Users\csjunxu\Desktop\PGPD_TIP\20images\';
fpath = fullfile(Original_image_dir, '*.png');
im_dir  = dir(fpath);
im_num = length(im_dir);

nSig = 40;
par.ps = 7; % patch size
par.step = 3; % the step of two neighbor patches
par.win = max(3*par.ps, 20);

par.outerIter = 8;
par.innerIter = 2;
par.WWIter = 100;
par.epsilon = 0.01;
par.model = 2;

for delta = 0.08
    par.delta = delta;
    for lambdasc = 0.1:0.1:0.5
        par.lambdasc = lambdasc;
        for lambdals = 0:0.1:0.5
            par.lambdals = lambdals;
            % record all the results in each iteration
            par.PSNR = zeros(par.outerIter, im_num, 'single');
            par.SSIM = zeros(par.outerIter, im_num, 'single');
            T512 = [];
            T256 = [];
            for i = 1:im_num
                par.nlsp = 90;  % number of non-local patches
                par.image = i;
                par.nSig = nSig/255;
                par.I =  single( imread(fullfile(Original_image_dir, im_dir(i).name)) )/255;
                S = regexp(im_dir(i).name, '\.', 'split');
                randn('seed',0);
                par.nim =   par.I + par.nSig*randn(size(par.I));
                %
                fprintf('%s :\n',im_dir(i).name);
                PSNR =   csnr( par.nim*255, par.I*255, 0, 0 );
                SSIM      =  cal_ssim( par.nim*255, par.I*255, 0, 0 );
                fprintf('The initial value of PSNR = %2.4f, SSIM = %2.4f \n', PSNR,SSIM);
                %
                time0 = clock;
                [im_out, par]  =  WLSWSC_RBF_1AG(par);
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
            name = sprintf(['WLSWSC_RBF_1AG_nSig' num2str(nSig) '_delta' num2str(delta) '_lsc' num2str(lambdasc) '_lls' num2str(lambdals) '_WIter' num2str(par.WWIter) '.mat']);
            save(name,'nSig','PSNR','SSIM','mPSNR','mSSIM','mT512','sT512','mT256','sT256');
        end
    end
end