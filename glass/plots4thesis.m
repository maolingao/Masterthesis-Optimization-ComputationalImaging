%% toy example single problem
% see mfd_test.m

%% motion blur with different boundary conditions
figPath = '/is/ei/mgao/Documents/thesis/article/figure/fig4thesis';

f = im2double(imread('Blurry1_5_kernel.png'));
f = f(:,:,1);

figure, imagesc(f)
axis image off; colormap gray;
filename = 'kernel';
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)

% conv
x = im2double(imread('cameraman.tif'));
shape = 'same';

figure, imagesc(x)
axis image off; colormap gray;
filename = 'gt';
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)

% conv2
y = conv2(x,f,shape);

figure, imagesc(y)
axis image off; colormap gray;
filename = strcat('conv_matlab_',shape);
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
% ECO
xs = size(x);
F = conv2MatOp(f,xs,shape);
y = F*x;

figure, imagesc(y)
axis image off; colormap gray;
filename = strcat('conv_ECO_',shape);
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)

xhat = F'*y;

figure, imagesc(xhat)
axis image off; colormap gray;
filename = strcat('corr_ECO_',shape);
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)

%% atmospheric blur + satellite
figPath = '/is/ei/mgao/Documents/thesis/article/figure/fig4thesis';

load('AtmosphericBlur30');
f = PSF;

figure, imagesc(f)
axis image off; colormap gray;
filename = 'kernel';
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)

% conv
load('satel.mat')
x  =  im2double(X);
shape = 'same';

figure, imagesc(x)
axis image off; colormap gray;
filename = 'gt';
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)

% conv2
y = conv2(x,f,shape);

figure, imagesc(y)
axis image off; colormap gray;
filename = strcat('conv_matlab_',shape);
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)
% ECO
xs = size(x);
F = conv2MatOp(f,xs,shape);
y = F*x;

figure, imagesc(y)
axis image off; colormap gray;
filename = strcat('conv_ECO_',shape);
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)

xhat = F'*y;

figure, imagesc(xhat)
axis image off; colormap gray;
filename = strcat('corr_ECO_',shape);
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)


%% gradient img(cameraman) + atmpsphericblur
figPath = '/is/ei/mgao/Documents/thesis/article/figure/fig4thesis';

load('AtmosphericBlur30');
f = PSF;

figure, imagesc(f)
axis image off; colormap gray;
filename = 'kernel';
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)

% conv
% x = im2double(imread('cameraman.tif'));
load('satel.mat')
x = im2double(X);
shape = 'same';

[xu,xv] = gradient(x);

figure, imagesc(xu)
axis image off; colormap gray;
filename = 'gt_u';
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)

figure, imagesc(xv)
axis image off; colormap gray;
filename = 'gt_v';
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)


% ECO
xf = size(f);
X = conv2gradMat(x,xf,shape);
y = X*f;

figure, imagesc(y{1})
axis image off; colormap gray;
filename = strcat('conv_ECO_u_',shape);
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)

figure, imagesc(y{2})
axis image off; colormap gray;
filename = strcat('conv_ECO_v_',shape);
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)

fhat = X'*y;

figure, imagesc(fhat)
axis image off; colormap gray;
filename = strcat('corr_ECO_',shape);
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename)

%% typical frame of different SNR
figPath = '/is/ei/mgao/Documents/thesis/article/figure/fig4thesis';

multiFilt    = betterImRead; % 100 speckle samples
for k = 1 : 10
f  =  multiFilt{round(1e2*rand(1))};
f  = f./sum(vec(f));

load('satel.mat')
x  =  im2double(X);
xs =  size(x);

shape = 'same';
F = conv2MatOp(f,xs,shape);

y = F*x;
SNR = 5*k; 
y =  addnoise(y, SNR, x);

figure, imagesc(f)
axis image off; colormap gray;
filename = sprintf('atmosphericPSF_SNR%ddB',SNR);
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename), close gcf

figure, imagesc(y)
axis image off; colormap gray;
filename = sprintf('conv_satel_SNR%ddB',SNR);
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename), close gcf
end

%% different regularization parameter for same SNR img
figPath = '/home/gao/Documents/MPI/thesis/article/figure/fig4thesis';
close all hidden
multiFilt    = betterImRead; % 100 speckle samples
f  =  multiFilt{round(1e2*rand(1))};
f  = f./sum(vec(f));

load('satel.mat')
x  =  im2double(X);
xs =  size(x);

shape = 'same';
F = conv2MatOp(f,xs,shape);

figure, imagesc(f)
axis image off; colormap gray;
filename = sprintf('atmosphericPSF');
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename), close gcf

y = F*x;
%##### settings ##### 
option.mode =  'eta'; % 'SNR'; %
etaVec = [0, 1e-4, 1e-3, 1e-2, 1e-1, 1e0];
SNRVec = 10:10:60;
%##### settings ##### 

iterN = 50;
tolN = 0;
option.version = 'FH';
option.figPath = figPath;
option.plotFlag = 1;


switch option.mode
    case 'SNR'
        trgVec = SNRVec;
    case 'eta'
        trgVec = etaVec;
    otherwise
        error('wrong target!')
end

for k = 1 :  length(trgVec) % length(etaVec)
    
%     eta = 10^(-k);
    trg = trgVec(k);
    
    clear HN
    HN = hessianMatrix(eye(size(f)));
    switch trg
        case trgVec(1)
            option.color = 'dre'; % set color and linestyle of hData
            option.LineStyle = '-'; % set color and linestyle of hData
            option.LineWidth = 2;
        case trgVec(2)
            option.color = 'ora'; % set color and linestyle of hData
            option.LineStyle = ':'; % set color and linestyle of hData
            option.LineWidth = 2;
        case trgVec(3)
            option.color = 'blu'; % set color and linestyle of hData
            option.LineStyle = '-.'; % set color and linestyle of hData
            option.LineWidth = 2;
        case trgVec(4)
            option.color = 'gra'; % set color and linestyle of hData
            option.LineStyle = '--'; % set color and linestyle of hData
            option.LineWidth = 2;
        case trgVec(5)
            option.color = 'mpg'; % set color and linestyle of hData
            option.LineStyle = '-.'; % set color and linestyle of hData
            option.LineWidth = 2;
        case trgVec(6)
            option.color = 'bla'; % set color and linestyle of hData
            option.LineStyle = '-'; % set color and linestyle of hData
            option.LineWidth = 2;
        otherwise
            error('[plot4thesis.m] : malformed eta.')
    end
        
    if strcmp(option.mode,'SNR') || k == 1  % mode SNR or first loop of eta
        y = F*x;
        SNR = SNRVec(k); 
        y =  addnoise(y, SNR, x);
        
        figure, imagesc(y)
        axis image off; colormap gray;
        filename = sprintf('conv_satel_SNR%ddB',SNR);
        filename = fullfile(figPath,filename);
        print(gcf, '-depsc2', filename), close gcf
        start = F'*y;
        
        eta = 0;
    else                % mode eta
        eta = trg;
    end
    
    
    [pncg_dI, HN, data] = deconv_pncg(F, y, x, HN, iterN, start, tolN, eta, option); % pncg
    
    data_snr{k} = data;
    %
    figure, imagesc(pncg_dI)
    axis image off; colormap gray;
    filename = sprintf('conv_satel_SNR%ddB_eta%d',SNR,eta);
    filename = fullfile(figPath,filename);
    print(gcf, '-depsc2', filename), close gcf
    
end
switch option.mode
    case 'SNR'
      save('data_snrDiff_eta0.mat','data_snr')
    case 'eta'
      save('data_snrFix_etaDiff.mat','data_snr')
end
%     saveErrCurve;
%% direct kernel estimation
figPath = '/is/ei/mgao/Documents/thesis/article/figure/fig4thesis';

option.win          =   'barthann';
load('AtmosphericBlur30.mat');  f = PSF;
figure, imagesc(f)
axis image off; colormap gray;
filename = sprintf('kernel_true');
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename), close gcf


natureI = im2double(imread('cameraman.tif'));
natureI = natureI./max(vec(natureI));

figure, imagesc(natureI)
axis image off; colormap gray;
filename = sprintf('natureI');
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename), close gcf


f1 = fspecial('gaussian',[2,2],1);
X1 = conv2MatOp(natureI,size(f1),'same');
natureIEsti   = X1 * f1;
SNR = 100;
natureIEsti = addnoise(natureIEsti,SNR,natureI);

figure, imagesc(natureIEsti)
axis image off; colormap gray;
filename = sprintf('natureIEsti');
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename), close gcf

X = conv2MatOp(natureI,size(f),'same');
frame   = X * f;
frame = addnoise(frame,SNR,natureI);

figure, imagesc(frame)
axis image off; colormap gray;
filename = sprintf('frame');
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename), close gcf

natureI4convmat     =   betterEdgeTaper(natureI,option);
natureIEsti4convmat =   betterEdgeTaper(natureIEsti,option);
frameEdgeTaperred   =   betterEdgeTaper(frame,  option);
kernel_nb_noiseFree =   fast_kernel_estimate(X, natureI4convmat, frameEdgeTaperred, option);
kernel_b_noiseFree  =   fast_kernel_estimate(X, natureIEsti4convmat, frameEdgeTaperred, option);

figure, imagesc(kernel_nb_noiseFree)
axis image off; colormap gray;
filename = sprintf('kernel_nb_noiseFree');
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename), close gcf
    
    
figure, imagesc(clip(kernel_b_noiseFree,1,0e-5))
axis image off; colormap gray;
filename = sprintf('kernel_b_noiseFree');
filename = fullfile(figPath,filename);
print(gcf, '-depsc2', filename), close gcf
    
%% illustration of our model
figPath = '/is/ei/mgao/Documents/thesis/article/figure/fig4thesis';

localsetup;
multiFilt    = betterImRead; % 100 speckle samples
numFrame     = option.numFrame;
multiFilt_ds = multiFilt(randperm(length(multiFilt),numFrame)); 
[multiFrame,multiFilt_ds,F,nature] = generateMultiFrame(numFrame, multiFilt_ds, option);

%% 
figPath = '/home/gao/Documents/MPI/thesis/article/figure/fig4thesis';
startup; 

fstp_err = figure;
fstp_rerr = figure;

%%% 4 legend
load('dataK_nb_cg_gradImg.mat')
stepLine = 1:length(dataK{1}.time);

figure(fstp_err),  
hData = loglog(stepLine,dataK{1}.errs,'Color', dre, 'linewidth',2); hold on
hData = loglog(stepLine,dataK{1}.errs,'Color', mpg, 'linewidth',2);
figure(fstp_rerr),  
hData = loglog(stepLine,dataK{1}.rerrs,'Color', dre, 'linewidth',2); hold on
hData = loglog(stepLine,dataK{1}.rerrs,'Color', mpg, 'linewidth',2);

load('dataK_nb_cg_gradImg.mat')
stepLine = 1:length(dataK{1}.time);

numset = length(dataK);
color = dre;

for k = 1 : numset
    
    figure(fstp_err), hData = loglog(stepLine, dataK{k}.errs,'Color',color,'lineWidth',2); thisFigure; hold on
    figure(fstp_rerr), hData = loglog(stepLine, dataK{k}.rerrs,'Color',color,'lineWidth',2);  thisFigure; hold on
    
end

load('dataK_nb_cg_normalImg.mat')

stepLine = 1:length(dataK{1}.time);

numset = length(dataK);
color = mpg;

for k = 1 : numset
    
    figure(fstp_err), hData = loglog(stepLine, dataK{k}.errs,'Color',color,'lineWidth',2);  thisFigure; hold on
    figure(fstp_rerr), hData = loglog(stepLine, dataK{k}.rerrs,'Color',color,'lineWidth',2);  thisFigure; hold on
    
end

%
figure(fstp_err)
hXLabel = xlabel('$\#\text{steps}$', 'Interpreter','Latex');
hYLabel = ylabel('$\|f * x - y\| / pixel$');
hLegend = legend('gradient image', 'normal image');
set(hLegend,'location','northeast')
thisFigure
figName = strcat('fstp_err_grad_normal_img','.tikz');
figname = fullfile(figPath,figName);
printTikz;
%
figure(fstp_rerr)
hXLabel = xlabel('$\#\text{steps}$', 'Interpreter','Latex');
hYLabel = ylabel('$\text{relative error}$');
hLegend = legend('gradient image', 'normal image');
set(hLegend,'location','northeast')
thisFigure
figName = strcat('fstp_rerr_grad_normal_img','.tikz');
figname = fullfile(figPath,figName);
printTikz;
