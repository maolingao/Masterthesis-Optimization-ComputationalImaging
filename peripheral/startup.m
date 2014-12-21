% Philipp's startup file, obtained from Michael S.

mpg = [0,0.4717,0.4604]; % color [0,125,122]
dre = [0.4906,0,0]; % color [130,0,0]
ora = [255,153,51] ./ 255;
blu = [0,0,0.509];
gra = 0.5 * ones(3,1);

lightmpg = [1,1,1] - 0.5 * ([1,1,1] - mpg);
lightdre = [1,1,1] - 0.5 * ([1,1,1] - dre);
lightblu = [1,1,1] - 0.5 * ([1,1,1] - blu);
lightora = [1,1,1] - 0.5 * ([1,1,1] - ora);

mpg2white = bsxfun(@minus,[1,1,1],bsxfun(@times,(linspace(0,0.6,2024)').^0.5,[1,1,1]-mpg));
dre2white = bsxfun(@minus,[1,1,1],bsxfun(@times,(linspace(0,0.6,2024)').^0.5,[1,1,1]-dre));
blu2white = bsxfun(@minus,[1,1,1],bsxfun(@times,(linspace(0,0.6,2024)').^0.5,[1,1,1]-blu));
ora2white = bsxfun(@minus,[1,1,1],bsxfun(@times,(linspace(0,0.6,2024)').^0.5,[1,1,1]-ora));

cya2black = bsxfun(@times,(linspace(0,0.6,2024)').^0.5,lightmpg);
red2black = bsxfun(@times,(linspace(0,0.6,2024)').^0.5,lightdre);
blu2black = bsxfun(@times,(linspace(0,0.6,2024)').^0.5,lightblu);
ora2black = bsxfun(@times,(linspace(0,0.6,2024)').^0.5,lightora);

GaussDensity = @(y,m,v) ...
    (bsxfun(@rdivide,exp(-0.5*bsxfun(@rdivide,bsxfun(@minus,y,m').^2,v'))./sqrt(2*pi),sqrt(v')));

vec = @(x)(x(:));

load handel;
beep = @() sound(y(1:2e4),Fs);
clear y Fs; % -- this is optional, if you don't want a cluttered workspace
