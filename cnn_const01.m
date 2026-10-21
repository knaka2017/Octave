function cnn_constellation_to_w_b_demo()
% Constellation density image -> b / w (LTI channel / equalizer)
% GNU Octave compatible for data generation
% CNN training requires MATLAB + Deep Learning Toolbox

rng(1);

%% ===== 設定 =====
targetType = "w";     % "w" / "b" / "both"
M = 16;               % 16QAM
Lsym = 4096;          % 判定点数（密度安定用）

Nb = 7;               % チャネル長
Nw = 21;              % 等化器長（Nbより長く）

imgH = 32; imgW = 32; % 画像サイズ
useLogDensity = true;

Ntrain = 2000;
Nval   = 300;

SNRdB_min = 10;
SNRdB_max = 30;

chanDecay = 0.7;

%% ===== データ生成 =====
fprintf("Generating training data...\n");
[XTrain, YTrain, outDim] = make_dataset( ...
    Ntrain, M, Lsym, Nb, Nw, imgH, imgW, ...
    SNRdB_min, SNRdB_max, chanDecay, useLogDensity, targetType);

fprintf("Generating validation data...\n");
[XVal, YVal] = make_dataset( ...
    Nval, M, Lsym, Nb, Nw, imgH, imgW, ...
    SNRdB_min, SNRdB_max, chanDecay, useLogDensity, targetType);

save cnn_dataset_octave.mat XTrain YTrain XVal YVal outDim
fprintf("Dataset saved: cnn_dataset_octave.mat\n");

%% ===== MATLAB only: CNN training =====
if exist("trainNetwork","file") ~= 2
    warning("Deep Learning Toolbox not found. Dataset generation only.");
    return;
end

layers = [
    imageInputLayer([imgH imgW 1], Normalization="none")

    convolution2dLayer(3,16,"Padding","same")
    reluLayer
    maxPooling2dLayer(2,"Stride",2)

    convolution2dLayer(3,32,"Padding","same")
    reluLayer
    maxPooling2dLayer(2,"Stride",2)

    convolution2dLayer(3,64,"Padding","same")
    reluLayer

    globalAveragePooling2dLayer
    fullyConnectedLayer(128)
    reluLayer
    fullyConnectedLayer(outDim)
    regressionLayer
];

opts = trainingOptions("adam", ...
    "MaxEpochs",10, ...
    "MiniBatchSize",64, ...
    "Shuffle","every-epoch", ...
    "ValidationData",{XVal,YVal}, ...
    "Verbose",true, ...
    "Plots","training-progress");

net = trainNetwork(XTrain,YTrain,layers,opts);
save cnn_trained_net.mat net
end

function [X4D, Y, outDim] = make_dataset(N, M, Lsym, Nb, Nw, imgH, imgW, ...
    SNRdB_min, SNRdB_max, chanDecay, useLogDensity, targetType)

X4D = zeros(imgH,imgW,1,N,"single");

if targetType=="b"
    outDim = 2*Nb;
elseif targetType=="w"
    outDim = 2*Nw;
else
    outDim = 2*Nb + 2*Nw;
end
Y = zeros(N,outDim,"single");

for i = 1:N
    x = qam_mod(M,Lsym);
    b = random_channel_fir(Nb,chanDecay);

    y0 = filter(b,1,x);
    SNRdB = SNRdB_min + (SNRdB_max-SNRdB_min)*rand();
    sigma = snr_to_sigma(y0,SNRdB);
    n = (sigma/sqrt(2))*(randn(size(y0))+1j*randn(size(y0)));
    y = y0 + n;

    img = constellation_density_image(y,imgH,imgW,useLogDensity);
    X4D(:,:,1,i) = single(img);

    if targetType=="b"
        Y(i,:) = [real(b).' imag(b).'];
    elseif targetType=="w"
        w = design_w_mmse(b,Nw,SNRdB);
        Y(i,:) = [real(w).' imag(w).'];
    else
        w = design_w_mmse(b,Nw,SNRdB);
        Y(i,:) = [real(b).' imag(b).' real(w).' imag(w).'];
    end
end
end

function img = constellation_density_image(y,H,W,useLogDensity)
I = real(y);
Q = imag(y);

p = 0.995;
lim = max([quantile(abs(I),p), quantile(abs(Q),p), 1e-6]);

edgesI = linspace(-lim,lim,W+1);
edgesQ = linspace(-lim,lim,H+1);

N = hist2d_octave(I,Q,edgesI,edgesQ);

if useLogDensity
    N = log1p(N);
end
N = N - min(N(:));
if max(N(:))>0
    N = N / max(N(:));
end
img = N;
end

function N = hist2d_octave(x,y,xedges,yedges)
Nx = length(xedges)-1;
Ny = length(yedges)-1;

[~,ix] = histc(x,[-inf xedges(2:end-1) inf]);
[~,iy] = histc(y,[-inf yedges(2:end-1) inf]);

valid = ix>=1 & ix<=Nx & iy>=1 & iy<=Ny;
ix = ix(valid);
iy = iy(valid);

N = accumarray([iy ix],1,[Ny Nx]);
end

function x = qam_mod(M,L)
m = sqrt(M);
re = randi([0 m-1],L,1);
im = randi([0 m-1],L,1);
x = (2*re-(m-1)) + 1j*(2*im-(m-1));
x = x / sqrt(mean(abs(x).^2));
end

function b = random_channel_fir(Nb,decay)
k = (0:Nb-1).';
b = decay.^k .* (randn(Nb,1)+1j*randn(Nb,1));
b = b / sqrt(sum(abs(b).^2)+eps);
end

function sigma = snr_to_sigma(y,SNRdB)
Ps = mean(abs(y).^2);
sigma = sqrt(Ps / 10^(SNRdB/10));
end

function w = design_w_mmse(b,Nw,SNRdB)
Nb = length(b);
[~,k] = max(abs(b));
d = k-1;

B = convmtx(b(:),Nw);
e = zeros(Nw+Nb-1,1);
e(d+1)=1;

lambda = 1/10^(SNRdB/10);
w = (B'*B + lambda*eye(Nw))\(B'*e);
end


