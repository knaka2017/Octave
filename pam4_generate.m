function [x,t,symbols] = pam4_generate(varargin)

% ============================================================
% PAM4 Generator for GNU Octave
%
% Data modes:
%   random
%   prbs7
%   prbs15
%   prbs31
%   markov
%
% Pulse shapes:
%   nrz
%   rz
%   rc
%   rrc
%   gauss
%
% Trigger frame:
%
%   0
%   0
%   +3
%   -3
%   0
%   0
%
% ============================================================

% ---------- Defaults ----------

fs = 2.5e9;
Rs = 625e6;
N  = 2^20;

pulse = "nrz";
data_mode = "prbs15";

levels = [-3 -1 1 3];

rolloff = 0.25;
span = 8;

max_run = 16;

trigger_enable = true;
trigger_scale  = 0.8;
trigger_period = 256;

% ---------- Parse ----------

for k = 1:2:length(varargin)

    switch lower(varargin{k})

        case "fs"
            fs = varargin{k+1};

        case "rs"
            Rs = varargin{k+1};

        case "n"
            N = varargin{k+1};

        case "pulse"
            pulse = lower(varargin{k+1});

        case "data_mode"
            data_mode = lower(varargin{k+1});

        case "rolloff"
            rolloff = varargin{k+1};

        case "span"
            span = varargin{k+1};

        case "max_run"
            max_run = varargin{k+1};

        case "trigger_enable"
            trigger_enable = varargin{k+1};

        case "trigger_scale"
            trigger_scale = varargin{k+1};

        case "trigger_period"
            trigger_period = varargin{k+1};

    endswitch

endfor

% ---------- Oversampling ----------

Or = round(fs/Rs);

if bitand(Or,Or-1)
    warning("Oversampling ratio is not power of two");
endif

% ---------- Symbol count ----------

Nsym = ceil(N/Or);

data_levels = trigger_scale * levels;

% ============================================================
% Generate data symbols
% ============================================================

switch data_mode

    case "random"

        idx = randi(4,Nsym,1);
        data_symbols = data_levels(idx);

    case "prbs7"

        data_symbols = ...
            generate_prbs_pam4(7,Nsym,data_levels);

    case "prbs15"

        data_symbols = ...
            generate_prbs_pam4(15,Nsym,data_levels);

    case "prbs31"

        data_symbols = ...
            generate_prbs_pam4(31,Nsym,data_levels);

    case "markov"

        data_symbols = ...
            generate_markov_pam4( ...
                Nsym,...
                data_levels,...
                max_run);

    otherwise

        error("Unknown data_mode");

endswitch

data_symbols = data_symbols(:);

% ============================================================
% Insert trigger frame
% ============================================================

if trigger_enable

    trigger_seq = [
         levels(3) ;
         levels(2) ;
         levels(3) ;
         levels(2) ;
         levels(3) ;
         levels(2) ;
         levels(4) ;   % +3
         levels(1) ;   % -3
         levels(3) ;
         levels(2) ;
         levels(3) ;
         levels(2) ;
         levels(3) ;
         levels(2)
    ];

    trigger_len = length(trigger_seq);

    payload_len = trigger_period - trigger_len;

    if payload_len <= 0
        error("trigger_period too small");
    endif

    symbols = [];

    idx = 1;

    while idx <= length(data_symbols)

        symbols = [
            symbols;
            trigger_seq
        ];

        ncopy = min( ...
            payload_len,...
            length(data_symbols)-idx+1);

        symbols = [
            symbols;
            data_symbols(idx:idx+ncopy-1)
        ];

        idx = idx + ncopy;

    endwhile

else

    symbols = data_symbols;

endif

symbols = symbols(:);

% ============================================================
% Pulse shaping
% ============================================================

switch pulse

    case "nrz"

        x = repelem(symbols,Or);

    case "rz"

        p = [
            ones(floor(Or/2),1);
            zeros(Or-floor(Or/2),1)
        ];

        x = upsample(symbols,Or);
        x = conv(x,p,"same");

    case {"rc","raisedcosine"}

        h = rc_filter(rolloff,span,Or);

        x = upsample(symbols,Or);
        x = conv(x,h,"same");

    case {"rrc","rootraisedcosine"}

        h = rrc_filter(rolloff,span,Or);

        x = upsample(symbols,Or);
        x = conv(x,h,"same");

    case {"gauss","gaussian"}

        h = gaussian_filter(0.5,span,Or);

        x = upsample(symbols,Or);
        x = conv(x,h,"same");

    otherwise

        error("Unknown pulse");

endswitch

% ---------- Trim ----------

if length(x) < N

    x(end+1:N) = 0;

else

    x = x(1:N);

endif

x = x(:);

t = (0:N-1)'/fs;

endfunction


% ============================================================
% PRBS -> PAM4 (Gray mapping)
% ============================================================

function symbols = generate_prbs_pam4(order,Nsym,levels)

bits = prbs_bits(order,2*Nsym+100);

symbols = zeros(Nsym,1);

for k=1:Nsym

    b1 = bits(2*k-1);
    b2 = bits(2*k);

    if b1==0 && b2==0

        symbols(k)=levels(1);

    elseif b1==0 && b2==1

        symbols(k)=levels(2);

    elseif b1==1 && b2==1

        symbols(k)=levels(3);

    else

        symbols(k)=levels(4);

    endif

endfor

endfunction


% ============================================================
% PRBS generator
% ============================================================

function bits = prbs_bits(order,N)

reg = ones(1,order);

bits = zeros(N,1);

for k=1:N

    bits(k)=reg(end);

    switch order

        case 7
            fb = xor(reg(7),reg(6));

        case 15
            fb = xor(reg(15),reg(14));

        case 31
            fb = xor(reg(31),reg(28));

        otherwise
            error("Unsupported PRBS order");

    endswitch

    reg(2:end)=reg(1:end-1);
    reg(1)=fb;

endfor

endfunction


% ============================================================
% Markov PAM4
% ============================================================

function symbols = generate_markov_pam4( ...
        Nsym,levels,max_run)

symbols = zeros(Nsym,1);

symbols(1)=levels(randi(4));

runlen = 1;

for k=2:Nsym

    candidates = levels;

    if runlen >= max_run
        candidates(candidates==symbols(k-1)) = [];
    endif

    symbols(k)= ...
        candidates(randi(length(candidates)));

    if symbols(k)==symbols(k-1)
        runlen++;
    else
        runlen=1;
    endif

endfor

endfunction


% ============================================================
% Raised Cosine
% ============================================================

function h = rc_filter(beta,span,sps)

t = (-span*sps/2:span*sps/2)/sps;

h = zeros(size(t));

for k=1:length(t)

    tk = t(k);

    if abs(tk) < 1e-12

        h(k)=1;

    elseif abs(abs(2*beta*tk)-1) < 1e-12

        h(k)=pi/4*sinc_oct(1/(2*beta));

    else

        h(k)= ...
            sinc_oct(tk).* ...
            cos(pi*beta*tk)./ ...
            (1-(2*beta*tk).^2);

    endif

endfor

h = h/sum(h);

endfunction


% ============================================================
% Root Raised Cosine
% ============================================================

function h = rrc_filter(beta,span,sps)

t = (-span*sps/2:span*sps/2)/sps;

h = zeros(size(t));

for k=1:length(t)

    tk = t(k);

    if abs(tk) < 1e-12

        h(k)=1-beta+4*beta/pi;

    elseif abs(abs(tk)-1/(4*beta)) < 1e-12

        h(k)=beta/sqrt(2)* ...
              ((1+2/pi)*sin(pi/(4*beta)) + ...
               (1-2/pi)*cos(pi/(4*beta)));

    else

        h(k)= ...
            (sin(pi*tk*(1-beta)) + ...
             4*beta*tk*cos(pi*tk*(1+beta))) ...
            / ...
            (pi*tk*(1-(4*beta*tk)^2));

    endif

endfor

h = h/sqrt(sum(h.^2));

endfunction


% ============================================================
% Gaussian filter
% ============================================================

function h = gaussian_filter(bt,span,sps)

n = -span*sps/2 : span*sps/2;

alpha = sqrt(log(2))/bt;

h = exp(-(pi*n/(alpha*sps)).^2);

h = h/sum(h);

endfunction


% ============================================================
% sinc
% ============================================================

function y = sinc_oct(x)

y = ones(size(x));

idx = abs(x) > 1e-12;

y(idx) = sin(pi*x(idx))./(pi*x(idx));

endfunction
