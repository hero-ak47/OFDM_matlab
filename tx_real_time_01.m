%% =====================================================================
%  OFDM TRANSMITTER (Streaming-friendly, Acoustic Robust)
%  Dua tren kien truc no-buffer:
%   - BPSK
%   - Chirp Preamble (chong nhiễu/multipath tot hon Schmidl-Cox)
%   - Baseband -> Passband (de de dang dua len FPGA)
%% =====================================================================
clear; clc; close all;

%% ==================== THAM SO HE THONG ==============================
Fs_bb   = 8000;          % Tan so lay mau baseband (Hz)
Fs_dac  = 48000;         % Tan so DAC / sound card (Hz)
Lup     = Fs_dac/Fs_bb;  % He so upsample = 6
N       = 256;           % So diem FFT
CP      = 96;            % Do dai cyclic prefix (baseband samples)
fc      = 13500;         % Tan so song mang (Hz)
df      = Fs_bb / N;     % Khoang cach subcarrier (31.25 Hz)

% Phan bo subcarrier (32 SC -> BW = 1000 Hz)
kUsed    = [-16:-1, 1:16];   
% Phan bo subcarrier (16 SC -> BW = 500 Hz)
%kUsed    = [-8:-1, 1:8];      
bin      = mod(kUsed, N) + 1; 
isP      = mod(kUsed, 2) == 0;        % Pilot moi 2 SC
nP       = sum(isP);
nD       = sum(~isP);
pilotVal = ones(nP, 1);               % Pilot BPSK (+1)
bitsPerSym = nD;                      % BPSK: 1 bit/symbol/SC

fprintf('=== OFDM Transmitter (Robust Acoustic) ===\n');
fprintf('Fs_bb = %d Hz, Fs_dac = %d Hz, fc = %d Hz\n', Fs_bb, Fs_dac, fc);
fprintf('N = %d, CP = %d, df = %.2f Hz\n', N, CP, df);
fprintf('Data SC = %d, Pilot SC = %d (BPSK)\n', nD, nP);

%% ==================== BAN TIN =======================================
msg = 'Ra di mang nang loi the, chua thang giac My chua ve Bach Khoa.';
fprintf('Ban tin: "%s"\n', msg);

bits     = reshape(dec2bin(double(msg), 8).' - '0', [], 1);
nSym     = ceil(numel(bits) / bitsPerSym);
txBits   = [bits; zeros(nSym*bitsPerSym - numel(bits), 1)];
fprintf('So OFDM symbol: %d\n', nSym);

%% ==================== PREAMBLE CHIRP (BASEBAND) =====================
% Tao chirp tai baseband tu -800 Hz den 800 Hz (quyet qua bang thong OFDM)
chirp_dur = 0.05; % 50ms
t_c = (0:round(chirp_dur*Fs_bb)-1)' / Fs_bb;
f0 = -800; f1 = 800;
% Baseband LFM chirp: exp(j*2*pi*(f0*t + (f1-f0)/(2*T)*t^2))
c_phase = 2*pi * (f0*t_c + (f1-f0)/(2*chirp_dur) * t_c.^2);
syncSym = exp(1j * c_phase);
% Windowing de giam side-lobe
win = hann(numel(syncSym));
syncSym = syncSym .* win;

%% ==================== TAO KHUNG OFDM ===============================
rng(7);   % Seed cho training

% --- Training symbol: tat ca SC, BPSK biet truoc ---
trainBits = randi([0 1], 2*numel(kUsed), 1);
Xtrain    = bpskMap(trainBits(1:numel(kUsed)));
trainSym  = ofdmMod(fillBins(Xtrain, bin, N), N, CP);

% --- Data symbols: pilot + du lieu BPSK ---
dBits    = reshape(txBits, bitsPerSym, nSym);
dataSyms = zeros((N + CP) * nSym, 1);
for s = 1:nSym
    X = zeros(N, 1);
    X(bin(isP))  = pilotVal;
    X(bin(~isP)) = bpskMap(dBits(:, s));
    dataSyms((s-1)*(N+CP)+1 : s*(N+CP)) = ofdmMod(X, N, CP);
end

% --- Ghep khung: guard + Chirp + gap + training + data + guard ---
guardLen = Fs_bb * 0.2;  % 0.2s
gapLen   = Fs_bb * 0.02; % 20ms giua chirp va OFDM
tx_bb = [ zeros(guardLen, 1);
          syncSym;             
          zeros(gapLen, 1);
          trainSym;            
          dataSyms;            
          zeros(guardLen, 1) ];

%% ==================== UPSAMPLE + UPCONVERT =========================
Nfilt    = 96;
h_interp = Lup * fir1(Nfilt - 1, 1/Lup);

tx_up = zeros(numel(tx_bb) * Lup, 1);
tx_up(1:Lup:end) = tx_bb;
tx_up = filter(h_interp, 1, tx_up);

n       = (0:numel(tx_up)-1).';
tx_pass = real(tx_up .* exp(1j * 2*pi * fc/Fs_dac * n));
tx_pass = 0.8 * tx_pass / max(abs(tx_pass));

%% ==================== LUU FILE =====================================
wavFile = fullfile(pwd, 'ofdm_tx_signal.wav');
audiowrite(wavFile, tx_pass, Fs_dac, 'BitsPerSample', 24);

refFile = fullfile(pwd, 'ofdm_tx_ref.mat');
save(refFile, 'Fs_bb', 'Fs_dac', 'Lup', 'N', 'CP', 'fc', 'df', ...
     'kUsed', 'bin', 'isP', 'nP', 'nD', 'pilotVal', 'bitsPerSym', ...
     'syncSym', 'gapLen', 'Xtrain', 'trainBits', 'txBits', 'bits', 'msg', ...
     'nSym', 'h_interp', 'Nfilt');
fprintf('Da luu ofdm_tx_signal.wav va ofdm_tx_ref.mat\n');

%% ==================== PHAT QUA LOA ==================================
fprintf('\n>>> Nhan Enter de phat qua loa...\n');
pause;
fprintf('Dang phat... (%.2f s)\n', numel(tx_pass)/Fs_dac);
sound(tx_pass, Fs_dac);
pause(numel(tx_pass)/Fs_dac + 0.5);
fprintf('Phat xong!\n');

%% ==================== HAM PHU ======================================
function s = ofdmMod(X, N, CP)
    x = ifft(X) * sqrt(N);
    s = [x(end-CP+1:end); x];
end
function X = fillBins(v, bin, N)
    X = zeros(N, 1);  X(bin) = v;
end
function s = bpskMap(b)
    s = 1 - 2*b(:); % 0->1, 1->-1
end
