% ofdm_receive.m
% OFDM Receiver - script độc lập, không GUI
% Tác giả gốc: Nguyen Quoc Khuong - HUST

clear; clc;

%% ========== THAM SỐ HỆ THỐNG (phải khớp với TX) ==========
fs          = 96000;
NFFT        = 256;
GI          = 64;
f1          = 3000;
f2          = 12000;
M_ary       = 4;
D_f         = 5;
Num_Sym     = 40;
sub_pwr     = 6;
FEC_enable  = false;
Send_Text   = true;

% Tham số Doppler
doppler_adjust  = 0;
max_doppler     = 60;   % Hz

%% ========== ĐỌC TÍN HIỆU (từ file hoặc thu trực tiếp) ==========
USE_FILE = true;   % true = đọc file, false = thu mic

if USE_FILE
    fprintf('Đọc file WAV...\n');
    [y, fs_file] = audioread('ofdm_tx.wav');
    y = y(:);
    if fs_file ~= fs
        y = resample(y, fs, fs_file);
    end
else
    fprintf('Thu âm %.1f giây...\n', Num_Sym * (2*NFFT+GI)/fs * 5);
    rec_time = Num_Sym * (2*NFFT+GI)/fs * 5 + 0.5;
    r = audiorecorder(fs, 16, 1);
    recordblocking(r, rec_time);
    y = getaudiodata(r, 'double');
end

fprintf('Độ dài tín hiệu thu: %d samples (%.2f s)\n', length(y), length(y)/fs);

%% ========== TÍNH TOÁN SUBCARRIER ==========
SubL = floor(2 * f1 * NFFT / fs);
SubL = SubL - mod(SubL, 2);
SubH = ceil(2 * f2 * NFFT / fs);
SubH = SubH + mod(SubH, 2);
N_D  = SubH - SubL + 1;

if mod(N_D, D_f) ~= 0
    SubH = SubH + D_f - mod(N_D, D_f);
    N_D  = SubH - SubL + 1;
end

K = 2*NFFT + GI + 1;   % Độ dài 1 OFDM symbol (có CP)

%% ========== PILOT SEQUENCE ==========
switch NFFT
    case 64
        Pilot_base = [1 1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 -1 1 1 -1 ...
                     -1 -1 1 1 -1 1 0 -1 1 -1 1 -1 1 1 1 -1 1 -1 -1 ...
                      1 -1 1 1 -1 -1 -1 1 -1 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    case 256
        Pilot_base = [ 1 1 -1 1 1 1 -1 -1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 -1 ...
            -1 -1 -1 -1 1 1 -1 1 1 1 1 -1 -1 1 1 -1 1 -1 -1 1 1 ...
             1 -1 1 1 1 1 1 1 -1 -1 1 -1 1 1 -1 -1 1 -1 1 -1 -1 ...
            -1 -1 -1 1 1 -1 1 -1 -1 -1 1 1 1 -1 1 -1 -1 1 1 1 1 ...
             1 -1 1 1 -1 1 1 -1 1 -1 -1 -1 1 -1 1 -1 -1 -1 1 -1 1 ...
             1 1 1 1 1 1 -1 1 -1 -1 -1 1 1 1 -1 -1 1 -1 -1 -1 -1 ...
             1 1 -1 -1 1 1 1 1 1 1 -1 1 -1 -1 -1 -1 -1 1 -1 1 -1 ...
            -1 1 1 1 -1 -1 1 -1 -1 1 -1 1 -1 1 -1 1 -1 1 1 1 -1 ...
             1 1 1 -1 1 -1 1 -1 1 -1 1 1 1 -1 -1 -1 1 -1 -1 -1 -1 ...
            -1 1 -1 1 1 -1 1 1 -1 -1 1 -1 1 -1 -1 1 -1 -1 -1 -1 1 ...
             1 1 -1 -1 1 -1 -1 -1 -1 -1 -1 -1 1 1 1 1 1 1 1 -1 -1 ...
            -1 -1 -1 1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 1 1 -1 1 1 -1 ...
             1 1 -1 1];
    otherwise
        error('NFFT=%d chưa hỗ trợ, cần load pilot từ file .mat', NFFT);
end

pilot_mask = [zeros(1, SubL-1), ones(1, N_D), zeros(1, NFFT - SubH)];
Pilot = sqrt(2*(M_ary-1)/3) * Pilot_base .* pilot_mask;
Pilot(SubL) = sub_pwr;  % pilot đồng bộ frame

%% ========== BƯỚC 1: LỌC BĂNG THÔNG ==========
fprintf('Lọc BPF [%.0fHz - %.0fHz]...\n', f1, f2);
[B, A] = butter(3, [1.8*f1/fs, 2.2*f2/fs]);
y_filt = filter(B, A, y);

%% ========== BƯỚC 2: ƯỚC LƯỢNG DOPPLER ==========
fprintf('Ước lượng Doppler...\n');
Y_full = abs(fft(y_filt));
N_y    = length(Y_full);

% Tìm peak trong vùng [f1, f2]
idx_lo = floor(f1 * N_y / fs) + 1;
idx_hi = floor(f2 * N_y / fs);
[~, rel_idx] = max(Y_full(idx_lo:idx_hi));
ftg_r = (idx_lo + rel_idx - 1) * fs / N_y;

% Tần số pilot tham chiếu (subcarrier đồng bộ)
ftg_0 = SubL * fs / (2 * NFFT);

delta_f  = (ftg_0 - ftg_r) * fs / ftg_0;
doppler  = ftg_0 - ftg_r;
fprintf('  ftg_ref=%.1fHz, ftg_rx=%.1fHz, Doppler=%.2fHz\n', ftg_0, ftg_r, doppler);

% Bù Doppler bằng resample (nếu lệch nhỏ)
if abs(doppler) < max_doppler
    fs_tx = round(fs - delta_f);
    if fs_tx ~= fs && abs(fs_tx - fs) < 200
        % Tìm tỉ lệ đơn giản nhất
        g = gcd(fs_tx, fs);
        y_filt = resample(y_filt, fs/g, fs_tx/g);
        fprintf('  Resample: %d → %d\n', fs_tx, fs);
    end
end

%% ========== BƯỚC 3: ĐỒNG BỘ THỜI GIAN (CP Correlation) ==========
fprintf('Đồng bộ thời gian...\n');

% Tìm điểm bắt đầu signal mạnh
ss = find(y_filt > max(y_filt)*0.5);
if isempty(ss)
    error('Không tìm thấy tín hiệu. Kiểm tra file WAV.');
end

start_search = max(1, ss(1) - K);
end_search   = min(length(y_filt), ss(1) + 3*K);
y_win = y_filt(start_search:end_search);
L_win = length(y_win);

% Tính CP correlation metric
syn_corr = zeros(1, L_win - 2*NFFT - GI - 2);
syn_eng  = zeros(1, L_win - 2*NFFT - GI - 2);

for i = 1:length(syn_corr)
    seg1 = y_win(i:i+GI);
    seg2 = y_win(i+K-GI:i+K);
    syn_corr(i) = max(syn_corr) - sum(abs(seg1 - seg2));  % diff nhỏ → tốt
    syn_eng(i)  = abs(seg1' * seg2);
end
% Kết hợp 2 metric
syn_metric = (max(syn_corr) - syn_corr) .* syn_eng;
syn_metric = max(syn_metric) - syn_metric;   % đảo: đỉnh = vị trí sync

[~, best_idx] = max(syn_metric);
start_point = start_search + best_idx - 1;
fprintf('  Start point: %d (%.3f s)\n', start_point, start_point/fs);

y_sync = y_filt(start_point:end);

%% ========== BƯỚC 4: TÁCH SYMBOL + XÓA CP ==========
fprintf('Tách %d symbol...\n', Num_Sym);

num_sym_avail = floor(length(y_sync) / K);
if num_sym_avail < Num_Sym
    warning('Chỉ có %d/%d symbol, thêm zero-pad.', num_sym_avail, Num_Sym);
    y_sync = [y_sync; zeros(Num_Sym*K - length(y_sync), 1)];
end

% Reshape → tách CP
data_GI  = reshape(y_sync(1:Num_Sym*K), K, Num_Sym)';   % Num_Sym x K
data_noCP = data_GI(:, GI+1:end);                        % Num_Sym x 2*NFFT+1

%% ========== BƯỚC 5: FFT ==========
fprintf('FFT...\n');
data_FFT = fft(data_noCP, [], 2);       % FFT theo từng hàng

% Lấy subcarrier [SubL+1 : SubH+1] (1-indexed, tương ứng freq dương)
data_sub = data_FFT(:, SubL+1:SubH+1);  % Num_Sym x N_D

%% ========== BƯỚC 6: CHANNEL ESTIMATION & EQUALIZATION ==========
fprintf('Channel estimation & equalization...\n');

Pilot_sub = Pilot(SubL:SubH);           % Pilot trong vùng subcarrier

data_eq = zeros(Num_Sym, floor(N_D/D_f)*(D_f-1));

for s = 1:Num_Sym
    % Ước lượng H tại vị trí pilot (mỗi D_f subcarrier)
    Y_pilot = data_sub(s, 1:D_f:end);
    X_pilot = Pilot_sub(1:D_f:end);
    H_est_sparse = Y_pilot ./ X_pilot;
    
    % Nội suy H cho tất cả subcarrier
    H_interp = interp(H_est_sparse, D_f);
    H_interp = H_interp(1:N_D);
    
    % Equalize từng data subcarrier (bỏ pilot)
    R = [];
    for k = 2:D_f
        idx_data = k:D_f:N_D;
        H_data   = H_interp(idx_data);
        Y_data   = data_sub(s, idx_data);
        R        = [R, Y_data ./ H_data];
    end
    data_eq(s, :) = R;
end

%% ========== BƯỚC 7: QAM DEMODULATION ==========
fprintf('QAM-%d demodulation...\n', M_ary);

rx_syms = data_eq.';
rx_syms = rx_syms(:);

%% ========== BƯỚC 8: GIẢI MÃ TEXT ==========
if Send_Text
    FEC_code = poly2trellis(3, [7 5]);
    
    % QAM demod → symbol index
    rx_idx = qamdemod(rx_syms, M_ary);
    
    % Symbol → bits
    rx_bits = de2bi(rx_idx, log2(M_ary), 'left-msb')';
    rx_bits = rx_bits(:)';
    
    % FEC decode (Viterbi)
    if FEC_enable
        tb = 11;
        rx_bits = vitdec([rx_bits, zeros(1, 2*tb)], FEC_code, tb, 'cont', 'hard');
        rx_bits = rx_bits(tb+1:end);
    end
    
    % Bits → bytes → ASCII
    n_bytes = floor(length(rx_bits) / 8);
    rx_bits = rx_bits(1:n_bytes*8);
    bit_mat = reshape(rx_bits, 8, n_bytes)';
    rx_text = char(bi2de(bit_mat, 'left-msb'))';
    
    % Lọc ký tự in được
    rx_text = rx_text(rx_text >= 32 & rx_text < 127);
    
    fprintf('\n====== VĂN BẢN GIẢI MÃ ======\n');
    disp(rx_text);
    fprintf('==============================\n');
else
    fprintf('Send_Text=false, bỏ qua giải mã text.\n');
end

%% ========== VẼ ĐỒ THỊ ==========
figure('Name', 'OFDM RX Analysis', 'Position', [100 100 1200 800]);

subplot(3,2,1);
plot((1:length(y))/fs, y);
title('Tín hiệu thu (gốc)'); xlabel('Thời gian (s)'); grid on;

subplot(3,2,2);
f_ax = (0:length(y)-1) * fs / length(y);
plot(f_ax(1:end/2)/1000, 20*log10(abs(fft(y_filt))(1:end/2) + 1e-10));
title('Phổ sau BPF (dB)'); xlabel('Tần số (kHz)'); grid on;
xlim([0 fs/2/1000]);

subplot(3,2,3);
plot(syn_metric);
title('CP Correlation (đỉnh = start symbol)'); xlabel('Sample'); grid on;

subplot(3,2,4);
plot(real(rx_syms(1:min(end,200))), imag(rx_syms(1:min(end,200))), '.');
axis equal; grid on;
title(sprintf('Constellation QAM-%d sau equalization', M_ary));
xlabel('I'); ylabel('Q');

subplot(3,2,5);
plot(abs(data_FFT(1, :)));
xline(SubL, 'r--', 'SubL'); xline(SubH, 'g--', 'SubH');
title('FFT symbol đầu tiên'); xlabel('Subcarrier index'); grid on;

subplot(3,2,6);
H_plot = data_sub(round(Num_Sym/2), 1:D_f:end) ./ Pilot_sub(1:D_f:end);
plot(abs(H_plot));
title('Channel estimate |H| (symbol giữa)');
xlabel('Pilot index'); ylabel('|H|'); grid on;

fprintf('\nHoàn thành! Kiểm tra constellation để đánh giá chất lượng.\n');
