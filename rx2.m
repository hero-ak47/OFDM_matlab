% ofdm_receive_optimized.m
% OFDM Receiver - Tối ưu hiệu năng (vectorized, pre-allocated)
% Tương thích hoàn toàn với ofdm_transmit_optimized.m
clear; clc; close all;

%% ========== THAM SỐ HỆ THỐNG ==========
fs       = 48000;
NFFT     = 256;
GI       = 128;
f1       = 9000;
M_ary    = 2;
D_f      = 4;
Num_Sym  = 180;
sub_pwr  = 2;
chirp_dur = 0.04;
pad_dur   = 0.25;

%% ========== CHỈ SỐ SUBCARRIER ==========
SubL = floor(2 * f1 * NFFT / fs);
SubL = SubL - mod(SubL, 2);
N_D  = 4;
SubH = SubL + N_D - 1;
N_data_per_sym = floor(N_D / D_f) * (D_f - 1);

fprintf('SubL=%d | SubH=%d | Subcarrier hoạt động=%d | Data/Symbol=%d\n', ...
        SubL, SubH, N_D, N_data_per_sym);

%% ========== THU QUA MICROPHONE ==========
duration_sec = 10;
fprintf('\n--- BẮT ĐẦU THU ---\nPhát tín hiệu OFDM ngay bây giờ (%d giây)...\n', duration_sec);
recObj = audiorecorder(fs, 16, 1);
recordblocking(recObj, duration_sec);
fprintf('Thu âm kết thúc. Đang xử lý...\n');

y = getaudiodata(recObj);
y = y(:);   % cột

figure('Name', 'Tín hiệu thu từ Mic');
plot(y); title('Tín hiệu thu (Time Domain)');
xlabel('Mẫu'); ylabel('Biên độ');

%% ========== PILOT SEQUENCE ==========
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

pilot_mask          = zeros(1, NFFT);
pilot_mask(SubL:SubH) = 1;
Pilot               = sqrt(2*(M_ary-1)/3 + 1e-6) * Pilot_base .* pilot_mask;
Pilot(SubL)         = sub_pwr;
Pilot_sub           = Pilot(SubL:SubH);   % 1 x N_D

%% ========== ĐỒNG BỘ THỜI GIAN (CHIRP MATCHED FILTER) ==========
fprintf('Tìm kiếm Preamble Chirp...\n');
f2_chirp  = ceil(SubH * fs / (2 * NFFT));
t_chirp   = (0 : round(chirp_dur * fs) - 1) / fs;
local_chirp = chirp(t_chirp, f1, t_chirp(end), f2_chirp, 'linear');
win_chirp   = raised_cosine_window(length(local_chirp));
local_chirp = local_chirp(:) .* win_chirp(:);

% Đảm bảo cả hai là cột double trước khi xcorr
y           = double(y(:));
local_chirp = double(local_chirp(:));

% Cross-correlation → tìm đỉnh
corr_result = xcorr(y, local_chirp);
corr_result = corr_result(length(y):end);   % chỉ lấy phần causal
[~, peak_idx] = max(abs(corr_result));

pad_samples = round(pad_dur * fs);
start_point = peak_idx + length(local_chirp) + pad_samples;
fprintf('Preamble tại mẫu %d → OFDM data bắt đầu tại mẫu %d\n', peak_idx, start_point);

%% ========== CẮT & KHỬ CP (VECTORIZED) ==========
K        = GI + 2*NFFT;
need_len = Num_Sym * K;

if length(y) < start_point + need_len - 1
    error('Tín hiệu thu quá ngắn. Kiểm tra thời gian thu hoặc đồng bộ.');
end

y_sync   = y(start_point : start_point + need_len - 1);
data_cp  = reshape(y_sync, K, Num_Sym).';          % Num_Sym x K
data_noCP = data_cp(:, GI+1:end);                  % Num_Sym x 2*NFFT

%% ========== FFT & TRÍCH XUẤT SUBCARRIER ==========
fprintf('FFT và cân bằng kênh...\n');
data_FFT = fft(data_noCP, 2*NFFT, 2);              % Num_Sym x 2*NFFT
% Lấy N_D subcarrier (1-indexed: SubL → SubH, offset +1 vì DC ở bin 1)
data_sub = data_FFT(:, SubL+1 : SubH+1);           % Num_Sym x N_D

%% ========== ƯỚC LƯỢNG & CÂN BẰNG KÊNH PHẲNG (VECTORIZED) ==========
% Vị trí pilot và data trong nhóm N_D (1-indexed)
pilot_rel = 1 : D_f : N_D;                         % [1] với N_D=4, D_f=4
data_rel  = setdiff(1:N_D, pilot_rel);              % [2,3,4]

% Ước lượng kênh: Y_pilot / X_pilot → H, mỗi symbol 1 scalar (kênh phẳng)
Y_pilot = data_sub(:, pilot_rel);                   % Num_Sym x n_pilot
X_pilot = Pilot_sub(pilot_rel);                     % 1 x n_pilot
H_est   = Y_pilot ./ X_pilot;                       % Num_Sym x n_pilot (= Num_Sym x 1)

% Zero-Forcing: áp dụng H cho toàn bộ subcarrier
X_est   = data_sub ./ repmat(H_est, 1, N_D);        % Num_Sym x N_D

% Tách data subcarrier
data_eq = X_est(:, data_rel);                       % Num_Sym x N_data_per_sym

%% ========== SẮP XẾP KÝ HIỆU ==========
rx_syms = data_eq.';
rx_syms = rx_syms(:);                               % cột (Num_Sym*N_data_per_sym) x 1

%% ========== SỬA XOAY PHA TOÀN CỤC (BPSK: dùng bình phương) ==========
if M_ary == 2
    % BPSK: pha tập trung tại 0 và π → bình phương để gom về 1 điểm
    phi     = angle(mean(rx_syms.^2)) / 2;
    rx_syms = rx_syms * exp(-1j * phi);
else
    phi     = angle(mean(rx_syms.^M_ary)) / M_ary;
    rx_syms = rx_syms * exp(-1j * phi);
end

%% ========== CONSTELLATION ==========
figure('Name', 'Constellation RX - BPSK');
plot(real(rx_syms), imag(rx_syms), 'r.', 'MarkerSize', 8);
grid on; hold on;
xl = xlim; yl = ylim;
plot(xl, [0 0], 'k--'); plot([0 0], yl, 'k--');
title('Chòm sao sau cân bằng (BPSK - 4 Subcarriers)');
xlabel('In-Phase (I)'); ylabel('Quadrature (Q)'); axis equal;

%% ========== GIẢI ĐIỀU CHẾ ==========
fprintf('Giải điều chế...\n');
if M_ary == 2
    rx_bits = double(real(rx_syms) > 0);
else
    rx_idx  = qamdemod(rx_syms, M_ary, 'gray');
    rx_bits = de2bi(rx_idx, log2(M_ary), 'left-msb').';
    rx_bits = rx_bits(:);
end

%% ========== CHUYỂN BITS → TEXT ==========
n_bytes          = floor(length(rx_bits) / 8);
rx_bits_trunc    = rx_bits(1 : n_bytes*8);
bit_mat          = reshape(rx_bits_trunc, 8, n_bytes).';
rx_bytes         = bi2de(bit_mat, 'left-msb');
rx_text          = char(rx_bytes).';
rx_text_clean    = rx_text(rx_text >= 32 & rx_text <= 126);

fprintf('\n========== KẾT QUẢ RX ==========\n');
if isempty(rx_text_clean)
    disp('[Trống] – Không trích xuất được text. Kiểm tra SNR / đồng bộ.');
else
    disp(rx_text_clean);
end
fprintf('=================================\n');

%% ========== HÀM PHỤ ==========
function w = raised_cosine_window(L)
    alpha = 0.1;
    n     = (0:L-1)';
    N     = L - 1;
    w     = ones(L, 1);
    rise  = n < alpha*N/2;
    fall  = n > N - alpha*N/2;
    w(rise) = 0.5 * (1 + cos(pi * (n(rise) - alpha*N/2) / (alpha*N/2)));
    w(fall) = 0.5 * (1 + cos(pi * (n(fall) - N + alpha*N/2) / (alpha*N/2)));
    % luôn trả về cột
    w = w(:);
end
