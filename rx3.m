% ofdm_receive_optimized.m
% RX khớp hoàn toàn với TX — giải toàn bộ super_frame
clear; clc; close all;

%% ========== THAM SỐ HỆ THỐNG — PHẢI KHỚP VỚI TX ==========
fs        = 48000;
NFFT      = 256;
GI        = 128;
f1        = 10000;
M_ary     = 2;
D_f       = 4;
Num_Sym   = 180;
sub_pwr   = 2;
chirp_dur = 0.04;
pad_dur   = 0.1;

N_D = 4;   % ← PHẢI BẰNG TX: 4 | 8 | 16 | 32

%% ========== TÍNH CHỈ SỐ SUBCARRIER ==========
SubL = floor(2 * f1 * NFFT / fs);
SubL = SubL - mod(SubL, 2);
SubH = SubL + N_D - 1;
N_data_per_sym = floor(N_D / D_f) * (D_f - 1);

fprintf('SubL=%d | SubH=%d | N_D=%d | Data/sym=%d\n', SubL, SubH, N_D, N_data_per_sym);

%% ========== PILOT (GIỐNG TX) ==========
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

pilot_pos = 1 : D_f : N_D;
data_pos  = setdiff(1:N_D, pilot_pos);

Pilot_ND            = zeros(1, N_D);
Pilot_ND(pilot_pos) = Pilot_base(1 : N_D/D_f) * sub_pwr;
X_pilot             = Pilot_ND(pilot_pos);   % giá trị pilot đã biết tại RX

%% ========== TÍNH THỜI GIAN THU ==========
K            = GI + 2*NFFT;
frame_len    = Num_Sym * K;
pad_len      = round(pad_dur * fs);
chirp_len    = round(chirp_dur * fs);
frame_stride = frame_len + pad_len;

duration_sec = 8;   % ← TX in ra số chính xác cần thiết khi chạy

fprintf('Bắt đầu thu %d giây — phát TX NGAY...\n', duration_sec);
recObj = audiorecorder(fs, 16, 1);
recordblocking(recObj, duration_sec);
fprintf('Thu xong. Đang xử lý...\n');

y = double(getaudiodata(recObj));
y = y(:);

figure('Name','Tín hiệu thu từ Mic');
plot((0:length(y)-1)/fs, y);
xlabel('Thời gian (s)'); ylabel('Biên độ'); title('Tín hiệu thu');

% [b,a] = butter(5, [f1-400 f1+400]/(fs/2), 'bandpass');
% y = y / max(abs(y));  % normalize
% y = filtfilt(b, a, y);
%% ========== ĐỒNG BỘ THỜI GIAN — CHIRP MATCHED FILTER ==========
fprintf('Tìm Preamble Chirp...\n');
f2_chirp    = ceil(SubH * fs / (2*NFFT));
t_chirp     = (0 : chirp_len-1) / fs;
local_chirp = chirp(t_chirp, f1, t_chirp(end), f2_chirp, 'linear');
win_c       = raised_cosine_window(length(local_chirp));
local_chirp = double(local_chirp(:)) .* double(win_c(:));

corr_out      = xcorr(y, local_chirp);
corr_out      = corr_out(length(y):end);
[~, peak_idx] = max(abs(corr_out));

start_data = peak_idx + chirp_len + pad_len;
fprintf('Preamble tại mẫu %d → Data bắt đầu mẫu %d\n', peak_idx, start_data);

%% ========== ĐẾM SỐ FRAME KHẢ DỤNG ==========
available   = length(y) - start_data + 1;
super_frame = floor((available + pad_len) / frame_stride);
super_frame = max(super_frame, 1);
fprintf('Số frame khả dụng: %d\n', super_frame);

%% ========== GIẢI MÃ TỪNG FRAME ==========
all_bits = [];
all_syms = [];

for fr = 1:super_frame
    fr_start = start_data + (fr-1) * frame_stride;
    fr_end   = fr_start + frame_len - 1;

    if fr_end > length(y)
        fprintf('Frame %d vượt tín hiệu — dừng tại frame %d.\n', fr, fr-1);
        break;
    end

    % Cắt & khử CP
    data_cp = reshape(y(fr_start:fr_end), K, Num_Sym).';
    data_td = data_cp(:, GI+1:end);               % Num_Sym x 2*NFFT

    % FFT & trích subcarrier
    data_FFT = fft(data_td, 2*NFFT, 2);
    data_sub = data_FFT(:, SubL+1 : SubH+1);      % Num_Sym x N_D

    % Ước lượng kênh & ZF (vectorized)
    Y_pilot = data_sub(:, pilot_pos);             % Num_Sym x N_pilot
    H_pilot = Y_pilot ./ X_pilot;                 % Num_Sym x N_pilot
    H_full  = repmat(H_pilot(:,1), 1, N_D);       % kênh phẳng: copy H pilot đầu

    X_est   = data_sub ./ H_full;
    data_eq = X_est(:, data_pos);                 % Num_Sym x N_data_per_sym

    % Sắp xếp & sửa pha
    rx_syms = data_eq.';
    rx_syms = rx_syms(:);
    phi     = angle(mean(rx_syms.^2)) / 2;
    rx_syms = rx_syms * exp(-1j * phi);

    % Gom
    all_syms = [all_syms; rx_syms];
    all_bits = [all_bits; double(real(rx_syms) > 0)];

    fprintf('Frame %d/%d giải xong (%d bits)\n', fr, super_frame, length(rx_syms));
end

%% ========== CONSTELLATION — TẤT CẢ FRAME ==========
figure('Name', sprintf('Constellation RX — %d SC, %d frames', N_D, super_frame));
plot(real(all_syms), imag(all_syms), 'r.', 'MarkerSize', 4);
grid on; hold on;
plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');
title(sprintf('Chòm sao BPSK — %d SC | %d frames | %d symbols', ...
      N_D, super_frame, length(all_syms)));
xlabel('I'); ylabel('Q'); axis equal;

%% ========== BITS → TEXT ==========
n_bytes       = floor(length(all_bits) / 8);
bit_mat       = reshape(all_bits(1:n_bytes*8), 8, n_bytes).';
rx_bytes      = bi2de(bit_mat, 'left-msb');
rx_text       = char(rx_bytes).';
rx_text_clean = rx_text(rx_text >= 32 & rx_text <= 126);

fprintf('\n========== KẾT QUẢ RX (%d frames, %d bytes) ==========\n', super_frame, n_bytes);
if isempty(rx_text_clean)
    disp('[Trống] — Kiểm tra đồng bộ hoặc SNR.');
else
    disp(rx_text_clean);
end
fprintf('======================================================\n');

%% ========== HÀM PHỤ ==========
function w = raised_cosine_window(L)
    alpha=0.1; n=(0:L-1)'; N=L-1; w=ones(L,1);
    r=n<alpha*N/2; f=n>N-alpha*N/2;
    w(r)=0.5*(1+cos(pi*(n(r)-alpha*N/2)/(alpha*N/2)));
    w(f)=0.5*(1+cos(pi*(n(f)-N+alpha*N/2)/(alpha*N/2)));
    w=w(:);
end
