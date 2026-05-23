% ofdm_transmit_optimized.m
% OFDM Transmitter - Tối ưu hiệu năng (vectorized, pre-allocated)
% Tương thích hoàn toàn với ofdm_receive_optimized.m
clear; clc; close all;

%% ========== THAM SỐ HỆ THỐNG ==========
fs          = 48000;
NFFT        = 256;
GI          = 128;
f1          = 9000;
M_ary       = 2;        % BPSK
D_f         = 4;        % 1 Pilot + 3 Data trên mỗi nhóm 4 subcarrier
Num_Sym     = 180;
OP          = 1.5;
sub_pwr     = 2;
chirp_dur   = 0.04;     % Độ dài xung chirp (s)
pad_dur     = 0.25;     % Khoảng im lặng trước/sau chirp (s)

input_text  = 'Ra di mang nang loi the, chua thang giac My chua ve Bach Khoa.Ra di mang nang loi the, chua thang giac My chua ve Bach Khoa.Ra di mang nang loi the, chua thang giac My chua ve Bach Khoa ';

%% ========== TÍNH CHỈ SỐ SUBCARRIER ==========
SubL = floor(2 * f1 * NFFT / fs);
SubL = SubL - mod(SubL, 2);
N_D  = 4;
SubH = SubL + N_D - 1;
N_data_per_sym = floor(N_D / D_f) * (D_f - 1);   % = 3

fprintf('SubL=%d | SubH=%d | Subcarrier hoạt động=%d | Data/Symbol=%d\n', ...
        SubL, SubH, N_D, N_data_per_sym);

%% ========== PILOT ==========
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
Pilot(SubL)         = sub_pwr;   % Ghi đè subcarrier đầu = Pilot cố định

%% ========== CHUẨN BỊ DỮ LIỆU ==========
y_bytes   = double(uint8(input_text))';
y_bits    = de2bi(y_bytes, 8, 'left-msb')';
Data_bits = y_bits(:);

NoBit_frame = Num_Sym * N_data_per_sym * log2(M_ary);
rem_bits    = mod(length(Data_bits), NoBit_frame);
if rem_bits ~= 0
    Data_bits = [Data_bits; zeros(NoBit_frame - rem_bits, 1)];
end
super_frame = length(Data_bits) / NoBit_frame;
fprintf('Tổng Frame: %d | Tổng bits (padded): %d\n', super_frame, length(Data_bits));

%% ========== PREAMBLE CHIRP ==========
f2_chirp      = ceil(SubH * fs / (2 * NFFT));
t_chirp       = (0 : round(chirp_dur * fs) - 1) / fs;
preamble_chirp = chirp(t_chirp, f1, t_chirp(end), f2_chirp, 'linear');
preamble_chirp = preamble_chirp .* raised_cosine_window(length(preamble_chirp)) * 0.07;

%% ========== OFDM MODULATION (VECTORIZED) ==========
pad_noise    = zeros(1, round(pad_dur * fs));
sym_len_cp   = GI + 2*NFFT;

% Pre-allocate toàn bộ tín hiệu TX
n_pad_total  = (super_frame + 2) * length(pad_noise);  % ước lượng sơ bộ
tx_parts     = cell(1, 3 + super_frame * 2);
tx_parts{1}  = pad_noise;
tx_parts{2}  = preamble_chirp;
tx_parts{3}  = pad_noise;
part_idx     = 4;

% Mask subcarrier vị trí pilot & data
pilot_idx  = SubL : D_f : SubH;            % [SubL, SubL+4, ...] → chỉ SubL vì N_D=4
data_rel   = setdiff(1:N_D, 1:D_f:N_D);   % vị trí data trong nhóm N_D (1-indexed)
data_abs   = SubL - 1 + data_rel;          % index tuyệt đối trong NFFT (1-indexed, 1-based)

for frame = 1:super_frame
    bits_frame = Data_bits((frame-1)*NoBit_frame + 1 : frame*NoBit_frame);
    
    % Điều chế BPSK
    sym_idx  = bi2de(reshape(bits_frame, log2(M_ary), []).', 'left-msb');
    qam_syms = qammod(sym_idx, M_ary);          % (Num_Sym*N_data_per_sym) x 1
    data_mat = reshape(qam_syms, N_data_per_sym, Num_Sym);  % N_data x Num_Sym

    % Xây dựng frame tần số: NFFT x Num_Sym (1-indexed)
    freq_frame              = zeros(NFFT, Num_Sym);
    freq_frame(SubL, :)     = sub_pwr;           % Pilot cố định
    freq_frame(data_abs, :) = data_mat;          % Data subcarriers

    % IFFT vectorized: đầu vào Hermitian để ra tín hiệu thực
    % Dạng: [DC, positive, mirror(conj(positive))]
    freq_pos  = freq_frame(1:NFFT, :);           % NFFT x Num_Sym
    freq_sym  = [zeros(1, Num_Sym);              % DC = 0
                 freq_pos;                        % positive freqs
                 flipud(conj(freq_pos(1:NFFT-1,:)))];  % negative freqs (mirror)
    % → (2*NFFT) x Num_Sym
    frame_td  = real(ifft(freq_sym, 2*NFFT, 1)).';  % Num_Sym x 2*NFFT

    % Thêm Cyclic Prefix (vectorized)
    cp_part   = frame_td(:, end-GI+1:end);
    frame_cp  = [cp_part, frame_td];             % Num_Sym x (GI+2*NFFT)

    frame_vec = frame_cp.';
    frame_vec = frame_vec(:).';                  % 1 x (Num_Sym*(GI+2*NFFT))

    tx_parts{part_idx}   = frame_vec;
    tx_parts{part_idx+1} = pad_noise;
    part_idx             = part_idx + 2;
end

tx_signal = [tx_parts{1:part_idx-1}];

%% ========== NORMALIZE & SCALE ==========
tx_signal = tx_signal / max(abs(tx_signal) + eps);
tx_signal = max(min(tx_signal * OP * 0.5, 1), -1);

%% ========== PHÁT ÂM THANH ==========
fprintf('Thời gian phát: %.2f s\n', length(tx_signal)/fs);
sound(tx_signal, fs);

%% ========== ĐỒ THỊ ==========
figure('Name', 'OFDM TX - 4 Subcarriers');
t = (0:length(tx_signal)-1) / fs;
subplot(2,1,1);
plot(t, tx_signal); xlabel('Thời gian (s)'); ylabel('Biên độ'); grid on;
title('Tín hiệu thời gian TX');

N_plot = 2^nextpow2(length(tx_signal));
Xf     = abs(fft(tx_signal, N_plot));
f_ax   = (0:N_plot/2-1) * fs / N_plot;
subplot(2,1,2);
plot(f_ax/1000, 20*log10(Xf(1:N_plot/2) + 1e-10));
xlabel('Tần số (kHz)'); ylabel('Biên độ (dB)'); xlim([0 fs/2000]); grid on;
title('Phổ tần số TX');

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
    w = w';
end
