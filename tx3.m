% ofdm_transmit_optimized.m
clear; clc; close all;

%% ========== THAM SỐ HỆ THỐNG ==========
fs          = 48000;
NFFT        = 256;
GI          = 128;
f1          = 10000;
M_ary       = 2;        % BPSK
D_f         = 4;        % 1 Pilot + 3 Data trên mỗi nhóm 4 subcarrier
Num_Sym     = 180;
OP          = 1.5;
sub_pwr     = 2;
chirp_dur   = 0.04;
pad_dur     = 0.1;

N_D = 4;   % ← SỐ SÓNG MANG (bội số của D_f=4): 4 | 8 | 16 | 32

input_text = 'Ra di mang nang loi the, chua thang giac My chua ve Bach Khoa.Ra di mang nang loi the, chua thang giac My chua ve Bach Khoa.';

%% ========== TÍNH CHỈ SỐ SUBCARRIER ==========
SubL = floor(2 * f1 * NFFT / fs);
SubL = SubL - mod(SubL, 2);
SubH = SubL + N_D - 1;
N_data_per_sym = floor(N_D / D_f) * (D_f - 1);   % = N_D * 3/4

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

% Vị trí pilot và data trong N_D subcarrier (1-indexed)
pilot_pos = 1 : D_f : N_D;                % [1,5,9,13] với N_D=16
data_pos  = setdiff(1:N_D, pilot_pos);    % [2,3,4, 6,7,8, ...]

% Pilot vector cho đúng N_D subcarrier
Pilot_ND            = zeros(1, N_D);
Pilot_ND(pilot_pos) = Pilot_base(1 : N_D/D_f) * sub_pwr;

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

frame_sec = Num_Sym * (GI + 2*NFFT) / fs;
fprintf('Tổng Frame: %d | Tổng bits (padded): %d\n', super_frame, length(Data_bits));
fprintf('>>> RX cần thu tối thiểu: %.1f s\n', ...
        pad_dur + chirp_dur + pad_dur + super_frame*(frame_sec + pad_dur) + 1);

%% ========== PREAMBLE CHIRP ==========
f2_chirp      = ceil(SubH * fs / (2 * NFFT));
t_chirp       = (0 : round(chirp_dur*fs)-1) / fs;
preamble_chirp = chirp(t_chirp, f1, t_chirp(end), f2_chirp, 'linear');
preamble_chirp = preamble_chirp .* raised_cosine_window(length(preamble_chirp)) * 0.07;

%% ========== OFDM MODULATION ==========
pad_noise  = zeros(1, round(pad_dur*fs));
tx_parts   = cell(1, 3 + super_frame*2);
tx_parts{1} = pad_noise;
tx_parts{2} = preamble_chirp;
tx_parts{3} = pad_noise;
part_idx = 4;

for frame = 1:super_frame
    bits_frame = Data_bits((frame-1)*NoBit_frame+1 : frame*NoBit_frame);

    sym_idx  = bi2de(reshape(bits_frame, log2(M_ary), []).', 'left-msb');
    qam_syms = qammod(sym_idx, M_ary);
    data_mat = reshape(qam_syms, N_data_per_sym, Num_Sym);   % N_data x Num_Sym

    % Ghép pilot + data vào block N_D x Num_Sym
    freq_ND             = repmat(Pilot_ND.', 1, Num_Sym);
    freq_ND(data_pos,:) = data_mat;

    % Nhúng vào FFT frame NFFT x Num_Sym
    freq_frame              = zeros(NFFT, Num_Sym);
    freq_frame(SubL:SubH,:) = freq_ND;

    % IFFT Hermitian → tín hiệu thực
    freq_sym = [zeros(1, Num_Sym);
                freq_frame;
                flipud(conj(freq_frame(1:NFFT-1,:)))];
    frame_td = real(ifft(freq_sym, 2*NFFT, 1)).';   % Num_Sym x 2*NFFT

    % CP
    cp_part   = frame_td(:, end-GI+1:end);
    frame_cp  = [cp_part, frame_td];
    frame_vec = frame_cp.';
    frame_vec = frame_vec(:).';

    tx_parts{part_idx}   = frame_vec;
    tx_parts{part_idx+1} = pad_noise;
    part_idx = part_idx + 2;
end

tx_signal = [tx_parts{1:part_idx-1}];
tx_signal = tx_signal / (max(abs(tx_signal)) + eps);
tx_signal = max(min(tx_signal * OP * 0.5, 1), -1);

fprintf('Thời gian phát: %.2f s\n', length(tx_signal)/fs);
sound(tx_signal, fs);

%% ========== ĐỒ THỊ ==========
figure('Name', sprintf('OFDM TX — %d Subcarriers', N_D));
t = (0:length(tx_signal)-1)/fs;
subplot(2,1,1); plot(t, tx_signal);
xlabel('Thời gian (s)'); ylabel('Biên độ'); grid on; title('Tín hiệu thời gian TX');

N_plot = 2^nextpow2(length(tx_signal));
Xf = abs(fft(tx_signal, N_plot));
f_ax = (0:N_plot/2-1)*fs/N_plot;
subplot(2,1,2); plot(f_ax/1000, 20*log10(Xf(1:N_plot/2)+1e-10));
xlabel('Tần số (kHz)'); ylabel('Biên độ (dB)'); xlim([0 fs/2000]); grid on;
title('Phổ tần số TX');

%% ========== HÀM PHỤ ==========
function w = raised_cosine_window(L)
    alpha=0.1; n=(0:L-1)'; N=L-1; w=ones(L,1);
    r=n<alpha*N/2; f=n>N-alpha*N/2;
    w(r)=0.5*(1+cos(pi*(n(r)-alpha*N/2)/(alpha*N/2)));
    w(f)=0.5*(1+cos(pi*(n(f)-N+alpha*N/2)/(alpha*N/2)));
    w=w';
end
