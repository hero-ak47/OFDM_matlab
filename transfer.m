% % ofdm_transmit_optimized.m
clear; clc; close all;

%% ========== THAM SỐ HỆ THỐNG (CẤU HÌNH CỐ ĐỊNH 4 SUBCARRIERS) ==========
fs          = 48000;    
NFFT        = 256;      
GI          = 128;      
f1          = 7000;     % Tần số trung tâm bắt đầu dải khống chế
M_ary       = 2;        % BPSK
D_f         = 4;        % Khớp với cấu hình 4 subcarriers (1 Pilot + 3 Data)
Num_Sym     = 180;       
OP          = 1.5;      
sub_pwr     = 2;        
FEC_enable  = false;    
sin_len     = 0.25;      

input_text  = 'Ra di mang nang loi the, chua thang giac My chua ve Bach Khoa';

%% ========== TÍNH TOÁN CỐ ĐỊNH CHÍNH XÁC 4 SUBCARRIER ==========
% Tính chỉ số bắt đầu dựa trên f1
SubL = floor(2 * f1 * NFFT / fs);
SubL = SubL - mod(SubL, 2);

% ÉP BUỘC ĐỘ RỘNG: Chỉ cho phép chạy đúng 4 subcarriers hoạt động
N_D  = 4; 
SubH = SubL + N_D - 1; 

% Số lượng subcarrier thực tế dùng cho data trên mỗi Symbol
N_data_per_sym = floor(N_D/D_f) * (D_f - 1); % = 3 Data subcarriers
fprintf('CẤU HÌNH SIÊU HẸP: SubL=%d, SubH=%d, Số subcarrier hoạt động=%d, Data/Symbol=%d\n', ...
        SubL, SubH, N_D, N_data_per_sym);

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
pilot_mask = [zeros(1, SubL-1), ones(1, N_D), zeros(1, NFFT - SubH)];
Pilot = sqrt(2*(M_ary-1)/3 + 1e-6) * Pilot_base .* pilot_mask; 

%% ========== CHUẨN BỊ DỮ LIỆU & ĐỒNG BỘ ==========
y_bytes = double(uint8(input_text))';
y_bits  = de2bi(y_bytes, 8, 'left-msb')';
Data_bits = y_bits(:);

NoS     = Num_Sym * N_data_per_sym;
NoBit_data = NoS * log2(M_ary);

rem_bits = mod(length(Data_bits), NoBit_data);
if rem_bits ~= 0
    Data_bits = [Data_bits; zeros(NoBit_data - rem_bits, 1)];
end
super_frame = length(Data_bits) / NoBit_data;
fprintf('Tổng số Frame phát: %d | Tổng số bits (đã padding): %d\n', super_frame, length(Data_bits));

%% TẠO XUNG ĐỒNG BỘ CHỚP (LFM PREAMBLE)
% Để quét tần số mượt mà, xung Chirp sẽ quét trong phạm vi thực tế của 4 subcarriers này
f_start_chirp = f1;
f_end_chirp   = ceil(SubH * fs / (2 * NFFT)); 
t_chirp = 0 : 1/fs : 0.04; 
preamble_chirp = chirp(t_chirp, f_start_chirp, t_chirp(end), f_end_chirp, 'linear');
win = raised_cosine_window(length(preamble_chirp));
preamble_chirp = preamble_chirp .* win * 0.07;

%% ========== OFDM MODULATION LOOP ==========
pad_noise = zeros(1, floor(sin_len * fs)); 

tx_signal = [pad_noise, preamble_chirp, pad_noise]; 
for frame = 1:super_frame
    bits_frame = Data_bits((frame-1)*NoBit_data + 1 : frame*NoBit_data);
    
    bits_mat = reshape(bits_frame, log2(M_ary), [])';
    sym_idx  = bi2de(bits_mat, 'left-msb');
    
    qam_syms = qammod(sym_idx, M_ary);
    data_mat = reshape(qam_syms, N_data_per_sym, Num_Sym);
    
    Pilot_temp = Pilot;
    Pilot_temp(SubL) = sub_pwr; 
    Pilot_frame = repmat(Pilot_temp', 1, Num_Sym);
    
    for i = 1:floor(N_D/D_f)
        idx_data = (D_f-1)*(i-1)+1 : (D_f-1)*i;
        idx_freq = SubL + D_f*(i-1)+1 : SubL + D_f*(i-1) + D_f-1;
        Pilot_frame(idx_freq, :) = data_mat(idx_data, :);
    end
    
    dataP = Pilot_frame'; 
    frame_td = zeros(Num_Sym, 2*NFFT);
    for s = 1:Num_Sym
        frame_td(s,:) = ifft([0, dataP(s,:), fliplr(conj(dataP(s, 1:NFFT-1)))]);
    end
    
    cp_part   = frame_td(:, end-GI+1:end);
    frame_cp  = [cp_part, frame_td]; 
    
    frame_vec = real(frame_cp');
    frame_vec = frame_vec(:)';
    
    tx_signal = [tx_signal, frame_vec, pad_noise];
end

%% ========== NORMALIZE & AMPLITUDE SCALE ==========
tx_signal = tx_signal / max(abs(tx_signal));
tx_signal = tx_signal * OP * 0.5; 
tx_signal = max(min(tx_signal, 1), -1); 

%% ========== HIỂN THỊ VÀ PHÁT AUDIO ==========
fprintf('Thời gian phát: %.2f giây\n', length(tx_signal)/fs);
sound(tx_signal, fs);

%% ========== ĐỒ THỊ KIỂM TRA ==========
figure('Name', 'Hệ Thống OFDM Phát Thực Tế - 4 Subcarriers');
subplot(2,1,1); t = (0:length(tx_signal)-1) / fs; plot(t, tx_signal);
xlabel('Thời gian (s)'); ylabel('Biên độ'); grid on;
subplot(2,1,2); N_fft_plot = 2^nextpow2(length(tx_signal)); Xf = abs(fft(tx_signal, N_fft_plot));
f_axis = (0:N_fft_plot/2-1) * fs / N_fft_plot; plot(f_axis/1000, 20*log10(Xf(1:N_fft_plot/2) + 1e-10));
xlabel('Tần số (kHz)'); ylabel('Biên độ (dB)'); xlim([0 fs/2/1000]); grid on;

function w = raised_cosine_window(L)
    % Tạo cửa sổ vát cạnh giúp tín hiệu không bật đột ngột gây sốc loa
    N = L - 1; w = zeros(L, 1); alpha = 0.1; 
    for idx = 0:N
        if idx < alpha*N/2
            w(idx+1) = 0.5 * (1 + cos(pi * (idx - alpha*N/2) / (alpha*N/2)));
        elseif idx > N - alpha*N/2
            w(idx+1) = 0.5 * (1 + cos(pi * (idx - N + alpha*N/2) / (alpha*N/2)));
        else
            w(idx+1) = 1;
        end
    end
    w = w';
end
