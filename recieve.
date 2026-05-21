% ofdm_transmit.m
% OFDM Transmitter - bỏ GUI, chạy trực tiếp từ script
% Tác giả gốc: Nguyen Quoc Khuong - HUST
% Tối ưu: script độc lập, có thể chỉnh tham số dễ dàng

clear; clc;

%% ========== THAM SỐ HỆ THỐNG ==========
fs          = 96000;    % Sample rate (Hz)
NFFT        = 256;      % Kích thước FFT
GI          = 64;       % Guard Interval (Cyclic Prefix)
f1          = 3000;     % Tần số subcarrier thấp nhất (Hz)
f2          = 12000;    % Tần số subcarrier cao nhất (Hz)
M_ary       = 4;        % Bậc QAM (4=QPSK, 16, 64...)
D_f         = 5;        % Khoảng cách pilot (1 pilot / D_f subcarrier)
Num_Sym     = 40;       % Số OFDM symbol mỗi frame
super_frame = 4;        % Số frame
OP          = 1.6;      % Hệ số khuếch đại
sub_pwr     = 6;        % Công suất pilot
FEC_enable  = false;    % Bật/tắt FEC (convolutional code rate 1/2)
sin_len     = 1;        % Độ dài noise đệm giữa các frame (x fs)

input_text  = 'Hello OFDM World! Testing transmission via sound card.';

%% ========== TÍNH TOÁN SUBCARRIER ==========
SubL = floor(2 * f1 * NFFT / fs);
SubL = SubL - mod(SubL, 2);
SubH = ceil(2 * f2 * NFFT / fs);
SubH = SubH + mod(SubH, 2);
N_D  = SubH - SubL + 1;

% Căn chỉnh N_D chia hết cho D_f
if mod(N_D, D_f) ~= 0
    SubH = SubH + D_f - mod(N_D, D_f);
    N_D  = SubH - SubL + 1;
end

fprintf('SubL=%d, SubH=%d, N_D=%d, Subcarriers dữ liệu/symbol=%d\n', ...
        SubL, SubH, N_D, floor(N_D/D_f)*(D_f-1));

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
        error('NFFT=%d: cần load pilot từ file .mat', NFFT);
end

% Tạo pilot mask và scale
pilot_mask = [zeros(1, SubL-1), ones(1, N_D), zeros(1, NFFT - SubH)];
Pilot = sqrt(2*(M_ary-1)/3) * Pilot_base .* pilot_mask;

%% ========== CHUẨN BỊ DỮ LIỆU ==========
% Đọc text → bytes → bits
y_bytes = double(uint8(input_text))';
y_bits  = de2bi(y_bytes, 8, 'left-msb')';
Data_bits = y_bits(:);

% FEC (tùy chọn)
FEC_code = poly2trellis(3, [7 5]);

% Số bit/symbol
NoS     = Num_Sym * floor(N_D/D_f) * (D_f-1);
NoBit   = NoS * log2(M_ary);
if FEC_enable
    NoBit_data = floor(NoBit / 2);  % rate 1/2
else
    NoBit_data = NoBit;
end

% Tính số frame thực tế từ dữ liệu
super_frame = floor(length(Data_bits) / NoBit_data);
if super_frame < 1
    % Nếu text ngắn hơn 1 frame: padding
    Data_bits = [Data_bits; zeros(NoBit_data - length(Data_bits), 1)];
    super_frame = 1;
end
fprintf('Số frame: %d, Bit/frame: %d\n', super_frame, NoBit_data);

%% ========== OFDM MODULATION ==========
% Noise padding giữa frame
pad_noise = randn(1, sin_len * fs) / 10000;

tx_signal = [];

for frame = 1:super_frame
    %--- Lấy bits của frame này ---
    bits_frame = Data_bits((frame-1)*NoBit_data + 1 : frame*NoBit_data);
    
    %--- FEC encode ---
    if FEC_enable
        bits_frame = convenc(bits_frame', FEC_code)';
    end
    
    %--- Bit → symbol index ---
    bits_mat = reshape(bits_frame, log2(M_ary), [])';
    sym_idx  = bi2de(bits_mat, 'left-msb');
    
    %--- QAM modulate ---
    qam_syms = qammod(sym_idx, M_ary);          % NoS x 1
    
    %--- Reshape: (N_data_per_sym x Num_Sym) ---
    N_data_per_sym = floor(N_D/D_f) * (D_f-1);
    data_mat = reshape(qam_syms, N_data_per_sym, Num_Sym);
    
    %--- Gán pilot frame (pilot xen kẽ data) ---
    Pilot_temp = Pilot;
    Pilot_temp(SubL) = sub_pwr;                  % pilot đồng bộ frame
    Pilot_frame = repmat(Pilot_temp', 1, Num_Sym);
    
    for i = 1:floor(N_D/D_f)
        idx_data = (D_f-1)*(i-1)+1 : (D_f-1)*i;
        idx_freq = SubL + D_f*(i-1)+1 : SubL + D_f*(i-1) + D_f-1;
        Pilot_frame(idx_freq, :) = data_mat(idx_data, :);
    end
    
    %--- IFFT (real-valued output trick) ---
    dataP = Pilot_frame';                        % Num_Sym x NFFT
    frame_td = zeros(Num_Sym, 2*NFFT);
    for s = 1:Num_Sym
        frame_td(s,:) = ifft([0, dataP(s,:), fliplr(conj(dataP(s,:)))]);
    end
    
    %--- Thêm Cyclic Prefix ---
    cp_part   = frame_td(:, end-GI+1:end);
    frame_cp  = [cp_part, frame_td];             % Num_Sym x (2*NFFT + GI)
    
    %--- Ghép thành vector 1D ---
    frame_vec = real(frame_cp');
    frame_vec = frame_vec(:)';
    
    tx_signal = [tx_signal, frame_vec, pad_noise];
end

%% ========== NORMALIZE & SCALE ==========
scale = OP * 2 * NFFT / N_D / sqrt(2*(M_ary-1)/3);
tx_signal = scale * tx_signal;

% Clipping để tránh clipping trên DAC
tx_signal = max(min(tx_signal, 1), -1);

%% ========== PHÁT QUA SOUND CARD ==========
fprintf('Thời gian phát: %.2f giây\n', length(tx_signal)/fs);
fprintf('Đang phát...\n');

% Lưu file wav (để kiểm tra hoặc phát lại)
audiowrite('ofdm_tx.wav', tx_signal, fs);

% Phát trực tiếp
sound(tx_signal, fs);

%% ========== VẼ ĐỒ THỊ KIỂM TRA ==========
figure('Name', 'OFDM TX Signal');

subplot(2,1,1);
t = (0:length(tx_signal)-1) / fs;
plot(t(1:min(end, 4*(2*NFFT+GI)*Num_Sym)), ...
     tx_signal(1:min(end, 4*(2*NFFT+GI)*Num_Sym)));
xlabel('Thời gian (s)'); ylabel('Biên độ');
title('Tín hiệu OFDM miền thời gian (4 frame đầu)');
grid on;

subplot(2,1,2);
N_fft_plot = 2^nextpow2(length(tx_signal));
Xf = abs(fft(tx_signal, N_fft_plot));
f_axis = (0:N_fft_plot/2-1) * fs / N_fft_plot;
plot(f_axis/1000, 20*log10(Xf(1:N_fft_plot/2) + 1e-10));
xlabel('Tần số (kHz)'); ylabel('Biên độ (dB)');
title(sprintf('Phổ tín hiệu OFDM | f1=%.0fHz, f2=%.0fHz', f1, f2));
xlim([0 fs/2/1000]); grid on;
