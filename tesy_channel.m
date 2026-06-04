% main_test_ofdm_unified.m
clear; clc; close all;

%% ========================================================================
%% ========== THÀNH PHẦN 1: THAM SỐ HỆ THỐNG & KÊNH THỦY ÂM ==============
%% ========================================================================
% --- Tham số hình học kênh truyền dưới nước (SISO) ---
load pdp_Halong_150m.mat;         % Tải mảng rho và tau từ dữ liệu mẫu
load x_opt_time_SISO_D_500m.mat;  % Tải mảng tọa độ phản xạ x
Ni_1 = 40;                        % Số tia mặt nước
Ni_2 = 40;                        % Số tia đáy biển
D = 150;                          % Khoảng cách Tx - Rx (m)
x_suf_1 = x(1:Ni_1);
x_bot_1 = x(Ni_1+1 : Ni_2+Ni_1);
y_1 = 0.5;                        % Khoảng cách từ sensor tới mặt nước
y_2 = 1.5;                        % Khoảng cách từ sensor tới đáy
fD  = 15;                          % Tần số Doppler cực đại (Hz)

% --- Tham số cấu hình OFDM ---
fs          = 96000;
NFFT        = 256;
GI          = 128;
f1          = 12000;
M_ary       = 2;        % BPSK
D_f         = 4;        % 1 Pilot + 3 Data trên mỗi nhóm 4 subcarrier
Num_Sym     = 100;
OP          = 1.5;
sub_pwr     = 2;
chirp_dur   = 0.04;
pad_dur     = 0.1;
N_D         = 4;        % <--- Số sóng mang hoạt động nhỏ (Sẽ kích hoạt cơ chế 1-Pilot)
B           = fs;       % Băng thông quy đổi đồng bộ tốc độ lấy mẫu kênh
SNR_dB      = 12;       % <--- THAY ĐỔI MỨC NHIỄU GIẢ LẬP TẠI ĐÂY (dB)

input_text = 'Ra di mang nang loi the, chua thang giac My chua ve Bach Khoa.Ra di mang nang loi the, chua thang giac My chua ve Bach KhoaRa di mang nang loi the, chua thang giac My chua ve Bach Khoa Ra di mang nang loi the, chua thang giac My chua ve Bach Khoa.';

%% ========== TÍNH CHỈ SỐ SUBCARRIER & KHỞI TẠO PILOT ==========
SubL = floor(2 * f1 * NFFT / fs);
SubL = SubL - mod(SubL, 2);
SubH = SubL + N_D - 1;
N_data_per_sym = floor(N_D / D_f) * (D_f - 1);
pilot_pos = 1 : D_f : N_D;
data_pos  = setdiff(1:N_D, pilot_pos);

% Bổ sung định nghĩa Pilot_base gốc chuẩn pha PN
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

Pilot_ND            = zeros(1, N_D);
Pilot_ND(pilot_pos) = Pilot_base(1:length(pilot_pos)) * sub_pwr;
X_pilot             = Pilot_ND(pilot_pos);

%% ========================================================================
%% ========== THÀNH PHẦN 2: KHỐI PHÁT (OFDM TRANSMITTER) =================
%% ========================================================================
fprintf('=== KHỐI PHÁT (TX) ===\n');
y_bytes   = double(uint8(input_text))';
y_bits    = de2bi(y_bytes, 8, 'left-msb')';
Data_bits = y_bits(:);
NoBit_frame = Num_Sym * N_data_per_sym * log2(M_ary);
rem_bits    = mod(length(Data_bits), NoBit_frame);
if rem_bits ~= 0
    Data_bits = [Data_bits; zeros(NoBit_frame - rem_bits, 1)];
end
super_frame = length(Data_bits) / NoBit_frame;

% Tạo Preamble Chirp định dạng chuẩn vector hàng
chirp_len      = round(chirp_dur * fs);
f2_chirp       = ceil(SubH * fs / (2 * NFFT));
t_chirp        = (0 : chirp_len-1) / fs;
preamble_chirp = chirp(t_chirp, f1, t_chirp(end), f2_chirp, 'linear');
preamble_chirp = preamble_chirp(:).' .* raised_cosine_window(length(preamble_chirp));
preamble_chirp = preamble_chirp * 0.07;
pad_noise   = zeros(1, round(pad_dur*fs));

tx_parts    = cell(1, 3 + super_frame*2);
tx_parts{1} = pad_noise;
tx_parts{2} = preamble_chirp;
tx_parts{3} = pad_noise;
part_idx = 4;

for frame = 1:super_frame
    bits_frame = Data_bits((frame-1)*NoBit_frame+1 : frame*NoBit_frame);
    sym_idx  = bi2de(reshape(bits_frame, log2(M_ary), []).', 'left-msb');
    qam_syms = qammod(sym_idx, M_ary);
    data_mat = reshape(qam_syms, N_data_per_sym, Num_Sym);
    
    freq_ND             = repmat(Pilot_ND.', 1, Num_Sym);
    freq_ND(data_pos,:) = data_mat;
    freq_frame              = zeros(NFFT, Num_Sym);
    freq_frame(SubL:SubH,:) = freq_ND;
    
    freq_sym = [zeros(1, Num_Sym);
                freq_frame;
                flipud(conj(freq_frame(1:NFFT-1,:)))];
    frame_td = real(ifft(freq_sym, 2*NFFT, 1)).';
    
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

%% ========================================================================
%% ========================================================================
%% ========================================================================
%% ========== THÀNH PHẦN 3: TÁCH LUỒNG TRUYỀN & CHẬP KÊNH (TỐI ƯU NHANH) ===
%% ========================================================================
fprintf('\n=== ĐANG XỬ LÝ QUA KÊNH TRUYỀN THỰC TẾ (VECTORIZED FAST DOPPLER) ===\n');

% --- Luồng 1: Không qua kênh truyền ---
rx_ideal_signal = tx_signal;

% --- Luồng 2: Có qua kênh truyền thủy âm đa đường & Nhiễu Doppler ---
teta = 2*pi*rand(1, Ni_1+Ni_2+1);

% 1. Tính toán đáp ứng xung gốc (Base Channel) tại t = 0 để lấy các hệ số biên độ và trễ
[h_base] = UnderWaterChannel_SISO_LOS(0, teta, D, fD, rho, tau, Ni_1, Ni_2, y_1, y_2, y_1, y_2, x_suf_1, x_bot_1, B);
L_h = length(h_base);

% 2. Tính toán tần số Doppler thực tế cho TỪNG ĐƯỜNG TRUYỀN (Path)
% (Đoạn này copy lại logic tính toán góc từ hàm của bạn để biết tần số Doppler fn của mỗi tap)
alphav = 2*pi*180/360;
phi_LOS = pi/2+acos((y_1-y_1)/sqrt(D^2+(y_1-y_1)^2));

phi_1_suf = pi/2 + atan((D - x_suf_1)./y_1); 
phi_1_bot = 3*pi/2 - atan((D - x_bot_1)./y_2);

f0_LOS = fD * cos(phi_LOS - alphav);
fn_suf = fD * cos(phi_1_suf - alphav);
fn_bot = fD * cos(phi_1_bot - alphav);

% Mảng tần số Doppler tổng hợp của tất cả các tia (Đúng thứ tự như trong hàm kênh)
fD_all_paths = [f0_LOS, fn_suf, fn_bot]; 

% 3. Khởi tạo mảng thời gian và tín hiệu đầu ra
Nsamp = length(tx_signal);
t_samples = (0 : Nsamp-1) / fs; 
fading_signal = zeros(1, Nsamp);

% 4. VÒNG LẶP THEO TIA (PATH) - Thay vì vòng lặp theo mẫu (Sample)
% Số lượng tia (1 + Ni_1 + Ni_2) nhỏ hơn rất nhiều so với Nsamp nên chạy cực kỳ nhanh!

% Đầu tiên, ta cần biết mỗi tia trong h_base nằm ở tap trễ (index) nào
% Để đơn giản và chính xác, ta tái tạo lại đóng góp của từng tia:
tau_LOS = sqrt(D^2)/1500; % ước lượng sơ bộ trễ LOS
D_T_suf = sqrt(y_1^2 + x_suf_1.^2); D_R_suf = sqrt(y_1^2 + (D-x_suf_1).^2);
tau_suf = (D_T_suf + D_R_suf)/1500;
D_T_bot = sqrt(y_2^2 + x_bot_1.^2); D_R_bot = sqrt(y_2^2 + (D-x_bot_1).^2);
tau_bot = (D_T_bot + D_R_bot)/1500;

tau_all = [tau_LOS, tau_suf, tau_bot];
tau_all_after = tau_all - min(tau_all);
tap_indices = round(tau_all_after / (1/B)) + 1; % Vị trí tap trong mảng h

% Lọc từng tia và áp đặt Doppler động
for p = 1:length(fD_all_paths)
    delay_tap = tap_indices(p);
    
    % Trích xuất biên độ phức của tia thứ p tại t = 0 từ hàm của bạn
    % (Vì trong hàm của bạn các tia trùng tap sẽ bị cộng dồn, ta lấy tỷ lệ tạm tính hoặc gọi trực tiếp)
    % Cách chuẩn nhất: Tính thành phần Doppler động trực tiếp lên tín hiệu dịch trễ:
    
    % Dịch trễ tín hiệu phát tương ứng với tia này
    tx_delayed = zeros(1, Nsamp);
    if delay_tap <= Nsamp
        tx_delayed(delay_tap:end) = tx_signal(1:end-delay_tap+1);
    end
    
    % Áp pha Doppler động theo thời gian cho riêng tia này
    % Nhân chập tín hiệu đã dịch trễ với biên độ tia tại t=0 và thành phần Doppler động
    % (Tìm biên độ chuẩn của tia p bằng cách tính toán nhanh từ công thức gốc của bạn)
    c_R = 4;
    if p == 1 % LOS
        h_p0 = sqrt(c_R/(1+c_R)) * exp(1j*teta(1));
    elseif p <= 1 + Ni_1 % Surface
        n = p - 1;
        c = interpolate_c(rho, tau, tau_suf(n)-tau_LOS);
        h_p0 = sqrt(c)*exp(1j*teta(n+1)) / sqrt((Ni_1+Ni_2)*(1+c_R));
    else % Bottom
        n = p - 1 - Ni_1;
        c = interpolate_c(rho, tau, tau_bot(n)-tau_LOS);
        h_p0 = sqrt(c)*exp(1j*teta(n+1+Ni_1)) / sqrt((Ni_1+Ni_2)*(1+c_R));
    end
    
    % Cộng dồn vào tín hiệu nhận tổng tổng thể
    fading_signal = fading_signal + h_p0 * tx_delayed .* exp(1j * 2 * pi * fD_all_paths(p) * t_samples);
end

% Thực hiện lấy phần thực (giả lập truyền dải cơ sở/dải thông thực tế)
fading_signal = real(fading_signal);

% Cộng thêm nhiễu Gauss trắng theo SNR cấu hình
rx_channel_signal = awgn(fading_signal, SNR_dB, 'measured');

fprintf('-> Đã mô phỏng xong kênh Doppler bằng thuật toán Vectorized!\n');
%% ========================================================================
%% ========== THÀNH PHẦN 4: GIẢI ĐIỀU CHẾ SONG SONG (DEMODULATION) =======
%% ========================================================================
% Chạy giải điều chế cho cả 2 trường hợp bằng cách gọi hàm phụ chung ở dưới
all_syms_ideal   = demodulate_ofdm_backend(rx_ideal_signal(:), super_frame, fs, chirp_len, pad_dur, f1, f2_chirp, Num_Sym, GI, NFFT, SubL, SubH, pilot_pos, data_pos, X_pilot, true);
all_syms_channel = demodulate_ofdm_backend(rx_channel_signal(:), super_frame, fs, chirp_len, pad_dur, f1, f2_chirp, Num_Sym, GI, NFFT, SubL, SubH, pilot_pos, data_pos, X_pilot, false);

% Giải mã bits luồng qua kênh truyền để in chuỗi văn bản kiểm tra hiệu năng
all_bits_channel = double(real(all_syms_channel) > 0);
n_bytes       = floor(length(all_bits_channel) / 8);
bit_mat       = reshape(all_bits_channel(1:n_bytes*8), 8, n_bytes).';
rx_bytes      = bi2de(bit_mat, 'left-msb');
rx_text       = char(rx_bytes).';
rx_text_clean = rx_text(rx_text >= 32 & rx_text <= 126);
fprintf('\n========== KẾT QUẢ GIẢI MÃ VĂN BẢN (CÓ QUA KÊNH) ==========\n');
disp(rx_text_clean);
fprintf('==========================================================\n');

% ĐÃ SỬA TIÊU ĐỀ: In chuỗi văn bản từ luồng LÝ TƯỞNG để đối chứng
all_bits_ideal = double(real(all_syms_ideal) > 0);
n_bytes       = floor(length(all_bits_ideal) / 8);
bit_mat       = reshape(all_bits_ideal(1:n_bytes*8), 8, n_bytes).';
rx_bytes      = bi2de(bit_mat, 'left-msb');
rx_text       = char(rx_bytes).';
rx_text_clean = rx_text(rx_text >= 32 & rx_text <= 126);
fprintf('\n========== KẾT QUẢ GIẢI MÃ VĂN BẢN (LÝ TƯỞNG) ============\n');
disp(rx_text_clean);
fprintf('==========================================================\n');

%% ========================================================================
%% ========== THÀNH PHẦN 5: VẼ ĐỒ THỊ SO SÁNH CHÒM SAO ===================
%% ========================================================================
figure('Name', 'So Sanh Chom Sao OFDM: Ly Tuong vs Qua Kenh Thuy Am', 'Position', [150, 150, 1000, 480]);
% --- Bên trái: Không qua kênh ---
subplot(1,2,1);
plot(real(all_syms_ideal), imag(all_syms_ideal), 'b.', 'MarkerSize', 8);
grid on; axis equal; xlim([-2 2]); ylim([-2 2]);
xline(0, 'k--'); yline(0, 'k--');
title('Chòm sao KHÔNG qua kênh (Lý tưởng)');
xlabel('In-phase (I)'); ylabel('Quadrature (Q)');

% --- Bên phải: Có qua kênh thủy âm ---
subplot(1,2,2);
plot(real(all_syms_channel), imag(all_syms_channel), 'r.', 'MarkerSize', 6);
grid on; axis equal; xlim([-2 2]); ylim([-2 2]);
xline(0, 'k--'); yline(0, 'k--');
title(sprintf('Chòm sao CÓ qua kênh thủy âm Halong (SNR = %d dB), (fD = %d Hz)', SNR_dB, fD));
xlabel('In-phase (I)'); ylabel('Quadrature (Q)');

%% ========================================================================
%% ========== CÁC HÀM PHỤ TRỢ (SUB-FUNCTIONS) ============================
%% ========================================================================
function all_syms = demodulate_ofdm_backend(y, super_frame, fs, chirp_len, pad_dur, f1, f2_chirp, Num_Sym, GI, NFFT, SubL, SubH, pilot_pos, data_pos, X_pilot, is_ideal)
    pad_len = round(pad_dur * fs);
    K = GI + 2*NFFT;
    frame_len = Num_Sym * K;
    frame_stride = frame_len + pad_len;
    
    if is_ideal
        % Luồng lý tưởng: Đồng bộ thời gian trực tiếp không lệch mẫu
        start_data = pad_len + chirp_len + pad_len + 1;
    else
        % Luồng qua kênh: Định vị khung dựa trên bộ lọc Matched Filter tìm đỉnh Chirp
        t_chirp = (0 : chirp_len-1) / fs;
        local_chirp = chirp(t_chirp, f1, t_chirp(end), f2_chirp, 'linear');
        
        % Thiết kế cửa sổ giảm nhiễu búp phụ
        alpha=0.1; n_w=(0:chirp_len-1)'; N_w=chirp_len-1; win=ones(chirp_len,1);
        r_w = n_w < alpha*N_w/2; f_w = n_w > N_w-alpha*N_w/2;
        win(r_w) = 0.5*(1+cos(pi*(n_w(r_w)-alpha*N_w/2)/(alpha*N_w/2)));
        win(f_w) = 0.5*(1+cos(pi*(n_w(f_w)-N_w+alpha*N_w/2)/(alpha*N_w/2)));
        local_chirp = local_chirp(:).' .* win(:).';
        
        corr_out = xcorr(double(y), double(local_chirp));
corr_out = corr_out(length(y):end);

corr_abs = abs(corr_out);

[max_corr, peak_idx] = max(corr_abs);

fprintf('\n===== DEBUG CHIRP =====\n');
fprintf('Peak index      = %d\n',peak_idx);
fprintf('Peak value      = %.3e\n',max_corr);
fprintf('Signal length   = %d\n',length(y));
fprintf('Chirp length    = %d\n',chirp_len);
fprintf('Pad length      = %d\n',pad_len);

figure;
plot(corr_abs);
hold on;
plot(peak_idx,max_corr,'ro');
title('Matched Filter Output');
grid on;

start_data = peak_idx + chirp_len + pad_len;

fprintf('start_data      = %d\n',start_data);
fprintf('=======================\n');
    end
    
    all_syms = [];
    N_D = length(pilot_pos) + length(data_pos);
    
    % Đảm bảo X_pilot là vector hàng để đồng bộ phép chia
    X_pilot = X_pilot(:).';

    fprintf('\n===== DEBUG FRAME =====\n');
fprintf('frame_len    = %d\n',frame_len);
fprintf('frame_stride = %d\n',frame_stride);
fprintf('signal_len   = %d\n',length(y));
fprintf('=======================\n');
    
    for fr = 1:super_frame
        fr_start = start_data + (fr-1) * frame_stride;
        fr_end   = fr_start + frame_len - 1;
        if fr_end > length(y), break; end

        fprintf('\nFrame %d\n',fr);
fprintf('fr_start = %d\n',fr_start);
fprintf('fr_end   = %d\n',fr_end);

if fr_end > length(y)

    fprintf('!!! FRAME OUT OF RANGE !!!\n');
    fprintf('fr_end=%d > signal=%d\n', ...
            fr_end,length(y));

    break;
end
        
        % Khử bỏ khoảng bảo vệ CP
        data_cp = reshape(y(fr_start:fr_end), K, Num_Sym).';
        data_td = data_cp(:, GI+1:end);
        
        % Chuyển đổi sang miền tần số bằng FFT
        data_FFT = fft(data_td, 2*NFFT, 2);
        data_sub = data_FFT(:, SubL+1 : SubH+1); 
        
        % --- KHỐI THU THÍCH ỨNG: NỘI SUY THEO SỐ LƯỢNG PILOT ---
        X_est = zeros(Num_Sym, N_D);
        for s_idx = 1:Num_Sym
            Y_sym = data_sub(s_idx, :); 
            H_pilot = Y_sym(pilot_pos) ./ X_pilot; 
            
            if length(pilot_pos) < 2
                % Nếu chỉ có 1 Pilot (Khi N_D = 4): Coi như kênh đáp ứng phẳng
                H_interpolated = ones(1, N_D) * H_pilot(1);
            else
                % Nếu từ 2 Pilot trở lên (Khi N_D >= 8): Nội suy tuyến tính và ép vector hàng
                H_interpolated = interp1(pilot_pos, H_pilot, 1:N_D, 'linear', 'extrap');
                H_interpolated = H_interpolated(:).'; 
            end
            X_est(s_idx, :) = Y_sym ./ H_interpolated;
        end
        
        data_eq = X_est(:, data_pos);
        rx_syms = data_eq.';
        rx_syms = rx_syms(:);
        
        % Khử sai pha tĩnh tổng thể (M-th power phase estimation cho BPSK)
        phi     = angle(mean(rx_syms.^2)) / 2;
        rx_syms = rx_syms * exp(-1j * phi);
        all_syms = [all_syms; rx_syms];
    end
end

function w = raised_cosine_window(L)
    alpha=0.1; n=(0:L-1)'; N=L-1; w=ones(L,1);
    r=n<alpha*N/2; f=n>N-alpha*N/2;
    w(r)=0.5*(1+cos(pi*(n(r)-alpha*N/2)/(alpha*N/2)));
    w(f)=0.5*(1+cos(pi*(n(f)-N+alpha*N/2)/(alpha*N/2)));
    w=w(:).';
end
