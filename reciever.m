% % ofdm_receive_optimized.m
% RX tương thích hoàn toàn với mã phát ổn định (BPSK + Chirp Sync)
% CẬP NHẬT: Cấu hình cố định chỉ thu 4 subcarriers
clear; clc; close all;
%% ================= SYSTEM PARAMETERS =================
fs = 48000; % Tần số lấy mẫu
NFFT = 256; % Kích thước FFT gốc
GI = 128; % Chiều dài CP
f1 = 7000; % Tần số subcarrier thấp nhất (Hz)
M_ary = 2; % BPSK (M=2)
D_f = 4; % Khớp với cấu hình 4 subcarriers (1 Pilot + 3 Data)
Num_Sym = 180; % Số OFDM symbol mỗi frame
sub_pwr = 2; % Năng lượng pilot tương ứng
%% ================= CONFIG & THU QUA MICROPHONE REAL-TIME =================
duration_sec = 6;
fprintf('--- CHUẨN BỊ THU THỰC TẾ QUA MIC ---\n');
fprintf('Mẹo: Bật script thu này TRƯỚC, sau đó nhấn PHÁT bên máy phát ngay lập tức.\n\n');
recObj = audiorecorder(fs, 16, 1);
fprintf('Đang mở Mic... Hãy phát tín hiệu OFDM ngay bây giờ (%d giây)...\n', duration_sec);
recordblocking(recObj, duration_sec);
fprintf('Thu âm kết thúc. Đang xử lý dữ liệu...\n');
y = getaudiodata(recObj);
y = y(:);
figure('Name', 'Tín hiệu thô thu từ Mic');
plot(y); title('Tín hiệu biên độ thu từ Mic (Time Domain)');
xlabel('Mẫu (Samples)'); ylabel('Biên độ');
%% ================= SUBCARRIER CALCULATIONS (CỐ ĐỊNH 4 SUBCARRIERS) =================
% Tính toán vị trí bắt đầu dựa trên f1
SubL = floor(2 * f1 * NFFT / fs);
SubL = SubL - mod(SubL, 2);
% ÉP CỐ ĐỊNH ĐỒNG BỘ HOÀN TOÀN VỚI PHÍA PHÁT
N_D = 4;
SubH = SubL + N_D - 1;
% Số subcarrier thực tế mang data trên 1 Symbol (= 3)
N_data_per_sym = floor(N_D/D_f) * (D_f - 1);
fprintf('CẤU HÌNH SIÊU HẸP: SubL=%d, SubH=%d, Số subcarrier hoạt động=%d, Data/Symbol=%d\n', ...
SubL, SubH, N_D, N_data_per_sym);
%% ================= CHUẨN BỊ PILOT MẪU =================
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
Pilot(SubL) = sub_pwr;
Pilot_sub = Pilot(SubL:SubH);
%% ================= ĐỒNG BỘ THỜI GIAN (CHIRP MATCHED FILTER) =================
fprintf('Đang tìm kiếm xung đồng bộ Preamble Chirp...\n');
% Tính toán f2 thực tế của 4 subcarrier để tái tạo chính xác xung phát
f2_actual = ceil(SubH * fs / (2 * NFFT));
t_chirp = 0 : 1/fs : 0.04;
local_chirp = chirp(t_chirp, f1, t_chirp(end), f2_actual, 'linear')';
% Áp cửa sổ giống phía phát
N_c = length(local_chirp);
alpha = 0.1; win = ones(N_c, 1);
for idx = 0:N_c-1
if idx < alpha*N_c/2
win(idx+1) = 0.5 * (1 + cos(pi * (idx - alpha*N_c/2) / (alpha*N_c/2)));
elseif idx > N_c - alpha*N_c/2
win(idx+1) = 0.5 * (1 + cos(pi * (idx - N_c + alpha*N_c/2) / (alpha*N_c/2)));
end
end
local_chirp = local_chirp .* win;
% Tương quan chéo
corr_result = xcorr(y, local_chirp);
corr_result = corr_result(length(y):end);
[~, peak_idx] = max(abs(corr_result));
pad_samples = floor(0.25 * fs); % Đồng bộ với sin_len = 0.5s từ TX
start_point = peak_idx + length(local_chirp) + pad_samples;
fprintf('Đã tìm thấy Preamble tại mẫu số: %d. Điểm bắt đầu OFDM data: %d\n', peak_idx, start_point);
%% ================= TÁCH CẮT VÀ KHỬ CP =================
K = 2*NFFT + GI;
need_len = Num_Sym * K;
if length(y) < (start_point + need_len - 1)
error('Tín hiệu thu bị ngắn hoặc mất mát dữ liệu frame, không thể giải mã.');
end
y_sync = y(start_point : start_point + need_len - 1);
data_cp = reshape(y_sync, K, Num_Sym).';
data_noCP = data_cp(:, GI+1:end);
%% ================= KHỐI FFT & TRÍCH XUẤT TẦN SỐ =================
fprintf('Thực hiện FFT và Cân bằng kênh...\n');
data_FFT = fft(data_noCP, 2*NFFT, 2);
data_sub = data_FFT(:, SubL+1 : SubH+1);
%% ================= ƯỚC LƯỢNG VÀ CÂN BẰNG KÊNH PHẲNG =================
data_eq = zeros(Num_Sym, N_data_per_sym);
for s = 1:Num_Sym
% 1. Trích xuất phiến Pilot nhận được từ kênh truyền (Chỉ có 1 phần tử)
Y_pilot = data_sub(s, 1:D_f:end);
X_pilot = Pilot_sub(1:D_f:end);
% 2. Ước lượng đáp ứng H kênh tại vị trí đặt Pilot đầu tiên
H_single = Y_pilot ./ X_pilot;
% 3. THAY THẾ INTERP1: Sao chép giá trị H này cho cả 4 subcarriers (Kênh phẳng)
H_interp = repmat(H_single, 1, N_D);
% 4. Cân bằng Zero-Forcing triệt nhiễu méo kênh
X_est = data_sub(s, :) ./ H_interp;
% 5. Lọc tách bỏ vị trí Pilot để thu hồi mảng Data mốc
data_only = [];
for k = 1:N_D
if mod(k-1, D_f) ~= 0
data_only = [data_only, X_est(k)];
end
end
data_eq(s, :) = data_only(1:N_data_per_sym);
end
%% ================= SẮP XẾP CHUỖI KÝ TỰ =================
rx_syms = data_eq.';
rx_syms = rx_syms(:);
%% ================= SỬA LỖI XOAY PHA TOÀN CỤC =================
if M_ary == 2
phi = angle(mean(rx_syms.^2)) / 2;
else
phi = angle(mean(rx_syms.^4)) / 4;
end
rx_syms = rx_syms * exp(-1j * phi);
%% ================= HIỂN THỊ CHÒM SAO KHÔNG GIAN =================
figure('Name', 'OFDM RX Constellation - 4 Subcarriers');
plot(real(rx_syms), imag(rx_syms), 'r.', 'MarkerSize', 8);
grid on; hold on;
xl = xlim; yl = ylim;
plot([xl(1) xl(2)], [0 0], 'k--'); plot([0 0], [yl(1) yl(2)], 'k--');
title('Chòm sao tín hiệu thu sau cân bằng (4 Subcarriers - BPSK)');
xlabel('In-Phase (I)'); ylabel('Quadrature (Q)');
axis equal;
%% ================= GIẢI ĐIỀU CHẾ VÀ CHUYỂN ĐỔI BITS =================
fprintf('Giải điều chế dữ liệu bit...\n');
if M_ary == 2
rx_bits = real(rx_syms) > 0;
else
rx_idx = qamdemod(rx_syms, M_ary, 'gray');
rx_bits = de2bi(rx_idx, log2(M_ary), 'left-msb').';
rx_bits = rx_bits(:);
end
%% ================= DỊCH BITS RA TEXT NGUYÊN BẢN =================
n_bytes = floor(length(rx_bits) / 8);
rx_bits_truncated = rx_bits(1 : n_bytes*8);
bit_mat = reshape(rx_bits_truncated, 8, n_bytes).';
rx_bytes = bi2de(bit_mat, 'left-msb');
rx_text = char(rx_bytes).';
rx_text_clean = rx_text(rx_text >= 32 & rx_text <= 126);
fprintf('\n========== KẾT QUẢ ĐẦU RA BỘ THU (RX TEXT) ==========\n');
if isempty(rx_text_clean)
disp('[Trống] - Không thể trích xuất văn bản hợp lệ. Kiểm tra SNR hoặc pha.');
else
disp(rx_text_clean);
end
fprintf('=====================================================\n');
