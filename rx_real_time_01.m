%% =====================================================================
%  OFDM RECEIVER (Streaming, Acoustic Robust)
%  Dua tren kien truc no-buffer:
%   - BPSK thay vi QPSK
%   - Dong bo bang Preamble Chirp (Matched Filter) thay vi Schmidl-Cox
%   - BPSK Blind Phase Tracking (giai quyet CFO/Phase offset)
%   - Khong can Farrow resampler (bo qua Doppler Gian thoi gian neu nho)
%% =====================================================================
clear; clc; close all;

%% ==================== CHE DO HOAT DONG ==============================
mode = 'record';           % 'record' / 'loopback' / 'sim' / 'file'
recDuration  = 10;
rxWavFile    = 'ofdm_rx_recorded.wav';

sim_snr      = 15;
sim_delays   = [0 15 35 60];  % Kênh âm thanh dội rất dài
sim_gains    = [1 .6 .4 .2];

fprintf('=== OFDM Receiver (Robust Acoustic) [mode: %s] ===\n', mode);

%% ==================== LOAD THAM SO TU TX ============================
if exist('ofdm_tx_ref.mat', 'file')
    load('ofdm_tx_ref.mat');
    fprintf('Da load tham so tu ofdm_tx_ref.mat\n');
else
    error('Khong tim thay ofdm_tx_ref.mat! Chay ofdm_tx.m truoc.');
end

L_sync = numel(syncSym);
E_sync = sum(abs(syncSym).^2);

%% ==================== KHOI TAO NGUON TIN HIEU =======================
isRealtime = strcmp(mode, 'record');
chunkSize_dac = 1024;

if ~isRealtime
    switch mode
        case 'loopback'
            [y_pass_all, ~] = audioread('ofdm_tx_signal.wav');
        case 'sim'
            [y_pass_all, ~] = audioread('ofdm_tx_signal.wav');
            rng(42);
            h_ch = zeros(max(sim_delays)+1, 1);
            h_ch(sim_delays+1) = sim_gains .* exp(1j*2*pi*rand(1, numel(sim_delays)));
            y_pass_all = filter(h_ch, 1, y_pass_all);
            sigP = mean(y_pass_all.^2);
            y_pass_all = y_pass_all + sqrt(sigP/10^(sim_snr/10)) * randn(size(y_pass_all));
        case 'file'
            [y_pass_all, ~] = audioread(rxWavFile);
    end
    if size(y_pass_all, 2) > 1, y_pass_all = y_pass_all(:,1); end
    passPtr = 1;
end

%% ==================== KHOI TAO FRONT-END ============================
nco_phase = 0;
nco_freq  = 2*pi * fc / Fs_dac;
h_aa      = 2 * fir1(Nfilt - 1, 1/Lup);
zf_aa     = zeros(Nfilt - 1, 1);
decim_phase = 0;

%% ==================== KHOI TAO RECEIVER STATE MACHINE ===============
state     = 0;        
buf       = complex([]);  
rp        = 1;        
symIdx    = 0;        
syncFound = false;
nSymDecoded = 0;      

rxBits  = zeros(bitsPerSym, nSym);
eqAll   = [];
snrEst  = [];
rxTextSoFar = '';     

pilotPos = find(isP);
dataPos  = find(~isP);

% Bo loc cho Chirp Matched Filter
h_match = flipud(conj(syncSym));
zf_match = zeros(numel(h_match)-1, 1);
h_energy = ones(L_sync, 1);
zf_energy = zeros(numel(h_energy)-1, 1);
sync_thr = 0.3; % Nguong tuong quan (0 -> 1)
metric_buf = [];

if isRealtime
    hasAudioTbx = ~isempty(ver('audio'));
    if hasAudioTbx
        reader = audioDeviceReader('SampleRate', Fs_dac, 'SamplesPerFrame', chunkSize_dac);
        cleanupObj = onCleanup(@() release(reader));
    else
        rec = audiorecorder(Fs_dac, 16, 1);
        record(rec);   
        lastReadIdx = 0;
        cleanupObj = onCleanup(@() stop(rec));
    end
    fprintf('\n>>> DANG THU (%ds)...\n', recDuration);
    tic;
end

%% ====================================================================
finished = false;
y_pass_log = [];    

while ~finished
    % DOC CHUNK
    if isRealtime
        if hasAudioTbx
            audioChunk = reader();
        else
            pause(chunkSize_dac / Fs_dac * 0.8);
            allData = getaudiodata(rec, 'double');
            nTotal  = numel(allData);
            if nTotal <= lastReadIdx, continue; end
            audioChunk  = allData(lastReadIdx+1 : nTotal);
            lastReadIdx = nTotal;
        end
        if toc > recDuration, finished = true; end
    else
        iEnd = min(passPtr + chunkSize_dac - 1, numel(y_pass_all));
        if passPtr > numel(y_pass_all), finished = true; continue; end
        audioChunk = y_pass_all(passPtr : iEnd);
        passPtr    = iEnd + 1;
    end

    nSamp = numel(audioChunk);
    if nSamp == 0, continue; end
    y_pass_log = [y_pass_log; audioChunk]; %#ok<AGROW>

    % DOWNCONVERT & LPF
    n_vec  = (0:nSamp-1).';
    y_iq   = audioChunk .* exp(-1j * (nco_phase + nco_freq * n_vec));
    nco_phase = mod(nco_phase + nco_freq * nSamp, 2*pi);
    [y_filt, zf_aa] = filter(h_aa, 1, y_iq, zf_aa);

    % DECIMATE
    firstOut = mod(Lup - decim_phase, Lup) + 1;
    if firstOut > nSamp
        decim_phase = mod(decim_phase + nSamp, Lup);
        bb_chunk = complex([]);
    else
        outIdx   = firstOut : Lup : nSamp;
        bb_chunk = y_filt(outIdx);
        decim_phase = mod(decim_phase + nSamp, Lup);
    end

    if isempty(bb_chunk), continue; end
    buf = [buf; bb_chunk]; %#ok<AGROW>

    % ---------- STATE 0: MATCHED FILTER SYNC ----------
    if state == 0
        % Chay matched filter tren chunk moi
        [corr_out, zf_match] = filter(h_match, 1, bb_chunk, zf_match);
        [energy_out, zf_energy] = filter(h_energy, 1, abs(bb_chunk).^2, zf_energy);
        
        M = abs(corr_out).^2 ./ (energy_out * E_sync + 1e-6);
        metric_buf = [metric_buf; M]; %#ok<AGROW>

        % Tim dinh (peak) neu vuot nguong
        [maxM, mIdx] = max(M);
        if maxM > sync_thr
            % Kiem tra xem dinh da "ha xuong" chua (de chac chan la dinh cuoi)
            if mIdx < numel(M) - 5
                % Tinh toan vi tri doc
                % mIdx la vi tri ket thuc cua syncSym trong bb_chunk nay.
                % Chuyen mIdx tu he toa do bb_chunk sang he toa do buf.
                % buf_len hien tai = numel(buf). bb_chunk co do dai numel(bb_chunk).
                % Vay vi tri tuong doi trong buf la:
                idx_in_buf = numel(buf) - numel(bb_chunk) + mIdx;
                
                rp = idx_in_buf + gapLen + 1; 
                symIdx = 0;
                state  = 1;
                syncFound = true;
                fprintf('[SYNC] Tim thay Preamble (Metric = %.2f)\n', maxM);
            end
        end
        
        % Giai phong buffer
        if state == 0 && numel(buf) > 3*(N+CP)
            buf(1:numel(buf) - 2*(N+CP)) = [];
        end
    end

    % ---------- STATE 1: BPSK OFDM DEMOD ----------
    while state == 1 && nSymDecoded < nSym
        needEnd = rp + N + CP - 1;
        if needEnd > numel(buf), break; end

        % Trich xuat 1 symbol
        seg = buf(rp : rp + N + CP - 1);
        
        % Khong dung Farrow resampler nua (gian thoi gian ngan, bo qua)
        % Khong bu CFO truoc FFT (xu ly sau FFT bang BPSK phase tracker)
        
        % Xoa CP & FFT
        seg_no_cp = seg(CP+1:end);
        Yu = fft(seg_no_cp) / sqrt(N);
        Yu = Yu(bin);

        if symIdx == 0
            % TRAINING
            H_est = Yu ./ Xtrain;
            fprintf('[TRAIN] |H| mean = %.3f\n', mean(abs(H_est)));
        else
            % DATA
            Hp = Yu(isP) ./ pilotVal;
            H_all = zeros(numel(kUsed), 1);
            H_all(pilotPos) = Hp;
            H_all(dataPos)  = interp1(pilotPos, Hp, dataPos, 'linear', 'extrap');
            
            Xe = Yu(~isP) ./ H_all(~isP);
            
            % BLIND PHASE TRACKER (BPSK specific)
            % Giai quyet xoay pha do CFO hoac offset
            phi = angle(mean(Xe.^2)) / 2;
            Xe_corr = Xe * exp(-1j * phi);
            
            % Demap BPSK
            rxBits(:, symIdx) = real(Xe_corr) < 0; % 1->-1, 0->1
            nSymDecoded = nSymDecoded + 1;
            
            eqAll = [eqAll; Xe_corr]; %#ok<AGROW>
            
            % In text
            symBits = rxBits(:, symIdx);
            nCharBits = floor(numel(symBits)/8) * 8;
            if nCharBits >= 8
                newChars = char(bin2dec(reshape(char(symBits(1:nCharBits)+'0'), 8, []).').');
                rxTextSoFar = [rxTextSoFar, newChars]; %#ok<AGROW>
                fprintf('[SYM %d] "%s"\n', nSymDecoded, newChars);
            end
        end

        symIdx = symIdx + 1;
        rp     = rp + N + CP;

        drop = rp - 10;
        if drop > 0
            buf(1:drop) = [];
            rp = rp - drop;
        end
    end

    if state == 1 && nSymDecoded >= nSym
        fprintf('\n>>> Da giai ma xong!\n');
        finished = true;
    end
end

%% ==================== KET QUA =======================================
if ~syncFound, error('LOI: Khong tim thay Preamble!'); end

rb = rxBits(:);
if ~isempty(bits)
    rb_msg = rb(1:min(numel(rb), numel(bits)));
    ber = mean(rb_msg ~= bits(1:numel(rb_msg)));
    fprintf('\n>>> BER = %.4g\n', ber);
end

fprintf('\nBan tin thu:\n  "%s"\n', rxTextSoFar);
if isRealtime && ~isempty(y_pass_log)
    audiowrite(rxWavFile, y_pass_log / max(abs(y_pass_log)+eps), Fs_dac, 'BitsPerSample', 24);
end

%% ==================== HIEN THI =====================================
figure('Name', 'Robust Acoustic RX', 'Position', [100 100 1200 600]);

subplot(1,3,1);
plot(metric_buf); grid on;
title('Chirp Sync Metric'); xlabel('Samples'); ylabel('M');

subplot(1,3,2);
plot(real(eqAll), imag(eqAll), '.', 'MarkerSize', 6);
hold on; plot([-1 1], [0 0], 'r+', 'MarkerSize', 12, 'LineWidth', 2);
axis equal; grid on; xlim([-2 2]); ylim([-2 2]);
title(sprintf('BPSK Constellation (BER=%.2g)', ber));
xlabel('I'); ylabel('Q');

subplot(1,3,3);
if exist('H_est', 'var')
    plot(abs(H_est), 'o-'); grid on;
    title('Dap ung kenh (Training)'); xlabel('Subcarrier'); ylabel('|H|');
end
