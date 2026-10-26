clc; clear; close all;

Nfft = 2048;                    % Размер FFT (по заданию)
M_data = 16;                    % Модуляция данных: 16-QAM
M_pilot = 4;                    % Модуляция пилотов: QPSK 
CP_len = round(Nfft/4);         % Длина CP = 1/4 * Nfft = 512
T_sym = 100e-6;                 % Длительность OFDM-символа 100 мкс (по заданию)
Fs = Nfft / T_sym;              % Частота дискретизации = 20.48 МГц
SNR_dB = 1:2:30;                % Диапазон SNR (дБ)
NumOfdmSymbols = 200;           % Количество OFDM символов

% ПРОФИЛИ КАНАЛОВ EPA, EVA, ETU 
% Задержки в наносекундах
EPA_delays_ns = [0, 30, 70, 80, 110, 190, 410];
EVA_delays_ns = [0, 30, 150, 310, 370, 710, 1090, 1730, 2510];
ETU_delays_ns = [0, 50, 120, 200, 230, 500, 1600, 2300, 5000];

% Средние относительные мощности в дБ
EPA_gains_dB = [0, -1, -2, -3, -8, -17.2, -20.8];
EVA_gains_dB = [0, -1.5, -1.4, -3.6, -0.6, -9.1, -7, -12, -16.9];
ETU_gains_dB = [-1, -1, -1, 0, 0, 0, -3, -5, -7];

% Перевод задержек в секунды
EPA_delays = EPA_delays_ns * 1e-9;
EVA_delays = EVA_delays_ns * 1e-9;
ETU_delays = ETU_delays_ns * 1e-9;

% Линейные мощности (нормированные)
EPA_gains_lin = 10.^(EPA_gains_dB/10);
EVA_gains_lin = 10.^(EVA_gains_dB/10);
ETU_gains_lin = 10.^(ETU_gains_dB/10);

EPA_gains_lin = EPA_gains_lin / sum(EPA_gains_lin);
EVA_gains_lin = EVA_gains_lin / sum(EVA_gains_lin);
ETU_gains_lin = ETU_gains_lin / sum(ETU_gains_lin);

channel_names = {'EPA', 'EVA', 'ETU'};
colors_ch = {'b', 'r', 'g'};

%% ДОП ЗАДАНИЕ 2.1: ПОСТРОЕНИЕ ИХ И ЧХ КАНАЛОВ EPA, EVA, ETU
fprintf('\nДОП ЗАДАНИЕ 2.1: ИХ и ЧХ каналов EPA, EVA, ETU\n');

channels_delays = {EPA_delays, EVA_delays, ETU_delays};
channels_gains = {EPA_gains_lin, EVA_gains_lin, ETU_gains_lin};

figure(14);

for ch_idx = 1:3
    tau = channels_delays{ch_idx};
    P_lin = channels_gains{ch_idx};

    tau_max_samples = round(tau(end) * Fs);
    h_len = max(tau_max_samples + 1, 10);
    h = zeros(1, h_len);

    for p = 1:length(tau)
        idx = round(tau(p) * Fs) + 1;
        if idx <= h_len
            h(idx) = sqrt(P_lin(p));
        end
    end

    Nfft_ch = 2^14;
    H = fftshift(fft(h, Nfft_ch));
    f_ch = (-Nfft_ch/2 : Nfft_ch/2-1) * (Fs / Nfft_ch);
    H_db = 20 * log10(abs(H) + eps);

    tau_avg = sum(P_lin .* tau) / sum(P_lin);
    tau_rms = sqrt(sum(P_lin .* (tau - tau_avg).^2) / sum(P_lin));
    Bc = 1 / tau_rms;

    subplot(3,1,ch_idx);
    stem((0:h_len-1)/Fs*1e9, abs(h), [colors_ch{ch_idx}, '.'], 'MarkerSize', 6, 'LineWidth', 1.5);
    title(sprintf('Импульсная характеристика канала %s: τ_rms=%.2f нс, B_c=%.2f кГц', ...
        channel_names{ch_idx}, tau_rms*1e9, Bc/1e3));
    xlabel('Задержка, нс'); 
    ylabel('|h(t)|');
    xlim([0, tau(end)*1e9*1.2]);
    grid on;

    fprintf('Канал %s: τ_rms = %.2f нс, B_c = %.2f кГц, h_len = %d отсчетов\n', ...
        channel_names{ch_idx}, tau_rms*1e9, Bc/1e3, h_len);
end

figure(15);
hold on; grid on;

for ch_idx = 1:3
    tau = channels_delays{ch_idx};
    P_lin = channels_gains{ch_idx};

    tau_max_samples = round(tau(end) * Fs);
    h_len = max(tau_max_samples + 1, 10);
    h = zeros(1, h_len);

    for p = 1:length(tau)
        idx = round(tau(p) * Fs) + 1;
        if idx <= h_len
            h(idx) = sqrt(P_lin(p));
        end
    end

    Nfft_ch = 2^14;
    H = fftshift(fft(h, Nfft_ch));
    f_ch = (-Nfft_ch/2 : Nfft_ch/2-1) * (Fs / Nfft_ch);
    H_db = 20 * log10(abs(H) + eps);

    tau_avg = sum(P_lin .* tau) / sum(P_lin);
    tau_rms = sqrt(sum(P_lin .* (tau - tau_avg).^2) / sum(P_lin));
    Bc = 1 / tau_rms;

    plot(f_ch/1e6, H_db, [colors_ch{ch_idx}, '-'], 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%s: τ_rms=%.1f нс, B_c=%.1f кГц', ...
        channel_names{ch_idx}, tau_rms*1e9, Bc/1e3));
end

xlabel('Частота, МГц');
ylabel('|H(f)|, дБ');
title('Частотные характеристики каналов EPA, EVA, ETU');
legend('Location', 'best');
xlim([-Fs/2/1e6, Fs/2/1e6]);
ylim([-60, 5]);
hold off;

%% ДОП ЗАДАНИЕ 2.2: АНАЛИТИЧЕСКАЯ ОЦЕНКА СЕЛЕКТИВНОСТИ
fprintf('\nДОП ЗАДАНИЕ 2.2: Анализ селективности\n');

B_signal = 1 / T_sym;
B_signal_total = Fs;

fprintf('Полоса OFDM сигнала (одна поднесущая): B = %.2f кГц\n', B_signal/1e3);
fprintf('Полная полоса OFDM сигнала: B_total = %.2f МГц\n', B_signal_total/1e6);
fprintf('\n');

for ch_idx = 1:3
    tau = channels_delays{ch_idx};
    P_lin = channels_gains{ch_idx};

    tau_avg = sum(P_lin .* tau) / sum(P_lin);
    tau_rms = sqrt(sum(P_lin .* (tau - tau_avg).^2) / sum(P_lin));
    Bc = 1 / tau_rms;

    fprintf('Канал %s:\n', channel_names{ch_idx});
    fprintf('  τ_rms = %.2f нс\n', tau_rms*1e9);
    fprintf('  B_c = %.2f кГц\n', Bc/1e3);

    if B_signal_total < Bc
        fprintf('  B_total < B_c => ПЛОСКИЕ замирания\n');
    else
        fprintf('  B_total > B_c => ЧАСТОТНО-СЕЛЕКТИВНЫЕ замирания\n');
    end

    if tau_rms < T_sym
        fprintf('  τ_rms < T_sym => МСИ устраняется CP\n');
    else
        fprintf('  τ_rms > T_sym => МСИ не устраняется полностью\n');
    end
    fprintf('\n');
end

%% ДОП ЗАДАНИЕ 2.3-2.4: BER В КАНАЛАХ С ЗАМИРАНИЯМИ (EPA, EVA, ETU)
fprintf('ДОП ЗАДАНИЕ 2.3-2.4: BER в каналах с замираниями\n');

NumSC_fad = 1200;
NumGuardSC_fad = round((Nfft - NumSC_fad)/2);

pilot_step_fad = 16;
pilot_idx_fad = 1:pilot_step_fad:NumSC_fad;
data_idx_fad = setdiff(1:NumSC_fad, pilot_idx_fad);

nDataSymbols_fad = length(data_idx_fad);
nPilots_fad = length(pilot_idx_fad);

rayleighchan_EPA = comm.RayleighChannel( ...
    'SampleRate', Fs, ...
    'PathDelays', EPA_delays, ...
    'AveragePathGains', EPA_gains_dB, ...
    'MaximumDopplerShift', 0, ...
    'NormalizePathGains', true);

rayleighchan_EVA = comm.RayleighChannel( ...
    'SampleRate', Fs, ...
    'PathDelays', EVA_delays, ...
    'AveragePathGains', EVA_gains_dB, ...
    'MaximumDopplerShift', 0, ...
    'NormalizePathGains', true);

rayleighchan_ETU = comm.RayleighChannel( ...
    'SampleRate', Fs, ...
    'PathDelays', ETU_delays, ...
    'AveragePathGains', ETU_gains_dB, ...
    'MaximumDopplerShift', 0, ...
    'NormalizePathGains', true);

channels_fad = {rayleighchan_EPA, rayleighchan_EVA, rayleighchan_ETU};

BER_fading = zeros(length(channels_fad), length(SNR_dB));

for ch_idx = 1:length(channels_fad)
    fprintf('Канал %s:\n', channel_names{ch_idx});

    for k = 1:length(SNR_dB)
        errorBits = 0;
        totalBits = NumOfdmSymbols * nDataSymbols_fad * log2(M_data);

        for sym = 1:NumOfdmSymbols
            data_bits = randi([0 1], nDataSymbols_fad * log2(M_data), 1);
            data_symbols = bi2de(reshape(data_bits, log2(M_data), []).', 'left-msb');
            mod_data = qammod(data_symbols, M_data, 'UnitAveragePower', true);

            pilot_bits = randi([0 1], nPilots_fad * log2(M_pilot), 1);
            pilot_symbols = bi2de(reshape(pilot_bits, log2(M_pilot), []).', 'left-msb');
            mod_pilots = qammod(pilot_symbols, M_pilot, 'UnitAveragePower', true);

            Sp = zeros(Nfft, 1);
            Sp(NumGuardSC_fad + data_idx_fad) = mod_data;
            Sp(NumGuardSC_fad + pilot_idx_fad) = mod_pilots;

            Sig = ifft(ifftshift(Sp));
            Sig_CP = [Sig(end-CP_len+1:end); Sig];

            release(channels_fad{ch_idx});
            Rx_Signal = channels_fad{ch_idx}(Sig_CP);

            Ps = mean(abs(Rx_Signal).^2);
            SNR_linear = 10^(SNR_dB(k)/10);
            N0 = Ps / SNR_linear;
            sigma_noise = sqrt(N0/2);
            noise = sigma_noise * (randn(size(Rx_Signal)) + 1j*randn(size(Rx_Signal)));
            Rx_Signal = Rx_Signal + noise;

            Sig_rec = Rx_Signal(CP_len+1:CP_len+Nfft);
            Sp_rec = fftshift(fft(Sig_rec));

            H_pilots = Sp_rec(NumGuardSC_fad + pilot_idx_fad) ./ mod_pilots;
            H_est = interp1(pilot_idx_fad, H_pilots, (1:NumSC_fad)', 'linear', 'extrap');

            Sp_eq = Sp_rec(NumGuardSC_fad + data_idx_fad) ./ H_est(data_idx_fad);

            demod_symbols = qamdemod(Sp_eq, M_data, 'UnitAveragePower', true);
            bits_rx = reshape(de2bi(demod_symbols, log2(M_data), 'left-msb').', [], 1);

            errorBits = errorBits + sum(data_bits ~= bits_rx);
        end

        BER_fading(ch_idx, k) = errorBits / totalBits;
    end

    fprintf('  BER при SNR=10 дБ = %e\n', BER_fading(ch_idx, 6));
end

figure(16);
hold on; grid on;
markers_ch = {'o-', 's-', 'd-'};

for ch_idx = 1:length(channels_fad)
    plot(SNR_dB, BER_fading(ch_idx,:), [colors_ch{ch_idx} markers_ch{ch_idx}], ...
        'LineWidth', 1.5, 'MarkerSize', 6, ...
        'DisplayName', sprintf('%s (с эквалайзером)', channel_names{ch_idx}));
end

plot(SNR_dB, berawgn(SNR_dB, 'qam', M_data), '--k', 'LineWidth', 2, ...
    'DisplayName', 'Теория AWGN (berawgn)');

xlabel('SNR (дБ)');
ylabel('BER');
title('Помехоустойчивость OFDM в каналах с замираниями (EPA, EVA, ETU)');
legend('Location', 'southwest');
set(gca, 'YScale', 'log');
ylim([1e-5 1]);
hold off;

%% ДОП ЗАДАНИЕ 2.5: СОЗВЕЗДИЕ ДО И ПОСЛЕ ЭКВАЛАЙЗЕРА
fprintf('\nСозвездие до и после эквалайзера (канал EPA, SNR=15 дБ)\n');

SNR_test = 15;
SNR_linear_test = 10^(SNR_test/10);

release(channels_fad{1});

constellation_before_eq = [];
constellation_after_eq = [];
data_all = [];
errorBits_const = 0;
totalBits_const = NumOfdmSymbols * nDataSymbols_fad * log2(M_data);

for sym = 1:NumOfdmSymbols
    data_bits = randi([0 1], nDataSymbols_fad * log2(M_data), 1);
    data_symbols = bi2de(reshape(data_bits, log2(M_data), []).', 'left-msb');
    mod_data = qammod(data_symbols, M_data, 'UnitAveragePower', true);

    pilot_bits = randi([0 1], nPilots_fad * log2(M_pilot), 1);
    pilot_symbols = bi2de(reshape(pilot_bits, log2(M_pilot), []).', 'left-msb');
    mod_pilots = qammod(pilot_symbols, M_pilot, 'UnitAveragePower', true);

    Sp = zeros(Nfft, 1);
    Sp(NumGuardSC_fad + data_idx_fad) = mod_data;
    Sp(NumGuardSC_fad + pilot_idx_fad) = mod_pilots;

    Sig = ifft(ifftshift(Sp));
    Sig_CP = [Sig(end-CP_len+1:end); Sig];

    Rx_Signal = channels_fad{1}(Sig_CP);

    Ps = mean(abs(Rx_Signal).^2);
    N0 = Ps / SNR_linear_test;
    sigma_noise = sqrt(N0/2);
    noise = sigma_noise * (randn(size(Rx_Signal)) + 1j*randn(size(Rx_Signal)));
    Rx_Signal = Rx_Signal + noise;

    Sig_rec = Rx_Signal(CP_len+1:CP_len+Nfft);
    Sp_rec = fftshift(fft(Sig_rec));

    if sym <= 5
        constellation_before_eq = [constellation_before_eq; Sp_rec(NumGuardSC_fad + data_idx_fad)];
    end

    H_pilots = Sp_rec(NumGuardSC_fad + pilot_idx_fad) ./ mod_pilots;
    H_est = interp1(pilot_idx_fad, H_pilots, (1:NumSC_fad)', 'linear', 'extrap');

    Sp_eq = Sp_rec(NumGuardSC_fad + data_idx_fad) ./ H_est(data_idx_fad);

    if sym <= 5
        constellation_after_eq = [constellation_after_eq; Sp_eq];
        data_all = [data_all; mod_data];
    end

    demod_symbols = qamdemod(Sp_eq, M_data, 'UnitAveragePower', true);
    bits_rx = reshape(de2bi(demod_symbols, log2(M_data), 'left-msb').', [], 1);
    errorBits_const = errorBits_const + sum(data_bits ~= bits_rx);
end

BER_const = errorBits_const / totalBits_const;
fprintf('BER при SNR=%d дБ в канале EPA = %e\n', SNR_test, BER_const);

figure(17);
subplot(1,2,1);
plot(real(constellation_before_eq), imag(constellation_before_eq), 'b.', 'MarkerSize', 3);
title('Созвездие ДО эквалайзера (канал EPA)');
xlabel('I'); ylabel('Q');
grid on; axis equal;
xlim([-3, 3]); ylim([-3, 3]);

subplot(1,2,2);
plot(real(constellation_after_eq), imag(constellation_after_eq), 'r.', 'MarkerSize', 3);
hold on;
ideal_const = qammod(0:15, 16, 'UnitAveragePower', true);
plot(real(ideal_const), imag(ideal_const), 'kx', 'MarkerSize', 10, 'LineWidth', 2);
title(sprintf('Созвездие ПОСЛЕ эквалайзера (BER=%.2e)', BER_const));
xlabel('I'); ylabel('Q');
grid on; axis equal;
xlim([-3, 3]); ylim([-3, 3]);
legend('Принятые символы', 'Идеальное созвездие');
hold off;

