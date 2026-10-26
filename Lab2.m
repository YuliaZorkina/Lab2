clc; clear; close all;

Nfft = 2048;                    % Размер FFT 
M_data = 16;                    % Модуляция данных: 16-QAM
M_pilot = 4;                    % Модуляция пилотов: QPSK 
CP_len = round(Nfft/4);         % Длина CP = 1/4 * Nfft = 512
T_sym = 100e-6;                 % Длительность OFDM-символа 100 мкс 
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
%% ДОП ЗАДАНИЕ 2.6: ДОПЛЕРОВСКИЕ СПЕКТРЫ EPA5, EVA70, ETU300
fprintf('\nДОП ЗАДАНИЕ 2.6: Доплеровские спектры\n');

fD_values_dop = [5 70 300];
N_symbols_doppler = [100000, 30000, 20000];

channel_names_dop = {'EPA5','EVA70','ETU300'};
delays_dop = {EPA_delays, EVA_delays, ETU_delays};
gains_dop = {EPA_gains_dB, EVA_gains_dB, ETU_gains_dB}; 
colors_dop = {'b','r','g'};

Fs_doppler = 1 / T_sym;  

Nfft_doppler = 2^18;
window_len_EPA = 8192;
window_len_default = 2048;
win_EPA = hamming(window_len_EPA);
win_default = hamming(window_len_default);
noverlap_EPA = window_len_EPA/2;
noverlap_default = window_len_default/2;

figure(18);
clf;
hold on;
grid on;

for idx = 1:length(fD_values_dop)
    fD = fD_values_dop(idx);
    N_sym = N_symbols_doppler(idx);
    
    if idx == 1
        win = win_EPA;
        noverlap = noverlap_EPA;
    else
        win = win_default;
        noverlap = noverlap_default;
    end

    chan = comm.RayleighChannel( ...
        'SampleRate', Fs_doppler, ...
        'PathDelays', delays_dop{idx}, ...
        'AveragePathGains', gains_dop{idx}, ...
        'MaximumDopplerShift', fD, ...
        'NormalizePathGains', true, ...
        'DopplerSpectrum', doppler('Jakes'), ...
        'PathGainsOutputPort', true);

    release(chan);
    [~, pathGains] = chan(ones(N_sym,1));

    h_dop = pathGains(:,1);

    [Pxx, f] = pwelch(h_dop, win, noverlap, ...
        Nfft_doppler, Fs_doppler, 'centered');

    Pxx_dB = 10*log10(Pxx/max(Pxx) + eps);

    mask = abs(f) <= min(Fs_doppler/2, 1.5*fD);

    plot(f(mask), Pxx_dB(mask), colors_dop{idx}, ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('%s (f_D=%d Гц)', ...
        channel_names_dop{idx}, fD));
end

xlabel('Доплеровская частота, Гц');
ylabel('Нормированная PSD, дБ');
title('Доплеровские спектры каналов (модель Джейкса)');
legend('Location','best');
xlim([-500 500]);
ylim([-50 5]);
hold off;

fprintf('Доплеровские спектры построены.\n');
%% ДОП ЗАДАНИЕ 2.7: ВРЕМЯ КОГЕРЕНТНОСТИ
fprintf('\nДОП ЗАДАНИЕ 2.7: Время когерентности\n');

fD_EPA5 = 5;
fD_EVA70 = 70;
fD_ETU300 = 300;

Tc_EPA5 = 9 / (16 * pi * fD_EPA5);
Tc_EVA70 = 9 / (16 * pi * fD_EVA70);
Tc_ETU300 = 9 / (16 * pi * fD_ETU300);

fprintf('Канал EPA5 (fD=5 Гц): Tc = %.2f мс\n', Tc_EPA5 * 1000);
fprintf('Канал EVA70 (fD=70 Гц): Tc = %.2f мс\n', Tc_EVA70 * 1000);
fprintf('Канал ETU300 (fD=300 Гц): Tc = %.2f мс\n', Tc_ETU300 * 1000);

fprintf('\nДлительность OFDM-символа: T_sym = %.2f мс\n', T_sym * 1000);
fprintf('Для EPA5: Tc/T_sym = %.1f => медленные замирания\n', Tc_EPA5/T_sym);
fprintf('Для EVA70: Tc/T_sym = %.2f\n', Tc_EVA70/T_sym);
fprintf('Для ETU300: Tc/T_sym = %.3f => быстрые замирания\n', Tc_ETU300/T_sym);

%% ДОП ЗАДАНИЕ 2.8: MSE ОЦЕНКИ КАНАЛА ОТ ОСШ
fprintf('\nДОП ЗАДАНИЕ 2.8: MSE оценки канала\n');

fD_test = 10;
NumSC_mse = 1200;
NumGuardSC_mse = round((Nfft - NumSC_mse)/2);
pilot_step_mse = 16;
pilot_idx_mse = 1:pilot_step_mse:NumSC_mse;
data_idx_mse = setdiff(1:NumSC_mse, pilot_idx_mse);  % информационные поднесущие
nPilots_mse = length(pilot_idx_mse);

rayleighchan_MSE = comm.RayleighChannel( ...
    'SampleRate', Fs, ...
    'PathDelays', EVA_delays, ...
    'AveragePathGains', EVA_gains_dB, ...
    'MaximumDopplerShift', fD_test, ...
    'NormalizePathGains', true);

SNR_MSE_dB = 0:2:30;
MSE_results = zeros(length(SNR_MSE_dB), 1);
N_MSE_symbols = 500; 

for k = 1:length(SNR_MSE_dB)
    MSE_sum = 0;

    for sym = 1:N_MSE_symbols
        % Генерируем пилоты
        pilot_bits = randi([0 1], nPilots_mse * log2(M_pilot), 1);
        pilot_symbols = bi2de(reshape(pilot_bits, log2(M_pilot), []).', 'left-msb');
        mod_pilots = qammod(pilot_symbols, M_pilot, 'UnitAveragePower', true);

        % Формируем OFDM символ только с пилотами
        Sp = zeros(Nfft, 1);
        Sp(NumGuardSC_mse + pilot_idx_mse) = mod_pilots;

        Sig = ifft(ifftshift(Sp));
        Sig_CP = [Sig(end-CP_len+1:end); Sig];

        % Передаем через канал БЕЗ шума 
        release(rayleighchan_MSE);
        Rx_ideal = rayleighchan_MSE(Sig_CP);

        % Добавляем шум
        Ps = mean(abs(Rx_ideal).^2);
        SNR_linear = 10^(SNR_MSE_dB(k)/10);
        N0 = Ps / SNR_linear;
        sigma_noise = sqrt(N0/2);
        noise = sigma_noise * (randn(size(Rx_ideal)) + 1j*randn(size(Rx_ideal)));
        Rx_noisy = Rx_ideal + noise;

        % Обрабатываем идеальный сигнал
        Sig_rec_ideal = Rx_ideal(CP_len+1:CP_len+Nfft);
        Sp_rec_ideal = fftshift(fft(Sig_rec_ideal));
        H_true_pilots = Sp_rec_ideal(NumGuardSC_mse + pilot_idx_mse) ./ mod_pilots;

        % Обрабатываем зашумленный сигнал - оценка канала
        Sig_rec_noisy = Rx_noisy(CP_len+1:CP_len+Nfft);
        Sp_rec_noisy = fftshift(fft(Sig_rec_noisy));
        H_est_pilots = Sp_rec_noisy(NumGuardSC_mse + pilot_idx_mse) ./ mod_pilots;

        % Интерполируем оценку и истинный канал на все активные поднесущие
        H_true_all = interp1(pilot_idx_mse, H_true_pilots, (1:NumSC_mse)', 'linear', 'extrap');
        H_est_all  = interp1(pilot_idx_mse, H_est_pilots,  (1:NumSC_mse)', 'linear', 'extrap');

        % MSE по информационным поднесущим
        MSE_sum = MSE_sum + mean(abs(H_est_all(data_idx_mse) - H_true_all(data_idx_mse)).^2);
    end

    MSE_results(k) = MSE_sum / N_MSE_symbols;
end

figure(19);
semilogy(SNR_MSE_dB, MSE_results, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('SNR (дБ)');
ylabel('MSE оценки канала');
title(sprintf('Зависимость MSE оценки канала от ОСШ (канал EVA, fD=%d Гц)', fD_test));
grid on;

fprintf('MSE оценки канала рассчитаны.\n');
