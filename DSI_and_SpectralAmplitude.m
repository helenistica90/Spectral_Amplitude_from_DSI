% Script to obtain DSI time series from Rayleigh wave ellipticity time series 
% and detect periodicities in different frequency bands.

% REQUIRED DATA:
% The .mat file must contain the following variables:
% 1. **dateX**: A time vector (datetime or numerical format) representing the timestamps of the data.
% 2. **MatFreq**: A matrix where each row corresponds to a specific frequency and each column corresponds to a timestamp, storing ellipticity values.
% 3. **Frecusmedia**: A vector containing the frequencies associated with the rows of MatFreq.

% Load the input .mat file. 
%  IMPORTANT: Replace 'FILE.mat' with your own .mat file containing:
% - dateX: Time vector (datetime or numerical format).
% - MatFreq: Matrix of ellipticity values.
% - Frecusmedia: Frequency vector.


clear all
load('FILE.mat')

% Define parameters
PASOENHORAS = 1;
station = 'TBT';
stackDSI = [];
DSI_semanal = [];
idctiempo = ~isnan(dateX);
dateC = dateX(idctiempo);
HVV = MatFreq(:, idctiempo);
clear MatFreq dateX
dateX = dateC;
MatFreq = HVV;
clear HVV dateC

% Detect sampling interval
intervalo_original = median(diff(dateX)) * 24;

% VENTANAESPECTRO defines the length of the time window for computing the fft, 5*PASOENHORAS, means that it will have a length of 5 days

VENTANAESPECTRO = 5 * PASOENHORAS;

% OVERLAP defines the overlap between consecutive time windows to compute the spectrogram, that will be saved in Espect vector
OVERLAP = round(VENTANAESPECTRO / 2);
tiempo_redondeado = dateshift(datetime(datestr(dateX)), 'start', 'hour');

% Define time range
hora_inicio = min(tiempo_redondeado);
hora_fin = max(tiempo_redondeado);
tiempo_resampleado = (hora_inicio:hours(PASOENHORAS):hora_fin)';

% Define frequency range
FrequencyPairs = [30,40];

for pairIdx = 1
    frecu1 = FrequencyPairs(pairIdx, 1);
    frecu2 = FrequencyPairs(pairIdx, 2);

    idx5 = (Frecusmedia >= frecu1);
    idx6 = (Frecusmedia <= frecu2);
    idx7 = idx5 + idx6;
    idx8 = (idx7 > 1.5);
    Frecus = Frecusmedia(idx8);
    HVZ = MatFreq(idx8, :);
    HV_resampleado = NaN(length(Frecus), length(tiempo_resampleado));

    % Resampling process
    for ij = 1:length(tiempo_resampleado)
        indices = (tiempo_redondeado >= tiempo_resampleado(ij) - hours(PASOENHORAS/2)) & ...
                  (tiempo_redondeado <= tiempo_resampleado(ij) + hours(PASOENHORAS/2));

        if any(indices)
            for kk = 1:length(Frecus)
                AG = HVZ(kk, indices);
                idc2 = AG < 20;
                idc3 = abs(AG) > 0.0;
                idc4 = idc2 + idc3;
                idc5 = idc4 > 1.5;
                if sum(idc5) > 1
                    HV_resampleado(kk, ij) = median(AG(idc5), 'omitnan');
                else
                    HV_resampleado(kk, ij) = NaN;
                end
            end
        end
    end

    % Compute DSI series
    idx_valid = ~isnan(HV_resampleado(1, :));
    hvm = HV_resampleado(:, idx_valid);
    hvref = [];
    
    for hh = 1:size(hvm, 1)
        indis = abs(hvm(hh, :)) > 0.01;
        hvref(hh) = median(hvm(hh, indis), 'omitnan');
    end
    
    idxhvref = ~isnan(hvref);
    DSI = NaN(1, size(HV_resampleado, 2));

    for jj = 1:size(HV_resampleado, 2)
        if sum(idxhvref) >= 2
            hv_data = HV_resampleado(idxhvref, jj);
            idx_valid = ~isnan(hv_data);
            if sum(idx_valid) > 1
                x = find(idx_valid);
                y = hv_data(idx_valid);
                xi = find(~idx_valid);
                hv_data(~idx_valid) = interp1(x, y, xi, 'linear');
            end
            DSI(jj) = 1 - xcorr(hvref(idxhvref), hv_data, 0, 'normalized');
        else
            DSI(jj) = 0;
        end
    end

    DSI = DSI - mean(DSI, 'omitnan');

    % Compute periodicity
    pasitosemana = round((VENTANAESPECTRO * 24) / PASOENHORAS);
    pasitosemana2 = round((24 * OVERLAP) / PASOENHORAS);
    hz = 0;

    for hhu = 1:pasitosemana2:(length(DSI) - pasitosemana)
        ZR = DSI(hhu:hhu + pasitosemana);
        ZRR = ZR - mean(ZR, 'omitnan');
        DSI_semanal = [DSI_semanal; (DSI(hhu:hhu + pasitosemana) - mean(ZR, 'omitnan'))];
        ZRR(isnan(ZRR)) = mean(ZRR, 'omitnan');

        if sum(isnan(HV_resampleado(1, hhu:hhu + pasitosemana))) < (VENTANAESPECTRO * 24 / 3)
            hz = hz + 1;
            [Z(hz, :), Espect(hz, :), Frekis, Ampli_24(hz), F24(hz), Ampli_12(hz), F12(hz), Pnormm(hz, :)] = ...
                estimatePeriodicityWithPeaks(ZRR, PASOENHORAS, 0, VENTANAESPECTRO);
            Tiempo(hz) = tiempo_resampleado(hhu + round(pasitosemana / 2));  % Store time reference
        end
    end
end

% Save results of the DSI time series
save('DSI_results.mat', 'DSI');

% If we are only interested on studying the 12 hour peak, we save only Ampli_12
save('Ampli_12_results.mat', 'Ampli_12');

function [Periodo, P11, f, Ampli24, F24, Ampli12, F12, Pnorm] = estimatePeriodicityWithPeaks(DSI, PASOENHORAS, KEY, VENTANAESPECTRO)
    if length(DSI) > (VENTANAESPECTRO / 1.2)
        N = length(DSI);
        Fs = 1 / (PASOENHORAS * 3600);
        Y = fft(DSI);
        P2 = abs(Y / N);
        P1 = P2(1:N/2+1);
        P1(2:end-1) = 2 * P1(2:end-1);
        f = Fs * (0:(N/2)) / N;
        valid_f = f > 0;
        f = f(valid_f);
        P1 = P1(valid_f);
        Periodo = (1 ./ f) / 3600;
        valid_idx = (Periodo >= 6.2) & (Periodo <= 30);
        Periodo = Periodo(valid_idx);
        P11 = P1(valid_idx);
        Pnorm = P11 / max(P11);

        idx_24 = (Periodo >= 20) & (Periodo <= 30);
        if any(idx_24)
            PP = Periodo(idx_24);
            PP11 = P11(idx_24);
            [~, idx_closest] = min(abs(PP - 24));
            Ampli24 = PP11(idx_closest);
            F24 = PP(idx_closest);
        else
            Ampli24 = NaN;
            F24 = NaN;
        end

        idx_12 = (Periodo >= 9) & (Periodo <= 15);
        if any(idx_12)
            PP = Periodo(idx_12);
            PP11 = P11(idx_12);
            [~, idx_closest] = min(abs(PP - 12));
            Ampli12 = PP11(idx_closest);
            F12 = PP(idx_closest);
        else
            Ampli12 = NaN;
            F12 = NaN;
        end
    else
        Periodo = NaN; P11 = NaN; f = NaN;
        Ampli24 = NaN; F24 = NaN; Ampli12 = NaN; F12 = NaN;
        Pnorm = NaN;
    end
end
