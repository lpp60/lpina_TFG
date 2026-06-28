% =========================================================================
% TFG - FASE 2: Obtención de Señal Cinemática Equivalente
% Análisis de un paciente representativo (Ejemplo para grado 4 - sujeto 4)
% =========================================================================
clear; clc; close all;

%% 1. DEFINIR RUTAS
% GRADO 4
% Archivo Muñeca
archivo_W = '';
% Archivo Dedo
archivo_F = '';

%% 2. IMPORTACIÓN Y TRATAMIENTO TEMPORAL
% Importamos la Muñeca
opts_W = detectImportOptions(archivo_W, 'Delimiter', ';');
datos_W = readtable(archivo_W, opts_W);
% Importamos el Dedo
opts_F = detectImportOptions(archivo_F, 'Delimiter', ';');
datos_F = readtable(archivo_F, opts_F);

% Vectores de tiempo absolutos en segundos
t_W = seconds(duration(datos_W.Time));
t_F = seconds(duration(datos_F.Time));

% Aislamos únicamente el Eje Vertical (Z)
acc_W_Z_raw = datos_W.Z;
acc_F_Z_raw = datos_F.Z;

%% 3. REMUESTREO (Dedo a 100 Hz con filtro anti-aliasing)
% Remuestreo con filtro pasabajos integrado (de 2500 Hz a 100 Hz)
temp_Z = resample(acc_F_Z_raw, 100, 2500);

% Ajustamos la longitud para que coincida exactamente con la señal de la muñeca
N_W = length(t_W);
acc_F_Z_100 = temp_Z(1:N_W);

% Establecemos el tiempo común y la Fs común para el resto del código
t = t_W;
acc_W_Z_100 = acc_W_Z_raw;
Fs = 100;

%% 4. PREPROCESAMIENTO: RECORTE Y FILTRADO
% Recorte de 3 segundos iniciales
idx = t >= 3;
t_valido = t(idx) - 3;
acc_W_Z_rec = acc_W_Z_100(idx);
acc_F_Z_rec = acc_F_Z_100(idx);

% Comprobación del sentido de ambos ejes
disp(['Gravedad media muñeca: ', num2str(mean(acc_W_Z_raw))]);
disp(['Gravedad media dedo: ', num2str(mean(acc_F_Z_raw))]);
disp(' '); 

% Filtrado Butterworth (3-15 Hz) para eliminar la gravedad y centrar en cero
[b, a] = butter(2, [3, 15] / (Fs / 2), 'bandpass');
acc_W_Z_filt = filtfilt(b, a, acc_W_Z_rec);
acc_F_Z_filt = filtfilt(b, a, acc_F_Z_rec);

%% 5. CÁLCULO DE LA ACELERACIÓN DIFERENCIAL
acc_rel_Z = acc_F_Z_filt - acc_W_Z_filt;

%% 6. COMPROBACIÓN TEMPORAL (Gráficas)
figure('Name', 'Dominio del Tiempo: Aceleraciones Verticales', 'Color', 'w', 'Position', [100, 100, 1000, 600]);

subplot(2,1,1);
plot(t_valido, acc_W_Z_filt, 'Color', [1.0 0.6 0.8], 'LineWidth', 1.2); hold on;
plot(t_valido, acc_F_Z_filt, 'Color', [0.95 0.85 0.35], 'LineWidth', 1.2);
title('Aceleraciones Absolutas Preprocesadas (Eje Z)');
ylabel('Aceleración (g)'); legend('Muñeca (a_{W,z})', 'Dedo (a_{F,z})'); grid on;

subplot(2,1,2);
plot(t_valido, acc_rel_Z, 'Color', [0.4 0.75 0.95], 'LineWidth', 1.2);
title('Aceleración Diferencial (a_{rel,z} = a_{F,z} - a_{W,z})');
xlabel('Tiempo (s)'); ylabel('Aceleración (g)'); grid on;

%% 7. EXTRACCIÓN DE MÉTRICAS
% Función para métricas temporales (RMS y Pico-Pico)
calc_temp = @(x) [rms(x), max(x) - min(x)];

met_W = calc_temp(acc_W_Z_filt);
met_F = calc_temp(acc_F_Z_filt);
met_Rel = calc_temp(acc_rel_Z);

% Configuración Welch (2.56 s)
Tw = 2.56; 
Nw = round(Tw * Fs); 
win = hamming(Nw); 
noverlap = round(0.5 * Nw); 
nfft = 2^(nextpow2(Nw) + 1);

% Cálculo PSD
[P_W, f] = pwelch(acc_W_Z_filt, win, noverlap, nfft, Fs);
[P_F, ~] = pwelch(acc_F_Z_filt, win, noverlap, nfft, Fs);
[P_Rel, ~] = pwelch(acc_rel_Z, win, noverlap, nfft, Fs);

% Aislar banda 3-15 Hz
idx_b = f >= 3 & f <= 15;
f_b = f(idx_b); P_W = P_W(idx_b); P_F = P_F(idx_b); P_Rel = P_Rel(idx_b);

% Extraer métricas espectrales
calc_frec = @(P, f_vec) [f_vec(P == max(P)), max(P), trapz(f_vec, P)]; % [Frec_Dom, Pot_Max, Area]

metF_W = calc_frec(P_W, f_b);
metF_F = calc_frec(P_F, f_b);
metF_Rel = calc_frec(P_Rel, f_b);

% Cálculo del ancho de banda
% Muñeca
[pks_W, ~, w_W] = findpeaks(P_W, f_b, 'WidthReference', 'halfheight');
bw_W = w_W(pks_W == max(pks_W)); if isempty(bw_W), bw_W = NaN; end

% Dedo
[pks_F, ~, w_F] = findpeaks(P_F, f_b, 'WidthReference', 'halfheight');
bw_F = w_F(pks_F == max(pks_F)); if isempty(bw_F), bw_F = NaN; end

% Relativa (Diferencial)
[pks_Rel, ~, w_Rel] = findpeaks(P_Rel, f_b, 'WidthReference', 'halfheight');
bw_Rel = w_Rel(pks_Rel == max(pks_Rel)); if isempty(bw_Rel), bw_Rel = NaN; end

% Tabla Resumen
Nombres = {'Muñeca', 'Dedo', 'Relativa'};
Tabla_2A = table([met_W(1); met_F(1); met_Rel(1)], ...
                 [met_W(2); met_F(2); met_Rel(2)], ...
                 [metF_W(1); metF_F(1); metF_Rel(1)], ...
                 [metF_W(2); metF_F(2); metF_Rel(2)], ...
                 [metF_W(3); metF_F(3); metF_Rel(3)], ...
                 [bw_W(1); bw_F(1); bw_Rel(1)], ... % (1) por seguridad matricial
    'VariableNames', {'RMS', 'Pico_a_Pico', 'Frec_Dom_Hz', 'Potencia_Max', 'Area_PSD', 'Ancho_Banda_Hz'}, ...
    'RowNames', Nombres);

disp('--- RESULTADOS SEÑAL DIFERENCIAL ---');
disp(' ');
disp(Tabla_2A);

%% COMPROBACIÓN METODOLÓGICA
% Cálculo de la densidad espectral cruzada (Cross-Spectral Density)
[P_FW, ~] = cpsd(acc_F_Z_filt, acc_W_Z_filt, win, noverlap, nfft, Fs);

% Aislar la banda de 3-15 Hz para el espectro cruzado
P_FW = P_FW(idx_b);

% Cálculo de la PSD relativa teórica mediante la fórmula de espectros cruzados:
% S_rel(f) = S_FF(f) + S_WW(f) - 2*Re[S_FW(f)]
S_rel_teorica = P_F + P_W - 2 * real(P_FW);

%% 8. COMPROBACIÓN FRECUENCIAL
figure('Name', 'Dominio de la Frecuencia: PSD (Welch)', 'Color', 'w', 'Position', [150, 150, 800, 500]);
plot(f_b, P_W, 'Color', [1.0 0.6 0.8], 'LineWidth', 1.5); hold on;
plot(f_b, P_F, 'Color', [0.95 0.85 0.35], 'LineWidth', 1.5);
plot(f_b, P_Rel, 'Color', [0.4 0.75 0.95], 'LineWidth', 1.5, 'LineStyle', '--');
title('Comparativa de Densidad Espectral de Potencia (Eje Z)');
xlabel('Frecuencia (Hz)'); ylabel('Potencia (g^2/Hz)');
legend('Muñeca (a_{W,z})', 'Dedo (a_{F,z})', 'Relativa (a_{rel,z})'); grid on;

%% 9. VALIDACIÓN METODOLÓGICA DE LA PSD RELATIVA
figure('Name', 'Validación Metodológica: PSD Relativa', 'Color', 'w', 'Position', [200, 200, 800, 500]);
plot(f_b, P_Rel, 'Color', [0.4 0.75 0.95], 'LineWidth', 3); hold on;
plot(f_b, S_rel_teorica, 'Color', [0.8 0.2 0.2], 'LineWidth', 1.5, 'LineStyle', '--');
title('Validación de la Densidad Espectral de Potencia Diferencial');
xlabel('Frecuencia (Hz)'); ylabel('Potencia (g^2/Hz)');
legend('PSD de la resta temporal (a_{rel,z})', 'PSD estimada por espectros cruzados');
grid on;

% Comprobación numérica del error medio entre ambas curvas
error_psd = mean(abs(P_Rel - S_rel_teorica));
disp(['Error medio entre métodos de estimación de PSD relativa: ', num2str(error_psd)]);

%% OBTENCIÓN DEL MOMENTO EXCITADOR EQUIVALENTE (M_tremor)

% 1. PARÁMETROS DEL MODELO (Valores promedio preliminares)
L = 0.15;   % Longitud efectiva muñeca-dedo (m)
I = 0.0025;  % Momento de inercia equivalente (kg*m^2)
c0 = 0.03;  % Amortiguamiento rotacional equivalente (N*m*s/rad)
k0 = 0.65;   % Rigidez rotacional equivalente (N*m/rad)

% 2. CONVERSIÓN A ACELERACIÓN ANGULAR
% Convertimos la aceleración relativa de 'g' a m/s^2
g_val = 9.81; 
a_rel_Z_ms2 = acc_rel_Z * g_val; 

% Calculamos la aceleración angular (rad/s^2)
theta_ddot = a_rel_Z_ms2 / L;

% 3. INTEGRACIÓN CONTROLADA EN LA BANDA DEL TEMBLOR (Dominio del Tiempo)
% Evita la deriva integrando con cumtrapz y aplicando el filtro paso banda 
% definido en la Fase 1 para centrar la señal en cero

% A) Velocidad angular
theta_dot_raw = cumtrapz(t_valido, theta_ddot);
theta_dot = filtfilt(b, a, theta_dot_raw); % Limpieza de la deriva

% B) Desplazamiento angular
theta_raw = cumtrapz(t_valido, theta_dot);
theta = filtfilt(b, a, theta_raw); % Limpieza de la deriva

% 4. CÁLCULO DEL MOMENTO EXCITADOR EQUIVALENTE
% Ecuación: M_tremor(t) = I*theta_ddot(t) + c0*theta_dot(t) + k0*theta(t)
M_tremor = I * theta_ddot + c0 * theta_dot + k0 * theta;

% 5. CÁLCULO DE LA PSD DEL MOMENTO EXCITADOR
[P_M, f_M] = pwelch(M_tremor, win, noverlap, nfft, Fs);

% Aislar banda 3-15 Hz para el momento
idx_b_M = f_M >= 3 & f_M <= 15;
f_b_M = f_M(idx_b_M); 
P_M = P_M(idx_b_M);

%% 6. REPRESENTACIÓN GRÁFICA
% Gráfica temporal
figure('Name', 'Fase 2B: Cinemática y Dinámica Equivalente', 'Color', 'w', 'Position', [100, 100, 1200, 800]);

subplot(4,1,1);
plot(t_valido, theta_ddot, 'Color', [0.95 0.85 0.35] , 'LineWidth', 1.2);
title('Aceleración Angular (\theta^{\prime\prime})');
ylabel('rad/s^2'); grid on;

subplot(4,1,2);
plot(t_valido, theta_dot, 'Color', [0.4 0.75 0.95], 'LineWidth', 1.2);
title('Velocidad Angular (\theta^{\prime})');
ylabel('rad/s'); grid on;

subplot(4,1,3);
plot(t_valido, theta, 'Color', [1.0 0.6 0.8], 'LineWidth', 1.2);
title('Desplazamiento Angular (\theta)');
ylabel('rad'); grid on;

subplot(4,1,4);
plot(t_valido, M_tremor, 'Color', [0.5 0.85 0.6], 'LineWidth', 1.5);
title('Momento Excitador Equivalente (M_{tremor})');
xlabel('Tiempo (s)'); ylabel('N\cdotm'); grid on;

% Gráfica espectral de comprobación
figure('Name', 'Fase 2B: PSD del Momento Excitador', 'Color', 'w', 'Position', [150, 150, 800, 400]);
plot(f_b_M, P_M, 'Color',[0.5 0.85 0.6], 'LineWidth', 2);
title('Densidad Espectral de Potencia del Momento Excitador Equivalente');
xlabel('Frecuencia (Hz)'); ylabel('Potencia ((N\cdotm)^2/Hz)');
grid on;

% Comprobación por consola
f_dom_M = f_b_M(P_M == max(P_M));
disp(' ');
disp(['Frecuencia dominante original (a_rel): ', num2str(metF_Rel(1)), ' Hz']);
disp(['Frecuencia dominante M_tremor: ', num2str(f_dom_M), ' Hz']);