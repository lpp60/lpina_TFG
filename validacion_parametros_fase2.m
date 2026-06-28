% =========================================================================
% FASE 2: VALIDACIÓN Y CLONADO DE SEÑAL SINTÉTICA
% =========================================================================
clear; clc; close all;

% 1. PARÁMETROS DEL MODELO BIOMECÁNICO
L  = 0.15;   % Longitud efectiva muñeca-dedo (m)
I  = 0.0025; % Momento de inercia equivalente (kg*m^2)
c0 = 0.03;   % Amortiguamiento rotacional (N*m*s/rad)
k0 = 0.65;   % Rigidez rotacional equivalente (N*m/rad)
Fs = 100;
g_val = 9.81;

% RUTAS DEL SUJETO 4 (Grado 4 Severo)
ruta_W = '';
ruta_F = '';

%% PARTE 1: OBTENER SEÑAL EXPERIMENTAL REAL
opts_W = detectImportOptions(ruta_W, 'Delimiter', ';'); datos_W = readtable(ruta_W, opts_W);
opts_F = detectImportOptions(ruta_F, 'Delimiter', ';'); datos_F = readtable(ruta_F, opts_F);

t = seconds(duration(datos_W.Time));
temp_Z = resample(datos_F.Z, 100, 2500);
acc_F_Z = temp_Z(1:length(t));
acc_W_Z = datos_W.Z;

idx = t >= 3; 
t_valido = t(idx) - 3; 
if isrow(t_valido), t_valido = t_valido'; end

[b, a] = butter(2, [3, 15] / (Fs / 2), 'bandpass');
acc_W_Z_filt = filtfilt(b, a, acc_W_Z(idx));
acc_F_Z_filt = filtfilt(b, a, acc_F_Z(idx));

% SEÑAL REAL A CLONAR
acc_rel_exp = acc_F_Z_filt - acc_W_Z_filt; 
if isrow(acc_rel_exp), acc_rel_exp = acc_rel_exp'; end

%% PARTE 2: EXTRAER PARÁMETROS DE LA SEÑAL EXPERIMENTAL
Tw = 2.56; Nw = round(Tw * Fs); win = hamming(Nw); noverlap = round(0.5 * Nw); nfft = 2^(nextpow2(Nw) + 1);
[P_exp, f_exp] = pwelch(acc_rel_exp, win, noverlap, nfft, Fs);

idx_b = f_exp >= 3 & f_exp <= 15;
f_b = f_exp(idx_b); 
P_exp_b = P_exp(idx_b);

f0_exp = f_b(P_exp_b == max(P_exp_b)); 
f0_exp = f0_exp(1); % Frecuencia dominante real
PotMax_exp = max(P_exp_b); % Potencia máxima real

% DETECCIÓN DE ARMÓNICO
f2_target = 2 * f0_exp;
idx_arm = find(f_b >= f2_target - 1.5 & f_b <= f2_target + 1.5); % Margen de +/- 1.5 Hz
[P_arm_max, idx_max_arm] = max(P_exp_b(idx_arm));

% Si el pico armónico tiene al menos el 5% de la potencia fundamental, lo incluimos
if ~isempty(P_arm_max) && (P_arm_max > 0.05 * PotMax_exp)
    f2_exp = f_b(idx_arm(idx_max_arm));
    incluir_armonico = true;
else
    f2_exp = 2 * f0_exp;
    incluir_armonico = false;
end

RMS_exp = rms(acc_rel_exp); 
PP_exp = max(acc_rel_exp) - min(acc_rel_exp); 
AUC_exp = trapz(f_b, P_exp_b); 

%% PARTE 3: MÍNIMOS CUADRADOS (El Clonado)
w0 = 2 * pi * f0_exp;

if incluir_armonico
    w2 = 2 * pi * f2_exp;
    X = [sin(w0*t_valido), cos(w0*t_valido), sin(w2*t_valido), cos(w2*t_valido)];
else
    X = [sin(w0*t_valido), cos(w0*t_valido)];
end

coefs = X \ acc_rel_exp; 

A1 = sqrt(coefs(1)^2 + coefs(2)^2);
if incluir_armonico
    A2 = sqrt(coefs(3)^2 + coefs(4)^2);
else
    A2 = 0;
end

acc_armonica = X * coefs; 

residuo_exp = acc_rel_exp - acc_armonica;
ruido_blanco = randn(size(t_valido));
ruido_filt = filtfilt(b, a, ruido_blanco);
ruido_esc = ruido_filt * (rms(residuo_exp) / rms(ruido_filt)); 

acc_sint_temp = acc_armonica + ruido_esc;
alpha = RMS_exp / rms(acc_sint_temp);
acc_rel_sint = alpha * acc_sint_temp; % Clon definitivo

% Ajustamos amplitudes finales para la tabla según el escalado alpha
A1 = A1 * alpha;
A2 = A2 * alpha;

%% PARTE 4: CALCULAR CINEMÁTICA Y DINÁMICA CON EL CLON
theta_ddot_sint = (acc_rel_sint * g_val) / L;

theta_dot_raw_s = cumtrapz(t_valido, theta_ddot_sint);
theta_dot_sint = filtfilt(b, a, theta_dot_raw_s);

theta_raw_s = cumtrapz(t_valido, theta_dot_sint);
theta_sint = filtfilt(b, a, theta_raw_s);

M_tremor_sint = I * theta_ddot_sint + c0 * theta_dot_sint + k0 * theta_sint;

%% PARTE 5: MÉTRICAS DE VALIDACIÓN
[P_sint, ~] = pwelch(acc_rel_sint, win, noverlap, nfft, Fs);
P_sint_b = P_sint(idx_b);

RMS_sint = rms(acc_rel_sint);
PP_sint = max(acc_rel_sint) - min(acc_rel_sint);
AUC_sint = trapz(f_b, P_sint_b);
f0_sint = f_b(P_sint_b == max(P_sint_b)); f0_sint = f0_sint(1);
PotMax_sint = max(P_sint_b);

% Calcular ancho de banda a -3dB
P_exp_b_norm = P_exp_b / max(P_exp_b);
idx_band_exp = P_exp_b_norm >= 0.5;
ancho_banda_exp = f_b(find(idx_band_exp, 1, 'last')) - f_b(find(idx_band_exp, 1, 'first'));

P_sint_b_norm = P_sint_b / max(P_sint_b);
idx_band_sint = P_sint_b_norm >= 0.5;
ancho_banda_sint = f_b(find(idx_band_sint, 1, 'last')) - f_b(find(idx_band_sint, 1, 'first'));

disp(' ');
disp(' TABLA DE VALIDACIÓN: EXPERIMENTAL VS SINTÉTICA');
fprintf('Frecuencia Dominante (Hz): Real = %.2f  | Sintética = %.2f\n', f0_exp, f0_sint);
fprintf('RMS (g):                   Real = %.4f  | Sintética = %.4f\n', RMS_exp, RMS_sint);
fprintf('Pico a Pico (g):           Real = %.4f  | Sintética = %.4f\n', PP_exp, PP_sint);
fprintf('Potencia Máxima (g^2/Hz):  Real = %.4f  | Sintética = %.4f\n', PotMax_exp, PotMax_sint);
fprintf('Área PSD (Potencia total): Real = %.4f  | Sintética = %.4f\n', AUC_exp, AUC_sint);
fprintf('Ancho de Banda (-3dB) (Hz): Real = %.4f | Sintética = %.4f\n', ancho_banda_exp, ancho_banda_sint);
if incluir_armonico
    fprintf('Ratio Amplitud (A2/A1):    Sintética = %.2f%%\n', (A2/A1)*100);
else
    fprintf('Armónico:                  No detectado con suficiente claridad.\n');
end

%% PARTE 6: GRÁFICAS COMPARATIVAS
figure('Name', 'Validación del Clon Sintético', 'Color', 'w', 'Position', [100 100 1000 600]);

subplot(2,1,1);
plot(t_valido, acc_rel_exp, 'Color', [0.5 0.85 0.6] , 'LineWidth', 1.2); hold on;
plot(t_valido, acc_rel_sint, '--', 'Color', [0 0 0], 'LineWidth', 1.5);
xlim([0 2]);
title('Comparación Temporal: Aceleración Relativa (Zoom de 2 segundos)');
xlabel('Tiempo (s)'); ylabel('Aceleración (g)');
legend('Señal Experimental', 'Señal Sintética', 'Location', 'best');
grid on;

subplot(2,1,2);
plot(f_b, P_exp_b, 'Color', [0.5 0.85 0.6], 'LineWidth', 1.2); hold on;
plot(f_b, P_sint_b, '--', 'Color', [0 0 0], 'LineWidth', 1.5);
title('Comparación Frecuencial: Densidad Espectral de Potencia (PSD)');
xlabel('Frecuencia (Hz)'); ylabel('Potencia (g^2/Hz)');
legend('PSD Experimental', 'PSD Sintética', 'Location', 'best');
grid on;

%% EXTRACCIÓN DE PARÁMETROS REALES POR GRADO PARA APP DESIGNER

% DEFINIR RUTAS DE TODOS LOS ARCHIVOS (9 pacientes en total, seguir el mismo orden para Muñeca y Dedo)
rutas_W = {
    '', ... % 1. Muñeca Grado 0
    '', ... % 2. Muñeca Grado 1 (Paciente A)
    '', ... % 3. Muñeca Grado 1 (Paciente B)
    '', ... % 4. Muñeca Grado 2 (Paciente A)
    '', ... % 5. Muñeca Grado 2 (Paciente B)
    '', ... % 6. Muñeca Grado 3 (Paciente A)
    '', ... % 7. Muñeca Grado 3 (Paciente B)
    '', ... % 8. Muñeca Grado 4 (Paciente A)
    ''  ... % 9. Muñeca Grado 4 (Paciente B)
};

rutas_F = {
    '', ... % 1. Dedo Grado 0
    '', ... % 2. Dedoa Grado 1 (Paciente A)
    '', ... % 3. Dedo Grado 1 (Paciente B)
    '', ... % 4. Dedo Grado 2 (Paciente A)
    '', ... % 5. Dedo Grado 2 (Paciente B)
    '', ... % 6. Dedo Grado 3 (Paciente A)
    '', ... % 7. Dedo Grado 3 (Paciente B)
    '', ... % 8. Dedo Grado 4 (Paciente A)
    ''  ... % 9. Dedo Grado 4 (Paciente B)
};

grados_pacientes = [0; 1; 1; 2; 2; 3; 3; 4; 4];
num_pacientes = length(grados_pacientes);

v_f0 = zeros(num_pacientes, 1); 
v_A1 = zeros(num_pacientes, 1); 
v_A2 = zeros(num_pacientes, 1);
v_ruido = zeros(num_pacientes, 1);

for i = 1:num_pacientes
    opts_W = detectImportOptions(rutas_W{i}, 'Delimiter', ';'); datos_W = readtable(rutas_W{i}, opts_W);
    opts_F = detectImportOptions(rutas_F{i}, 'Delimiter', ';'); datos_F = readtable(rutas_F{i}, opts_F);
    
    t = seconds(duration(datos_W.Time));
    temp_Z = resample(datos_F.Z, 100, 2500); acc_F_Z = temp_Z(1:length(t)); acc_W_Z = datos_W.Z;
    
    idx = t >= 3; t_valido = t(idx) - 3;
    if isrow(t_valido), t_valido = t_valido'; end
    
    acc_W_Z_filt = filtfilt(b, a, acc_W_Z(idx)); acc_F_Z_filt = filtfilt(b, a, acc_F_Z(idx));
    
    acc_rel_exp = acc_F_Z_filt - acc_W_Z_filt;
    if isrow(acc_rel_exp), acc_rel_exp = acc_rel_exp'; end
    
    [P_exp, f_exp] = pwelch(acc_rel_exp, win, noverlap, nfft, Fs);
    idx_b = f_exp >= 3 & f_exp <= 15; f_b = f_exp(idx_b); P_exp_b = P_exp(idx_b);
    
    f0_exp = f_b(P_exp_b == max(P_exp_b)); f0_exp = f0_exp(1);
    v_f0(i) = f0_exp;
    
    % Evaluar si existe armónico
    f2_target = 2 * f0_exp;
    idx_arm = find(f_b >= f2_target - 1.5 & f_b <= f2_target + 1.5);
    [P_arm_max, idx_max_arm] = max(P_exp_b(idx_arm));
    
    if ~isempty(P_arm_max) && (P_arm_max > 0.05 * max(P_exp_b))
        f2_exp = f_b(idx_arm(idx_max_arm));
        incluir_arm = true;
    else
        incluir_arm = false;
    end
    
    w0 = 2 * pi * f0_exp; 
    if incluir_arm
        w2 = 2 * pi * f2_exp;
        X = [sin(w0*t_valido), cos(w0*t_valido), sin(w2*t_valido), cos(w2*t_valido)];
    else
        X = [sin(w0*t_valido), cos(w0*t_valido)];
    end
    
    coefs = X \ acc_rel_exp;
    acc_armonica = X * coefs;
    
    % Ruido residual
    residuo_exp = acc_rel_exp - acc_armonica;
    ruido_blanco = randn(size(t_valido));
    ruido_filt = filtfilt(b, a, ruido_blanco);
    ruido_esc = ruido_filt * (rms(residuo_exp) / rms(ruido_filt));
    
    % Factor alpha con la señal sintética COMPLETA
    acc_sint_temp = acc_armonica + ruido_esc;
    alpha = rms(acc_rel_exp) / rms(acc_sint_temp); 
    
    % Guardamos amplitudes y ruido final ya escalados
    v_A1(i) = sqrt(coefs(1)^2 + coefs(2)^2) * alpha;
    if incluir_arm
        v_A2(i) = sqrt(coefs(3)^2 + coefs(4)^2) * alpha;
    else
        v_A2(i) = 0;
    end
    v_ruido(i) = rms(residuo_exp) * alpha;
end

Tabla_Individual = table(grados_pacientes, v_f0, v_A1, v_A2, v_ruido, ...
    'VariableNames', {'Grado', 'f0_Hz', 'A1_g', 'A2_g', 'RMS_Ruido_g'});

Tabla_App = groupsummary(Tabla_Individual, 'Grado', 'mean', ...
    {'f0_Hz', 'A1_g', 'A2_g', 'RMS_Ruido_g'});

Tabla_App.Properties.VariableNames = {'Grado', 'Num_Pacientes', 'Frecuencia_f0_Hz', ...
                                      'Amplitud_A1_g', 'Amplitud_A2_g', 'RMS_Ruido_g'};
disp(' ');
disp('   PARÁMETROS POR GRADO ');
disp(Tabla_App);