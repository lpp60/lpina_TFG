% =========================================================================
% FASE 2B: ANÁLISIS AUTOMATIZADO Y PROMEDIO POR GRADOS (0 a 4)
% =========================================================================
clear; clc; close all;

% 1. PARÁMETROS DEL MODELO BIOMECÁNICO
L  = 0.15;   % Longitud efectiva muñeca-dedo (m)
I  = 0.0025; % Momento de inercia equivalente (kg*m^2)
c0 = 0.03;   % Amortiguamiento rotacional (N*m*s/rad)
k0 = 0.65;   % Rigidez rotacional equivalente (N*m/rad)

% 2. DEFINIR RUTAS DE TODOS LOS ARCHIVOS (9 pacientes en total, seguir el mismo orden para Muñeca y Dedo)
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

% Vector que indica a qué grado pertenece cada uno de los 9 archivos
grados_pacientes = [0; 1; 1; 2; 2; 3; 3; 4; 4];

% Inicializar vectores para guardar resultados de los 9 pacientes
num_pacientes = length(grados_pacientes);
res_RMS = zeros(num_pacientes, 1);
res_Max = zeros(num_pacientes, 1);
res_Frec = zeros(num_pacientes, 1);
res_Pot = zeros(num_pacientes, 1);

Fs = 100;
g_val = 9.81;


% 3. BUCLE PARA PROCESAR CADA PACIENTE
for i = 1:num_pacientes
    % IMPORTACIÓN
    opts_W = detectImportOptions(rutas_W{i}, 'Delimiter', ';');
    datos_W = readtable(rutas_W{i}, opts_W);
    opts_F = detectImportOptions(rutas_F{i}, 'Delimiter', ';');
    datos_F = readtable(rutas_F{i}, opts_F);
    
    t = seconds(duration(datos_W.Time));
    
    % REMUESTREO Y RECORTE
    temp_Z = resample(datos_F.Z, 100, 2500);
    acc_F_Z = temp_Z(1:length(t));
    acc_W_Z = datos_W.Z;
    
    idx = t >= 3; 
    t_valido = t(idx) - 3;
    acc_W_Z_rec = acc_W_Z(idx);
    acc_F_Z_rec = acc_F_Z(idx);
    
    % FILTRADO
    [b, a] = butter(2, [3, 15] / (Fs / 2), 'bandpass');
    acc_W_Z_filt = filtfilt(b, a, acc_W_Z_rec);
    acc_F_Z_filt = filtfilt(b, a, acc_F_Z_rec);
    
    acc_rel_Z = acc_F_Z_filt - acc_W_Z_filt;
    
    % DINÁMICA EQUIVALENTE
    theta_ddot = (acc_rel_Z * g_val) / L;
    
    theta_dot_raw = cumtrapz(t_valido, theta_ddot);
    theta_dot = filtfilt(b, a, theta_dot_raw);
    
    theta_raw = cumtrapz(t_valido, theta_dot);
    theta = filtfilt(b, a, theta_raw);
    
    M_tremor = I * theta_ddot + c0 * theta_dot + k0 * theta;
    
    % EXTRACCIÓN DE PARÁMETROS
    res_RMS(i) = rms(M_tremor);
    res_Max(i) = max(M_tremor) - min(M_tremor); 
    
    Tw = 2.56; Nw = round(Tw * Fs); win = hamming(Nw); noverlap = round(0.5 * Nw); nfft = 2^(nextpow2(Nw) + 1);
    [P_M, f_M] = pwelch(M_tremor, win, noverlap, nfft, Fs);
    
    idx_b = f_M >= 3 & f_M <= 15;
    f_b = f_M(idx_b); P_M_b = P_M(idx_b);
    
    frec_dom = f_b(P_M_b == max(P_M_b));
    res_Frec(i) = frec_dom(1);
    res_Pot(i) = max(P_M_b);
end

%% 4. AGRUPACIÓN Y CÁLCULO DE MEDIAS POR GRADO
% Creamos una tabla con todos los resultados individuales
Tabla_Individual = table(grados_pacientes, res_RMS, res_Max, res_Frec, res_Pot, ...
    'VariableNames', {'Grado', 'RMS_Nm', 'PicoAPico_Nm', 'Frecuencia_Hz', 'Potencia_M'});

% Agrupamos por grado y calculamos la media (ignorando los posibles NaN)
Tabla_Resumen = groupsummary(Tabla_Individual, 'Grado', 'mean', ...
    {'RMS_Nm', 'PicoAPico_Nm', 'Potencia_M'});

% Para la frecuencia dominante media
Tabla_Frec = groupsummary(Tabla_Individual, 'Grado', 'mean', 'Frecuencia_Hz');
Tabla_Resumen.mean_Frecuencia_Hz = Tabla_Frec.mean_Frecuencia_Hz;
Tabla_Resumen.Properties.VariableNames = {'Grado', 'Num_Pacientes', 'Momento_RMS_Medio', ...
                                          'Momento_PicoAPico_Medio', 'Potencia_Max_Media', 'Frecuencia_Media_Hz'};

disp(' ');
disp('   RESULTADOS MOMENTO EXCITADOR (M_tremor) MEDIO POR GRADO     ');
disp(Tabla_Resumen);

%% 5. GRÁFICA COMPARATIVA FINAL
figure('Name', 'Evolución del Momento Excitador', 'Color', 'w', 'Position', [200 200 800 500]);
bar(Tabla_Resumen.Grado, Tabla_Resumen.Momento_RMS_Medio, 'FaceColor', [0.5 0.85 0.6]);
title('Evolución de la Fuerza del Temblor (Momento Excitador RMS)');
xlabel('Grado de Severidad Clínica (MDS-UPDRS)');
ylabel('Momento Excitador Medio (N\cdotm)');
grid on;