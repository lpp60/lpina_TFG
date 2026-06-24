% =========================================================================
% TFG - Análisis de Temblor Esencial Postural
% FASE 1: análisis dinámico cruzado (pareja dedo - muñeca)
% =========================================================================
clear; clc; close all;

%% 1 DEFINIR RUTAS DE LOS ARCHIVOS .csv 
% Archivo Muñeca
archivo_wrist = '';
% Archivo Dedo
archivo_finger = '';

%% 2 IMPORTAR LOS DATOS Y REMUESTREAR (SINCRONIZACIÓN DE Fs)
% Importamos la Muñeca
opts_W = detectImportOptions(archivo_wrist, 'Delimiter', ';');
datos_W = readtable(archivo_wrist, opts_W);

% Importamos el Dedo
opts_F = detectImportOptions(archivo_finger, 'Delimiter', ';');
datos_F = readtable(archivo_finger, opts_F);

% Extraemos los vectores de tiempo absolutos en segundos
t_W = seconds(duration(datos_W.Time));
t_F = seconds(duration(datos_F.Time));

% Calculamos la Frecuencia de Muestreo (Fs) de cada uno
Fs_W = 1 / mean(diff(t_W));
Fs_F = 1 / mean(diff(t_F));

fprintf('Longitud Muñeca: %d muestras | Fs = %.2f Hz | Duración: %.2f s\n', length(t_W), Fs_W, t_W(end));
fprintf('Longitud Dedo: %d muestras | Fs = %.2f Hz | Duración: %.2f s\n', length(t_F), Fs_F, t_F(end));

% Remuestreo con filtro anti-aliasing integrado (de 2500 Hz a 100 Hz)
temp_X = resample(datos_F.X, 100, 2500);
temp_Y = resample(datos_F.Y, 100, 2500);
temp_Z = resample(datos_F.Z, 100, 2500);

% Ajustamos la longitud para que coincida exactamente con la señal de la muñeca
N_W = length(t_W);
acc_F_X = temp_X(1:N_W);
acc_F_Y = temp_Y(1:N_W);
acc_F_Z = temp_Z(1:N_W);


% Ahora ambas señales tienen exactamente la misma longitud y la misma Fs
acc_W = [datos_W.X, datos_W.Y, datos_W.Z];
acc_F = [acc_F_X, acc_F_Y, acc_F_Z];

% Establecemos el tiempo y la Fs unificadas para el resto del programa
t = t_W; 
Fs = Fs_W;

fprintf('\nSeñales sincronizadas por interpolación a Fs = %.2f Hz\n\n', Fs);

%% 3 PREPROCESADO: RECORTAR TRAMO INICIAL Y FILTRAR
t_recorte = 3; 
idx = t >= t_recorte;
t_valido = t(idx) - t_recorte;
acc_W_val = acc_W(idx, :);
acc_F_val = acc_F(idx, :);

% Filtro Butterworth (3-15 Hz)
[b, a] = butter(2, [3, 15] / (Fs / 2), 'bandpass');

acc_W_filt = zeros(size(acc_W_val));
acc_F_filt = zeros(size(acc_F_val));

for eje = 1:3
    acc_W_filt(:,eje) = filtfilt(b, a, acc_W_val(:,eje));
    acc_F_filt(:,eje) = filtfilt(b, a, acc_F_val(:,eje));
end

%% 4 ANÁLISIS ESPECTRAL, MÉTRICAS CRUZADAS Y VISUALIZACIÓN
muestras_base = 256;
Nw = round((muestras_base / 100) * Fs); 
win = hamming(Nw); % Función ventana
noverlap = round(0.5 * Nw); % Solapamiento del 50%
nfft = 2^(nextpow2(Nw) + 1); % FFT optimizada y suavizada

Nombres_Ejes = {'X', 'Y', 'Z'};
Resultados_Cruzados = table();
figure('Name', 'Relación Dinámica Muñeca-Dedo (3 Ejes)', 'Color', 'w', 'Position', [50, 50, 1400, 800]);

for eje = 1:3
    % Cálculo matemático (Welch, Espectro Cruzado y Coherencia)
    [Pxx_W, f] = pwelch(acc_W_filt(:,eje), win, noverlap, nfft, Fs);
    [Pxx_F, ~] = pwelch(acc_F_filt(:,eje), win, noverlap, nfft, Fs);
    [Pxy, ~]   = cpsd(acc_W_filt(:,eje), acc_F_filt(:,eje), win, noverlap, nfft, Fs);
    [Cxy, ~]   = mscohere(acc_W_filt(:,eje), acc_F_filt(:,eje), win, noverlap, nfft, Fs);
    
    % Filtramos la banda 3-15 Hz para la extracción de métricas y gráficas
    idx_b = f >= 3 & f <= 15;
    f_b = f(idx_b);
    P_W = Pxx_W(idx_b); 
    P_F = Pxx_F(idx_b);
    Pxy_b = Pxy(idx_b);
    Cxy_b = Cxy(idx_b);
    
    [Pot_Max_W, i_W] = max(P_W); Frec_Dom_W = f_b(i_W);
    [Pot_Max_F, i_F] = max(P_F); Frec_Dom_F = f_b(i_F);
    
    idx_ref = i_W; % Frecuencia dominante de la muñeca como referencia
    Ratio_Amplitud = sqrt(P_F(idx_ref) / P_W(idx_ref));
    Fase_Grados = angle(Pxy_b(idx_ref)) * (180 / pi);
    Coherencia_Dom = Cxy_b(idx_ref);
    
    nueva_fila = table(Nombres_Ejes(eje), Frec_Dom_W, Pot_Max_W, Frec_Dom_F, Pot_Max_F, ...
                       Ratio_Amplitud, Coherencia_Dom, Fase_Grados, ...
        'VariableNames', {'Eje', 'Frec_W_Hz', 'Pot_Max_W', 'Frec_F_Hz', 'Pot_Max_F', ...
                          'Ratio_Amp_F_W', 'Coherencia', 'Fase_Grados'});
    Resultados_Cruzados = [Resultados_Cruzados; nueva_fila];

    % Visualización
    % Fila 1: PSDs Superpuestas
    subplot(3, 3, eje);
    plot(f_b, P_W, 'Color', [1.0 0.6 0.8], 'LineWidth', 1.5); hold on;
    plot(f_b, P_F, 'Color', [0.95 0.85 0.35], 'LineWidth', 1.5);
    title(['PSD - Eje ', Nombres_Ejes{eje}]);
    if eje == 1, ylabel('PSD (g^2/Hz)'); legend('Muñeca', 'Dedo', 'Location', 'best'); end
    grid on;
    
    % Fila 2: Coherencia
    subplot(3, 3, eje + 3);
    plot(f_b, Cxy_b, 'Color', [0.5 0.85 0.6], 'LineWidth', 1.2);
    title(['Coherencia - Eje ', Nombres_Ejes{eje}]);
    if eje == 1, ylabel('Coherencia (0 a 1)'); end
    ylim([0 1.1]); grid on;
    
    % Fila 3: Fase
    subplot(3, 3, eje + 6);
    fase_plot = angle(Pxy_b) * (180/pi);
    plot(f_b, fase_plot, 'Color',[0.4 0.75 0.95], 'LineWidth', 1.2);
    title(['Fase - Eje ', Nombres_Ejes{eje}]);
    xlabel('Frecuencia (Hz)'); 
    if eje == 1, ylabel('Fase (Grados)'); end
    ylim([-180 180]); yticks(-180:90:180); grid on;
end

%% 6 GENERACIÓN DE LA TABLA RESUMEN POR REGISTRO

% Separar la ruta del nombre del archivo
[ruta_carpeta, nombre_archivo, ~] = fileparts(archivo_wrist);

% Extraer la Severidad Clínica (es el nombre de la última carpeta)
carpetas = split(ruta_carpeta, filesep);
Severidad_Clinica = str2double(carpetas{end});

% Extraer el Sujeto y el Lado desde el nombre del archivo
% El nombre es tipo: PosturalTremor_Measure317_Subject17_left
partes_nombre = split(nombre_archivo, '_');
Sujeto_ID = partes_nombre{3}; % Se queda con 'Subject17'
Lado = partes_nombre{4};      % Se queda con 'left'

% Preasignamos variables para crear columnas del mismo tamaño (3 ejes)
Num_Ejes = height(Resultados_Cruzados);
Array_Sujeto = repmat({Sujeto_ID}, Num_Ejes, 1);
Array_Lado = repmat({Lado}, Num_Ejes, 1);
Array_Severidad = repmat(Severidad_Clinica, Num_Ejes, 1);
Interpretacion = cell(Num_Ejes, 1);

% Evaluamos la lógica de interpretación para cada eje (X, Y, Z) en un bucle
for i = 1:Num_Ejes
    Ratio = Resultados_Cruzados.Ratio_Amp_F_W(i);
    Coher = Resultados_Cruzados.Coherencia(i);
    Fase  = Resultados_Cruzados.Fase_Grados(i);
    Frec_W = Resultados_Cruzados.Frec_W_Hz(i);
    Frec_F = Resultados_Cruzados.Frec_F_Hz(i);
    
    % Escenario 1: Movimiento común / Sólido rígido
    % (Frecuencia igual, Ratio pequeño, sin apenas desfase, alta coherencia)
    if (Frec_W == Frec_F) && (Ratio <= 1.5) && (abs(Fase) < 20) && (Coher >= 0.7)
        Interpretacion{i} = 'Movimiento comun (solido rigido)';

    % Escenario 2: Amplificación distal con movimiento relativo
    % (Frecuencia igual, y o hay mucho ratio o hay mucho desfase)
    elseif (Frec_W == Frec_F) && (Ratio > 1.2 || abs(Fase) >= 20) && (Coher > 0.5)
        Interpretacion{i} = 'Amplificacion distal con mto. relativo';

    % Escenario 3: Dinámica local no coherente o posible armónico
    elseif (Frec_W ~= Frec_F) || (Coher <= 0.5)
        Interpretacion{i} = 'Posible componente distal local o armonico';

    % Escenario 4: No concluyente
    else
        Interpretacion{i} = 'Resultado no concluyente';
    end
end

% Construimos la tabla final
Tabla_Resumen_Registro = table(Array_Sujeto, Array_Lado, Array_Severidad, ...
    Resultados_Cruzados.Eje, ...
    Resultados_Cruzados.Frec_W_Hz, Resultados_Cruzados.Pot_Max_W, ...
    Resultados_Cruzados.Frec_F_Hz, Resultados_Cruzados.Pot_Max_F, ...
    Resultados_Cruzados.Ratio_Amp_F_W, Resultados_Cruzados.Fase_Grados, ...
    Resultados_Cruzados.Coherencia, Interpretacion, ...
    'VariableNames', {'ID_Sujeto', 'Lado', 'Grado', 'Eje', ...
                      'Frec_W_Hz', 'Pot_W', 'Frec_F_Hz', 'Pot_F', ...
                      'Ratio_Amp', 'Fase_Grados', 'Coherencia', 'Interpretacion'});
disp(' ');
disp('      TABLA RESUMEN');
disp(' ');
disp(Tabla_Resumen_Registro);