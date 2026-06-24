% =========================================================================
% TFG - Análisis de Temblor Esencial Postural
% FASE 1: caracterización dinámica
% =========================================================================
clear; clc; close all;

%% 1 DEFINICIÓN DE RUTAS Y VARIABLES

ruta_principal = '';

localizaciones = {'Finger', 'Wrist'};
severidades = {'0', '1', '2', '3', '4'};

% Tabla vacía donde iremos apilando los resultados de todos los archivos
Dataset_Final = table(); 

%% 2 BUCLE PARA RECORRER CARPETAS
for loc = 1:length(localizaciones)
    loc_actual = localizaciones{loc};
    
    for sev = 1:length(severidades)
        sev_actual = severidades{sev};
        
        % Ruta a la carpeta específica (ej. .../Senyales/Finger/3)
        ruta_carpeta = fullfile(ruta_principal, loc_actual, sev_actual);
        
        % Buscamos todos los archivos .csv en esa carpeta
        archivos_csv = dir(fullfile(ruta_carpeta, '*.csv'));
        
        % Si la carpeta está vacía, saltamos a la siguiente
        if isempty(archivos_csv)
            continue;
        end
        
        % Procesamos cada archivo encontrado en esa carpeta
        for i = 1:length(archivos_csv)
            nombre_archivo = archivos_csv(i).name;
            ruta_completa = fullfile(ruta_carpeta, nombre_archivo);
            
            % Importar Datos
            opts = detectImportOptions(ruta_completa, 'Delimiter', ';');
            datos = readtable(ruta_completa, opts);
            acc_X = datos.X; acc_Y = datos.Y; acc_Z = datos.Z;
            t = seconds(duration(datos.Time));
            
            % Calcular Fs_real
            Fs_real = 1 / mean(diff(t));
            
            % Preprocesado (Recorte 3s)
            idx_validos = t >= 3;
            acc_X = acc_X(idx_validos); acc_Y = acc_Y(idx_validos); acc_Z = acc_Z(idx_validos);
            
            % Filtro Butterworth (3-15 Hz)
            [b, a] = butter(2, [3, 15] / (Fs_real / 2), 'bandpass');
            acc_X = filtfilt(b, a, acc_X); acc_Y = filtfilt(b, a, acc_Y); acc_Z = filtfilt(b, a, acc_Z);
            
            % Variables temporales (con norma euclídea)
            acc_Norma = sqrt(acc_X.^2 + acc_Y.^2 + acc_Z.^2);
            
            RMS_val = rms(acc_Norma);
            Pico_Pico_val = max(acc_Norma) - min(acc_Norma);
            MAV_val = mean(abs(acc_Norma));
            Std_val = std(acc_Norma);
            
            % Variables frecuenciales (PSD por eje y suma)
            muestras_base = 256;
            Nw = round((muestras_base / 100) * Fs_real); 
            win = hamming(Nw); % Función ventana
            noverlap = round(0.5 * Nw); % Solapamiento del 50%
            nfft = 2^(nextpow2(Nw) + 1); % FFT optimizada y suavizada
            
            % 1. Calculamos la PSD de cada eje por separado
            [Pxx_X, f_welch] = pwelch(acc_X, win, noverlap, nfft, Fs_real);
            [Pxx_Y, ~]       = pwelch(acc_Y, win, noverlap, nfft, Fs_real);
            [Pxx_Z, ~]       = pwelch(acc_Z, win, noverlap, nfft, Fs_real);
            
            % 2. Sumamos las PSDs linealmente para obtener la potencia total sin distorsión
            Pxx_Total = Pxx_X + Pxx_Y + Pxx_Z;
            
            % Filtramos para quedarnos con la banda de interés (3 a 15 Hz)
            idx_banda = f_welch >= 3 & f_welch <= 15;
            f_banda = f_welch(idx_banda);
            Pxx_banda = Pxx_Total(idx_banda);
            
            % 3. Extraemos las métricas de la PSD Total
            [Potencia_Max, idx_max] = max(Pxx_banda);
            Frec_Dom = f_banda(idx_max);
            Area_PSD = bandpower(Pxx_Total, f_welch, 'psd', [3 15]);
            
            % Creamos una fila con las etiquetas y las variables calculadas
            nueva_fila = table({loc_actual}, str2double(sev_actual), {nombre_archivo}, ...
                RMS_val, Pico_Pico_val, MAV_val, Std_val, Frec_Dom, Potencia_Max, Area_PSD, ...
                'VariableNames', {'Localizacion', 'Severidad', 'Archivo', ...
                'RMS', 'Pico_Pico', 'MAV', 'Std', 'Frec_Dominante', 'Potencia_Max', 'Area_PSD'});
            
            Dataset_Final = [Dataset_Final; nueva_fila];
        end
    end
end

%% RESUMEN Y EXPORTACIÓN
disp('======================================================');
disp('Mostrando las primeras 5 filas del Dataset:');
disp(head(Dataset_Final, 5));

% Guardar el dataset en un archivo CSV:
writetable(Dataset_Final, 'Dataset_Temblor.csv');
disp('Dataset guardado como "Dataset_Temblor.csv"');

%% 4 RESULTADOS
% Extraemos los datos separando Muñeca (Wrist) y Dedo (Finger)
severidades_plot = 0:4;
rms_wrist = zeros(1, 5); rms_finger = zeros(1, 5);
pot_wrist = zeros(1, 5); pot_finger = zeros(1, 5);

for i = 1:length(severidades_plot)
    sev = severidades_plot(i);
    
    % Índices lógicos para filtrar la tabla
    idx_w = strcmp(Dataset_Final.Localizacion, 'Wrist') & Dataset_Final.Severidad == sev;
    idx_f = strcmp(Dataset_Final.Localizacion, 'Finger') & Dataset_Final.Severidad == sev;
    
    % Calculamos medias omitiendo posibles NaNs (valores nulos si faltara algún dato)
    rms_wrist(i) = mean(Dataset_Final.RMS(idx_w), 'omitnan');
    rms_finger(i) = mean(Dataset_Final.RMS(idx_f), 'omitnan');
    
    pot_wrist(i) = mean(Dataset_Final.Potencia_Max(idx_w), 'omitnan');
    pot_finger(i) = mean(Dataset_Final.Potencia_Max(idx_f), 'omitnan');
end

% Gráficas de Barras Agrupadas
figure('Name', 'Comparativa Clínica: Muñeca vs Dedo', 'Color', 'w', 'Position', [100, 100, 1000, 450]);

% Colores para la gráfica
color_wrist = [1.00 0.65 0.80]; % Rosa (Muñeca)
color_finger = [1.00 0.93 0.65]; % Amarillo (Dedo)

% Gráfica 1: RMS (Escala Logarítmica)
subplot(1,2,1);
bar_data_rms = [rms_wrist', rms_finger']; 
b1 = bar(severidades_plot, bar_data_rms);
b1(1).FaceColor = color_wrist; 
b1(2).FaceColor = color_finger; 
title('Evolución de la Amplitud (RMS) vs Severidad');
xlabel('Grado de Severidad Clínica (MDS-UPDRS)');
ylabel('RMS de la Aceleración (g) - Escala Logarítmica');
legend('Muñeca (Wrist)', 'Dedo (Finger)', 'Location', 'northwest');
set(gca, 'YScale', 'log'); % <-- Transforma el eje Y a logarítmico
grid on;

% Gráfica 2: Potencia Máxima (Escala Logarítmica)
subplot(1,2,2);
bar_data_pot = [pot_wrist', pot_finger'];
b2 = bar(severidades_plot, bar_data_pot);
b2(1).FaceColor = color_wrist;
b2(2).FaceColor = color_finger;
title('Evolución de la Potencia Máxima vs Severidad');
xlabel('Grado de Severidad Clínica (MDS-UPDRS)');
ylabel('Potencia Espectral Máxima (g^2/Hz) - Escala Logarítmica');
legend('Muñeca (Wrist)', 'Dedo (Finger)', 'Location', 'northwest');
set(gca, 'YScale', 'log'); % <-- Transforma el eje Y a logarítmico
grid on;

sgtitle('Análisis del Temblor Postural: Comparativa de Localización');

%% 5 GENERACIÓN DE TABLA RESUMEN
% Agrupamos por Severidad Y Localizacion para sacar las medias de todo
tabla_resumen_completa = groupsummary(Dataset_Final, {'Severidad', 'Localizacion'}, 'mean', ...
    {'RMS', 'Potencia_Max'});

% Tabla con todas las combinaciones posibles en Command Window:
disp('Tabla Resumen:');
disp(tabla_resumen_completa);

% Guardar el dataset en un archivo CSV:
writetable(tabla_resumen_completa, 'Tabla_Resultados_Resumen.xlsx');