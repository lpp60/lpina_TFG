% =========================================================================
% Análisis de Temblor Esencial Postural
% FASE 1: análisis individual
% =========================================================================
%clear; clc; close all;

%% 1 DEFINIR LA RUTA DEL ARCHIVO .csv

archivo_csv = '';

%% 2 IMPORTAR LOS DATOS
% Leemos el .csv delimitado por punto y coma
opts = detectImportOptions(archivo_csv, 'Delimiter', ';');
datos = readtable(archivo_csv, opts);

% Extraemos los datos de aceleración en cada eje
acc_X = datos.X;
acc_Y = datos.Y;
acc_Z = datos.Z;

%% 3 TRATAMIENTO DEL VECTOR TIEMPO Y CÁLCULO DE Fs
% Convertimos la columna 'Time' a un formato de duración para poder pasarlo a segundos.
tiempo_duracion = duration(datos.Time); 
t = seconds(tiempo_duracion); % Vector de tiempo relativo en segundos

% Cálculo de la frecuencia de muestreo real (Fs)
dt = diff(t);                 % Diferencia de tiempo entre muestras
dt_medio = mean(dt);          % Periodo de muestreo medio
Fs_real = 1 / dt_medio;       % Frecuencia de muestreo (Hz)

% Mostramos el resultado en la Command Window
fprintf('Archivo analizado: %s\n', archivo_csv);
fprintf('Número total de muestras: %d\n', length(t));
fprintf('Frecuencia de Muestreo (Fs) calculada: %.2f Hz\n\n', Fs_real);

%% 4 VISUALIZACIÓN DE LAS ACELERACIONES (3 EJES)
figure('Name', 'Aceleraciones Triaxiales', 'Color', 'w', 'Position', [100, 100, 800, 600]);

% Colores para las gráficas
color_rosa  = [1.0 0.6 0.8];
color_verde = [0.5 0.85 0.6];
color_azul  = [0.4 0.75 0.95];
color_gris = [0.7 0.7 0.7];

% Gráfica Eje X
subplot(3,1,1);
plot(t, acc_X, 'Color',color_verde, 'LineWidth', 1);
title('Aceleración Eje X');
ylabel('Acel. (g o m/s^2)');
grid on;

% Gráfica Eje Y
subplot(3,1,2);
plot(t, acc_Y, 'Color',color_rosa, 'LineWidth', 1);
title('Aceleración Eje Y');
ylabel('Acel. (g o m/s^2)');
grid on;

% Gráfica Eje Z
subplot(3,1,3);
plot(t, acc_Z, 'Color',color_azul, 'LineWidth', 1);
title('Aceleración Eje Z');
xlabel('Tiempo (s)');
ylabel('Acel. (g o m/s^2)');
grid on;

sgtitle('Señales en bruto de temblor postural');

%% 5 PREPROCESADO: RECORTAR TRAMO INICIAL
t_recorte = 3; % Segundos a eliminar al principio
idx_validos = t >= t_recorte; % Índices lógicos de los datos válidos

% Aplicamos el recorte a los vectores
t_recortado = t(idx_validos) - t_recorte; % Restamos t_recorte para que empiece en 0s
acc_X_recortado = acc_X(idx_validos);
acc_Y_recortado = acc_Y(idx_validos);
acc_Z_recortado = acc_Z(idx_validos);

%% 6 PREPROCESADO: FILTRADO DIGITAL
% Frecuencias de corte del filtro paso banda (3 Hz a 15 Hz)
f_corte_inf = 3;
f_corte_sup = 15;

% Diseño del filtro Butterworth de 4º orden (2º orden bidireccional se vuelve 4º)
orden_filtro = 2; 
% Frecuencia de Nyquist es Fs/2. Normalizamos las frecuencias de corte:
Wn = [f_corte_inf, f_corte_sup] / (Fs_real / 2); 

% Obtenemos los coeficientes del filtro
[b, a] = butter(orden_filtro, Wn, 'bandpass');

% Aplicamos filtfilt (filtrado de fase cero) a cada eje
acc_X_filt = filtfilt(b, a, acc_X_recortado);
acc_Y_filt = filtfilt(b, a, acc_Y_recortado);
acc_Z_filt = filtfilt(b, a, acc_Z_recortado);

%% 7 COMPARACIÓN VISUAL (ANTES Y DESPUÉS DEL FILTRADO)
figure('Name', 'Comparativa Filtrado', 'Color', 'w', 'Position', [150, 150, 900, 700]);

% Comparación Eje X
subplot(3,1,1);
plot(t_recortado, acc_X_recortado, 'Color', color_gris, 'LineWidth', 1); hold on; % Señal original en gris
plot(t_recortado, acc_X_filt, 'Color',color_verde, 'LineWidth', 1);
title('Eje X: Señal Original vs Filtrada (3-15 Hz)');
ylabel('Acel. (g)');
legend('Original (recortada)', 'Filtrada', 'Location', 'eastoutside');
grid on;

% Comparación Eje Y
subplot(3,1,2);
plot(t_recortado, acc_Y_recortado, 'Color', color_gris, 'LineWidth', 1); hold on;
plot(t_recortado, acc_Y_filt, 'Color',color_rosa, 'LineWidth', 1);
title('Eje Y: Señal Original vs Filtrada (3-15 Hz)');
ylabel('Acel. (g)');
legend('Original (recortada)', 'Filtrada', 'Location', 'eastoutside');
grid on;

% Comparación Eje Z
subplot(3,1,3);
plot(t_recortado, acc_Z_recortado, 'Color', color_gris,'LineWidth', 1); hold on;
plot(t_recortado, acc_Z_filt, 'Color',color_azul, 'LineWidth', 1);
title('Eje Z: Señal Original vs Filtrada (3-15 Hz)');
xlabel('Tiempo (s)');
ylabel('Acel. (g)');
legend('Original (recortada)', 'Filtrada', 'Location', 'eastoutside');
grid on;

sgtitle('Efecto del filtrado sobre la señal de aceleración (centrado y limpieza)');

%% 8 EXTRACCIÓN DE CARACTERÍSTICAS: DOMINIO DEL TIEMPO
% Calculamos todas las métricas para cada eje usando las funciones matemáticas de Matlab
% Creamos una función para calcular los 6 parámetros para un vector dado (x)
calc_metricas_tiempo = @(x) [ ...
    rms(x), ...               % 1 Valor Cuadrático Medio (RMS)
    max(x), ...               % 2 Máximo
    min(x), ...               % 3 Mínimo
    max(x) - min(x), ...      % 4 Amplitud Pico a Pico
    mean(abs(x)), ...         % 5 Valor Medio Absoluto (MAV)
    std(x) ...                % 6 Desviación Estándar
];

% Aplicamos la función a las señales filtradas de los 3 ejes
metricas_X = calc_metricas_tiempo(acc_X_filt);
metricas_Y = calc_metricas_tiempo(acc_Y_filt);
metricas_Z = calc_metricas_tiempo(acc_Z_filt);

% Creamos una tabla resumen para mostrar los resultados en consola
Nombres_Variables = {'RMS', 'Maximo', 'Minimo', 'Pico_a_Pico', 'MAV', 'Desviacion_Std'};
Nombres_Filas = {'Eje_X', 'Eje_Y', 'Eje_Z'};

Tabla_Variables_Tiempo = array2table([metricas_X; metricas_Y; metricas_Z], ...
    'VariableNames', Nombres_Variables, 'RowNames', Nombres_Filas);

disp('--- RESULTADOS: VARIABLES EN EL DOMINIO DEL TIEMPO ---');
disp(Tabla_Variables_Tiempo);

%% 9 EXTRACCIÓN DE CARACTERÍSTICAS: DOMINIO DE LA FRECUENCIA (Welch)
% Configuración del método de Welch
muestras_base = 256;

Nw = round((muestras_base / 100) * Fs_real); 
win = hamming(Nw); % Función ventana
noverlap = round(0.5 * Nw); % Solapamiento del 50%
nfft = 2^(nextpow2(Nw) + 1); % FFT optimizada y suavizada

% Cálculo de la Densidad Espectral de Potencia (PSD) para los 3 ejes
[Pxx_X, f_welch] = pwelch(acc_X_filt, win, noverlap, nfft, Fs_real);
[Pxx_Y, ~]       = pwelch(acc_Y_filt, win, noverlap, nfft, Fs_real);
[Pxx_Z, ~]       = pwelch(acc_Z_filt, win, noverlap, nfft, Fs_real);

% Filtramos los resultados para quedarnos solo con la banda de interés (3 a 15 Hz)
idx_banda = f_welch >= 3 & f_welch <= 15;
f_banda = f_welch(idx_banda);

Pxx_X_banda = Pxx_X(idx_banda);
Pxx_Y_banda = Pxx_Y(idx_banda);
Pxx_Z_banda = Pxx_Z(idx_banda);

% --- CÁLCULO DE MÉTRICAS FRECUENCIALES POR EJE ---
% Eje X
[Potencia_Max_X, idx_max_X] = max(Pxx_X_banda);
Frecuencia_Dom_X = f_banda(idx_max_X);
Area_PSD_X = bandpower(Pxx_X, f_welch, 'psd', [3 15]); % Integral de la potencia en la banda

% Eje Y
[Potencia_Max_Y, idx_max_Y] = max(Pxx_Y_banda);
Frecuencia_Dom_Y = f_banda(idx_max_Y);
Area_PSD_Y = bandpower(Pxx_Y, f_welch, 'psd', [3 15]);

% Eje Z
[Potencia_Max_Z, idx_max_Z] = max(Pxx_Z_banda);
Frecuencia_Dom_Z = f_banda(idx_max_Z);
Area_PSD_Z = bandpower(Pxx_Z, f_welch, 'psd', [3 15]);

% Mostrar resultados en la Command Window
Nombres_Var_Frec = {'Frecuencia_Dominante_Hz', 'Potencia_Maxima', 'Area_PSD'};
Tabla_Variables_Frec = array2table([Frecuencia_Dom_X, Potencia_Max_X, Area_PSD_X; ...
                                    Frecuencia_Dom_Y, Potencia_Max_Y, Area_PSD_Y; ...
                                    Frecuencia_Dom_Z, Potencia_Max_Z, Area_PSD_Z], ...
    'VariableNames', Nombres_Var_Frec, 'RowNames', Nombres_Filas);

disp('--- RESULTADOS: VARIABLES EN EL DOMINIO DE LA FRECUENCIA ---');
disp(Tabla_Variables_Frec);

%% 10 VISUALIZACIÓN DE LA DENSIDAD ESPECTRAL DE POTENCIA (PSD)
figure('Name', 'Espectro de Potencia (Welch)', 'Color', 'w', 'Position', [200, 200, 800, 500]);
plot(f_banda, Pxx_X_banda, 'Color',color_verde, 'LineWidth', 1.5); hold on;
plot(f_banda, Pxx_Y_banda, 'Color',color_rosa, 'LineWidth', 1.5);
plot(f_banda, Pxx_Z_banda, 'Color',color_azul, 'LineWidth', 1.5);
title('Densidad Espectral de Potencia (Método de Welch) - Banda 3-15 Hz');
xlabel('Frecuencia (Hz)');
ylabel('Potencia / Frecuencia (g^2/Hz)');
legend('Eje X', 'Eje Y', 'Eje Z');
grid on;
xlim([3 15]);

%% 11. VALIDACIÓN DEL FILTRADO (Comparativa 0.5-40 Hz vs 3-15 Hz en 3 ejes)
% Se evalúa el efecto del filtro en los tres ejes ortogonales

% 1. Filtrado amplio (0.5 - 40 Hz) de las señales recortadas (acc_X, acc_Y, acc_Z)
fc_low_wide = 0.5; 
fc_high_wide = 40;
[b_wide, a_wide] = butter(4, [fc_low_wide, fc_high_wide]/(Fs_real/2), 'bandpass');

acc_X_wide = filtfilt(b_wide, a_wide, acc_X);
acc_Y_wide = filtfilt(b_wide, a_wide, acc_Y);
acc_Z_wide = filtfilt(b_wide, a_wide, acc_Z);

% 2. PSD del filtro amplio (Aprovechamos win, noverlap y nfft calculados en la sección de Welch)
[Pxx_X_wide, f_wide] = pwelch(acc_X_wide, win, noverlap, nfft, Fs_real);
[Pxx_Y_wide, ~]      = pwelch(acc_Y_wide, win, noverlap, nfft, Fs_real);
[Pxx_Z_wide, ~]      = pwelch(acc_Z_wide, win, noverlap, nfft, Fs_real);

% 3. Extracción de frecuencias dominantes para comparativa numérica
[max_P_X_wide, idx_X_wide] = max(Pxx_X_wide); f_dom_X_wide = f_wide(idx_X_wide);
[max_P_Y_wide, idx_Y_wide] = max(Pxx_Y_wide); f_dom_Y_wide = f_wide(idx_Y_wide);
[max_P_Z_wide, idx_Z_wide] = max(Pxx_Z_wide); f_dom_Z_wide = f_wide(idx_Z_wide);

% Frecuencias dominantes del filtro estricto (Pxx_X, Pxx_Y, Pxx_Z calculadas previamente)
[max_P_X_strict, idx_X_strict] = max(Pxx_X); f_dom_X_strict = f_welch(idx_X_strict);
[max_P_Y_strict, idx_Y_strict] = max(Pxx_Y); f_dom_Y_strict = f_welch(idx_Y_strict);
[max_P_Z_strict, idx_Z_strict] = max(Pxx_Z); f_dom_Z_strict = f_welch(idx_Z_strict);

% Imprimir por consola para comprobar
fprintf('\n--- VALIDACIÓN DEL FILTRO (AMPLIO vs ESTRICTO) ---\n');
fprintf('Eje X -> Amplio (0.5 - 40 Hz): %.4f Hz | Seleccionado (3-15 Hz): %.4f Hz\n', f_dom_X_wide, f_dom_X_strict);
fprintf('Eje Y -> Amplio (0.5 - 40 Hz): %.4f Hz | Seleccionado (3-15 Hz): %.4f Hz\n', f_dom_Y_wide, f_dom_Y_strict);
fprintf('Eje Z -> Amplio (0.5 - 40 Hz): %.4f Hz | Seleccionado (3-15 Hz): %.4f Hz\n', f_dom_Z_wide, f_dom_Z_strict);

% 4. Representación Gráfica Comparativa
figure('Name', 'Validación del Filtrado (3 Ejes)', 'Color', 'w', 'Position', [100, 100, 1200, 400]);

% Gráfica Eje X
subplot(3,1,1);
plot(f_wide, Pxx_X_wide, 'k', 'LineWidth', 1.5, 'DisplayName', 'Amplio (0.5-40 Hz)'); hold on;
plot(f_welch, Pxx_X, 'Color', color_verde, 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Seleccionado (3-15 Hz)');
xlim([0 20]); xlabel('Frecuencia (Hz)'); ylabel('Potencia (g^2/Hz)'); title('Eje X'); grid on; legend('Location','northeast');

% Gráfica Eje Y
subplot(3,1,2);
plot(f_wide, Pxx_Y_wide, 'k', 'LineWidth', 1.5, 'DisplayName', 'Amplio (0.5-40 Hz)'); hold on;
plot(f_welch, Pxx_Y, 'Color', color_rosa, 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Seleccionado (3-15 Hz)');
xlim([0 20]); xlabel('Frecuencia (Hz)'); ylabel('Potencia (g^2/Hz)');  title('Eje Y'); grid on; legend('Location','northeast');

% Gráfica Eje Z
subplot(3,1,3);
plot(f_wide, Pxx_Z_wide, 'k', 'LineWidth', 1.5, 'DisplayName', 'Amplio (0.5-40 Hz)'); hold on;
plot(f_welch, Pxx_Z, 'Color', color_azul, 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Seleccionado (3-15 Hz)');
xlim([0 20]); xlabel('Frecuencia (Hz)'); ylabel('Potencia (g^2/Hz)');  title('Eje Z'); grid on; legend('Location','northeast');
