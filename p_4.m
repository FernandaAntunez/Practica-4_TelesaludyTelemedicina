clear; 
close all; 
clc;

% === 0.1 Leer audio (ajusta la ruta a tus archivos) ===
% Ejemplos del repositorio de Michigan (cambia por tus rutas):
 % [x, fs] = audioread('09_apex_holo_sys_mur_supine_bell.mp3');
 % [x, fs] = audioread('02_apex_split_s1_supine_bell.mp3');
  [x, fs] = audioread('07_apex_mid_sys_mur_supine_bell.mp3');

% === 0.2 Si es estéreo, convertir a mono ===
if size(x,2) > 1
    x = mean(x, 2);
end

% === 0.3 Normalizar a [-1, 1] ===
a = -1; 
b = 1;
x = (x - min(x)) / max(eps, (max(x) - min(x)));  % [0,1]
x = x * (b - a) + a;                              % [-1,1]

% === 0.4 Ventana de trabajo (2 s para empezar) ===
t = (0:length(x)-1)/fs;
max_time = 10;                         % segundos
sel = t <= max_time;
x = x(sel);  
t = t(sel);

% === 1.1 Probabilidad normalizada p[n] ===
p = abs(x);
p = p ./ max(eps, max(p));

% === 1.2 Entropía/energía de Shannon puntual (log10) ===
E = -p .* log10(p + eps);     % variante equivalente a -|x|log|x|

% === 1.3 Estandarización y normalización a [0,1] ===
E_z  = (E - mean(E)) / (std(E) + eps);
Env0 = (E_z - min(E_z)) / max(eps, (max(E_z) - min(E_z)));   % envolvente base [0,1]

% === 2.1 LPF sobre la envolvente (no sobre x) ===
fc = 10;                                         % corte ~10 Hz
[b,a] = butter(4, fc/(fs/2), 'low');
Env = filtfilt(b, a, Env0);

% === 2.2 Normalización final a [0,1] ===
Env = (Env - min(Env)) / max(eps, (max(Env) - min(Env)));

% Visual rápido
figure('Name','Señal & Envolvente (Shannon LPF)');
subplot(2,1,1); 
plot(t, x, 'k'); 
grid on;
title('Señal PCG (normalizada)'); 
xlabel('Tiempo (s)'); ylabel('Amplitud');
subplot(2,1,2); 
plot(t, Env, 'b'); 
grid on;
title('Envolvente de Shannon (LPF)'); 
xlabel('Tiempo (s)'); 
ylabel('Amplitud norm.');



% === 3.1 Derivada discreta y detección de cambios de signo ===
d = diff(Env);
idx_ext = []; 
tipo = [];  % tipo: 1 = min, 2 = max

for i = 1:length(d)-1
    if d(i) < 0 && d(i+1) > 0
        idx_ext(end+1) = i+1; tipo(end+1) = 1;  % mínimo
    elseif d(i) > 0 && d(i+1) < 0
        idx_ext(end+1) = i+1; tipo(end+1) = 2;  % máximo
    end
end

% Visual
figure('Name','Envolvente con min/máx');
plot(t, Env, 'b'); 
hold on; 
grid on;
plot(t(idx_ext(tipo==1)), Env(idx_ext(tipo==1)), 'go', 'DisplayName','Mínimos');
plot(t(idx_ext(tipo==2)), Env(idx_ext(tipo==2)), 'ro', 'DisplayName','Máximos');
legend; 
xlabel('Tiempo (s)'); ylabel('Amplitud norm.');
title('Extremos locales sobre la envolvente');
% === 4.1 Construcción de tripletes mín–máx–mín ===
tri_samp = []; 
tri_time = []; 
tri_amp = [];
for k = 1:length(tipo)-2
    if tipo(k)==1 && tipo(k+1)==2 && tipo(k+2)==1
        i1 = idx_ext(k); i2 = idx_ext(k+1); i3 = idx_ext(k+2);
        tri_samp(end+1,:) = [i1,i2,i3];                  %#ok<SAGROW>
        tri_time(end+1,:) = [t(i1), t(i2), t(i3)];       %#ok<SAGROW>
        tri_amp(end+1,:)  = [Env(i1), Env(i2), Env(i3)]; %#ok<SAGROW>
    end
end

% === 4.2 Área de cada triángulo ===
areas = zeros(size(tri_time,1),1);
for i = 1:size(tri_time,1)
    x1=tri_time(i,1); y1=tri_amp(i,1);
    x2=tri_time(i,2); y2=tri_amp(i,2);
    x3=tri_time(i,3); y3=tri_amp(i,3);
    areas(i) = 0.5 * abs( x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2) );
end

% Visual
figure('Name','Triángulos mín–máx–mín');
plot(t, Env, 'g', 'LineWidth', 1.2); hold on; grid on;
for i = 1:size(tri_time,1)
    tx = [tri_time(i,:), tri_time(i,1)];
    ty = [tri_amp(i,:),  tri_amp(i,1)];
    plot(tx, ty, 'r-', 'LineWidth', 1.0);
end
plot(t(idx_ext), Env(idx_ext), 'ko', 'MarkerSize', 4);
xlabel('Tiempo (s)'); ylabel('Amplitud norm.');
title('Envolvente y triángulos candidatos');



% === 5.1 Selección por área (umbral: promedio) ===
Amed = mean(areas);
mask_big = areas > Amed;

% === 5.2 Propuesta de ciclos: [min1, min2] de cada triángulo grande ===
ciclos_idx = tri_samp(mask_big, [1 3]);  % índices
ciclos_idx = sortrows(ciclos_idx, 1);

% === 5.3 Limpieza por duración fisiológica y no solape ===
minRR = round(0.30*fs);   % 0.30 s
maxRR = round(1.50*fs);   % 1.50 s
ciclos_ref = [];
for i = 1:size(ciclos_idx,1)
    L = ciclos_idx(i,2) - ciclos_idx(i,1);
    if L >= minRR && L <= maxRR
        if isempty(ciclos_ref) || ciclos_idx(i,1) > ciclos_ref(end,2)
            ciclos_ref = [ciclos_ref; ciclos_idx(i,:)]; %#ok<AGROW>
        end
    end
end



% === 5.4 Visual: sombrear ciclos en la envolvente ===
figure('Name','Ciclos cardiacos detectados');
plot(t, Env, 'k'); 
grid on; 
hold on;
for i = 1:size(ciclos_ref,1)
    i1 = ciclos_ref(i,1); i2 = ciclos_ref(i,2);
    fill([t(i1) t(i2) t(i2) t(i1)], [0 0 1 1], ...
         'c', 'FaceAlpha', 0.18, 'EdgeColor','none');
end
legend('Envolvente','Ciclos'); xlabel('Tiempo (s)'); ylabel('Amplitud norm.');
title('Ciclos cardiacos (propuestos)');

% Normal y envolvente
figure('Name','Original & Envolvente (Sobrepuesta)'); 
plot(t, x, 'k'); 
hold on;
plot(t, Env, 'b');
grid on;
title('Señal PCG (Original y Envolvente)'); 
xlabel('Tiempo (s)'); ylabel('Amplitud');

% Original, envolvente y triángulos máximos con S1 y S2
figure('Name','Señal Original, envolvente y triángulos máximos');
plot(t, x, 'k'); 
hold on;
grid on;
plot(t, Env, 'r', 'LineWidth', 1.2); 
hold on; 
grid on;

% Identificación S1 y S2 basada en tiempo entre crestas
colors = {'c','r'};   % S1=cian, S2=rojo
labels = {'S1','S2'};

big_idx = find(mask_big);                     % índices de triángulos grandes
peak_times = mean(tri_time(big_idx,:), 2);    % tiempo promedio de cada triángulo
[peak_times, sort_idx] = sort(peak_times);    % ordenar por tiempo
big_idx_sorted = big_idx(sort_idx);

prev_time = -Inf;
cycle_threshold = 0.4;   % ajustar según tu señal (s)
for j = 1:length(big_idx_sorted)
    idx = big_idx_sorted(j);
    
    % Determinar si es S1 o S2
    if (peak_times(j) - prev_time) > cycle_threshold
        color = colors{1};  % S1
        label = labels{1};
    else
        color = colors{2};  % S2
        label = labels{2};
    end
    prev_time = peak_times(j);
    
    % Graficar el triángulo
    tx = [tri_time(idx,:), tri_time(idx,1)];
    ty = [tri_amp(idx,:), tri_amp(idx,1)];
    plot(tx, ty, color, 'LineWidth', 1.5);
    
    % Etiqueta en la cresta
    text(peak_times(j), max(tri_amp(idx,:)), label, 'Color', color, ...
         'FontWeight','bold', 'HorizontalAlignment','center');
end

% Graficar puntos de máximos externos si los tienes
plot(t(idx_ext), Env(idx_ext), 'ko', 'MarkerSize', 4);

xlabel('Tiempo (s)'); ylabel('Amplitud norm.');
title('Envolvente y triángulos candidatos con S1 y S2');




% Calcular s1 y s2
% === 6.1 Calcular altura de cada triángulo ===
alturas = zeros(size(tri_amp,1),1);
for i = 1:size(tri_amp,1)
    y_min = min(tri_amp(i,[1 3]));
    y_max = tri_amp(i,2);
    alturas(i) = y_max - y_min;
end

% === 6.2 Obtener los 4 triángulos de mayor altura ===
[~, idx_sorted] = sort(alturas, 'descend');
top4_idx = idx_sorted(1:4);

% === 6.3 Ordenar esos 4 triángulos por tiempo (para que sea S1, S2, S1, S2) ===
[~, order_by_time] = sort(tri_time(top4_idx,2));  % ordenar por el tiempo del pico
top4_idx = top4_idx(order_by_time);

% === 6.4 Graficar señal original, envolvente y triángulos ===
figure('Name','Señal Original, Envolvente y S1/S2');
plot(t, x, 'k'); 
hold on; 
grid on;
plot(t, Env, 'b', 'LineWidth', 1.2);

% Dibujar todos los triángulos "grandes"
for i = find(mask_big)'
    tx = [tri_time(i,:), tri_time(i,1)];
    ty = [tri_amp(i,:),  tri_amp(i,1)];
    plot(tx, ty, 'b', 'LineWidth', 0.8);
end


% === Ordenar triángulos grandes por tiempo ===
big_idx = find(mask_big);
peak_times = tri_time(big_idx,2);  % tiempo del pico
[peak_times, sort_idx] = sort(peak_times);
big_idx_sorted = big_idx(sort_idx);

% === Inicializar variables para inicios de ciclo ===
S1_inicios = [];
prev_time = -Inf;
cycle_threshold = 0.4;  % tiempo típico S1→S2 (ajustar según señal)
colors = {'c','r'};
labels = {'S1','S2'};

% === Graficar y asignar S1/S2 según tiempo entre picos ===
figure('Name','Señal Original, Envolvente y Triángulos S1/S2');
plot(t, x, 'k'); hold on; grid on;
plot(t, Env, 'b', 'LineWidth', 1.2);

for j = 1:length(big_idx_sorted)
    idx = big_idx_sorted(j);
    
    % Determinar S1 o S2 según distancia al pico anterior
    if (peak_times(j) - prev_time) > cycle_threshold
        color = colors{1}; label = labels{1};  % S1
        S1_inicios(end+1) = tri_time(idx,1);  % inicio de ciclo
    else
        color = colors{2}; label = labels{2};  % S2
    end
    prev_time = peak_times(j);
    
    % Graficar triángulo
    tx = [tri_time(idx,:), tri_time(idx,1)];
    ty = [tri_amp(idx,:), tri_amp(idx,1)];
    plot(tx, ty, color, 'LineWidth', 2);
    
    % Pico con etiqueta
    peak_t = tri_time(idx,2);
    peak_y = tri_amp(idx,2);
    plot(peak_t, peak_y, 'o', 'Color', color, 'MarkerFaceColor', color);
    text(peak_t, peak_y + 0.05, label, 'Color', color, 'FontSize', 11, ...
         'FontWeight','bold','HorizontalAlignment','center');
end

% === Marcar inicios de ciclo (S1) ===
for k = 1:length(S1_inicios)
    xline(S1_inicios(k), '--k', 'LineWidth', 1);
    plot(S1_inicios(k), Env(big_idx_sorted(1)), 'ko', 'MarkerFaceColor','y', 'MarkerSize',8);
end



% === 6.7 Marcar líneas de fin de ciclo (inicio del siguiente S1) ===
for k = 2:length(S1_inicios)
    xline(S1_inicios(k), ':r', 'LineWidth', 1.2, ...
        'DisplayName', 'Fin de ciclo');
end
