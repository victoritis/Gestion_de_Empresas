clear;

% Nombre del archivo de texto que deseas abrir
nombreArchivo = 'ejem_clase1.txt';

% Intenta abrir el archivo en modo de lectura
fid = dlmread(nombreArchivo);

% Extrae las columnas impares (índices 1, 3, 5, etc.)
Dij = fid(:, 2:2:end);

% Elimina la primera fila
Dij = Dij(2:end, :);

% Muestra la matriz resultante
disp(Dij);

orden = [4 2 5 1 3];
disp(orden);

% Calcula la duración total del proceso para el orden dado
n = length(orden); % Número de tareas
m = size(Dij, 2);   % Número de máquinas

% Crea una matriz para almacenar los tiempos por tarea para cada máquina
tiempos_maquinas = zeros(n, m);

for i = 1:n
    tarea = orden(i);
    for j = 1:m
        if j == 1 && i == 1
            tiempos_maquinas(i, j) = Dij(tarea, j);
        elseif j == 1
            tiempos_maquinas(i, j) = tiempos_maquinas(i - 1,j) + Dij(tarea, j);
        elseif i == 1
            tiempos_maquinas(i, j) = tiempos_maquinas(i, j - 1) + Dij(tarea, j);
        else
            tiempos_maquinas(i, j) = max(tiempos_maquinas(i, j - 1), tiempos_maquinas(i, j)) + Dij(tarea, j);
        end
    end
end

% Muestra los tiempos por tarea para cada máquina
for j = 1:m
    disp(['Tiempos de la máquina ', num2str(j), ' para cada tarea:']);
    disp(tiempos_maquinas(j,:));
end

% Tiempo total es el último elemento de la última máquina
tiempo_total = max(tiempos_maquinas(:));
disp(['El tiempo total es: ', num2str(tiempo_total)]);

