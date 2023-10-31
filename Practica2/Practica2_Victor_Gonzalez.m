clear;

% Nombre del archivo de texto que deseas abrir
nombreArchivo = 'Doc1.txt';

% Intenta abrir el archivo en modo de lectura
fid = dlmread(nombreArchivo);


% Extrae las columnas impares (índices 1, 3, 5, etc.)
Dij = fid(:, 2:2:end);

% Eliminar la primera fila
Dij = Dij(2:end, :);

% Muestra la matriz resultante
disp(Dij);



orden = [9    4    1    2    8    3   10    7    5   11    6];

disp(orden);

%%%%%%%%%EN LA formula de teoria, en el primer caso, dentro del max, LA
%%%%%%%%%PARTE DE ABAJO NO EXISTE



% Calcula la duración total del proceso para el orden dado
n = length(orden); % Número de tareas
m = size(Dij,2);    %numero de maquinas
tiempos_maquinas = zeros(1, m);

for i = 1:n
    tarea = orden(i);
    for j = 1:m
        if j == 1
            tiempos_maquinas(j) = tiempos_maquinas(j) + Dij(tarea,j);
            disp(['Tiempo de la maquina ', num2str(j), ' es: ', num2str(tiempos_maquinas(j))]);
        else
            tiempos_maquinas(j) = max(tiempos_maquinas(j), tiempos_maquinas(j-1));
            tiempos_maquinas(j) = tiempos_maquinas(j) + Dij(tarea,j);
            disp(['Tiempo de la maquina ', num2str(j), ' es: ', num2str(tiempos_maquinas(j))]);
        end
    end
end

tiempo_total = max(tiempos_maquinas);
disp(['El tiempo total es : ', num2str(tiempo_total)]);




