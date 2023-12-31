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


% Calcula la duración total del proceso para el orden dado
n = length(orden); % Número de tareas
m = size(Dij,2);    %numero de maquinas
tiempos_maquinas = zeros(n, m);

for i = 1:n  %i numero de tarea
    tarea = orden(i);  
    for j = 1:m  %j numero de maquina
        if j == 1 && i == 1 %en caso de ser la primera tarea y primera maquina
            tiempos_maquinas(i,j) = Dij(tarea,j);
            disp(['Tiempo de la maquina ', num2str(j),' en tarea ' ,num2str(i) , ' es: ', num2str(tiempos_maquinas(i,j))]);
        elseif i == 1   %si es primera tarea solo
            tiempos_maquinas(i,j) = max(tiempos_maquinas(i,j), tiempos_maquinas(i,j-1));
            tiempos_maquinas(i,j) = tiempos_maquinas(i,j) + Dij(tarea,j);
            disp(['Tiempo de la maquina ', num2str(j),' en tarea ' ,num2str(i) , ' es: ', num2str(tiempos_maquinas(i,j))]);

        elseif j == 1   %en caso de ser la primera maquina solo
            tiempos_maquinas(i,j) = tiempos_maquinas(i-1,j) + Dij(tarea,j);
            disp(['Tiempo de la maquina ', num2str(j),' en tarea ' ,num2str(i) , ' es: ', num2str(tiempos_maquinas(i,j))]);
        else
            tiempos_maquinas(i,j) = max(tiempos_maquinas(i-1,j), tiempos_maquinas(i,j-1));
            tiempos_maquinas(i,j) = tiempos_maquinas(i,j) + Dij(tarea,j);
            disp(['Tiempo de la maquina ', num2str(j),' en tarea ' ,num2str(i) , ' es: ', num2str(tiempos_maquinas(i,j))]);
        end
    end
end

tiempo_total = max(tiempos_maquinas(:));
disp(['El tiempo total es : ', num2str(tiempo_total)]);
