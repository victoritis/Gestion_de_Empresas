% Nombre del archivo de texto que deseas abrir
nombreArchivo = 'ejem_clase2.txt';

% Intenta abrir el archivo en modo de lectura
fid = fopen(nombreArchivo, 'r');

% Verifica si el archivo se abrió correctamente
if fid == -1
    error('No se pudo abrir el archivo.');
end
% Lee la primera línea y omítela (no la almacenamos)
fgetl(fid);

% Inicializa una matriz vacía para almacenar los datos deseados
Dij = [];

% Lee el contenido del archivo y extrae las columnas impares
while ~feof(fid)
    % Lee una línea
    linea = fgetl(fid);
    
    % Convierte la línea en un vector numérico
    valores = str2double(strsplit(linea, ' '));
    
    % Extrae las columnas impares (índices 1, 3, 5, etc.)
    columnas_impares = valores(2:2:end);
    
    % Agrega las columnas impares a la matriz de datos
    Dij = [Dij; columnas_impares];
end



% Muestra la matriz resultante
disp(Dij);

orden = [4 2 5 1 3];

disp(orden);

%%%%%%%%%EN LA formula de teoria, en el primer caso, dentro del max, LA
%%%%%%%%%PARTE DE ABAJO NO EXISTE

% Cierra el archivo después de usarlo
fclose(fid);


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




