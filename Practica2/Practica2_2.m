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

% Inicializa variables para almacenar la mejor permutación y su makespan
mejor_permutacion = [];
mejor_makespan = Inf;

% Genera todas las permutaciones posibles de las tareas
permutaciones = perms(orden);

% Itera a través de las permutaciones y calcula el makespan para cada una
for i = 1:size(permutaciones, 1)
    permutacion_actual = permutaciones(i, :);
    makespan_actual = calcularMakespan(Dij, permutacion_actual);
    
    % Compara con el mejor makespan encontrado hasta ahora
    if makespan_actual < mejor_makespan
        mejor_makespan = makespan_actual;
        mejor_permutacion = permutacion_actual;
    end
end

% Muestra la mejor permutación y su makespan
disp('Mejor Permutación:');
disp(mejor_permutacion);
disp(['Mejor Makespan: ', num2str(mejor_makespan)]);

% Función para calcular el makespan dado una matriz de tiempos y una permutación
function makespan = calcularMakespan(Dij, permutacion)
    n = size(Dij, 1);  % Número de tareas
    m = size(Dij, 2);  % Número de máquinas
    tiempo_acumulado = zeros(1, m);
    
    for i = 1:n
        tarea = permutacion(i);
        tiempo_maquinas = Dij(tarea, :);
        tiempo_acumulado = tiempo_acumulado + tiempo_maquinas;
    end
    
    makespan = max(tiempo_acumulado);
end
