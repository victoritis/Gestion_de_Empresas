% Definir las dimensiones del problema
numTareas = 75;
numMaquinas = 20;

% Cargar los datos del archivo
% Asegúrate de que el archivo 'datos.txt' tenga el formato correcto
Dij = dlmread('Doc1.txt');

% Verificar que los datos tienen el tamaño correcto
assert(size(Dij, 1) == numTareas && size(Dij, 2) == numMaquinas, 'Datos con dimensiones incorrectas');

% Definir parámetros de la Búsqueda Tabú
tamanoListaTabu = 10; % Este valor puede necesitar ajustes
maxIteraciones = 1000; % Ajustar según sea necesario

% Inicializar estadísticas de rendimiento
resultados = zeros(1, 5); % Para almacenar los resultados de las 5 ejecuciones





function fmed = Evaluar(orden, Dij, numTareas, numMaquinas, tiemposMaquinas)
    tiempoFinalizacion = 0;
    for i = 1:numTareas
        tarea = orden(i);
        tiemposMaquinas(1) = tiemposMaquinas(1) + Dij(tarea, 1);

        for j = 2:numMaquinas
            if tiemposMaquinas(j - 1) > tiemposMaquinas(j)
                tiemposMaquinas(j) = tiemposMaquinas(j - 1);
            end
            tiemposMaquinas(j) = tiemposMaquinas(j) + Dij(tarea, j);
        end

        tiempoFinalizacion = tiempoFinalizacion + tiemposMaquinas(numMaquinas);
    end
    
    fmed = tiempoFinalizacion / numTareas;
end
