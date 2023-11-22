% Nombre del archivo de texto
nombreArchivo = 'ejem_clase1.txt';

% Leer el archivo y almacenar los datos en 'Dij'
Dij = dlmread(nombreArchivo, ' ', 1, 0);  % Ignorar la primera línea

% Extraer el número de tareas y máquinas
[n, m] = size(Dij);

% Generar un orden aleatorio de las tareas
orden = randperm(n);

% Inicializar variables antes de llamar a Evaluar
tiempos_maquinas = zeros(1, m);
tiempo_finalizacion = 0;

% Llamar a Evaluar con todas las variables necesarias
fmed = Evaluar(orden, Dij, n, m, tiempos_maquinas, tiempo_finalizacion);


% Mostrar el resultado
disp(['Fmed para el orden aleatorio es: ', num2str(fmed)]);

function fmed = Evaluar(orden, Dij, n, m, tiempos_maquinas, tiempo_finalizacion)
    for i = 1:n
        tarea = orden(i);
        tiempos_maquinas(1) = tiempos_maquinas(1) + Dij(tarea, 1);

        for j = 2:m
            % Reemplazar la llamada a max por una operación condicional
            if tiempos_maquinas(j - 1) > tiempos_maquinas(j)
                tiempos_maquinas(j) = tiempos_maquinas(j - 1);
            end
            tiempos_maquinas(j) = tiempos_maquinas(j) + Dij(tarea, j);
        end

        tiempo_finalizacion = tiempo_finalizacion + tiempos_maquinas(m);
    end
    
    fmed = tiempo_finalizacion / n;
end



