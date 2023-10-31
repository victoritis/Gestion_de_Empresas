clear;


%%%%%%%%%%%%%%%%%%%   EL MEJOR VECINO   %%%%%%%%%%%%%%%%%%%%%%

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

%Numero de tareas, para generar el vector orden
num_tareas = size(Dij, 1);

%Generar permutacion de n a num_tareas
%orden = randperm(num_tareas);
orden = [6     3    11     7     8     5     1     2     4     9    10]
%Generamos la primera iteracion
mejor_solucion = Evaluar(orden,Dij);
mejor_orden = orden;
disp('EJECUTANDO ALGORITMO DEL MEJOR VECINO');
disp(['Orden generado aleatoriamente : ', mat2str(orden),' con tiempo : ', mat2str(mejor_solucion)]);

Vecino_mejor_Ecnontrado = true;
contador = 0;
contador_vecinos = 0;
n = length(orden);


while Vecino_mejor_Ecnontrado == true
    % Calcular y guardar los pares de vecinos 2 a 2 en cada iteración
    Vecino_mejor_Ecnontrado = false;
    orden_a_cambiar = mejor_orden;
    for i = 1:n-1
        for j = i+1:n
            orden = orden_a_cambiar;
            contador_vecinos = contador_vecinos + 1;
            A = orden(i);
            B = orden(j);
            orden(i) = B;
            orden(j) = A;
            solucion = Evaluar(orden,Dij);
            if solucion < mejor_solucion
                mejor_solucion = solucion;
                mejor_orden = orden;
                Vecino_mejor_Ecnontrado = true;
            end
        end
    end
    contador = contador + 1;
end

disp(['La mejor solucion encontrada es el orden : ',mat2str(mejor_orden),' con tiempo : ' ,num2str(mejor_solucion)]);
disp(['Se han explorado ',mat2str(contador),' veces los vecinos', '(Un total de ',mat2str(contador_vecinos), ' vecinos)' ]);



%%%%%%%%%%%%%%%%%%%%% FUNCION EVALUAR %%%%%%%%%%%%%%%%%%%%%%%%

function max_tiempo = Evaluar(orden, Dij)
    n = length(orden); % Número de tareas
    m = size(Dij, 2);  % Número de máquinas
    tiempos_maquinas = zeros(1, m);

    for i = 1:n
        tarea = orden(i);
        for j = 1:m
            if j == 1
                tiempos_maquinas(j) = tiempos_maquinas(j) + Dij(tarea, j);                
            else
                tiempos_maquinas(j) = max(tiempos_maquinas(j), tiempos_maquinas(j - 1));
                tiempos_maquinas(j) = tiempos_maquinas(j) + Dij(tarea, j);
            end
        end
    end
    
    max_tiempo = max(tiempos_maquinas);
end





