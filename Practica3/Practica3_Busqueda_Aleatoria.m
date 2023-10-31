clear;

% Nombre del archivo de texto que deseas abrir
nombreArchivo = 'Doc2.txt';

% Intenta abrir el archivo en modo de lectura
fid = dlmread(nombreArchivo);


% Extrae las columnas impares (índices 1, 3, 5, etc.)
Dij = fid(:, 2:2:end);

% Eliminar la primera fila
Dij = Dij(2:end, :);

% Muestra la matriz resultante
disp(Dij);

%para la cuenta de iteraciones
iteracion = 1;

%%%%%%%%%%%%%%%%%%%% PEDIR ITERACIONES USUARIO %%%%%%%%%%%%%%%%%%%%%

% Pregunta al usuario cuántas iteraciones desea en la búsqueda aleatoria
num_iteraciones = input('Por favor, ingrese el número de iteraciones deseadas: ');

% Asegúrate de que num_iteraciones sea un número positivo
if num_iteraciones <= 0
    error('El número de iteraciones debe ser un valor positivo.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Numero de tareas, para generar el vector orden
num_tareas = size(Dij, 1);

%Generar permutacion de n a num_tareas
orden = randperm(num_tareas);

%Generamos la primera iteracion
mejor_solucion = Evaluar(orden,Dij);
mejor_orden = orden;
disp(['El tiempo total de la iteracion ', num2str(iteracion) ,' es : ', num2str(mejor_solucion), ' con orden : ', mat2str(orden)]);


for x = 1:num_iteraciones - 1
    iteracion = iteracion + 1;
    %Generar permutacion de n a num_tareas
    orden = randperm(num_tareas);
    %Generamos la primera iteracion
    solucion = Evaluar(orden,Dij);
    disp(['El tiempo total de la iteracion ', num2str(iteracion) ,' es : ', num2str(solucion), ' con orden : ', mat2str(orden)]);
    if solucion < mejor_solucion
        mejor_solucion = solucion;
        mejor_orden = orden;
    end
    
end


disp(['La mejor solucion es el orden : ',mat2str(mejor_orden),' con tiempo : ' ,num2str(mejor_solucion)]);




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





