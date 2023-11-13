clear;


%%%%%%%%%%%%%%%%%%%   ALGORITMO GENETICO   %%%%%%%%%%%%%%%%%%%%%%

disp('EJECUTANDO ALGORITMO GENETICO');

% Nombre del archivo de texto que deseas abrir
nombreArchivo = 'Doc3.txt';

% Intenta abrir el archivo en modo de lectura
fid = dlmread(nombreArchivo);


% Extrae las columnas impares (índices 1, 3, 5, etc.)
Dij = fid(:, 2:2:end);

% Eliminar la primera fila
Dij = Dij(2:end, :);

%Numero de tareas, para generar la poblacion
num_tareas = size(Dij, 1);

% Definir el tamaño de la población
tam_poblacion = 50; % Puedes ajustar este valor según sea necesario

% Generar la población inicial
poblacion_inicial = CrearPoblacionInicial(num_tareas, tam_poblacion);

% Calcular el fitness para cada individuo en la población inicial
fitness_poblacion = CalcularFitnessPoblacion(poblacion_inicial, Dij);

% Definir el número de generaciones
num_generaciones = 100; % Puedes ajustar este valor

% Parámetros de la mutación
prob_mutacion = 0.05; % Ajusta esta probabilidad según sea necesario

% Número de individuos élite
num_elite = 5;

poblacion = poblacion_inicial;


% Bucle principal del algoritmo genético
for gen = 1:num_generaciones

    % Calcular el fitness inicial
    fitness_poblacion = CalcularFitnessPoblacion(poblacion, Dij);

    % Ordenar la población según el fitness y seleccionar los mejores
    [fitness_ordenado, indices] = sort(fitness_poblacion);
    elite = poblacion(indices(1:num_elite), :);

    % Selección
    num_padres = tam_poblacion; % Número de padres a seleccionar, ajustar según sea necesario
    padres_seleccionados = SeleccionarPadres(num_tareas, fitness_poblacion, num_padres);

    % Inicializar una nueva población para los hijos
    nueva_poblacion = zeros(size(poblacion));

    % Cruce
    for i = 1:2:size(padres_seleccionados, 1)
        % Selecciona dos padres al azar de los seleccionados
        indices = randperm(size(padres_seleccionados, 1), 2);
        padre1 = poblacion(padres_seleccionados(indices(1)), :);
        padre2 = poblacion(padres_seleccionados(indices(2)), :);

        % Realiza el cruce
        [hijo1, hijo2] = CruceDeOrden(padre1, padre2);

        % Añade los hijos a la nueva población
        nueva_poblacion(i, :) = hijo1;
        if i+1 <= size(poblacion, 1)
            nueva_poblacion(i+1, :) = hijo2;
        end
    end

    % Actualizar la población con los nuevos hijos
    poblacion = nueva_poblacion;

    % Mutación
    % Aplicar mutación a la población actualizada
    for i = 1:size(poblacion, 1)
        poblacion(i, :) = Mutar(poblacion(i, :), prob_mutacion);
    end

    % Reemplazar los peores individuos con los élite
    poblacion(end-num_elite+1:end, :) = elite;

    % Calcular el fitness para la población actualizada
    fitness_poblacion = CalcularFitnessPoblacion(poblacion, Dij);

    % (Opcional) Guardar/Imprimir resultados parciales o de seguimiento
    % ... (puedes imprimir o guardar la mejor solución encontrada hasta ahora) ...

end

% Calcular el fitness final de la población
fitness_poblacion = CalcularFitnessPoblacion(poblacion, Dij);

% Encontrar el mejor individuo
[mejor_fitness, indice_mejor] = min(fitness_poblacion);
mejor_individuo = poblacion(indice_mejor, :);

% Mostrar el mejor individuo y su fitness
disp(['La mejor solucion encontrada es el orden : ',mat2str(mejor_individuo),' con tiempo : ' ,num2str(Evaluar(mejor_individuo,Dij))]);
disp(['  ']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 % APLICAMOS BUSQUEDA LOCAL DEL MEJOR VECINO %
%Generamos la primera iteracion
mejor_solucion = Evaluar(mejor_individuo,Dij);
mejor_orden = mejor_individuo;
disp('EJECUTANDO ALGORITMO DEL MEJOR VECINO SOBRE SOLUCION ANTERIOR');
disp(['Orden calculado en algoritmo genetico : ', mat2str(mejor_individuo),' con tiempo : ', mat2str(Evaluar(mejor_individuo,Dij))]);

Vecino_mejor_Ecnontrado = true;
contador = 0;
contador_vecinos = 0;
n = length(mejor_orden);


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
disp(['  ']);
if mejor_orden == mejor_individuo
    disp(['El algoritmo de busqueda local no implementa ninguna mejora sobre la solucion del algoritmo genetico para este caso']);
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%%%%%%%%%%%%%%%%%%% FUNCION EVALUAR FITNESS %%%%%%%%%%%%%%%%%%%

% Calcular el fitness de toda la población
function fitness = CalcularFitnessPoblacion(poblacion, Dij)
    num_soluciones = size(poblacion, 1);
    fitness = zeros(1, num_soluciones);
    for i = 1:num_soluciones
        solucion = poblacion(i, :);
        fitness(i) = Evaluar(solucion, Dij);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%% CREAR POBLACION INICIAL %%%%%%%%%%%%%%%%%%%%%%%%

% Crear una población inicial
function poblacion = CrearPoblacionInicial(num_tareas, tam_poblacion)
    poblacion = zeros(tam_poblacion, num_tareas);
    for i = 1:tam_poblacion
        poblacion(i, :) = randperm(num_tareas);
    end
end


%%%%%%%%%%%%%%%%%% SELECCION POR RULETA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Función de selección tipo ruleta
function padres = SeleccionarPadres(num_tareas, fitness_poblacion, num_padres)
    % Calcular la probabilidad de selección para cada individuo
    probabilidad = 1 ./ fitness_poblacion; % Aquí, menor fitness significa mayor probabilidad
    probabilidad = probabilidad / sum(probabilidad);

    % Seleccionar padres
    padres = zeros(num_padres, num_tareas);
    for i = 1:num_padres
        acumulado = cumsum(probabilidad);
        r = rand();
        for j = 1:length(acumulado)
            if r <= acumulado(j)
                padres(i, :) = j;
                break;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%% FUNCION CRUCE %%%%%%%%%%%%%%%%%%%%%%%%

% Función de cruce de orden (Order Crossover - OX)
function [hijo1, hijo2] = CruceDeOrden(padre1, padre2)
    num_tareas = length(padre1);
    
    % Seleccionar dos puntos de corte aleatorios
    puntos_corte = sort(randperm(num_tareas, 2));
    inicio = puntos_corte(1);
    fin = puntos_corte(2);

    % Crear hijos parcialmente vacíos
    hijo1 = zeros(1, num_tareas);
    hijo2 = zeros(1, num_tareas);

    % Copiar la subsecuencia del primer padre a cada hijo
    hijo1(inicio:fin) = padre1(inicio:fin);
    hijo2(inicio:fin) = padre2(inicio:fin);

    % Completar los hijos con los elementos restantes del otro padre
    hijo1 = CompletarHijo(hijo1, padre2, inicio, fin, num_tareas);
    hijo2 = CompletarHijo(hijo2, padre1, inicio, fin, num_tareas);
end

% Función para completar un hijo en el cruce de orden
function hijo = CompletarHijo(hijo, padre, inicio, fin, num_tareas)
    posicion_hijo = 1;
    for i = 1:num_tareas
        if posicion_hijo == inicio
            posicion_hijo = fin + 1;
        end
        if ~ismember(padre(i), hijo)
            hijo(posicion_hijo) = padre(i);
            posicion_hijo = posicion_hijo + 1;
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%% FUNCION MUTACION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Función de mutación
function mutado = Mutar(individuo, prob_mutacion)
    mutado = individuo;
    if rand() < prob_mutacion
        % Selecciona dos puntos al azar para intercambiar
        puntos = randperm(length(individuo), 2);
        temp = mutado(puntos(1));
        mutado(puntos(1)) = mutado(puntos(2));
        mutado(puntos(2)) = temp;
    end
end
