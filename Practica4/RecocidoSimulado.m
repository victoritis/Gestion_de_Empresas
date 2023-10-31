clear;


%%%%%%%%%%%%%%%%%%%   RECOCIDO SIMULADO   %%%%%%%%%%%%%%%%%%%%%%

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
orden = [6     3    11     7     8     5     1     2     4     9    10]
%Generamos la primera iteracion
mejor_solucion = Evaluar(orden,Dij);
mejor_orden = orden;
disp('EJECUTANDO ALGORITMO RECOCIDO SIMULADO');
disp(['Orden generado aleatoriamente : ', mat2str(orden),' con tiempo : ', mat2str(mejor_solucion)]);


T0 = 1000;   %temperatura inicial
Tf = 0.00001; %temperatura final
reduccion = 0.99;  %reduccion de alfa y
n = length(orden); % Número de tareas
sol_act = orden;
e = exp(1);  %para la funcion probabilidad de aceptación de solucion peor   %%%%QUE DEBERIA VALER???????????????????????????????????????????????????????????????????????????????????'
L = 100;
T = T0;

while T >= Tf    %mientras no llegue a la T final
    
    i = L;   %calculamos numero iteraciones con L
    for j = 1:i
        
        sol_cand = sol_act;
        p1 = randi(n); % Genera un número aleatorio entre 1 y 'n' n esel numero de tareas
        p2 = randi(n); % Genera un número aleatorio entre 1 y 'n' n esel numero de tareas
        while p1 == p2
            p1 = randi(n); % Genera un número aleatorio entre 1 y 'n' n esel numero de tareas
            p2 = randi(n); % Genera un número aleatorio entre 1 y 'n' n esel numero de tareas
        end
        A = sol_cand(p1);
        B = sol_cand(p2);
        sol_cand(p1) = B;  %cambiamos de posicion 
        sol_cand(p2) = A;  %cambiamos de posicion
        
        delta = Evaluar(sol_cand, Dij) - Evaluar(sol_act, Dij);

        if delta < 0
            sol_act = sol_cand;

        elseif rand < e^(-delta/T) 
            sol_act = sol_cand;
        else
            sol_act = sol_act;
        end

    end

T = alfa(T, reduccion);

end

mejor_solucion = Evaluar(sol_act, Dij);
disp(['La mejor solucion encontrada es el orden : ',mat2str(sol_act),' con tiempo : ' ,num2str(mejor_solucion)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 % APLICAMOS BUSQUEDA LOCAL DEL MEJOR VECINO %
%Generamos la primera iteracion
mejor_solucion = Evaluar(sol_act,Dij);
mejor_orden = sol_act;
disp('EJECUTANDO ALGORITMO DEL MEJOR VECINO SOBRE SOLUCION ANTERIOR');
disp(['Orden calculado en recocido : ', mat2str(mejor_orden),' con tiempo : ', mat2str(mejor_solucion)]);

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


%%%%%%%%%%%%%%    FUNCION ALFA   %%%%%%%%%%%%%%%%%%%%%%%%%%


function resultado = alfa(temperatura, reduccion)
    resultado = temperatura * reduccion;
end


