#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <random>
#include <filesystem>
#include <omp.h>

std::vector<int> generarPermutacionInicial(const std::vector<std::vector<int>>& Dij) {
    size_t n = Dij.size();
    std::vector<int> permutacion(n);
    std::iota(permutacion.begin(), permutacion.end(), 1);

    // Pre-calculando los tiempos totales
    std::vector<int> tiemposTotales(n, 0);
    for (size_t i = 0; i < n; ++i) {
        tiemposTotales[i] = std::accumulate(Dij[i].begin(), Dij[i].end(), 0);
    }

    std::sort(permutacion.begin(), permutacion.end(), [&tiemposTotales](int a, int b) {
        return tiemposTotales[a - 1] < tiemposTotales[b - 1];
    });

    return permutacion;
}



double EvaluarFmed(const std::vector<int>& orden, const std::vector<std::vector<int>>& Dij) {
    size_t n = orden.size(); // Número de tareas
    size_t m = Dij[0].size(); // Número de máquinas
    std::vector<int> tiempos_maquinas(m, 0);

    for (size_t i = 0; i < n; ++i) {
        int tarea = orden[i] - 1; // Asumiendo que orden[i] está 1-indexado
        for (size_t j = 0; j < m; ++j) {
            tiempos_maquinas[j] = std::max(j == 0 ? 0 : tiempos_maquinas[j - 1], tiempos_maquinas[j]) + Dij[tarea][j];
        }
    }

    double fmed = std::accumulate(tiempos_maquinas.begin(), tiempos_maquinas.end(), 0.0) / m;
    return fmed;
}

std::vector<std::vector<int>> leerYProcesarArchivo(const std::string& nombreArchivo) {
    std::vector<std::vector<int>> datos;
    std::ifstream archivo(nombreArchivo);
    std::string linea;

    // Omitir la primera línea
    std::getline(archivo, linea);

    while (std::getline(archivo, linea)) {
        std::istringstream iss(linea);
        std::vector<int> fila;
        int valor;
        int indice = 0;

        while (iss >> valor) {
            if (indice % 2 != 0) {
                fila.push_back(valor);
            }
            ++indice;
        }

        datos.push_back(fila);
    }

    return datos;
}



int main() {

    int numHilosMaxSistema = omp_get_max_threads(); // Obtener el máximo de hilos del sistema
    int numHilosDeseados = 8; // Tu límite deseado de hilos
    int numHilos = std::min(numHilosMaxSistema, numHilosDeseados); // Usa el menor de estos dos valores

    std::string nombreArchivo = "C:/Users/Victor/Desktop/UNI/UBU/Gestion de Empresas/Practica_Competitiva/Doc2.txt";

    auto Dij = leerYProcesarArchivo(nombreArchivo);

    if (Dij.empty() || Dij[0].empty()) {
        std::cerr << "Error al leer el archivo o archivo vacío." << std::endl;
        return 1;
    }

    int num_tareas = Dij.size();



    const int tamanoPoblacion = 80; // 8 hilos x 10 individuos por hilo
    const int numGeneraciones = 1000;
    const int individuosPorHilo = tamanoPoblacion / numHilos; // Individuos por hilo

    std::vector<std::vector<int>> secuencias(tamanoPoblacion);
    std::vector<double> aptitudes(tamanoPoblacion);

    // Inicializar la población con soluciones basadas en la heurística SPT
    for (int i = 0; i < tamanoPoblacion; ++i) {
        secuencias[i] = generarPermutacionInicial(Dij);
    }

    // Bucle principal del algoritmo genético
    for (int gen = 0; gen < numGeneraciones; ++gen) {

        #pragma omp parallel num_threads(numHilos)
        {
            int idHilo = omp_get_thread_num(); // Obtiene el ID del hilo actual
            int inicio = idHilo * individuosPorHilo; // Índice de inicio para este hilo
            int fin = inicio + individuosPorHilo; // Índice de fin para este hilo

            for (int i = inicio; i < fin; ++i) {
                // Procesamiento para cada individuo
                aptitudes[i] = EvaluarFmed(secuencias[i], Dij);
            }
        }

        // Selección y cruce (realizar estos pasos de forma adecuada)
        // ...

        // Mutación (puede ser paralelizada también si es apropiado)
        // ...

        // Reemplazo de la población
        // ...

        // (Opcional) Guardar/Imprimir resultados parciales o de seguimiento
        // ...
    }

    // (Opcional) Imprimir la mejor solución encontrada
    // ...



    return 0;
}
