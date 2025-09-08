#include "DataTypes.h"
#include <iostream>
#include <fstream>
#include <sstream>

TimeSeriesData loadParticleData(const std::string& filepath) {
    TimeSeriesData data;
    std::ifstream file(filepath);

    if (!file.is_open()) {
        std::cerr << "ERRO: Nao foi possivel abrir o arquivo: " << filepath << std::endl;
        return data;
    }

    std::cout << "Lendo arquivo de dados: " << filepath << std::endl;

    std::string line;
    // std::getline(file, line); // Descomente se tiver cabeçalho

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        std::stringstream ss(line);
        int id;
        float x, y, vel, time;
        
        if (!(ss >> id >> x >> y >> vel >> time)) {
            std::cerr << "AVISO: Linha mal formatada, pulando: " << line << std::endl;
            continue;
        }

        // MODIFICADO: Multiplicando as coordenadas por 1000
        x *= 1000.0f;
        y *= 1000.0f;

        data[time].push_back({id, x, y, vel});
    }
    
    file.close();

    if (data.empty()) {
        std::cerr << "AVISO: Nenhum dado foi carregado." << std::endl;
    } else {
        std::cout << "Dados carregados com sucesso. " << data.size() << " passos de tempo encontrados." << std::endl;
    }
    
    return data;
}