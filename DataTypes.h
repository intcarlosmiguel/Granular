#include <vector>
#include <map>
#include <string>

// Representa o estado de uma única partícula em um dado momento.
// Usamos vec2 do GLM para posições, é uma boa prática.
struct ParticleState {
    int id;
    float x, y;
    float velocity;
};

// Mapeia um instante de tempo (float) para um vetor de todas as partículas naquele instante.
using TimeSeriesData = std::map<float, std::vector<ParticleState>>;

// Função para carregar os dados do arquivo.
TimeSeriesData loadParticleData(const std::string& filepath);