#include <iostream>
#include <vector>
#include <string>
#include <cstdio>
#include <cmath>       // Para a função tan()
#include <algorithm>   // Para std::max

// GLEW (precisa vir antes do GLFW)
#include <GL/glew.h>
// GLFW
#include <GLFW/glfw3.h>
// GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

// STB Image Write
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

// Nossas próprias estruturas e o parser
#include "DataTypes.h"

// --- Constantes ---
const int SCREEN_WIDTH = 1920;
const int SCREEN_HEIGHT = 1080;
const float PI = 3.14159265359f;

// --- Shaders para as Partículas (circulares) ---
const char* particleVertexShaderSource = R"(
    #version 330 core
    layout (location = 0) in vec2 aPos;
    uniform mat4 projection;
    void main() {
        gl_Position = projection * vec4(aPos, 0.0, 1.0);
        gl_PointSize = 15.0; // Diametro de 15px (raio 7.5)
    }
)";

const char* particleFragmentShaderSource = R"(
    #version 330 core
    out vec4 FragColor;
    uniform vec3 particleColor;
    void main() {
        vec2 coord = gl_PointCoord - vec2(0.5);
        if (dot(coord, coord) > 0.25) {
            discard;
        }
        FragColor = vec4(particleColor, 1.0);
    }
)";

// --- Shaders para as Retas (simples) ---
const char* lineVertexShaderSource = R"(
    #version 330 core
    layout (location = 0) in vec2 aPos;
    uniform mat4 projection;
    void main() {
        gl_Position = projection * vec4(aPos, 0.0, 1.0);
    }
)";

const char* lineFragmentShaderSource = R"(
    #version 330 core
    out vec4 FragColor;
    uniform vec3 lineColor;
    void main() {
        FragColor = vec4(lineColor, 1.0);
    }
)";


// --- Funções Auxiliares ---
void saveFrame(const char* filename, int width, int height);
unsigned int createShaderProgram(const char* vertexSource, const char* fragmentSource);


int main(int argc, char* argv[]) {
    // --- 1. Validação dos Argumentos de Entrada ---
    if (argc != 4) {
        std::cerr << "Uso: " << argv[0] << " <arquivo.txt> <abertura> <angulo>" << std::endl;
        return 1;
    }
    
    std::string dataFilepath = argv[1];
    float abertura, angulo;

    try {
        abertura = std::stof(argv[2]);
        angulo = std::stof(argv[3]);
    } catch (const std::exception& e) {
        std::cerr << "ERRO: 'abertura' e 'angulo' devem ser valores numericos." << std::endl;
        std::cerr << "Recebido: abertura=" << argv[2] << ", angulo=" << argv[3] << std::endl;
        return 1;
    }
    
    // --- 2. Inicialização (GLFW, GLEW) ---
    if (!glfwInit()) {
        std::cerr << "Falha ao inicializar GLFW" << std::endl;
        return -1;
    }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    GLFWwindow* window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Visualizacao Cientifica de Particulas", NULL, NULL);
    if (!window) {
        std::cerr << "Falha ao criar a janela GLFW" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    if (glewInit() != GLEW_OK) {
        std::cerr << "Falha ao inicializar GLEW" << std::endl;
        glfwTerminate();
        return -1;
    }

    // --- 3. Carregar Dados da Simulação ---
    TimeSeriesData simulationData = loadParticleData(dataFilepath);
    if (simulationData.empty()) {
        glfwTerminate();
        return -1;
    }
    
    std::vector<std::vector<ParticleState>> timeSteps;
    float minX = 1e9, maxX = -1e9, minY = 1e9, maxY = -1e9;
    for (const auto& pair : simulationData) {
        timeSteps.push_back(pair.second);
    }
    for (const auto& p : timeSteps[0]) {
        if (p.x < minX) minX = p.x;
        if (p.x > maxX) maxX = p.x;
        if (p.y < minY) minY = p.y;
        if (p.y > maxY) maxY = p.y;
    }
    
    // --- 4. Configurar Geometria das Retas ---
    std::cout << "Configurando geometria com abertura=" << abertura << ", angulo=" << angulo << std::endl;
    const float L1 = abertura * 7.5f * 2.0f;
    const float y0 = 98.0f + 154.0f * tan(angulo * PI / 180.0f);

    std::vector<float> lineVertices = {
        0.0f, 0.0f,      0.0f, 98.0f,        // Reta 1
        L1,   0.0f,      L1,   98.0f,        // Reta 2
        0.0f, 98.0f,    -154.0f, y0,         // Reta 3
       -154.0f, y0,     -154.0f, 910.0f,     // Reta 4
        L1,   98.0f,     L1 + 154.0f, y0,     // Reta 5
        L1 + 154.0f, y0,  L1 + 154.0f, 910.0f  // Reta 6
    };

    // --- 5. Configurar Shaders e Buffers OpenGL ---
    unsigned int particleShaderProgram = createShaderProgram(particleVertexShaderSource, particleFragmentShaderSource);
    unsigned int lineShaderProgram = createShaderProgram(lineVertexShaderSource, lineFragmentShaderSource);
    
    unsigned int particleVAO, particleVBO;
    glGenVertexArrays(1, &particleVAO);
    glGenBuffers(1, &particleVBO);
    glBindVertexArray(particleVAO);
    glBindBuffer(GL_ARRAY_BUFFER, particleVBO);
    glBufferData(GL_ARRAY_BUFFER, timeSteps[0].size() * 2 * sizeof(float), NULL, GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    unsigned int lineVAO, lineVBO;
    glGenVertexArrays(1, &lineVAO);
    glGenBuffers(1, &lineVBO);
    glBindVertexArray(lineVAO);
    glBindBuffer(GL_ARRAY_BUFFER, lineVBO);
    glBufferData(GL_ARRAY_BUFFER, lineVertices.size() * sizeof(float), lineVertices.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    // --- 6. Configurar a Câmera / Projeção ---
    // Expande os limites para garantir que as retas também fiquem visíveis
    for (size_t i = 0; i < lineVertices.size(); i += 2) {
        float vx = lineVertices[i];
        float vy = lineVertices[i+1];
        if (vx < minX) minX = vx;
        if (vx > maxX) maxX = vx;
        if (vy < minY) minY = vy;
        if (vy > maxY) maxY = vy;
    }

    // ATUALIZADO: Correção da proporção de aspecto (aspect ratio)
    float worldWidth = maxX - minX;
    float worldHeight = maxY - minY;
    float worldCenterX = minX + worldWidth / 2.0f;
    float worldCenterY = minY + worldHeight / 2.0f;

    float screenAspect = (float)SCREEN_WIDTH / (float)SCREEN_HEIGHT;

    float viewWidth, viewHeight;

    // Determina se a largura ou a altura do mundo deve ser a dimensão de referência
    if (worldWidth / worldHeight >= screenAspect) {
        // O mundo é mais "largo" que a tela, então a largura domina
        viewWidth = worldWidth;
        viewHeight = worldWidth / screenAspect;
    } else {
        // O mundo é mais "alto" que a tela, então a altura domina
        viewHeight = worldHeight;
        viewWidth = worldHeight * screenAspect;
    }

    // Adiciona uma margem de 10%
    viewWidth *= 1.1f;
    viewHeight *= 1.1f;

    float orthoLeft   = worldCenterX - viewWidth / 2.0f;
    float orthoRight  = worldCenterX + viewWidth / 2.0f;
    float orthoBottom = worldCenterY - viewHeight / 2.0f;
    float orthoTop    = worldCenterY + viewHeight / 2.0f;

    glm::mat4 projection = glm::ortho(orthoLeft, orthoRight, orthoBottom, orthoTop, -1.0f, 1.0f);


    glEnable(GL_PROGRAM_POINT_SIZE);

    // --- 7. Loop de Renderização e Captura ---
    int currentFrameIndex = 0;
    int savedFrameCounter = 0;
    while (currentFrameIndex < timeSteps.size() && !glfwWindowShouldClose(window)) {
        glClearColor(0.05f, 0.05f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        // Desenhar as Retas (fundo)
        glUseProgram(lineShaderProgram);
        glUniformMatrix4fv(glGetUniformLocation(lineShaderProgram, "projection"), 1, GL_FALSE, &projection[0][0]);
        glUniform3f(glGetUniformLocation(lineShaderProgram, "lineColor"), 0.8f, 0.7f, 0.2f);
        glLineWidth(2.0f);
        glBindVertexArray(lineVAO);
        glDrawArrays(GL_LINES, 0, lineVertices.size() / 2);

        // Desenhar as Partículas (frente)
        const auto& currentParticles = timeSteps[currentFrameIndex];
        std::vector<glm::vec2> positions;
        positions.reserve(currentParticles.size());
        for (const auto& p : currentParticles) {
            positions.emplace_back(p.x, p.y);
        }
        glBindBuffer(GL_ARRAY_BUFFER, particleVBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, positions.size() * sizeof(glm::vec2), positions.data());
        
        glUseProgram(particleShaderProgram);
        glUniformMatrix4fv(glGetUniformLocation(particleShaderProgram, "projection"), 1, GL_FALSE, &projection[0][0]);
        glUniform3f(glGetUniformLocation(particleShaderProgram, "particleColor"), 0.4f, 0.8f, 1.0f);
        
        glBindVertexArray(particleVAO);
        glDrawArrays(GL_POINTS, 0, positions.size());

        // Salva o frame
        char filename[256];
        snprintf(filename, sizeof(filename), "output/frame_%05d.png", savedFrameCounter++);
        saveFrame(filename, SCREEN_WIDTH, SCREEN_HEIGHT);

        currentFrameIndex++;
        glfwPollEvents();
        glfwSwapBuffers(window);
    }

    std::cout << "Renderizacao e captura de todos os frames concluidas!" << std::endl;

    // --- 8. Limpeza ---
    glDeleteVertexArrays(1, &particleVAO);
    glDeleteBuffers(1, &particleVBO);
    glDeleteProgram(particleShaderProgram);
    glDeleteVertexArrays(1, &lineVAO);
    glDeleteBuffers(1, &lineVBO);
    glDeleteProgram(lineShaderProgram);

    glfwTerminate();
    return 0;
}

// --- Implementação das Funções Auxiliares ---

void saveFrame(const char* filename, int width, int height) {
    unsigned char* pixels = new unsigned char[width * height * 3];
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels);
    stbi_flip_vertically_on_write(true);
    if (!stbi_write_png(filename, width, height, 3, pixels, width * 3)) {
        std::cerr << "ERRO: Falha ao salvar o frame " << filename << std::endl;
    } else {
        std::cout << "Frame salvo: " << filename << std::endl;
    }
    delete[] pixels;
}

unsigned int createShaderProgram(const char* vertexSource, const char* fragmentSource) {
    int success;
    char infoLog[512];

    unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexSource, NULL);
    glCompileShader(vertexShader);
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        std::cerr << "ERRO::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
    }

    unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentSource, NULL);
    glCompileShader(fragmentShader);
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
        std::cerr << "ERRO::SHADER::FRAGMENT::COMPILATION_FAILED\n" << infoLog << std::endl;
    }

    unsigned int shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
        std::cerr << "ERRO::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
    }
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    return shaderProgram;
}