# Makefile pour un projet C++ avec CUDA et OpenMP

# Compilateur CUDA
NVCC := nvcc

# Compilateur C++
CXX := clang++

# Options de compilation
CXXFLAGS := -std=c++20 -Wall -fopenmp=libomp -O3 -funroll-loops
CUDAFLAGS := -arch=sm_52

# Options de l'édition des liens
LDFLAGS := -lm -lcudart

# Nom de l'exécutable
TARGET := spinner_code

# Liste des fichiers sources C++
SRCS_CPP := main.cpp spinner_ppm.cpp
OBJS_CPP := $(SRCS_CPP:.cpp=.o)

# Liste des fichiers sources CUDA
SRCS_CUDA := spinner_CUDA.cu
OBJS_CUDA := $(SRCS_CUDA:.cu=.o)

# Liste des fichiers d'en-tête
HEADERS := spinner_ppm.h spinner_CUDA.cuh

# Règle par défaut (première règle dans le fichier)
all: $(TARGET)

# Règle pour la cible principale
$(TARGET): $(OBJS_CPP) $(OBJS_CUDA)
	$(NVCC) $(CUDAFLAGS) -c $(SRCS_CUDA) -o $(OBJS_CUDA)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS_CPP) $(OBJS_CUDA) $(LDFLAGS)

# Règle pour les fichiers sources C++
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Règle pour les fichiers sources CUDA
%.o: %.cu $(HEADERS)
	$(NVCC) $(CUDAFLAGS) -c $< -o $@

# Règle pour nettoyer les fichiers objets
clean:
	rm -f $(OBJS_CPP) $(OBJS_CUDA) $(TARGET)

