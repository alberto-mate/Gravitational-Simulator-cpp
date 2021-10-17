/* Librerias */
#include <iostream>
#include <math.h>
#include <string>
#include <fstream>
#include <random>
#include <chrono>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

/* Constantes */
const double GRAVITY = 6.674 * pow(10, -11); // Constante gravedad universal
const double M = pow(10, 21);                // Media (distribución normal)
const double SDM = pow(10, 15);              // Desviación (distribución normal)

/* Estructura objeto */
struct object {
    double pos_x;
    double pos_y;
    double pos_z;
    double speed_x;
    double speed_y;
    double speed_z;
    double mass;
};

int main(int argc, char const *argv[])
{

    /* Comprobación inicial variables */
    if (argc != 6) {
        cerr << "Número de argumentos incorrecto\n";
        return -1;
    }

    /* Comprobación de valores iniciales de variables */
    if ((atoi(argv[1]) <= 0 || atoi(argv[2]) <= 0 || atoi(argv[3]) <= 0 || atof(argv[4]) <= 0.0 || atof(argv[5]) <= 0.0) ||
        (atof(argv[1]) != atoi(argv[1]) || atof(argv[2]) != atoi(argv[2]) || atof(argv[3]) != atoi(argv[3]))) {
        cerr << "Datos erróneos de los argumentos\n";
        return -2;
    }
    
    /* Almacenamiento de los argumentos en sus respectivas variables */
    int num_objects = atoi(argv[1]);
    int num_iterations = atoi(argv[2]);
    int random_seed = atoi(argv[3]);
    float size_enclosure = atof(argv[4]);
    float time_step = atof(argv[5]);

    /* Impresión por pantalla de las variables */
    cout << "num_objects: " << num_objects << "\n"; // Número de objetos a simular (>0 entero)
    cout << "num_iterations: " << num_iterations << "\n"; // Número de iteraciones a simular (>0 entero)
    cout << "random_seed: " << random_seed << "\n"; // Semilla para distribuciones aleatorias
    cout << "size_enclosure: " << size_enclosure << "\n"; // Tamaño del recinto (>0 real)
    cout << "time_step: " << time_step << "\n"; // Incremento de tiempo en cada iteración (>0 real)

    /* Coordenadas aleatorias */
    mt19937 gen(seed);
    uniform_real_distribution<double> position_dist{0.0, size_enclosure};
    normal_distribution<double> mass_dist{M, SDM};

    // AOS
    // Array of structs
    vector<object> objects(num_objects);

    /* Fichero de configuracion inicial */
    ofstream file_init;
    file_init.open("init_conf.txt");
    file_init << fixed << setprecision(3) << size_enclosure << " " << time_step << " " << num_objects << endl;

    // Crear los objetos
    for (int i = 0; i < num_objects; i++){
        // Igual meter variables para ver si mejora rendimiento
        objects[i].pos_x = position_dist(gen);
        objects[i].pos_y = position_dist(gen);
        objects[i].pos_z = position_dist(gen);
        objects[i].speed_x = 0;
        objects[i].speed_y = 0;
        objects[i].speed_z = 0;
        objects[i].mass = mass_dist(gen);

        file_init << fixed << setprecision(3) << objects[i].pos_x << " " << objects[i].pos_y << " " << objects[i].pos_z << " " << objects[i].speed_x << " " << objects[i].speed_y << " " << objects[i].speed_z << " " << objects[i].mass << endl;
    }

}
