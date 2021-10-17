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

/* CONSTANTES */
const double GRAVITY_CONST = 6.674 * pow(10, -11); // Constante gravedad universal
const double M = pow(10, 21);                      // Media (distribución normal)
const double SDM = pow(10, 15);                    // Desviación (distribución normal)

/* ESTRUCTURA OBJETO*/
struct object
{
    double pos_x;
    double pos_y;
    double pos_z;
    double speed_x;
    double speed_y;
    double speed_z;
    double mass;
};

/* MAIN */
int main(int argc, char const *argv[])
{

    /* Comprobación inicial variables */
    if (argc != 6)
    {
        cerr << "Número de argumentos incorrecto\n";
        return -1;
    }

    /* Comprobación de valores iniciales de variables */
    if ((atoi(argv[1]) <= 0 || atoi(argv[2]) <= 0 || atoi(argv[3]) <= 0 || atof(argv[4]) <= 0.0 || atof(argv[5]) <= 0.0) ||
        (atof(argv[1]) != atoi(argv[1]) || atof(argv[2]) != atoi(argv[2]) || atof(argv[3]) != atoi(argv[3])))
    {
        cerr << "Datos erróneos de los argumentos\n";
        return -2;
    }

    /* Almacenamiento de los argumentos en sus respectivas variables */
    int num_objects = atoi(argv[1]);      // Número de objetos a simular (>0 entero)
    int num_iterations = atoi(argv[2]);   // Número de iteraciones a simular (>0 entero)
    int random_seed = atoi(argv[3]);      // Semilla para distribuciones aleatorias
    float size_enclosure = atof(argv[4]); // Tamaño del recinto (>0 real)
    float time_step = atof(argv[5]);      // Incremento de tiempo en cada iteración (>0 real)

    /* Impresión por pantalla de las variables */
    cout << "num_objects: " << num_objects << "\n";
    cout << "num_iterations: " << num_iterations << "\n";
    cout << "random_seed: " << random_seed << "\n";
    cout << "size_enclosure: " << size_enclosure << "\n";
    cout << "time_step: " << time_step << "\n";

    /* Coordenadas aleatorias */
    mt19937_64 gen(random_seed);
    uniform_real_distribution<double> position_dist(0.0, nextafter(size_enclosure, std::numeric_limits<double>::max()));
    normal_distribution<double> mass_dist(M, SDM);

    /* AOS - Array of Structs */
    vector<object> objects(num_objects);

    /* Fichero de configuracion inicial */
    ofstream file_init;
    file_init.open("init_config.txt");
    file_init << fixed << setprecision(3) << size_enclosure << " " << time_step << " " << num_objects << endl;

    /* Creación de objetos */
    for (int i = 0; i < num_objects; i++)
    {
        // Igual meter variables para ver si mejora rendimiento
        objects[i].pos_x = position_dist(gen);
        objects[i].pos_y = position_dist(gen);
        objects[i].pos_z = position_dist(gen);
        objects[i].speed_x = 0;
        objects[i].speed_y = 0;
        objects[i].speed_z = 0;
        objects[i].mass = mass_dist(gen);

        // Ponemos la precisión a 3 decimales. Imprimimos el objeto
        file_init << fixed << setprecision(3) << objects[i].pos_x << " " << objects[i].pos_y << " " << objects[i].pos_z << " " << objects[i].speed_x << " " << objects[i].speed_y << " " << objects[i].speed_z << " " << objects[i].mass << endl;
    }

    file_init.close(); // Cerramos el fichero

    /* Matriz de fuerzas */
    static double force_matrix[num_objects - 1];
    for (int i = 0; i < (num_objects - 1); i++){
        static double force_vector[num_objects - 1 - i];
        force_matrix[i] = (force_vector);
    }  
    /*

    */
} 

/* FUNCIONES */
/* Distancia euclídea entre dos objetos */
double euclidean_norm(object object_1, object object_2)
{
    return sqrt(pow(object_1.pos_x - object_2.pos_x, 2) + pow(object_1.pos_y - object_2.pos_y, 2) + pow(object_1.pos_z - object_2.pos_z, 2));
}

/* Fuerza gravitatoria entre dos objetos */
double *vector_gravitational_force(object object_1, object object_2)
{
    double force_module = (GRAVITY_CONST * object_1.mass * object_2.mass) / (pow(euclidean_norm(object_1, object_2), 3));
    static double result[3] = {force_module * (object_1.pos_x - object_2.pos_x), force_module * (object_1.pos_y - object_2.pos_y), force_module * (object_1.pos_y - object_2.pos_y)};

    return result;
}

/* Vector de aceleración */
double *vector_acceleration(object object_1, )
{
    vector<double> suma(3, 0.0)
    for (int i = 0; i < num_objects; i++){
        if (object_array[i] != object_1){
            suma += gravitational_force(object[i], object_1);
        }
    }
    return suma;
    
}

/* Vector velocidad */
double *vector_speed(object object_1, time_step)
{
    /* Calculamos la aceleración del objeto para obtener la velocidad */
    static double new_acceleration[3] = vector_acceleration(object_1);

    /* Cálculo del vector velocidad */
    static double new_speed[3] = {
        object_1.speed_x + (new_acceleration[0] * time_step), 
        object_1.speed_y + (new_acceleration[1] * time_step),
        object_1.speed_z + (new_acceleration[2] * time_step)};

    /* Devolvemos el vector velocidad */
    return new_speed;
}

/* Vector de posicion */
double *vector_posicion(object object_1, time_step)
{
    /* Calculamos la velocidad del objeto para obtener la posición */
    static double new_speed[3] = vector_speed(object_1);

    /* Cálculo del vector posición */
    static double new_position[3] = {object_1.pos_x + (new_speed[0] * time_step), object_1.pos_y + (new_speed[1] * time_step), object_1.pos_z + (new_speed[2] * time_step)};

    /* Devolvemos el vector posición */
    return new_posicion;
}

/* Función para recolocar al objeto si traspasa los límites */
void check_border(object object_1, size_enclosure){

    /* Checks posición x */
    if (object_1.pos_x <= 0){
        object_1.pos_x = 0;
        object_1.speed_x = -(object_1.speed_x);
    } else if (object_1.pos_x >= size_enclosure){
        object_1.pos_x = size_enclosure;
        object_1.speed_x = -(object_1.speed_x);
    }

    /* Checks posición y */
    if (object_1.pos_y <= 0){
        object_1.pos_y = 0;
        object_1.speed_y = -(object_1.speed_y);
    } else if (object_1.pos_y >= size_enclosure){
        object_1.pos_y = size_enclosure;
        object_1.speed_y = -(object_1.speed_y);
    }
    
    /* Checks posición z */
    if (object_1.pos_z <= 0){
        object_1.pos_z = 0;
        object_1.speed_z = -(object_1.speed_z);
    } else if (object_1.pos_z >= size_enclosure){
        object_1.pos_z = size_enclosure;
        object_1.speed_z = -(object_1.speed_z);
    }
    
}