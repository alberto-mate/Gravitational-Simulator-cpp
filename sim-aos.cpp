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

struct vector_elem
{
    double x;
    double y;
    double z;
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
    uniform_real_distribution<double> position_dist(0.0, nextafter(size_enclosure, numeric_limits<double>::max()));
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

    /* Declaración matriz de fuerzas */
    vector_elem *force_matrix[num_objects - 1];
    for (int i = 0; i < (num_objects - 1); i++)
    {
        vector_elem force_vector[num_objects - 1 - i];
        force_matrix[i] = force_vector;
    }

    /* Comprobar que no hay colisiones antes de las iteraciones */

    /* Iteraciones */
    for (int iteration = 0; iteration < num_iterations; iteration++)
    {
        /* Bucle rellenar matriz de fuerzas*/

        for (int i = 0; i < sizeof(force_matrix); i++)
        {
            for (int j = 0; j < sizeof(force_matrix[i]); j++)
            {
                force_matrix[i][j] = vector_gravitational_force(objeto[i], objeto[i + j + 1]);
            }
        }

        // PYTHON
        // for i in range(len(force_matrix)):
        //     for j in range(len(force_matrix[i])):
        //         force_matrix[i][j] = vector_gravitational_force(objeto[i],objeto[i+j+1])

        // Calcular todas las aceleraciones
        static vector_elem acceleration[num_objects];
        acceleration = all_vector_acceleration(objects, force_matrix, num_objects)

            /* Bucle obtener aceleración, velocidad y posición de los objetos en la iteración. Comprueba también que no se haya pasado de los límites.*/
            for (int i = 0; i < num_objects; i++)
        {

            // Vector velocidades
            vector_speed(objects[i], time_step);
            // Vector posiciones
            vector_position(objects[i], time_step);
            // Comprobar bordes
            check_border(object[i], size_enclosure);
        }

        /* Bucle anidado para comprobar colisiones entre objetos */
        bool collision = false for (int i = 0; i < num_objects; i++)
        {
            for (int j = 0; j < num_objects - i - 1; j++)
            { // Se resta i y 1 para evitar comprobaciones dobles
                // Comprobar colisiones
                if (check_collision(objects[i], objects[j])
                {
                    collision = true num_objects -= 1
                }
            }
        }
        if (collision == true)
        {

            // Borrar la matriz de fuerzas
            // Volver a crearla
            delete force_matrix;
            vector_elem *force_matrix[num_objects - 1];
            for (int i = 0; i < (num_objects - 1); i++)
            {
                vector_elem force_vector[num_objects - 1 - i];
                force_matrix[i] = force_vector;
            }
        }
    }
}

/* FUNCIONES */
/* Distancia euclídea entre dos objetos */
double euclidean_norm(object object_1, object object_2)
{
    return sqrt(pow(object_1.pos_x - object_2.pos_x, 2) + pow(object_1.pos_y - object_2.pos_y, 2) + pow(object_1.pos_z - object_2.pos_z, 2));
}

/* Fuerza gravitatoria entre dos objetos */
vector_elem *vector_gravitational_force(object object_1, object object_2)
{
    double force_module = (GRAVITY_CONST * object_1.mass * object_2.mass) / (pow(euclidean_norm(object_1, object_2), 3));
    vector_elem result;
    result.x = force_module * (object_1.pos_x - object_2.pos_x);
    result.y = force_module * (object_1.pos_y - object_2.pos_y);
    result.z = force_module * (object_1.pos_z - object_2.pos_z);

    return result;
}

/* Vector de aceleración */
vector_elem *all_vector_acceleration(*object objects, *double force_matrix, int num_objects)
{
    /* PYTHON (calculate all forces)
    matriz = [[1,2,3], [4,5], [6]]
    objeto = [0,0,0,0]
    for i in range(len(matriz)):
        for j in range(len(matriz[i])):
            objeto[i] += matriz[i][j]
            objeto[i+j+1] -= matriz[i][j]
    */

    static vector_elem new_acceleration[num_objects];
    for (int i = 0; i < num_objects; i++)
    {
        for (int j = 0; j < num_objects - i - 1; j++)
        {
            new_acceleration[i] = operate_forces(new_acceleration[i], force_matrix[i][j], 1);
            new_acceleration[i + j + 1] = operate_forces(new_acceleration[i + j + 1], force_matrix[i][j], -1);
        }
    }

    for (int i = 0; i < num_objects; i++)
    {
        new_acceleration[i] = divide_mass(new_acceleration[i], objects[i].mass);
    }
    return new_acceleration;
}

vector_elem *operate_forces(vector_elem force_1, vector_elem force_2, int type)
{
    // type = 1 suma y type = -1 resta
    vector_elem new_force;
    new_force.x = force_1.x + (type * force_2.x);
    new_force.y = force_1.y + (type * force_2.y);
    new_force.z = force_1.z + (type * force_2.z);
    return new_force;
}

vector_elem *divide_mass(vector_elem force, double mass)
{
    vector_elem new_acceleration;
    new_acceleration.x = force.x / mass;
    new_acceleration.y = force.y / mass;
    new_acceleration.z = force.z / mass;
    return new_acceleration;
}

/* Vector velocidad */
double *vector_speed(object object_1, int time_step)
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
double *vector_position(object object_1, int time_step)
{
    /* Calculamos la velocidad del objeto para obtener la posición */
    static double new_speed[3] = vector_speed(object_1);

    /* Cálculo del vector posición */
    static double new_position[3] = {object_1.pos_x + (new_speed[0] * time_step), object_1.pos_y + (new_speed[1] * time_step), object_1.pos_z + (new_speed[2] * time_step)};

    /* Devolvemos el vector posición */
    return new_position;
}

/* Función para recolocar al objeto si traspasa los límites */
void check_border(object object_1, int size_enclosure)
{

    /* Checks posición x */
    if (object_1.pos_x <= 0)
    {
        object_1.pos_x = 0;
        object_1.speed_x = -(object_1.speed_x);
    }
    else if (object_1.pos_x >= size_enclosure)
    {
        object_1.pos_x = size_enclosure;
        object_1.speed_x = -(object_1.speed_x);
    }

    /* Checks posición y */
    if (object_1.pos_y <= 0)
    {
        object_1.pos_y = 0;
        object_1.speed_y = -(object_1.speed_y);
    }
    else if (object_1.pos_y >= size_enclosure)
    {
        object_1.pos_y = size_enclosure;
        object_1.speed_y = -(object_1.speed_y);
    }

    /* Checks posición z */
    if (object_1.pos_z <= 0)
    {
        object_1.pos_z = 0;
        object_1.speed_z = -(object_1.speed_z);
    }
    else if (object_1.pos_z >= size_enclosure)
    {
        object_1.pos_z = size_enclosure;
        object_1.speed_z = -(object_1.speed_z);
    }
}

bool check_collision(object object_1, object object_2)
{

    /* Comprobar si colisionan (distancia euclídea entre 1 y 2 menor que 1) */
    if (euclidean_norm(object_1, object_2) < 1)
    {
        /* Reemplazar objeto1 por uno nuevo de suma de masas y velocidades con la misma posición que el objeto1*/
        object_1.mass += object_2.mass;
        object_1.speed_x += object_2.speed_x;
        object_1.speed_y += object_2.speed_y;
        object_1.speed_z += object_2.speed_z;

        /* Eliminar objeto2 */
        objects.remove(objects.begin(), objects.end(), object_2);
        /*Devulve si hay colision devulve una confirmación*/
        return True
    }
    return False
}