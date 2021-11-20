/* Librerias */
#include <iostream>
#include <math.h>
#include <fstream>
#include <random>
#include <vector>
#include <iomanip>
#include <chrono>
#include <omp.h>

using namespace std;

/* CONSTANTES */
const double GRAVITY_CONST = 6.674 * pow(10, -11); // Constante gravedad universal
const double M = pow(10, 21);                      // Media (distribución normal)
const double SDM = pow(10, 15);                    // Desviación (distribución normal)

/* ESTRUCTURAS */
/* Estructura objeto */
struct object {
    vector<double> pos_x;
    vector<double> pos_y;
    vector<double> pos_z;
    vector<double> speed_x;
    vector<double> speed_y;
    vector<double> speed_z;
    vector<double> mass;
};

/* Estructura vector_elem */
struct vector_elem{
    double x;
    double y;
    double z;
};

/* DECLARACIÓN PREVIA DE FUNCIONES */
double euclidean_norm(object objects, int index_1, int index_2);
void vector_gravitational_force(object objects, int index_1, int index_2, double* forces);
void calc_gravitational(int num_objects, int k, object objects, double* forces);
void vector_acceleration(object objects, int i, double* forces, vector_elem* acceleration);
void vector_speed(object *objects, int i, vector_elem *acceleration, double time_step);
void vector_position(object *objects, int i, double time_step);
void check_border(object *objects, int i, double size_enclosure);
bool check_collision(object objects, int i, int j);

/* MAIN */
int main(int argc, char const *argv[])
{
    /* Medición del tiempo */
    double t1 = omp_get_wtime();


    /* Comprobación inicial argumentos */
    if (argc != 6)
    {
        cerr << "Número de argumentos incorrecto\n";
        return -1;
    }

    /* Comprobación de valores iniciales de argumentos */
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

    /* SOA - Structure of Arrays */
    object objects; 
    objects.mass.resize(num_objects);
    objects.pos_x.resize(num_objects);
    objects.pos_y.resize(num_objects);
    objects.pos_z.resize(num_objects);
    objects.speed_x.resize(num_objects);
    objects.speed_y.resize(num_objects);
    objects.speed_z.resize(num_objects);

    /* Coordenadas y masas pseudoaleatorias */
    mt19937_64 gen(random_seed);
    uniform_real_distribution<double> position_dist(0.0, nextafter(size_enclosure, numeric_limits<double>::max()));
    normal_distribution<double> mass_dist(M, SDM);

    /* Fichero de configuracion inicial */
    ofstream file_init;
    file_init.open("init_config.txt");
    file_init << fixed << setprecision(3) << size_enclosure << " " << time_step << " " << num_objects << endl;

    /* Creación de objetos */
    for (int i = 0; i < num_objects; i++)
    {
        objects.pos_x[i] = position_dist(gen); // Posicion x, y, z
        objects.pos_y[i] = position_dist(gen);
        objects.pos_z[i] = position_dist(gen);
        objects.mass[i] = mass_dist(gen); // Masa

        // Ponemos la precisión a 3 decimales. Imprimimos el objeto
        file_init << fixed << setprecision(3) << objects.pos_x[i] << " " << objects.pos_y[i] << " " << objects.pos_z[i] << " " << objects.speed_x[i] << " " << objects.speed_y[i] << " " << objects.speed_z[i] << " " << objects.mass[i] << endl;
    }

    file_init.close(); // Cerramos el fichero "init_config.txt"


    /* Bucle anidado para comprobar colisiones entre objetos previas a las iteraciones */
    for (long unsigned int i = 0; i < objects.mass.size(); i++)
    {
        for (long unsigned int j = i + 1; j < objects.mass.size(); j++)
        {
            // Comprobar colisiones
            // Colision entre objetos diferentes que no hayan sido eliminados con anterioridad
            if (check_collision(objects, i, j))
            {   // Comprobar colisión
                // Actualización de la masa y velocidades del primer objeto que colisiona generando uno nuevo
                // Actualización de la masa y velocidades del primer objeto que colisiona generando uno nuevo
                objects.mass[i] += objects.mass[j];
                objects.speed_x[i] += objects.speed_x[j];
                objects.speed_y[i] += objects.speed_y[j];
                objects.speed_z[i] += objects.speed_z[j];
                
                // Lo borramos de los vectores
                objects.pos_x.erase(objects.pos_x.begin() + j);
                objects.pos_y.erase(objects.pos_y.begin() + j);
                objects.pos_z.erase(objects.pos_z.begin() + j);
                objects.speed_x.erase(objects.speed_x.begin() + j);
                objects.speed_y.erase(objects.speed_y.begin() + j);
                objects.speed_z.erase(objects.speed_z.begin() + j);
                objects.mass.erase(objects.mass.begin() + j);                    
                j--;
            }
        }
    }

    // Actualizamos el número de objetos en el vector
    num_objects = objects.mass.size();

    /* Iteraciones */
    for (int iteration = 0; iteration < num_iterations; iteration++)
    {
        /* Bucle para obtener nuevas propiedades de los objetos en la iteración (fuerzas, aceleración y velocidad) */
        for (int i = 0; i < num_objects; i++)
        {
            // Solo entrarán en el condicional objetos que no se han eliminado
            // Cálculo de la fuerza gravitatoria
            double forces[3] = {0.0, 0.0, 0.0};
            calc_gravitational(num_objects, i, objects, forces);
            // Cálculo del vector aceleración
            vector_elem *acceleration = (vector_elem *)malloc(sizeof(vector_elem));
            vector_acceleration(objects, i, forces, acceleration);
            //  Cálculo del vector velocidad
            vector_speed(&objects, i, acceleration, time_step);
        }


        /* Bucle para calcular posiciones y comprobar bordes */
        for (int i = 0; i < num_objects; i++)
        {
            // Cálculo del vector posiciones
            vector_position(&objects, i, time_step);
            //  Comprobar bordes
            check_border(&objects, i, size_enclosure);
        }

        /* Bucle anidado para comprobar colisiones entre objetos */
        for (long unsigned int i = 0; i < objects.mass.size(); i++)
        {
            for (long unsigned int j = i + 1; j < objects.mass.size(); j++)
            {
                // Comprobar colisiones
                // Colision entre objetos diferentes que no hayan sido eliminados con anterioridad
                if (check_collision(objects, i, j))
                {   // Comprobar colisión
                    // Actualización de la masa y velocidades del primer objeto que colisiona generando uno nuevo
                    objects.mass[i] += objects.mass[j];
                    objects.speed_x[i] += objects.speed_x[j];
                    objects.speed_y[i] += objects.speed_y[j];
                    objects.speed_z[i] += objects.speed_z[j];

                    // Lo borramos de los vectores
                    objects.pos_x.erase(objects.pos_x.begin() + j);
                    objects.pos_y.erase(objects.pos_y.begin() + j);
                    objects.pos_z.erase(objects.pos_z.begin() + j);
                    objects.speed_x.erase(objects.speed_x.begin() + j);
                    objects.speed_y.erase(objects.speed_y.begin() + j);
                    objects.speed_z.erase(objects.speed_z.begin() + j);
                    objects.mass.erase(objects.mass.begin() + j);                    
                    j--;
                }
            }
        }

        // Actualizamos el número de objetos en el vector
        num_objects = objects.mass.size();
    }

    /* Escribimos en el archivo "final_config.txt" los parámetros finales */
    ofstream file_final;
    file_final.open("final_config.txt");
    file_final << fixed << setprecision(3) << size_enclosure << " " << time_step << " " << num_objects << endl;

    for (int i = 0; i < num_objects; i++)
    {
        file_final << fixed << setprecision(3) << objects.pos_x[i] << " " << objects.pos_y[i] << " " << objects.pos_z[i] << " " << objects.speed_x[i] << " " << objects.speed_y[i] << " " << objects.speed_z[i] << " " << objects.mass[i] << endl;
    }

    file_init.close(); // Cerramos el fichero "final_config.txt"

    /* Medición del tiempo */
    double t2 = omp_get_wtime();
    double diff = t2 - t1;

}
/* FUNCIONES */
/* Distancia euclídea entre dos objetos */
double euclidean_norm(object objects, int i, int j)
{
    return std::sqrt((objects.pos_x[i]- objects.pos_x[j]) * (objects.pos_x[i]- objects.pos_x[j]) + (objects.pos_y[i]- objects.pos_y[j]) * (objects.pos_y[i]- objects.pos_y[j]) + (objects.pos_z[i]- objects.pos_z[j]) * (objects.pos_z[i]- objects.pos_z[j]));
}

/* Fuerza gravitatoria entre dos objetos */
void vector_gravitational_force(object objects, int i, int j, double *forces)
{
    double dist = euclidean_norm(objects, i, j);
    double Fg = GRAVITY_CONST * objects.mass[i] * objects.mass[j]/ (dist*dist*dist);

    /*
    #pragma omp parallel
    {
        #pragma omp sections
        {
            #pragma omp section
            forces[0] += (Fg * (objects.pos_x[i] - objects.pos_x[j]));
            #pragma omp section
            forces[1] += (Fg * (objects.pos_y[i] - objects.pos_y[j]));
            #pragma omp section
            forces[2] += (Fg * (objects.pos_z[i] - objects.pos_z[j]));

        }
    }
    */
    
    forces[0] += (Fg * (objects.pos_x[i] - objects.pos_x[j]));
    forces[1] += (Fg * (objects.pos_y[i] - objects.pos_y[j]));
    forces[2] += (Fg * (objects.pos_z[i] - objects.pos_z[j]));
}

/* Fuerza gravitatoria que ejerce un objeto */
void calc_gravitational(int num_objects, int i, object objects, double *forces){
    for (int j = 0; j < num_objects; j++){
        if (j != i){
            vector_gravitational_force(objects, j, i, forces);
        }
    }
}

/* Vector aceleración */
void vector_acceleration(object objects, int i, double *forces, vector_elem *acceleration)
{
    /* Cálculo del vector aceleración */
    acceleration->x = forces[0] / objects.mass[i];
    acceleration->y = forces[1] / objects.mass[i];
    acceleration->z = forces[2] / objects.mass[i];
}

/* Vector velocidad */
void vector_speed(object *objects, int i, vector_elem *acceleration, double time_step)
{
    /* Cálculo del vector velocidad */
    objects->speed_x[i] += (acceleration->x * time_step);
    objects->speed_z[i] += (acceleration->z * time_step);
    objects->speed_y[i] += (acceleration->y * time_step);
}

/* Vector de posicion */
void vector_position(object *objects, int i, double time_step)
{
    /* Cálculo del vector posición */
    objects->pos_x[i] += (objects->speed_x[i] * time_step);
    objects->pos_y[i] += (objects->speed_y[i] * time_step);
    objects->pos_z[i] += (objects->speed_z[i] * time_step);
}

/* Función para recolocar al objeto si traspasa los límites */
void check_border(object *objects, int i, double size_enclosure)
{

    /*#pragma omp parallel
    {
        #pragma omp sections
        {
            #pragma omp section
            // Checks posición x
            if (objects->pos_x[i] <= 0)
            {
                objects->pos_x[i] = 0;
                objects->speed_x[i] = -1 * (objects->speed_x[i]);
            }
            else if (objects->pos_x[i] >= size_enclosure)
            {
                objects->pos_x[i] = size_enclosure;
                objects->speed_x[i] = -1 * (objects->speed_x[i]);
            }

            #pragma omp section
            // Checks posición y
            if (objects->pos_y[i] <= 0)
            {
                objects->pos_y[i] = 0;
                objects->speed_y[i] = -1 * (objects->speed_y[i]);
            }
            else if (objects->pos_y[i] >= size_enclosure)
            {
                objects->pos_y[i] = size_enclosure;
                objects->speed_y[i] = -1 * (objects->speed_y[i]);
            }

            #pragma omp section
            // Checks posición z
            if (objects->pos_z[i] <= 0)
            {
                objects->pos_z[i] = 0;
                objects->speed_z[i] = -1 * (objects->speed_z[i]);
            }
            else if (objects->pos_z[i] >= size_enclosure)
            {
                objects->pos_z[i] = size_enclosure;
                objects->speed_z[i] = -1 * (objects->speed_z[i]);
            }    
        }
    }
    */

    // Checks posición x
    if (objects->pos_x[i] <= 0)
    {
        objects->pos_x[i] = 0;
        objects->speed_x[i] = -1 * (objects->speed_x[i]);
    }
    else if (objects->pos_x[i] >= size_enclosure)
    {
        objects->pos_x[i] = size_enclosure;
        objects->speed_x[i] = -1 * (objects->speed_x[i]);
    }

    // Checks posición y
    if (objects->pos_y[i] <= 0)
    {
        objects->pos_y[i] = 0;
        objects->speed_y[i] = -1 * (objects->speed_y[i]);
    }
    else if (objects->pos_y[i] >= size_enclosure)
    {
        objects->pos_y[i] = size_enclosure;
        objects->speed_y[i] = -1 * (objects->speed_y[i]);
    }

    // Checks posición z
    if (objects->pos_z[i] <= 0)
    {
        objects->pos_z[i] = 0;
        objects->speed_z[i] = -1 * (objects->speed_z[i]);
    }
    else if (objects->pos_z[i] >= size_enclosure)
    {
        objects->pos_z[i] = size_enclosure;
        objects->speed_z[i] = -1 * (objects->speed_z[i]);
    }
}

/* Comprobar colisión entre dos objetos (distancia euclídea entre objetos menor que 1) */
bool check_collision(object objects, int i, int j)
{
    if (euclidean_norm(objects, i, j) < 1)
    {
        return true;
    }
    return false;
}