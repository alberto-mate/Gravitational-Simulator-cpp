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

/* ESTRUCTURAS */
/* Estructura objeto */
struct object
{
    double *pos_x;
    double *pos_y;
    double *pos_z;
    double *speed_x;
    double *speed_y;
    double *speed_z;
    double *mass;
};

/* Estructura vector_elem */
struct vector_elem
{
    double x;
    double y;
    double z;
};

/* DECLARACIÓN PREVIA DE FUNCIONES */
double euclidean_norm(object objects, int index_1, int index_2);
void vector_gravitational_force(object objects, int index_1, int index_2, double* forces);
void calc_gravitational(int num_objects, int k, object objects, double* forces);
void vector_acceleration(object objects, int i, double* forces, vector_elem* acceleration);
void vector_speed(object *objects, int i, vector_elem *acceleration, float time_step);
void vector_position(object *objects, int i, float time_step);
void check_border(object *objects, int i, float size_enclosure);
bool check_collision(object objects, int i, int j);

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

    /* SOA - Structure of Arrays */
    object objects; 
    objects.mass = (double*) malloc(sizeof(double)*num_objects);
    objects.pos_x = (double*) malloc(sizeof(double)*num_objects);
    objects.pos_y = (double*) malloc(sizeof(double)*num_objects);
    objects.pos_z = (double*) malloc(sizeof(double)*num_objects);
    objects.speed_x = (double*) malloc(sizeof(double)*num_objects);
    objects.speed_y = (double*) malloc(sizeof(double)*num_objects);
    objects.speed_z = (double*) malloc(sizeof(double)*num_objects);

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
        objects.speed_x[i] = 0; // Velocidad x, y, z
        objects.speed_y[i] = 0;
        objects.speed_z[i] = 0;
        objects.mass[i] = mass_dist(gen); // Masa

        // Ponemos la precisión a 3 decimales. Imprimimos el objeto
        file_init << fixed << setprecision(3) << objects.pos_x[i] << " " << objects.pos_y[i] << " " << objects.pos_z[i] << " " << objects.speed_x[i] << " " << objects.speed_y[i] << " " << objects.speed_z[i] << " " << objects.mass[i] << endl;
    }

    file_init.close(); // Cerramos el fichero "init_config.txt"

    /* TODO Comprobar que no hay colisiones antes de las iteraciones */

    /* Iteraciones */
    for (int iteration = 0; iteration < num_iterations; iteration++)
    {
        /* Bucle para obtener nuevas propiedades de los objetos en la iteración. Comprueba también que no se haya pasado de los límites.*/
        for (int i = 0; i < num_objects; i++){
            if (objects.mass[i]!=0.0){ // Solo entrarán en el condicional objetos que no se han eliminado

                // Cálculo de la fuerza gravitatoria
                double forces[3] = {0,0,0};
                calc_gravitational(num_objects, i, objects, forces);

                // Cálculo del vector aceleración
                vector_elem *acceleration = (vector_elem*) malloc(sizeof(vector_elem));
                vector_acceleration(objects, i, forces, acceleration);
                
                // Cálculo del vector velocidad
                vector_speed(&objects, i, acceleration, time_step);

                // Cálculo del vector posiciones
                vector_position(&objects, i, time_step);

                // Comprobar bordes
                check_border(&objects, i, size_enclosure);
            }
        } 
        //cout<<"Nuevas posiciones calculadas \n";

        /* Bucle anidado para comprobar colisiones entre objetos */
        int num_objects_after_delete = num_objects;
        for (int i = 0; i < num_objects; i++){
            for (int j = 0; j < num_objects; j++){ // TODO AÑADIR OPTI num_objetos -i-1 objects[j+i+1]
                // Se resta i y 1 para evitar comprobaciones dobles
                // Comprobar colisiones
                if (i!=j && objects.mass[i]!=0.0 && objects.mass[i]!=0.0){ // Colision entre objetos diferentes que no hayan sido eliminados con anterioridad
                    if (check_collision(objects, i, j)){ // Comprobar colisión

                        // Actualización de la masa y velocidades del nuevo objeto resultante
                        objects.mass[i] += objects.mass[j];
                        objects.speed_x[i] += objects.speed_x[j];
                        objects.speed_y[i] += objects.speed_y[j];
                        objects.speed_z[i] += objects.speed_z[j];

                        // Eliminación (masa = 0.0) del segundo objeto que ha colisionado
                        num_objects_after_delete -= 1;
                        objects.mass[j] = 0.0;
                    }
                }
            }
        }
        cout<<num_objects_after_delete<<"\n";

        // Actualizar array objects debido a colisiones
        
        int counter = 0; // Contador que recorrerá los objetos de nuestro array hasta borrar los colisionados
        while (counter<num_objects_after_delete){
            if (objects.mass[counter]==0.0){ // Si el objeto tiene masa=0.0 -> ha sido eliminado y hay que actualizar el array
                for (int i = 0; i<num_objects-counter-1; i++){
                    objects.pos_x[i+counter]=objects.pos_x[i+counter+1]; // Mueve el resto de objetos para sustituir el eliminado
                    objects.pos_y[i+counter]=objects.pos_y[i+counter+1];
                    objects.pos_z[i+counter]=objects.pos_z[i+counter+1];
                    objects.speed_x[i+counter]=objects.speed_x[i+counter+1];
                    objects.speed_y[i+counter]=objects.speed_y[i+counter+1];
                    objects.speed_z[i+counter]=objects.speed_z[i+counter+1];
                    objects.mass[i+counter]=objects.mass[i+counter+1];
                }

            }else{ // Si el objeto tiene masa!=0.0 -> no hay que actualizar el array
                counter++; // Aumentamos el contador
            }
        }
        num_objects = num_objects_after_delete;

        cout<<"Fin iteración: "<<iteration<<" Num objetos:"<<num_objects<<"\n";
    }

}

/* FUNCIONES */
/* Distancia euclídea entre dos objetos */
double euclidean_norm(object objects, int index_1, int index_2)
{   
    return std::sqrt(pow(objects.pos_x[index_1] - objects.pos_x[index_2], 2) + pow(objects.pos_y[index_1] - objects.pos_y[index_2], 2) + pow(objects.pos_z[index_1] - objects.pos_z[index_2], 2));
}

/* Fuerza gravitatoria entre dos objetos */
void vector_gravitational_force(object objects, int index_1, int index_2, double* forces)
{
    double force_module = (GRAVITY_CONST * objects.mass[index_1] * objects.mass[index_2]) / (pow(euclidean_norm(objects, index_1, index_2), 3));
    forces[0] += force_module * (objects.pos_x[index_1] - objects.pos_x[index_2]);
    forces[1] += force_module * (objects.pos_y[index_1] - objects.pos_y[index_2]);
    forces[2] += force_module * (objects.pos_z[index_1] - objects.pos_z[index_2]);
}

/* Fuerza gravitatoria que ejerce un objeto */
void calc_gravitational(int num_objects, int k, object objects, double* forces){
    for (int j=0; j<num_objects; j++){
        if (j!=k){
            vector_gravitational_force(objects, j,  k , forces);
        }
    }
}

/* Vector aceleración */
void vector_acceleration(object objects, int i, double* forces, vector_elem* acceleration)
{
    /* Cálculo del vector aceleración */
    acceleration->x= forces[0]/objects.mass[i];
    acceleration->y= forces[1]/objects.mass[i];
    acceleration->z= forces[2]/objects.mass[i];
}

/* Vector velocidad */
void vector_speed(object *objects, int i, vector_elem *acceleration, float time_step)
{      
    /* Cálculo del vector velocidad */
    objects->speed_x[i] = objects->speed_x[i] + (acceleration->x * time_step);
    objects->speed_y[i] = objects->speed_y[i] + (acceleration->y * time_step);
    objects->speed_z[i] = objects->speed_z[i] + (acceleration->z * time_step);
}

/* Vector de posicion */
void vector_position(object *objects, int i, float time_step)
{
    /* Cálculo del vector posición */
    objects->pos_x[i] = objects->pos_x[i] + (objects->speed_x[i] * time_step);
    objects->pos_y[i] = objects->pos_y[i] + (objects->speed_y[i] * time_step);
    objects->pos_z[i] = objects->pos_z[i] + (objects->speed_z[i] * time_step);

}

/* Función para recolocar al objeto si traspasa los límites */
void check_border(object *objects, int i, float size_enclosure)
{
    /* Checks posición x */
    if (objects->pos_x[i] <= 0)
    {
        objects->pos_x[i] = 0;
        objects->speed_x[i] = -1 * objects->speed_x[i];
    }
    else if (objects->pos_x[i] >= size_enclosure)
    {
        objects->pos_x[i] = size_enclosure;
        objects->speed_x[i] = -1 * (objects->speed_x[i]);
    }

    /* Checks posición y */
    if (objects->pos_y[i] <= 0)
    {
        objects->pos_y[i] = 0;
        objects->speed_y[i] = -1 * (objects->speed_y)[i];
    }
    else if (objects->pos_y[i] >= size_enclosure)
    {
        objects->pos_y[i] = size_enclosure;
        objects->speed_y[i] = -1 * (objects->speed_y)[i];
    }

    /* Checks posición z */
    if (objects->pos_z[i] <= 0)
    {
        objects->pos_z[i] = 0;
        objects->speed_z[i] = -1 * objects->speed_z[i];
    }
    else if (objects->pos_z[i] >= size_enclosure)
    {
        objects->pos_z[i] = size_enclosure;
        objects->speed_z[i] =  -1 * objects->speed_z[i];
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