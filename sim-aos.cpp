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
const double GRAVITY_CONST = 6.674 * 1E-11; // Constante gravedad universal
const double M = 1e+21;                     // Media (distribución normal)
const double SDM = 1e+15;                   // Desviación (distribución normal)

/* ESTRUCTURAS */
/* Estructura objeto */
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

/* Estructura vector_elem */
struct vector_elem
{
    double x;
    double y;
    double z;
};

/* DECLARACIÓN PREVIA DE FUNCIONES */
void vector_gravitational_force(object object_1, object object_2, double *forces);
void calc_gravitational(int num_objects, int index_1, vector<object> &objects, double *forces);
void vector_acceleration(object object_1, double *forces, vector_elem *acceleration);
void vector_speed(object *object_1, vector_elem *acceleration, double time_step);
void vector_position(object *object_1, double time_step);
void check_border(object *object_1, double size_enclosure);
bool check_collision(object object_1, object object_2);

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
    int num_objects = atoi(argv[1]);       // Número de objetos a simular (>0 entero)
    int num_iterations = atoi(argv[2]);    // Número de iteraciones a simular (>0 entero)
    int random_seed = atoi(argv[3]);       // Semilla para distribuciones aleatorias
    double size_enclosure = stod(argv[4]); // Tamaño del recinto (>0 real)
    double time_step = stod(argv[5]);      // Incremento de tiempo en cada iteración (>0 real)

    /* Impresión por pantalla de las variables */
    cout << "num_objects: " << num_objects << "\n";
    cout << "num_iterations: " << num_iterations << "\n";
    cout << "random_seed: " << random_seed << "\n";
    cout << "size_enclosure: " << size_enclosure << "\n";
    cout << "time_step: " << time_step << "\n";

    /* Coordenadas aleatorias */
    mt19937_64 gen(random_seed);
    uniform_real_distribution<> position_dist(0.0, size_enclosure);
    normal_distribution<> mass_dist{1E21, 1E15};

    /* AOS - Array of Structs */
    // object *objects = (object*) malloc(sizeof(object)*num_objects);
    vector<object> objects(num_objects);

    /* Fichero de configuracion inicial */
    ofstream file_init;
    file_init.open("init_config.txt");
    file_init << fixed << setprecision(3) << size_enclosure << " " << time_step << " " << num_objects << endl;

    /* Creación de objetos */
    for (int i = 0; i < num_objects; i++)
    {
        objects[i].pos_x = position_dist(gen); // Posicion x, y, z
        objects[i].pos_y = position_dist(gen);
        objects[i].pos_z = position_dist(gen);
        objects[i].speed_x = 0; // Velocidad x, y, z
        objects[i].speed_y = 0;
        objects[i].speed_z = 0;
        objects[i].mass = mass_dist(gen); // Masa
        // cout<<(objects[i].mass)<<"\n";

        // Ponemos la precisión a 3 decimales. Imprimimos el objeto
        file_init << fixed << setprecision(3) << objects[i].pos_x << " " << objects[i].pos_y << " " << objects[i].pos_z << " " << objects[i].speed_x << " " << objects[i].speed_y << " " << objects[i].speed_z << " " << objects[i].mass << endl;
    }

    file_init.close(); // Cerramos el fichero "init_config.txt"

    /* TODO Comprobar que no hay colisiones antes de las iteraciones */

    /* Iteraciones */
    for (int iteration = 0; iteration < num_iterations; iteration++)
    {
        cout << "\nIteracion " << iteration << "\n";
        /* Bucle para obtener nuevas propiedades de los objetos en la iteración. Comprueba también que no se haya pasado de los límites.*/
        for (int i = 0; i < num_objects; i++)
        {
            if (objects[i].mass != 0.0)
            { // Solo entrarán en el condicional objetos que no se han eliminado
                // Cálculo de la fuerza gravitatoria
                double forces[3] = {0.0, 0.0, 0.0};
                calc_gravitational(num_objects, i, objects, forces);
                // cout<< "Obj: " << i<<" posx: "<<objects[i].pos_x<<" posy: "<< objects[i].pos_y<<" posz: "<<objects[i].pos_z<<" speedx: "<<objects[i].speed_x<<" speedy: "<< objects[i].speed_y<<" speedz: "<<objects[i].speed_z<<" mass: "<<objects[i].mass<<"\n";
                cout << "Forces " << i << " ax: " << forces[0] << " ay: " << forces[1] << " az: " << forces[2] << "\n";
                // Cálculo del vector aceleración
                vector_elem *acceleration = (vector_elem *)malloc(sizeof(vector_elem));
                vector_acceleration(objects[i], forces, acceleration);
                // cout<< "Obj: " << i<<" posx: "<<objects[i].pos_x<<" posy: "<< objects[i].pos_y<<" posz: "<<objects[i].pos_z<<" speedx: "<<objects[i].speed_x<<" speedy: "<< objects[i].speed_y<<" speedz: "<<objects[i].speed_z<<"\n";
                //  Cálculo del vector velocidad
                vector_speed(&objects[i], acceleration, time_step);
                // cout<< "Obj: " << i<<" posx: "<<objects[i].pos_x<<" posy: "<< objects[i].pos_y<<" posz: "<<objects[i].pos_z<<" speedx: "<<objects[i].speed_x<<" speedy: "<< objects[i].speed_y<<" speedz: "<<objects[i].speed_z<<"\n";
            }
        }

        for (int i = 0; i < num_objects; i++)
        {
            // Cálculo del vector posiciones
            vector_position(&objects[i], time_step);
            // cout<< "Obj: " << i<<" posx: "<<objects[i].pos_x<<" posy: "<< objects[i].pos_y<<" posz: "<<objects[i].pos_z<<" speedx: "<<objects[i].speed_x<<" speedy: "<< objects[i].speed_y<<" speedz: "<<objects[i].speed_z<<"\n";
            //  Comprobar bordes
            check_border(&objects[i], size_enclosure);
            // cout<< "Obj: " << i<<" posx: "<<objects[i].pos_x<<" posy: "<< objects[i].pos_y<<" posz: "<<objects[i].pos_z<<" speedx: "<<objects[i].speed_x<<" speedy: "<< objects[i].speed_y<<" speedz: "<<objects[i].speed_z<<"\n";
        }
        // cout<<"Nuevas posiciones calculadas \n";
        int deletes[num_objects];
        int num_objects_after_delete = num_objects;
        /* Bucle anidado para comprobar colisiones entre objetos */
        for (int i = 0; i < num_objects; i++)
        {
            for (int j = i + 1; j < num_objects; j++)
            { // TODO AÑADIR OPTI num_objetos -i-1 objects[j+i+1]
                // Se resta i y 1 para evitar comprobaciones dobles
                // Comprobar colisiones
                if (i != j)
                { // Colision entre objetos diferentes que no hayan sido eliminados con anterioridad
                    if (check_collision(objects[i], objects[j]))
                    { // Comprobar colisión
                        if (deletes[j] != 1 && deletes[i] != 1)
                        {
                            // Actualización de la masa y velocidades del primer objeto que colisiona generando uno nuevo
                            objects[i].mass += objects[j].mass;
                            objects[i].speed_x += objects[j].speed_x;
                            objects[i].speed_y += objects[j].speed_y;
                            objects[i].speed_z += objects[j].speed_z;

                            deletes[j] = 1;
                        }
                    }
                }
            }
        }
        int contador = 0;
        while (contador < num_objects)
        {
            if (deletes[contador] == 1)
            {
                objects.erase(objects.begin() + contador);
                num_objects -= 1;
            }
            contador++;
        }

        cout << "Fin iteración: " << iteration << " Num objetos:" << num_objects << "\n";

        // Eliminación del segundo objeto que colisiona del sistema

        /*
        // Actualizar array objects debido a colisiones

        int counter = 0; // Contador que recorrerá los objetos de nuestro array hasta borrar los colisionados
        while (counter<num_objects_after_delete){
            if (objects[counter].mass==0.0){ // Si el objeto tiene masa=0.0 -> ha sido eliminado y hay que actualizar el array
                for (int i = 0; i<num_objects-counter-1; i++){
                    objects[i+counter]=objects[i+counter+1]; // Mueve el resto de objetos para sustituir el eliminado
                    // num_objects = 4; counter = 0 i<3
                    // i=0 objects[0]=objects[1]
                    // i=1 objects[1]=objects[2]
                    // i=2 objects[2]=objects[3]

                    // num_objects = 4; counter = 1 i<2
                    // i=0 objects[1]=objects[2]
                    // i=1 objects[2]=objects[3]

                    // num_objects = 4; counter = 2 i<1
                    // i=0 objects[2]=objects[3]

                    // num_objects = 4; counter = 3 i<0
                }

            }else{ // Si el objeto tiene masa!=0.0 -> no hay que actualizar el array
                counter++; // Aumentamos el contador
            }
        }
        num_objects = num_objects_after_delete;
        */
    }
    ofstream file_final;
    file_final.open("final_config.txt");
    file_final << fixed << setprecision(3) << size_enclosure << " " << time_step << " " << num_objects << endl;

    for (int i = 0; i < num_objects; i++)
    {
        file_final << fixed << setprecision(3) << objects[i].pos_x << " " << objects[i].pos_y << " " << objects[i].pos_z << " " << objects[i].speed_x << " " << objects[i].speed_y << " " << objects[i].speed_z << " " << objects[i].mass << endl;
    }

    file_init.close(); // Cerramos el fichero "final_config.txt"
}

/* FUNCIONES */
/* Distancia euclídea entre dos objetos */
double euclidean_norm(object object_1, object object_2)
{
    // cout<<"Pos 2 "<<" px: "<<(object_1.pos_x)<<" "<< (object_2.pos_x)<<" "<<(object_1.pos_x - object_2.pos_x)<<" py: "<<(object_1.pos_y)<<" "<< (object_2.pos_y)<<" "<<(object_1.pos_y - object_2.pos_y)<<" pz: "<<(object_1.pos_z)<<" "<< (object_2.pos_z)<<" "<<(object_1.pos_z - object_2.pos_z)<<"\n";
    return std::sqrt((object_1.pos_x - object_2.pos_x) * (object_1.pos_x - object_2.pos_x) + (object_1.pos_y - object_2.pos_y) * (object_1.pos_y - object_2.pos_y) + (object_1.pos_z - object_2.pos_z) * (object_1.pos_z - object_2.pos_z));
}

/* Fuerza gravitatoria entre dos objetos */
void vector_gravitational_force(object object_1, object object_2, double *forces)
{
    // cout<<"Pos 1 "<<" px: "<<(object_1.pos_x)<<" "<< (object_2.pos_x)<<" "<<(object_1.pos_x - object_2.pos_x)<<" py: "<<(object_1.pos_y)<<" "<< (object_2.pos_y)<<" "<<(object_1.pos_y - object_2.pos_y)<<" pz: "<<(object_1.pos_z)<<" "<< (object_2.pos_z)<<" "<<(object_1.pos_z - object_2.pos_z)<<"\n";
    // cout<<"Forces 1 "<<" ax: "<<forces[0]<<" ay: "<<forces[1]<<" az: "<<forces[2]<<"\n";
    double x = (euclidean_norm(object_1, object_2));
    if (x != 0)
    {
        double ttt = (GRAVITY_CONST * object_1.mass * object_2.mass * (object_1.pos_x - object_2.pos_x)) / (x * x * x);
        //cout << "Vector x de la fuerza: " << ttt << endl;
        //cout << "Euclidean: " << x << endl;
        forces[0] += (GRAVITY_CONST * object_1.mass * object_2.mass * (object_1.pos_x - object_2.pos_x)) / (x * x * x);
        forces[1] += (GRAVITY_CONST * object_1.mass * object_2.mass * (object_1.pos_y - object_2.pos_y)) / (x * x * x);
        forces[2] += (GRAVITY_CONST * object_1.mass * object_2.mass * (object_1.pos_z - object_2.pos_z)) / (x * x * x);
    }
}

/* Fuerza gravitatoria que ejerce un objeto */
void calc_gravitational(int num_objects, int i, std::vector<object> &objects, double *forces)
{
    for (int j = 0; j < num_objects; j++)
    {
        if (j != i)
        {
            //cout << "I: " << i << "  J: " << j << endl;
            //cout << "Forces Antes " << j << " ax: " << forces[0] << " ay: " << forces[1] << " az: " << forces[2] << "\n";
            vector_gravitational_force(objects[j], objects[i], forces);
            //cout << "Forces Despues " << j << " ax: " << forces[0] << " ay: " << forces[1] << " az: " << forces[2] << "\n";
        }
    }
}

/* Vector aceleración */
void vector_acceleration(object object_1, double *forces, vector_elem *acceleration)
{
    if (object_1.mass != 0)
    {
        // cout<<"ax: "<<forces[0]<<" ay: "<<forces[1]<<" az: "<<forces[2]<<"\n";
        /* Cálculo del vector aceleración */
        acceleration->x = forces[0] / object_1.mass;
        acceleration->y = forces[1] / object_1.mass;
        acceleration->z = forces[2] / object_1.mass;
        // cout<<"ax: "<<forces[0]<<" ay: "<<forces[1]<<" az: "<<forces[2]<<"\n";
    }
    else
    {
        acceleration->x = 0;
        acceleration->y = 0;
        acceleration->z = 0;
    }
}

/* Vector velocidad */
void vector_speed(object *object_1, vector_elem *acceleration, double time_step)
{
    // cout<<" posx: "<<object_1->pos_x<<" posy: "<< object_1->pos_y<<" posz: "<<object_1->pos_z<<" speedx: "<<object_1->speed_x<<" speedy: "<< object_1->speed_y<<" speedz: "<<object_1->speed_z<<"\n";
    /* Cálculo del vector velocidad */
    object_1->speed_x = object_1->speed_x + (acceleration->x * time_step);
    object_1->speed_y = object_1->speed_y + (acceleration->y * time_step);
    object_1->speed_z = object_1->speed_z + (acceleration->z * time_step);
    // cout<<" posx: "<<object_1->pos_x<<" posy: "<< object_1->pos_y<<" posz: "<<object_1->pos_z<<" speedx: "<<object_1->speed_x<<" speedy: "<< object_1->speed_y<<" speedz: "<<object_1->speed_z<<"\n";
}

/* Vector de posicion */
void vector_position(object *object_1, double time_step)
{
    /* Cálculo del vector posición */
    object_1->pos_x = object_1->pos_x + (object_1->speed_x * time_step);
    object_1->pos_y = object_1->pos_y + (object_1->speed_y * time_step);
    object_1->pos_z = object_1->pos_z + (object_1->speed_z * time_step);
}

/* Función para recolocar al objeto si traspasa los límites */
void check_border(object *object_1, double size_enclosure)
{
    // Checks posición x
    if (object_1->pos_x <= 0)
    {
        object_1->pos_x = 0;
        object_1->speed_x = -1 * (object_1->speed_x);
        // cout << "Toca el borde x 0"<<object_1->pos_x<<" "<<object_1->speed_x<<"\n";
    }
    else if (object_1->pos_x >= size_enclosure)
    {
        object_1->pos_x = size_enclosure;
        object_1->speed_x = -1 * (object_1->speed_x);
        // cout << "Toca el borde x grande"<<object_1->pos_x<<"\n";
    }

    // Checks posición y
    if (object_1->pos_y <= 0)
    {
        object_1->pos_y = 0;
        object_1->speed_y = -1 * (object_1->speed_y);
        // cout << "Toca el borde y 0"<<object_1->pos_y<<"\n";
    }
    else if (object_1->pos_y >= size_enclosure)
    {
        object_1->pos_y = size_enclosure;
        object_1->speed_y = -1 * (object_1->speed_y);
        // cout << "Toca el borde y grande"<<object_1->pos_y<<"\n";
    }

    // Checks posición z
    if (object_1->pos_z <= 0)
    {
        object_1->pos_z = 0;
        object_1->speed_z = -1 * (object_1->speed_z);
        // cout << "Toca el borde z 0"<<object_1->pos_z<<"\n";
    }
    else if (object_1->pos_z >= size_enclosure)
    {
        object_1->pos_z = size_enclosure;
        object_1->speed_z = -1 * (object_1->speed_z);
        // cout << "Toca el borde z grande"<<object_1->pos_z<<"\n";
    }
}

/* Comprobar colisión entre dos objetos (distancia euclídea entre objetos menor que 1) */
bool check_collision(object object_1, object object_2)
{
    // cout<<"H: "<<euclidean_norm(object_1,object_2)<<"\n";
    if (euclidean_norm(object_1, object_2) < 1)
    {
        return true;
    }
    return false;
}