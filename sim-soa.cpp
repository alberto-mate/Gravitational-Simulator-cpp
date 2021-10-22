/* Librerias */
#include <iostream>
#include <math.h>
#include <fstream>
#include <random>
#include <vector>
#include <iomanip>

using namespace std;

/* CONSTANTES */
const double GRAVITY_CONST = 6.674 * pow(10, -11); // Constante gravedad universal
const double M = pow(10, 21);                      // Media (distribución normal)
const double SDM = pow(10, 15);                    // Desviación (distribución normal)

/* ESTRUCTURAS */
/* Estructura objeto */
struct object
{
    vector<double> pos_x;
    vector<double> pos_y;
    vector<double> pos_z;
    vector<double> speed_x;
    vector<double> speed_y;
    vector<double> speed_z;
    vector<double> mass;
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
void vector_speed(object *objects, int i, vector_elem *acceleration, double time_step);
void vector_position(object *objects, int i, double time_step);
void check_border(object *objects, int i, double size_enclosure);
bool check_collision(object objects, int i, int j);

/* MAIN */
int main(int argc, char const *argv[])
{

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

    /* Impresión por pantalla de las variables */
    cout << "num_objects: " << num_objects << "\n";
    cout << "num_iterations: " << num_iterations << "\n";
    cout << "random_seed: " << random_seed << "\n";
    cout << "size_enclosure: " << size_enclosure << "\n";
    cout << "time_step: " << time_step << "\n";

    /* Coordenadas y masas pseudoaleatorias */
    mt19937_64 gen(random_seed);
    uniform_real_distribution<double> position_dist(0.0, nextafter(size_enclosure, numeric_limits<double>::max()));
    normal_distribution<double> mass_dist{M, SDM};

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


    /* Bucle anidado para comprobar colisiones entre objetos previas a las iteraciones */
    for (long unsigned int i = 0; i < objects.mass.size(); i++)
    {
        for (long unsigned int j = i + 1; j < objects.mass.size(); j++)
        {
            // Comprobar colisiones
            if (i != j)
            { // Colision entre objetos diferentes que no hayan sido eliminados con anterioridad
                if (check_collision(objects, i, j))
                {   // Comprobar colisión
                    // Actualización de la masa y velocidades del primer objeto que colisiona generando uno nuevo

                    //cout << "Collision objects " << i << " and " << j << endl;
                    //cout << "Obj: " << i << " posx: " << objects.pos_x[i] << " posy: " << objects.pos_y[i] << " posz: " << objects.pos_z[i] << " speedx: " << objects.speed_x[i] << " speedy: " << objects.speed_y[i] << " speedz: " << objects.speed_z[i] << " mass " << objects.mass[i] << "\n";
                    //cout << "Obj: " << i << " posx: " << objects.pos_x[j] << " posy: " << objects.pos_y[j] << " posz: " << objects.pos_z[j] << " speedx: " << objects.speed_x[j] << " speedy: " << objects.speed_y[j] << " speedz: " << objects.speed_z[j] << " mass " << objects.mass[j] << "\n";
                    //cout << "Body " << j << " removed" << endl;

                    // Actualización de la masa y velocidades del primer objeto que colisiona generando uno nuevo
                    objects.mass[i] += objects.mass[j];
                    objects.speed_x[i] += objects.speed_x[j];
                    objects.speed_y[i] += objects.speed_y[j];
                    objects.speed_z[i] += objects.speed_z[j];

                    //cout << "Object " << i << " after collapse" << endl;
                    //cout << "Obj: " << i << " posx: " << objects.pos_x[i] << " posy: " << objects.pos_y[i] << " posz: " << objects.pos_z[i] << " speedx: " << objects.speed_x[i] << " speedy: " << objects.speed_y[i] << " speedz: " << objects.speed_z[i] << " mass " << objects.mass[i] << "\n";
                    
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
    }

    // Actualizamos el número de objetos en el vector
    num_objects = objects.mass.size();

    /* Iteraciones */
    for (int iteration = 0; iteration < num_iterations; iteration++)
    {
        cout << "\nIteracion " << iteration << "\n";
        /* Bucle para obtener nuevas propiedades de los objetos en la iteración (fuerzas, aceleración y velocidad) */
        for (int i = 0; i < num_objects; i++)
        {
            if (objects.mass[i] != 0.0)
            { // Solo entrarán en el condicional objetos que no se han eliminado
                // Cálculo de la fuerza gravitatoria
                double forces[3] = {0.0, 0.0, 0.0};
                calc_gravitational(num_objects, i, objects, forces);
                // cout<< "Obj: " << i<<" posx: "<<objects[i].pos_x<<" posy: "<< objects[i].pos_y<<" posz: "<<objects[i].pos_z<<" speedx: "<<objects[i].speed_x<<" speedy: "<< objects[i].speed_y<<" speedz: "<<objects[i].speed_z<<" mass: "<<objects[i].mass<<"\n";
                cout << "Forces " << i << " ax: " << forces[0] << " ay: " << forces[1] << " az: " << forces[2] << "\n";
                // Cálculo del vector aceleración
                vector_elem *acceleration = (vector_elem *)malloc(sizeof(vector_elem));
                vector_acceleration(objects, i, forces, acceleration);
                // cout<< "Obj: " << i<<" posx: "<<objects[i].pos_x<<" posy: "<< objects[i].pos_y<<" posz: "<<objects[i].pos_z<<" speedx: "<<objects[i].speed_x<<" speedy: "<< objects[i].speed_y<<" speedz: "<<objects[i].speed_z<<"\n";
                //  Cálculo del vector velocidad
                vector_speed(&objects, i, acceleration, time_step);
                // cout<< "Obj: " << i<<" posx: "<<objects[i].pos_x<<" posy: "<< objects[i].pos_y<<" posz: "<<objects[i].pos_z<<" speedx: "<<objects[i].speed_x<<" speedy: "<< objects[i].speed_y<<" speedz: "<<objects[i].speed_z<<"\n";
            }
        }


        /* Bucle para calcular posiciones y comprobar bordes */
        for (int i = 0; i < num_objects; i++)
        {
            // Cálculo del vector posiciones
            vector_position(&objects, i, time_step);
            // cout<< "Obj: " << i<<" posx: "<<objects[i].pos_x<<" posy: "<< objects[i].pos_y<<" posz: "<<objects[i].pos_z<<" speedx: "<<objects[i].speed_x<<" speedy: "<< objects[i].speed_y<<" speedz: "<<objects[i].speed_z<<"\n";
            //  Comprobar bordes
            check_border(&objects, i, size_enclosure);
            // cout<< "Obj: " << i<<" posx: "<<objects[i].pos_x<<" posy: "<< objects[i].pos_y<<" posz: "<<objects[i].pos_z<<" speedx: "<<objects[i].speed_x<<" speedy: "<< objects[i].speed_y<<" speedz: "<<objects[i].speed_z<<"\n";
        }
        // cout<<"Nuevas posiciones calculadas \n";

        /* Bucle anidado para comprobar colisiones entre objetos */
        for (long unsigned int i = 0; i < objects.mass.size(); i++)
        {
            for (long unsigned int j = i + 1; j < objects.mass.size(); j++)
            {
                // Comprobar colisiones
                if (i != j)
                { // Colision entre objetos diferentes que no hayan sido eliminados con anterioridad
                    if (check_collision(objects, i, j))
                    {   // Comprobar colisión
                        // Actualización de la masa y velocidades del primer objeto que colisiona generando uno nuevo

                        cout << "Collision objects " << i << " and " << j << endl;
                        cout << "Obj: " << i << " posx: " << objects.pos_x[i] << " posy: " << objects.pos_y[i] << " posz: " << objects.pos_z[i] << " speedx: " << objects.speed_x[i] << " speedy: " << objects.speed_y[i] << " speedz: " << objects.speed_z[i] << " mass " << objects.mass[i] << "\n";
                        cout << "Obj: " << i << " posx: " << objects.pos_x[j] << " posy: " << objects.pos_y[j] << " posz: " << objects.pos_z[j] << " speedx: " << objects.speed_x[j] << " speedy: " << objects.speed_y[j] << " speedz: " << objects.speed_z[j] << " mass " << objects.mass[j] << "\n";
                        cout << "Body " << j << " removed" << endl;

                        objects.mass[i] += objects.mass[j];
                        objects.speed_x[i] += objects.speed_x[j];
                        objects.speed_y[i] += objects.speed_y[j];
                        objects.speed_z[i] += objects.speed_z[j];

                        cout << "Object " << i << " after collapse" << endl;
                        cout << "Obj: " << i << " posx: " << objects.pos_x[i] << " posy: " << objects.pos_y[i] << " posz: " << objects.pos_z[i] << " speedx: " << objects.speed_x[i] << " speedy: " << objects.speed_y[i] << " speedz: " << objects.speed_z[i] << " mass " << objects.mass[i] << "\n";
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
        }

        // Actualizamos el número de objetos en el vector
        num_objects = objects.mass.size();
        cout << "Fin iteración: " << iteration << " Num objetos:" << num_objects << "\n";
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
}
/* FUNCIONES */
/* Distancia euclídea entre dos objetos */
double euclidean_norm(object objects, int i, int j)
{
    // cout<<"Pos 2 "<<" px: "<<(object_1.pos_x)<<" "<< (object_2.pos_x)<<" "<<(object_1.pos_x - object_2.pos_x)<<" py: "<<(object_1.pos_y)<<" "<< (object_2.pos_y)<<" "<<(object_1.pos_y - object_2.pos_y)<<" pz: "<<(object_1.pos_z)<<" "<< (object_2.pos_z)<<" "<<(object_1.pos_z - object_2.pos_z)<<"\n";
    return std::sqrt((objects.pos_x[i]- objects.pos_x[j]) * (objects.pos_x[i]- objects.pos_x[j]) + (objects.pos_y[i]- objects.pos_y[j]) * (objects.pos_y[i]- objects.pos_y[j]) + (objects.pos_z[i]- objects.pos_z[j]) * (objects.pos_z[i]- objects.pos_z[j]));
}

/* Fuerza gravitatoria entre dos objetos */
void vector_gravitational_force(object objects, int i, int j, double *forces)
{
    // cout<<"Pos 1 "<<" px: "<<(object_1.pos_x)<<" "<< (object_2.pos_x)<<" "<<(object_1.pos_x - object_2.pos_x)<<" py: "<<(object_1.pos_y)<<" "<< (object_2.pos_y)<<" "<<(object_1.pos_y - object_2.pos_y)<<" pz: "<<(object_1.pos_z)<<" "<< (object_2.pos_z)<<" "<<(object_1.pos_z - object_2.pos_z)<<"\n";
    // cout<<"Forces 1 "<<" ax: "<<forces[0]<<" ay: "<<forces[1]<<" az: "<<forces[2]<<"\n";
    double x = (euclidean_norm(objects, i, j));
    if (x != 0)
    {
        //cout << "Vector x de la fuerza: " << ttt << endl;
        //cout << "Euclidean: " << x << endl;
        forces[0] += (GRAVITY_CONST * objects.mass[i] * objects.mass[j] * (objects.pos_x[i] - objects.pos_x[j])) / (x * x * x);
        forces[1] += (GRAVITY_CONST * objects.mass[i] * objects.mass[j] * (objects.pos_y[i] - objects.pos_y[j])) / (x * x * x);
        forces[2] += (GRAVITY_CONST * objects.mass[i] * objects.mass[j] * (objects.pos_z[i] - objects.pos_z[j])) / (x * x * x);
    }
}

/* Fuerza gravitatoria que ejerce un objeto */
void calc_gravitational(int num_objects, int i, object objects, double *forces)
{
    for (int j = 0; j < num_objects; j++)
    {
        if (j != i)
        {
            //cout << "I: " << i << "  J: " << j << endl;
            //cout << "Forces Antes " << j << " ax: " << forces[0] << " ay: " << forces[1] << " az: " << forces[2] << "\n";
            vector_gravitational_force(objects, j, i, forces);
            //cout << "Forces Despues " << j << " ax: " << forces[0] << " ay: " << forces[1] << " az: " << forces[2] << "\n";
        }
    }
}

/* Vector aceleración */
void vector_acceleration(object objects, int i, double *forces, vector_elem *acceleration)
{
    if (objects.mass[i] != 0)
    {
        // cout<<"ax: "<<forces[0]<<" ay: "<<forces[1]<<" az: "<<forces[2]<<"\n";
        /* Cálculo del vector aceleración */
        acceleration->x = forces[0] / objects.mass[i];
        acceleration->y = forces[1] / objects.mass[i];
        acceleration->z = forces[2] / objects.mass[i];
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
void vector_speed(object *objects, int i, vector_elem *acceleration, double time_step)
{
    // cout<<" posx: "<<object_1->pos_x<<" posy: "<< object_1->pos_y<<" posz: "<<object_1->pos_z<<" speedx: "<<object_1->speed_x<<" speedy: "<< object_1->speed_y<<" speedz: "<<object_1->speed_z<<"\n";
    /* Cálculo del vector velocidad */
    objects->speed_x[i] = objects->speed_x[i] + (acceleration->x * time_step);
    objects->speed_y[i] = objects->speed_y[i] + (acceleration->y * time_step);
    objects->speed_z[i] = objects->speed_z[i] + (acceleration->z * time_step);
    // cout<<" posx: "<<object_1->pos_x<<" posy: "<< object_1->pos_y<<" posz: "<<object_1->pos_z<<" speedx: "<<object_1->speed_x<<" speedy: "<< object_1->speed_y<<" speedz: "<<object_1->speed_z<<"\n";
}

/* Vector de posicion */
void vector_position(object *objects, int i, double time_step)
{
    /* Cálculo del vector posición */
    objects->pos_x[i] = objects->pos_x[i] + (objects->speed_x[i] * time_step);
    objects->pos_y[i] = objects->pos_y[i] + (objects->speed_y[i] * time_step);
    objects->pos_z[i] = objects->pos_z[i] + (objects->speed_z[i] * time_step);
}

/* Función para recolocar al objeto si traspasa los límites */
void check_border(object *objects, int i, double size_enclosure)
{
    // Checks posición x
    if (objects->pos_x[i] <= 0)
    {
        objects->pos_x[i] = 0;
        objects->speed_x[i] = -1 * (objects->speed_x[i]);
        // cout << "Toca el borde x 0"<<object_1->pos_x<<" "<<object_1->speed_x<<"\n";
    }
    else if (objects->pos_x[i] >= size_enclosure)
    {
        objects->pos_x[i] = size_enclosure;
        objects->speed_x[i] = -1 * (objects->speed_x[i]);
        // cout << "Toca el borde x grande"<<object_1->pos_x<<"\n";
    }

    // Checks posición y
    if (objects->pos_y[i] <= 0)
    {
        objects->pos_y[i] = 0;
        objects->speed_y[i] = -1 * (objects->speed_y[i]);
        // cout << "Toca el borde y 0"<<object_1->pos_y<<"\n";
    }
    else if (objects->pos_y[i] >= size_enclosure)
    {
        objects->pos_y[i] = size_enclosure;
        objects->speed_y[i] = -1 * (objects->speed_y[i]);
        // cout << "Toca el borde y grande"<<object_1->pos_y<<"\n";
    }

    // Checks posición z
    if (objects->pos_z[i] <= 0)
    {
        objects->pos_z[i] = 0;
        objects->speed_z[i] = -1 * (objects->speed_z[i]);
        // cout << "Toca el borde z 0"<<object_1->pos_z<<"\n";
    }
    else if (objects->pos_z[i] >= size_enclosure)
    {
        objects->pos_z[i] = size_enclosure;
        objects->speed_z[i] = -1 * (objects->speed_z[i]);
        // cout << "Toca el borde z grande"<<object_1->pos_z<<"\n";
    }
}

/* Comprobar colisión entre dos objetos (distancia euclídea entre objetos menor que 1) */
bool check_collision(object objects, int i, int j)
{
    // cout<<"H: "<<euclidean_norm(object_1,object_2)<<"\n";
    if (euclidean_norm(objects, i, j) < 1)
    {
        return true;
    }
    return false;
}
