// Librerias que hemos incluido
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

int main( int argc, char const *argv[]) {
    
    /*Comprobación inicial variables*/
    if (argc!=6){
        cerr <<"Wrong number of arguments\n";
        return -1;
    }

    if ( (atoi(argv[1])<=0 || atoi(argv[2])<=0 || atoi(argv[3])<=0 || atof(argv[4])<=0.0 || atof(argv[5])<=0.0) 
        ||
        (atof(argv[1])!=atoi(argv[1]) || atof(argv[2])!=atoi(argv[2]) || atof(argv[3])!=atoi(argv[3])) ){
        cerr <<"Values of the arguments are incorrect\n";
        return -2;
    }
    int num_objects = atoi(argv[1]);
    int num_iterations = atoi(argv[2]);
    int random_seed = atoi(argv[3]);
    float size_enclosure = atof(argv[4]);
    float time_step = atof(argv[5]);

    cout<<"num_objects: "<<num_objects<<"\n";
    cout<<"num_iterations: "<<num_iterations<<"\n";
    cout<<"random_seed: "<<random_seed<<"\n";
    cout<<"size_enclosure: "<<size_enclosure<<"\n";
    cout<<"time_step: "<<time_step<<"\n";
}

