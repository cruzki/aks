#include <cstdlib>
#include <iostream>
#include <NTL/ZZ.h> //libreria con numeros ZZ
#include "imp_C_O(12).h" //libreria con el AKS
#include <fstream> //escritura de archivos

NTL_CLIENT

int main(int argc, char *argv[])
//Programa que realiza un test de esfuerzo bruto. Preparado para
//realizar computación GRID "manual". Testea la primalidad (mediante
//el algoritmo AKS y una implementación de Miller-Rabin) de los
//numeros dados desde rango.txt e imprime los resultados en 
//salida.txt.
//Con la variable opcion, tarcer parámetro leído de rango.txt, 
//controlamos que versión del algoritmo AKS deseamos ejecutar.
{
    ZZ a; //inicio
    ZZ i; //contador
    ZZ b; //final
    int opcion; //0 = algoritmo estándar, !=0 algoritmo modificado 
    unsigned long j; //contador de primos
    short aux; //auxiliar
    
    ifstream entrada("rango.txt", ios::in); //abro ficero con rango 
    entrada >> a >> b >> opcion; //leo limites
    entrada.close(); //cierro fichero
    
    ofstream salida;  // objeto de la clase ofstream
    salida.open("salida.txt",ios::out); //abro archivo

    clock_t total; //tiempo total
    clock_t t; //variable que controla el tiempo

    salida << "Tabla de tiempos AKS" << endl;

    j=0; //inicializo contador de primos
    total=clock(); //cronometro total
        
for (i=a;i <= b; i++)  //criba AKS
{ 
    t=clock(); //cronometro
    aux=AKS(to_ZZ(i),opcion); //auxiliar para no contar la suma de los primos
    salida << i << " " << aux << " " << (clock()-t)/(double)CLOCKS_PER_SEC << endl;
    j+=aux; //cuento el numero de primos
}   
    salida << endl << "Tiempo total empleado = " << (clock()-total)/(double)CLOCKS_PER_SEC << endl;
    salida << endl << "La cantidad de números primos en el intervalo ["<<a<<","<<b<<"] es " << j << endl;

    salida << endl << "Tabla de tiempos Miller-Rabin test" <<endl;
    
    j=0; //inicializo contador de primos
    total=clock(); //cronometro total
    
for (i=1; i <=b; i++)  //criba probabilistica
{  
    t=clock(); //cronometro
    aux=ProbPrime(to_ZZ(i), 1);
    salida << i << " " << aux << " " <<(clock()-t)/(double)CLOCKS_PER_SEC << endl;
    j+=aux;
} 
    salida << endl << "Tiempo total empleado = " << (clock()-total)/(double)CLOCKS_PER_SEC << endl << endl;
    salida << endl << "La cantidad de números primos en el intervalo ["<<a<<","<<b<<"] es " << j << endl;
    salida.close(); //cierro el fichero
} // main
