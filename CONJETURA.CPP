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
    unsigned long jA; //contador de primos AKS
    unsigned long jC; //contador de primos Miller-Rabin
    short auxA; //auxiliar AKS
    short auxC; //auxiliar Miller-Rabin
    
    ifstream entrada("rango.txt", ios::in); //abro fichero con rango 
    entrada >> a >> b; //leo limites
    entrada.close(); //cierro fichero
    
    ofstream salida;  // objeto de la clase ofstream
    salida.open("salida.txt",ios::out); //abro archivo

    salida << "Tabla de tiempos AKS" << endl;

    jA=0; //inicializo contador de primos
    jC=0; //inicializo contador de primos
        
for (i=a;i <= b; i++)  //criba AKS
{ 
    auxC=AKS(to_ZZ(i),1); 
    auxA=AKS(to_ZZ(i),0);
    if (auxA != auxC) 
    {
       salida << i << " " << auxA << " " << auxC << " " << //
       ProbPrime(to_ZZ(i), 1) << endl; //muestro quien falla
    }
    jA+=auxA; //cuento el numero de primos
    jC+=auxC; //cuento el numero de primos
}   
    salida << endl << "La cantidad de números primos en el intervalo ["//
    <<a<<","<<b<<"] es " << jA << "," << jC << endl;

    salida.close(); //cierro el fichero
} // main
