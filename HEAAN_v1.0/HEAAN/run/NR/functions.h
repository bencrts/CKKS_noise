#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <complex.h>
#include <string>
#include <fstream>

using namespace std;

typedef vector<double> dVec;
typedef complex<double>* cVec;

void print_array(int slots, double *array,int entries)
{
    cout << "(";
    for (int i = 0; i < entries; i++)
    {
        cout << array[i] << ",";
    }
    cout << array[entries - 1] << ")\n";
}

void print_vec(dVec vec){
    cout << "(";
    for(int i = 0; i<vec.size()-1;i++){
        cout << vec[i]<<",";
    }
    cout << vec[vec.size()-1]<<")\n";
}

void scale_vec(dVec &vec,double scale){
    for(int i = 0; i<vec.size();i++){
        vec[i]/=scale;
    }
}
double *random_array(int slots, double a, double b)
{
    //generate a random array of doubles in the interval [a,b] of dimension slots
    srand(time(NULL));
    double *array = new double[slots];
    for (int i = 0; i < slots; i++)
    {
        array[i] = (b - a) * rand() / RAND_MAX;
        array[i] += a;
    }
    return array;
}

double *fz(int slots, double *array)
{
    //calculate 1/z component wise
    double *fz = new double[slots];
    for (int i = 0; i < slots; i++)
    {
        fz[i] = 1 / array[i];
    }
    return fz;
}

double max_diff(int slots, double *array1, double *array2)
{
    //find infinity norm of array1 - array2. Must be real!
    double max = 0;
    for(int i = 0; i<slots;i++){
        if(abs(array1[i]-array2[i])>max)max = abs(array1[i]-array2[i]);
    }
    return max;
}

double max_log_diff(int slots,double *array1,double *array2){
    return log2(max_diff(slots,array1,array2));
}

double *real_part(int slots, cVec array){
    double *real = new double[slots];
    for(int i = 0; i< slots;i++){
        real[i] = array[i].real();
    }
    return real;
}

void update_max_avg(dVec &max,dVec &avg,int iter,double value){
    if(value > max[iter])max[iter]=value;
    avg[iter]+=value;
    //do you want to add logs here? avg log or log avg
}

void print_vec_toFile(string filename,dVec vec){
    std::ofstream outFile(filename, std::ios_base::app | std::ios_base::out);
    for(int i = 0; i<vec.size()-1;i++){
        outFile << vec[i]<<",";
    }
    outFile << vec[vec.size()-1]<<"\n";
}