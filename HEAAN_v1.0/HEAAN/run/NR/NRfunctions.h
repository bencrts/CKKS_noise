using namespace std;

double T1(double a,double b){
    return -8/(a*a+6*a*b+b*b);
}

double T0(double a,double b){
    return 8*(a+b)/(a*a+6*a*b+b*b);
}
