#ifndef MAINH
#define MAINH
//To store the input parameters for voltage. 
//V=V.bias+V.sin*sin(2*M_PI*f*T);
struct signal{
    double bias; //Voltage bias
    double sin; //Amplitude of sine
    double f; //Frequency of sine
};
#endif
