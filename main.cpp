#include "normal.hpp"

int main(int argc,const char *argv[])
{
    int run = atoi(argv[1]);//just run
    int iteration = atoi(argv[2]);//just iteration
    int DIM = atoi(argv[3]); //Function Dimension
    int pop = atoi(argv[4]); //Function Dimension
    int OUTPUT_NODE_QUANTITY =atoi(argv[5]);
    const char *F = argv[6];
    double CR = atof(argv[7]);
    double Factor = atof(argv[8]);


    if( argc > 1 )
    {
        
        DE de;  
        de.ALL(run,iteration, pop,DIM,F,OUTPUT_NODE_QUANTITY,CR,Factor);
        
    }
    
    return 0;
}