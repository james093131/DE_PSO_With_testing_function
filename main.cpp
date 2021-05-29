#include "normal.hpp"

int main(int argc,const char *argv[])
{
    int run = atoi(argv[1]);//just run
    int iteration = atoi(argv[2]);//just iteration
    int DIM = atoi(argv[3]); //Function Dimension
    int pop = atoi(argv[4]); //Function Dimension
    int OUTPUT_NODE_QUANTITY =atoi(argv[5]);
    int F = atoi(argv[6]);
    double CR = atof(argv[7]);
    double Factor = atof(argv[8]);
    const char *Algorithm = argv[9];

    if( argc > 1 )
    {
        if(Algorithm == std::string("D"))
        {
            DE de;  
            de.ALL(run,iteration, pop,DIM,F,OUTPUT_NODE_QUANTITY,CR,Factor,Algorithm);
        }
        else if((Algorithm == std::string("P")))
        {
            PSO pso;
            pso.ALL(run,iteration,pop,DIM,F,OUTPUT_NODE_QUANTITY,Algorithm);
        }
       
        
    }
    
    return 0;
}