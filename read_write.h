//
//// Created by David on 7/7/2018.
////
//
//#ifndef DAVIDFOLLY_READ_WRITE_H
//#define DAVIDFOLLY_READ_WRITE_H
//
//#include <stdio.h>
//#include <stdlib.h>
//
//// read_file
////
//// from a newly opend file (f) read the input parameters advancing the file to the end of input
void read_file(FILE* in);
//
//
//// write_file
////
// concatinate to file out the resultArray
#ifndef DAVIDFOLLY_READ_WRITE_H
#define DAVIDFOLLY_READ_WRITE_H

void write_file(double* resultArray, int size, FILE* out);



#endif //DAVIDFOLLY_READ_WRITE_H


