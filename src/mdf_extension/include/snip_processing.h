/*
qLD - High performance computation of Linkage disequilibrium
Copyright (C) 2020  C. Theodoris, N. Alachiotis

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
#ifndef SP_H
#define SP_H

#include "header.h"
#include "structs.h"
#include "mem.h"
#include "io.h"

int isValidVCFBase(char input);

int scanStateVector(char * stateVector, char X);

int getStatesREF(char * string, char * stateVector, int line);

int getStatesNum(char * stateVector);

int getStatesALT(char * string, char * stateVector, int line);

float getValueAF(char * string, int line);

int checkVTisSNP(char * string);

int getGTpos(char * string);

int getGTfield(char * string, int GTpos);

void dataShuffleKnuth(char * data, int startIndex, int endIndex);

void getGTdata(char * string, 
               char * stateVector, 
               int statesTotal, 
               char * sampleData);

void processSampleVCF(helper_t *helperData,
                      int GTpos,
                      char *stateVector,
                      int statesALT);

int mapCharToInt(char a);

int getDataType(int * states);

int determineStates(int snipSize, char * line, int *states);

int removeNonPolymorphicSiteBIN(char* line, int snipSize, int filterOut);

int removeNonPolymorphicSiteDNA(char* line, int snipSize, int filterOut);

int removeNonPolymorphicSite(char* line, int snipSize,int states, int filterOut);

void switchValues(int * input, char * sortedStates, int index0, int index1);

void sort(int * freqs, char * sortedStates);

int minorToMajor(int* inp, char* outp);

void countStates(int *numStates, int *gap, int *states);

char getCharBIN(char input, char state0, char state1);

void convertAlignmentBIN(int snipSize,char* line);

void binaryDeduction(int statesIn,int snipSize,char* line);

void compressSnip_x64(table_x64 *tableData,
                     char * snipChar,
                     int size_snp);

#endif
