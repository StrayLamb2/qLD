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
#include "../include/mem.h"

void *realloc_buff(void *buffer, unsigned int newSize, char type[32])
{
    void *ptr=NULL;
    if(!strcmp(type, "char"))
    {
        ptr=(char *)realloc((char *)buffer, newSize*sizeof(char));
    }
    else if(!strcmp(type, "char*"))
    {
        ptr=(char **)realloc((char **)buffer, newSize*sizeof(char *));
    }
    else if(!strcmp(type, "unsigned int"))
    {
        ptr=(unsigned int *)realloc((unsigned int *)buffer,
                newSize*sizeof(unsigned int));
    }
    else if(!strcmp(type, "int"))
    {
        ptr=(int *)realloc((int *)buffer, newSize*sizeof(int));
    }

    else if(!strcmp(type, "ResultDataType"))
    {
        ptr=(ResultDataType *)realloc((ResultDataType *)buffer,
                newSize*sizeof(ResultDataType));
    }
    else if(!strcmp(type, "inputDataType_x64"))
    {
        ptr=(inputDataType_x64 *)realloc((inputDataType_x64 *)buffer,
                newSize*sizeof(inputDataType_x64));
    }
    else
    {
#ifdef VERBOSE
        printf("ERROR: Realloc with type \"%s\" failed\n",type);
#else
        ;
#endif
    }
    assert(ptr);
    return ptr;
}

void print_alloc_size(char table[32], unsigned int size)
{
    printf("Size of %s is %u\n", table, size);
}

table_x64* create_table_x64(unsigned int tableSize)
{
    table_x64 *newTable=(table_x64 *)malloc(sizeof(table_x64));

    unsigned int *POStable=(unsigned int*)malloc(tableSize*sizeof(unsigned int));
    assert(POStable);
    newTable->POStable=POStable;

    char **IDtable=(char **)malloc(tableSize*sizeof(char *));
    assert(IDtable);
    for(unsigned int i=0; i < tableSize; i++)
    {
        IDtable[i]=(char *)malloc(IDLENGTH*sizeof(char));
        assert(IDtable[i]);
    }
    newTable->IDtable=IDtable;

    unsigned int *BCtable=(unsigned int*)malloc(tableSize*sizeof(unsigned int));
    assert(BCtable);
    newTable->BCtable=BCtable;

    inputDataType_x64 *SNPtable=(inputDataType_x64 *)malloc(sizeof(inputDataType_x64));
    assert(SNPtable);
    newTable->SNPtable=SNPtable;

    newTable->tableSize=tableSize;
    newTable->tableIndex=0;
    newTable->compSize=0;
    newTable->compIndex=0;

    return newTable;
}

void free_table_x64(table_x64 *tableData)
{
    free(tableData->SNPtable);
    free(tableData->POStable);
    for(unsigned int i=0; i < tableData->tableSize; i++)
    {
        free(tableData->IDtable[i]);
    }
    free(tableData->IDtable);
    free(tableData->BCtable);

    free(tableData);
}

helper_t* create_helper_t(int ploidy)
{
    helper_t *newHelper=(helper_t *)malloc(sizeof(helper_t));

    newHelper->ploidy=ploidy;

    char *snipChar=(char *)malloc(STRINGLENGTH*sizeof(char));
    assert(snipChar);
    newHelper->snipChar=snipChar;

    newHelper->snipCharLength=STRINGLENGTH;

    newHelper->snipCharIndex=0;

    char *snipPart=(char *)malloc((ploidy+1)*sizeof(char));
    assert(snipPart);
    newHelper->snipPart=snipPart;

    char *line=(char *)malloc(STRINGLENGTH*sizeof(char));
    assert(line);
    newHelper->line=line;

    newHelper->lineLength=STRINGLENGTH;

    char *word=(char *)malloc(STRINGLENGTH*sizeof(char));
    assert(word);
    newHelper->word=word;

    newHelper->wordLength=STRINGLENGTH;

    return newHelper;
}

void free_helper_t(helper_t *helperData)
{
    free(helperData->snipChar);
    free(helperData->snipPart);
    free(helperData->line);
    free(helperData->word);

    free(helperData);
}
