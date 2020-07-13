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
#include "../include/header.h"
#include "../include/preprocess_data.h"
#include "../include/io.h"

pre_t *init_preprocess_struct(void)
{
    pre_t *preData=(pre_t *)malloc(sizeof(pre_t));

    preData->input=(char *)malloc(INFILENAMESIZE*sizeof(char));
    assert(preData->input);
    preData->output=(char *)malloc(INFILENAMESIZE*sizeof(char));
    assert(preData->output);
    preData->sampleList=(char **)malloc(sizeof(char *));
    assert(preData->sampleList);
    preData->sampleList[0]=(char *)malloc(sizeof(char));
    assert(preData->sampleList[0]);

    strcpy(preData->input,"");
    strcpy(preData->output,"");
    strcpy(preData->sampleList[0],"");
    preData->ploidy=2;
    preData->sampleListSize=0;

    return preData;
}

void print_preprocess_struct(pre_t *preData)
{
    printf("Input: %s\n",preData->input);
    printf("Output: %s\n",preData->output);
    printf("%d Samples:\n",preData->sampleListSize);
    for(int i=0; i < preData->sampleListSize; i++)
    {
        printf("\t%s\n", preData->sampleList[i]);
    }
}

void free_preprocess_struct(pre_t *preData)
{
    free(preData->input);
    free(preData->output);
    // minor workaround to changing lots of stuff, cause cannot set list to null in
    // the struct. Let it flow
    if(!preData->sampleListSize)
        preData->sampleListSize++;
    for(int i=0; i < preData->sampleListSize; i++)
    {
        free(preData->sampleList[i]);
    }
    free(preData->sampleList);

    free(preData);
}

void preprocess_data(pre_t *preData, header_t *headerData)
{
    //we read the header file to get the two important header lines and
    //(last line of header) the snips per file, the snip size in bits, the total snips,
    //and the minimum and maximum positions inside the file
    readHeaderFile(preData, headerData);
    //this check is needed if the parser is given a large enough file size that the
    //total snips are less than the calculated snip size
    if(headerData->snipsPerFile < 0 ||
       headerData->snipsPerFile > headerData->totalSnips)
    {
        headerData->snipsPerFile=headerData->totalSnips;
    }

    //with the given window(s) or input list we find the files we need to
    //open to read the snips we need
    findFiles(preData, headerData);
}
