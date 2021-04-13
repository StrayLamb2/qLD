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
#include "../include/input_list.h"
#include "../include/correlate_MEM.h"
#include "../include/correlate_IO.h"
#include "../include/correlate_IOx64.h"
#include "../include/snip_processing.h"
#include "../include/read_file.h"

// Fill tables from file in MDF mode
int writeMDF_x64(table_x64 *tableData,
                 helper_t *helperData)
{
    int eol=0, index=0;
    unsigned int i, BC=0;
    int elementIndex=-1;
    int position=-1;
    char *strtoul_ptr;
    int wordLength=STRINGLENGTH;
    char ID[IDLENGTH]="";

    // Parse initial info from line and save it
    while(getWordFromString(helperData->line,
                            &(helperData->word),
                            &eol,
                            &wordLength,
                            &index) == 1)
    {
        elementIndex++;
        switch(elementIndex)
        {
            case 0:
                // CHROM
                //if(strcmp(VCF_alignment_name,*word)!=0)
                //{
                //strncpy(VCF_alignment_name, *word, MAX_CHROM_NAME_VCF);
                //assert(strlen(VCF_alignment_name)!=0);
                //}
                break;
            case 1: // POS
                position=atoi(helperData->word);
                break;
            case 2: // ID
                strcpy(ID,helperData->word);
                break;
            case 3: // BC
                BC=atoi(helperData->word);
                break;
            default:
                break;
        }
        if(elementIndex == 3 || eol == 1)
            break;
    }
    if (position == -1)
    {
        return -1;
    }
    helperData->snipChar[0]='\0';

    // Reallocate space if not enough
    if(tableData->tableIndex >= tableData->tableSize-1)
    {
        tableData->tableSize+=BUFFER_INCR;
        tableData->SNPtable=(inputDataType_x64 *)realloc_buff(
                                    (void *)tableData->SNPtable,
                                    tableData->tableSize*tableData->compSize,
                                    "inputDataType_x64");
        tableData->POStable=(unsigned int *)realloc_buff(
                                    (void *)tableData->POStable,
                                    tableData->tableSize,
                                    "unsigned int");
        tableData->BCtable=(unsigned int *)realloc_buff(
                                (void *)tableData->BCtable,
                                tableData->tableSize,
                                "unsigned int");
        tableData->IDtable=(char **)realloc_buff(
                                (void *)tableData->IDtable,
                                tableData->tableSize,
                                "char*");
        for(i=tableData->tableSize-BUFFER_INCR; i < tableData->tableSize; i++)
            tableData->IDtable[i]=(char *)malloc(IDLENGTH*sizeof(char));
    }
    tableData->compIndex=0;
    // Fill tables
    while(getWordFromString(helperData->line,
                            &(helperData->word),
                            &eol,
                            &wordLength,
                            &index) == 1)
    {
        tableData->SNPtable[tableData->tableIndex*
                            tableData->compSize+
                            tableData->compIndex++]=strtoul(helperData->word,
                                                            &strtoul_ptr,
                                                            10);
        if(eol == 1)
            break;
    }
    tableData->POStable[tableData->tableIndex]=position;
    tableData->BCtable[tableData->tableIndex]=BC;
    strcpy(tableData->IDtable[tableData->tableIndex],ID);
    tableData->tableIndex++;

    return 1;
}

// Fill tables from file in VCF MODE
int writeToTable_x64(task_t *nodeData,
                     table_x64 *tableData,
                     helper_t *helperData,
                     int pos,
                     int *skippedLines,
                     int *counter)
{
    int eol=0, index=0;
    static int procVCF_count=0;
    unsigned int i;
    int elementIndex=-1, lineSkipped=0, s, rnps_flag=0;
    char stateVector[MAX_STATES_VCF];
    int position=-1, statesREF=0, statesALT=0, VTisSNP=0, GTpos=-1;
    int sampleIndex=-1;
    float AF=0.5;

    int statesNum=0, states[STATESALL];
    for(int i=0; i < STATESALL; ++i)
        states[i]=0;

    int wordLength=STRINGLENGTH;
    char ID[IDLENGTH];

    // Parse initial info from line and save it
    while(getWordFromString(helperData->line,
                            &(helperData->word),
                            &eol,
                            &wordLength,
                            &index) == 1)
    {
        elementIndex++;
        switch(elementIndex)
        {
            case 0:
                // CHROM
                //if(strcmp(VCF_alignment_name,*word)!=0)
                //{
                //strncpy(VCF_alignment_name, *word, MAX_CHROM_NAME_VCF);
                //assert(strlen(VCF_alignment_name)!=0);
                //}
                break;

            case 1: // POS
                position=atoi(helperData->word);
                break;

            case 2: // ID
                if(helperData->word)
                {
                    strcpy(ID,helperData->word);
                    if(ID[strlen(ID)] != '\0')
                        ID[strlen(ID)]='\0';
                }
                else
                {
                    fprintf(stderr,"ERROR: word not found at writeToTable_x64\n");
                    exit(1);
                }
                break;

            case 3: // REF
                if(strlen(helperData->word) > 1)
                {
                    lineSkipped=1;
                    break;
                }
                getStatesREF(helperData->word, stateVector, pos);
                statesREF=getStatesNum(stateVector);

                if(statesREF == 0)
                    lineSkipped=1;
                break;

            case 4: // ALT
                if(strlen(helperData->word) > 1)
                {
                    lineSkipped=1;
                    break;
                }

                getStatesALT(helperData->word, stateVector, pos);
                statesALT=getStatesNum(stateVector);

                if(statesALT != 0 && stateVector[statesALT-1] == '.')
                    for( s = statesALT - 1; s < statesALT-1 + statesREF; ++s)
                        stateVector[s]=stateVector[s-statesREF];
                if(statesALT == 0)
                    lineSkipped=1;
                break;

            case 5: // QUAL
                break;

            case 6: // FILTER
                if(strcmp("PASS", helperData->word) != 0 &&
                   strcmp(".", helperData->word) != 0)
                {
                    lineSkipped=1;
                }
                break;

            case 7: // INFO
                if(!strcmp(".", helperData->word))
                {
                    AF=.9;
                    break;
                }
                AF=getValueAF((helperData->word), pos); // returns -1.0 if AF not found
                if(AF == 0.0 || AF == 1.0)
                {
                    lineSkipped=1;
                    break;
                }
                VTisSNP=checkVTisSNP(helperData->word); // returns -1.0 if VT not found
                if(VTisSNP == 0)
                    lineSkipped=1;
                break;

            case 8: // FORMAT
                GTpos=getGTpos(helperData->word); // returns -1 if GT not found
                if(GTpos == -1)
                    lineSkipped=1;
                break;

            default:
                break;

        }
        if(elementIndex == VCF_HLENGTH-1 || eol == 1 || lineSkipped == 1)
            break;
    }
    if (position == -1)
    {
        return -1;
    }
    if(lineSkipped != 1)
    {
        sampleIndex=0;
        helperData->snipChar[0]='\0';
        while(getWordFromString(helperData->line,
                                &(helperData->word),
                                &eol,
                                &wordLength,
                                &index) == 1)
        {
            assert(sampleIndex < nodeData->snipSize);
            // If valid sample, process it
            if(nodeData->valid_mask[sampleIndex] == 1)
            {
                processSampleVCF(helperData,
                                 GTpos,
                                 stateVector,
                                 statesALT);
                procVCF_count++;
            }
            sampleIndex++;

            if(eol == 1)
                break;
        }
        helperData->snipCharIndex=0;
        if(strlen(helperData->snipChar) != (size_t)nodeData->valid_count)
        {
            fprintf(stderr, "\n\n ERROR: There are %lu nucleotides in position %d. "
                            "Expected %d.\n\n", strlen(helperData->snipChar),
                            position, nodeData->valid_count);
            assert(strlen(helperData->snipChar) == (size_t)nodeData->valid_count);
        }
        // Determine the state of the site
        statesNum=determineStates(nodeData->valid_count, helperData->snipChar, states);
        if(statesNum == -1)
        {
            fprintf(stderr, "\n\n ERROR: Empty alignment (no states found).\n\n");
            exit(1);
        }
        // Remove non-polymorphic sites
        rnps_flag=removeNonPolymorphicSite(helperData->snipChar,
                                           nodeData->valid_count,
                                           statesNum,
                                           0);
        if(!rnps_flag)
        {
            (*counter)++;
            // Re-allocate space if not enough
            if(tableData->tableIndex >= tableData->tableSize-1)
            {
                tableData->tableSize+=BUFFER_INCR;
                tableData->SNPtable=(inputDataType_x64 *)realloc_buff(
                                            (void *)tableData->SNPtable,
                                            tableData->tableSize*tableData->compSize,
                                            "inputDataType_x64");
                tableData->POStable=(unsigned int *)realloc_buff(
                                            (void *)tableData->POStable,
                                            tableData->tableSize,
                                            "unsigned int");
                tableData->BCtable=(unsigned int *)realloc_buff(
                                        (void *)tableData->BCtable,
                                        tableData->tableSize,
                                        "unsigned int");
                tableData->IDtable=(char **)realloc_buff(
                                        (void *)tableData->IDtable,
                                        tableData->tableSize,
                                        "char*");
                for(i=tableData->tableSize-BUFFER_INCR; i < tableData->tableSize; i++)
                    tableData->IDtable[i]=(char *)malloc(IDLENGTH*sizeof(char));
            }
            // Process and save SNP 
            tableData->POStable[tableData->tableIndex]=position;
            tableData->BCtable[tableData->tableIndex]=0;
            strcpy(tableData->IDtable[tableData->tableIndex],ID);
            binaryDeduction(statesNum, nodeData->valid_count, (helperData->snipChar));
            compressSnip_x64(tableData, helperData->snipChar, nodeData->valid_count);
            assert(tableData->BCtable[tableData->tableIndex] <=
                   (unsigned int)nodeData->valid_count);

            tableData->tableIndex++;
        }
        else
            (*skippedLines)++;
    }
    else
        (*skippedLines)++;
    return 1;
}

// Function that reads lines from the files
void readTable_x64(threadData_t *threadData,
                   task_t *nodeData,
                   table_x64 *tableData,
                   helper_t *helperData)
{
    int pos, eol, eof, status, counter=0, index;
    int skippedLines=0;
    inFileType fpIn;
    
    for(int i=0; i < nodeData->filesListNum; i++)
    {

        eol=0;
        eof=0;
        fpIn = FOPEN(nodeData->filesList[i], "r");
        
        if(fpIn == NULL)
        {
            fprintf(stderr,"ERROR: Failed to open file %s\n", nodeData->filesList[i]);
            exit(1);
        }
        assert(fpIn);
        status=getNextLine(fpIn, &(helperData->line), &eol, &eof, &(helperData->lineLength));

        if(status == 1 && eol == 1)
        {
            if(strcmp(helperData->line, nodeData->headerLine1))
            {
                if(!threadData->mdf)
                {
                    fprintf(stderr,"ERROR: Header 1 different\n");
                    exit(1);
                }
            }
            status=getNextLine(fpIn,
                               &(helperData->line),
                               &eol,
                               &eof,
                               &(helperData->lineLength));
            if(status == 1 && eol==1)
            {
                if(strcmp(helperData->line, nodeData->headerLine2))
                {
                    if(!threadData->mdf)
                    {
                        fprintf(stderr,"Header 2 different\n");
                        exit(1);
                    }
                }
                do
                {
                    status=getNextLine(fpIn,
                                       &(helperData->line),
                                       &eol,
                                       &eof,
                                       &(helperData->lineLength));
                    if(status == 1 && eol == 1)
                    {
                        index=0;
                        status=getWordFromString(helperData->line,
                                                 &(helperData->word),
                                                 &eol,
                                                 &(helperData->wordLength),
                                                 &index); //chrom
                        assert(status);
                        status=getWordFromString(helperData->line,
                                                 &(helperData->word),
                                                 &eol,
                                                 &(helperData->wordLength),
                                                 &index); //pos
                        assert(status);
                        pos=atoi(helperData->word);
                        if(pos >= nodeData->posWmin &&
                           pos <= nodeData->posWmax)
                        {
                            if(threadData->mdf)
                            {
                                writeMDF_x64(tableData,
                                             helperData);
                            }
                            else
                            {
                                writeToTable_x64(nodeData,
                                                 tableData,
                                                 helperData,
                                                 pos,
                                                 &skippedLines,
                                                 &counter);
                            } 
                        }
                    }
                    else
                        pos=0;
                }while(!eof && nodeData->posWset && pos < nodeData->posWmax);
            }
        }
        else
        {
            fprintf(stderr,"ERROR: File %s empty\n", nodeData->filesList[i]);
            exit(1);
        }
        FCLOZE(fpIn);
    }

#ifdef VERBOSE
    fprintf(threadData[0].threadLog, "Table %d Depth:%d, Width:%d\n", nodeData->id,
                                            tableData->tableSize, tableData->compSize);
#else
    threadData[0].threadID=threadData[0].threadID;
#endif
}
