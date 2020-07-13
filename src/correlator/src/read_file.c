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
#include "../include/correlate_IO.h"
#include "../include/correlate_MEM.h"
#include "../include/read_file.h"

/*
returns current time as double in seconds
*/
double gettime(void)
{
    struct timeval ttime;
    gettimeofday(&ttime, NULL);
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
}

int getNextLine(inFileType fpIn, char **line, int *readEOL, int *readEOF, int *lineLength)
{
    *readEOL=*readEOF=0;
    if(gzgets(fpIn, *line, *lineLength) != NULL)
    {
        if((*line)[strlen(*line)-1] != '\n')
        {
            while((*line)[strlen(*line)-1] != '\n')
            {
                char *tmpLine=malloc((*lineLength)*sizeof(char));
                if(gzgets(fpIn,tmpLine,(*lineLength)>>1) != NULL)
                {
                    (*lineLength)=(*lineLength) << 1;
                    (*line)=realloc_buff((void *)(*line),(*lineLength),"char");
                    strcat((*line), tmpLine);
                }
                else
                {
                    (*line)[strlen(*line)]='\0';
                    *readEOF=1;
                    if(strlen(*line) >= 3)
                    {
                        *readEOL=1;
                        free(tmpLine);
                        return 1;
                    }
                    free(tmpLine);
                    return 0;
                }
                free(tmpLine);
            }
        }
        *readEOL=1;
        (*line)[strlen(*line)-1]='\0';
        if(strlen(*line) >= 3)
            return 1;
        else
            return 0;
    }
    else
    {
        *readEOF=1;
        return 0;
    }
    return 1;
}

int getWordFromString(char* line, char ** word, int *readEOL, int *wordLength, int *index)
{
    *readEOL=0;
    int i=0;
    char ent=line[(*index)++];
    while(ent == ' '|| ent == 9) // horizontal tab
        ent=line[(*index)++];
    if(ent == '\0')
    {
        *readEOL=1;
        (*word)[0]='\0';
        return 0;
    }

    while(!(ent == 9 || ent == 32) && ent != '\0')
    {
        if(i+1 >= (*wordLength))
        {
            (*wordLength)=(*wordLength) << 1;

            char* tmp=realloc((*word), (*wordLength)*sizeof(char));
            if(tmp)
            {
                *word=tmp;
            }
            else
            {
                printf("ERROR: word not found at getWordFromString\n");
                exit(1);
            }
        }
        (*word)[i++]=ent;
        ent=line[(*index)++];
    }
    (*word)[i]='\0';
    if(ent == 10 || ent == 13 || ent == '\0')
        *readEOL=1;
    return 1;
}
