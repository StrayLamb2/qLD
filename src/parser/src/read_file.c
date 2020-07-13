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
#include "header.h"

double gettime(void)
{
    struct timeval ttime;
    gettimeofday(&ttime , NULL);
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
}

int getNextWord(gzFile fp, char * word, int *readEOL, int *readEOF, int *wordLength)
{
    *readEOL = *readEOF = 0;
    
    int i=0;
    word[i] = gzgetc(fp);
    
    if(word[i] == EOF)
    {
        *readEOF = 1;
        word[i] = '\0';
        return 0;
    }
    
    if(word[i] == 10 || word[i] == 13)
    {
        *readEOL = 1;
        word[i] = '\0';
        return 0;
    }
    while(!(word[i] == 9 || word[i] == 32) && 
          !(word[i] == 10 || word[i] == 13)  && 
          word[i] != EOF)
    {
        if(i+1 >= (*wordLength) )
        {
            (*wordLength) = (*wordLength) << 1; 
            
            word = (char*)realloc( word, (*wordLength) * sizeof(char) );       
        }
        
        i++;
        word[i] = gzgetc(fp);
    }
    
    if(word[i] == 9 || word[i] == 32) 
    {
        while(word[i] == 9 || word[i] == 32)
        {
            if(i+1 >= (*wordLength) )
            {
                (*wordLength) = (*wordLength) << 1; 
                
                word = (char*)realloc( word, (*wordLength) * sizeof(char) );       
            }
            
            i++;
            word[i] = gzgetc(fp);                
        }
    }
    
    if(word[i] == 10 || word[i] == 13)
        *readEOL = 1;
    else
        gzungetc(word[i],fp);  
    
    word[i] = '\0';   
    return 1;    
}

int getNextLine(gzFile fpIn, char ** line, int *readEOL, int *readEOF, int *lineLength)
{
    *readEOL = *readEOF = 0;

    if(gzgets(fpIn, *line, *lineLength) != NULL) 
    {
        if((*line)[strlen(*line)-1] != '\n') 
        {
            while((*line)[strlen(*line)-1] != '\n')
            {
                char * tmpLine = (char*)malloc( (*lineLength) * sizeof(char));
                (*lineLength) = (*lineLength) << 1; 
                (*line) = (char*)realloc( (*line), (*lineLength) * sizeof(char) );
                if(gzgets(fpIn, tmpLine, (*lineLength)>>1) != NULL) 
                {
                    strcat((*line), tmpLine);
                }
                else 
                {
                    (*line)[strlen(*line)] = '\0';
                    *readEOF = 1;
                    if(strlen(*line) >= 3) 
                    {
                        *readEOL = 1;
                        return 1;
                    }
                    free(tmpLine);
                    return 0;
                }
                free(tmpLine);
            }
        }      
        *readEOL = 1;
        (*line)[strlen(*line)-1] = '\0';
        if(strlen(*line) >= 3) 
            return 1;
        else
            return 0;
    }
    else 
    {
        *readEOF = 1; 
        return 0; 
    }
    return 1;
}

void writeLine(gzFile fpOut,char **line)
{
    if(gzputs(fpOut,*line) < 0) 
    {
        printf("error writing to file");
        assert(0);
    }
    gzputc(fpOut,'\n');
}
//not really needed
int countLines(gzFile fpIn, char **line, int* lineLength)
{
  	int lines=0, notHeader = 0;
  	while(gzgets(fpIn, *line, *lineLength) != NULL)
	{
		if((*line)[strlen(*line)-1] != '\n') 
        {
			while((*line)[strlen(*line)-1] != '\n')
            {
	            char * tmpLine = (char*)malloc( (*lineLength) * sizeof(char));
                (*lineLength) = (*lineLength) << 1; 
                (*line) = (char*)realloc( (*line), (*lineLength) * sizeof(char) );
                if(gzgets(fpIn, tmpLine, (*lineLength)>>1) != NULL) 
                {
                    strcat((*line), tmpLine);
                }
                else 
                {
                    int tmp = strlen(*line);
                    (*line)[tmp] = '\n';
                    (*line)[tmp+1] = '\0';
                }
	            free(tmpLine);
	        }
		}
        if(strlen(*line) >=6 ) 
        {
            if(notHeader ==1)
                lines++;

            if((*line)[0] == '#' && 
                             (*line)[1] == 'C' && 
                             (*line)[2] == 'H' && 
                             (*line)[3] == 'R' && 
                             (*line)[4] == 'O' && 
                             (*line)[5] == 'M')
                notHeader = 1;
        }
	}
	return lines;
}

int getWordFromString(char* line, char ** word, int *readEOL, int *wordLength, int *index)
{
    *readEOL = 0;
    
    int i=0;
    char ent = line[(*index)++];
    while(ent==' '|| ent == 9) // horizontal tab
        ent = line[(*index)++]; 
    if(ent == '\0')
    {
        *readEOL = 1;
        (*word)[0] = '\0';
        return 0;
    }
    
    while(!(ent == 9 || ent == 32) && ent != '\0')
    {
        
        if(i+1 >= (*wordLength) )
        {
            (*wordLength) = (*wordLength) << 1; 
            
            (*word) = (char*)realloc( (*word), (*wordLength) * sizeof(char) );       
        }
        
        (*word)[i++] = ent;
        
        
        ent = line[(*index)++];
    }
    
    (*word)[i] = '\0';
    
    if(ent == 10 || ent == 13 || ent =='\0')
        *readEOL = 1;
    
    return 1;
}

int cmpfunc (const void * a, const void * b)
{
   return ( *(int*)a - *(int*)b );
}

int sortList(char* inputListName,int** list, char** line, int *lineLength)
{
    gzFile fpIn=NULL;
    fpIn = gzopen(inputListName,"r");
    int status,eol=0,eof=0,pos,listItems = 0;
    while(eof == 0) 
    {
        status =  getNextLine(fpIn, line, &eol, &eof, lineLength);
        if(status == 1 && eol==1) 
        {    
            pos = atoi(*line);
            listItems++;
            *list = (int*)realloc(*list,sizeof(int)*listItems);
            (*list)[listItems-1] = pos;
        }
    }
    qsort(*list, listItems, sizeof(int), cmpfunc);
    return listItems;
}
