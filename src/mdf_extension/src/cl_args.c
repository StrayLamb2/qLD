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
#include "../include/cl_args.h"
#include "../include/io.h"
#include "../include/help_menu.h"

int arg_parser(int argc, char **argv,  pre_t *preData)
{
    int i;
    for(i=0; i < argc; ++i)
    {
        if(!strcmp(argv[i], "-input"))
        {
            if(++i < argc)
            {
                getDir(preData->input, argv[i]);
            }
            else
            {
                fprintf(stderr, "ERROR: Insert input directory path after \"-input\" flag\n");
                exit(1);
            }
            DIR* directory=opendir(preData->input);

            if(!directory)
            {
                fprintf(stderr, "ERROR: Input directory \"%s\" does not exist\n", preData->input);
                exit(1);
            }

            closedir(directory);
        }
        else if(!strcmp(argv[i], "-output"))
        {
            struct stat st = {0};
 
            if(++i < argc)
            {
                if (stat(argv[i], &st) == -1) 
                {
                    mkdir(argv[i], 0700);
                }
                getDir(preData->output,argv[i]);
            }
            continue;
        }
        else if(!strcmp(argv[i], "-seed"))
        {
            srand(atoi(argv[i+1]));
        }
        else if(!strcmp(argv[i], "-impute"))
        {
            preData->impute=1;
        }
        else if(!strcmp(argv[i], "-sampleList"))
        {
            if(++i < argc)
            {
                getFilter(argv[i], &(preData->sampleList), &(preData->sampleListSize));
            }
            else
            {
                fprintf(stderr, "ERROR: Insert sample list file path after \"-sampleList\" flag\n");
                exit(1);
            }
        }
        else if(!strcmp(argv[i], "-ploidy"))
        {
            if(++i < argc)
            {
                if(!strcmp(argv[i],"haploid"))
                    preData->ploidy=1;
                else if(!strcmp(argv[i],"phased_diploid"))
                    preData->ploidy=2;
                else if(!strcmp(argv[i],"unphased_diploid"))
                    preData->ploidy=2;
                else
                {

                }
            }
            else
            {
                fprintf(stderr, "ERROR: Insert valid ploidy after \"-ploidy\" flag\n");
                fprintf(stderr, "Valid ploidies: [haploid, phased_diploid, unphased_diploid]\n");
                exit(1);
            }
        }
        else if(!strcmp(argv[i], "-help") ||
                !strcmp(argv[i], "-h") ||
                !strcmp(argv[i], "--help"))
        {
            printHelp();
            exit(0);
        }
    }

    if(!strcmp(preData->input,""))
    {
        fprintf(stderr, "ERROR: Please specify input file with \"-input\"\n");
        exit(1);
    }
    if(!strcmp(preData->output,""))
    {
        fprintf(stderr, "ERROR: Please specify output file with \"-output\"\n");
        exit(1);
    }
    return 0;
}
