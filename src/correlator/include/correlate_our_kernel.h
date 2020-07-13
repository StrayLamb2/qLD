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
#ifndef COK_H
#define COK_H

long get_corecount(void);

int get_blis(void);

void printHelp(void);

void commandLineParser(int argc,
        char** argv,
        char * inPath,
        char * inPath2,
        int * inPath2set,
        char * outfile,
        char **sampleListName,
        char **sampleListName2,
        int *ploidy,
        int * posWmin1,
        int * posWmax1,
        int *posWset1,
        int * posWmin2,
        int * posWmax2,
        int *posWset2,
        char * inList,
        int * inListSet,
        float * r2limit,
        int *gpu,
        int *blis,
        int *table,
        int *threads,
        int *sorted,
        int *mdf);
#endif
