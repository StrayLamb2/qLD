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
#ifndef RF_H
#define RF_H

#include "header.h"
#include "pthreads.h"

double gettime(void);

int	getNextLine(inFileType fpIn,
                char ** line,
                int *readEOL,
                int *readEOF,
                int *lineLength);

int	getWordFromString(char* line,
                      char ** word,
                      int *readEOL,
                      int *wordLength,
                      int *index);
#endif
