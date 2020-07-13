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
#ifndef HEADR
#define HEADR

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <inttypes.h>
#include <zlib.h>
#include <dirent.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/sysinfo.h>
#include <libgen.h>

typedef uint64_t inputDataType_x64;
typedef float ResultDataType;

typedef gzFile inFileType;
#define FOPEN gzopen
#define FCLOZE gzclose

#ifdef GZIP
typedef gzFile outFileType;
#define FPRINT gzprintf
#define FOPEN_OUT gzopen
#define FCLOZE_OUT gzclose
#else
typedef FILE* outFileType;
#define FPRINT fprintf
#define FOPEN_OUT fopen
#define FCLOZE_OUT fclose
#endif

enum various_sizes_e
{
    FILEBUFFERSIZE=1024*1024*16,
    GENOMESIZE_D=2,
    BUFFERSIZE=1024,
    BUFFER_INCR=128,
    INFILENAMESIZE=1024,
    STRINGLENGTH=1024,
    IDLENGTH=256,
    NOFLINES=100
};

enum vcf_sizes_e
{
    VCF_HLENGTH=9, // number of fields in the VCF header line
    MAX_CHROM_NAME_VCF=100,
    MAX_STATES_VCF=5,
    MAXAFLENGTH=9
};

enum formats_e
{
    OTHER_FORMAT=-1,
    MS_FORMAT=0,
    FASTA_FORMAT=1,
    MACS_FORMAT=2,
    VCF_FORMAT=3
};

#define STATESALL 8
enum states_e
{
    ZERO='0',
    ONE='1',
    GAP='-',
    AD='A',
    CY='C',
    GU='G',
    TH='T',
    UN='N',
    ad='a',
    cy='c',
    gu='g',
    th='t'
};

#define min(a,b)     ( (a) > (b) ? (b) : (a) )
#define max(a,b)     ( (a) <= (b) ? (b) : (a) )

#endif
