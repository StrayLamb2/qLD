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
#include "../include/structs.h"
#include "../include/cl_args.h"
#include "../include/io.h"
#include "../include/mem.h"
#include "../include/preprocess_data.h"
#include "../include/process_data.h"

int main(int argc, char** argv)
{   
    pre_t *preData=init_preprocess_struct();
    header_t *headerData=init_header_struct();

    arg_parser(argc, argv, preData);
    preprocess_data(preData, headerData);

    //print_preprocess_struct(preData);
    //print_header_struct(headerData);

    unsigned int tableSize1_s=BUFFERSIZE;

    table_x64 *tableData=create_table_x64(tableSize1_s);

    helper_t *helperData=create_helper_t(preData->ploidy);

    unsigned int compSize=(headerData->valid_count/(sizeof(inputDataType_x64)*8)+
            (headerData->valid_count%(sizeof(inputDataType_x64)*8)!=0?1:0));

    tableData->SNPtable=(inputDataType_x64*)realloc_buff(
                                (void *)tableData->SNPtable,
                                tableSize1_s*compSize,
                                "inputDataType_x64");
    assert(tableData->SNPtable);

    tableData->tableIndex=0;
    tableData->compSize=compSize;
    tableData->compIndex=0;

    helperData->snipCharIndex=0;

    readTable_x64(preData,
                  headerData,
                  tableData,
                  helperData);

    //print_header_struct(headerData);

    free_table_x64(tableData);
    free_helper_t(helperData);
    free_preprocess_struct(preData);
    free_header_struct(headerData);
    return 0;
}
