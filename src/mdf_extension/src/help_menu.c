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
#include "../include/help_menu.h"
#include "../include/header.h"

void printHelp(void)
{
printf("qLD-parse-2MDF manual\n");
printf("---------------------\n");
    printf("\t-input input_Directory\n");
    printf("\t-output output_Directory\n");
    printf("\t-sampleList input_File\n");
    printf("\nDescription:\n");
    printf("\t-input <STRING>\t\tSpecifies the path of the input alignment parsed files\n");
    printf("\t-output <STRING>\tSpecifies the path of the output alignment directory.\n");
    printf("\t-sampleList <STRING>\ttxt file with Format:\n"
            "\t\t\t\t\"sample1\n"
            "\t\t\t\t sample2\n"
            "\t\t\t\t sample3\n"
            "\t\t\t\t   ...\n"
            "\t\t\t\t sampleN\"\n"
            "\t\t\t\tSpecifies the name of the file that includes\n"
            "\t\t\t\ta list of valid samples from the input\n"
            "\t\t\t\tthat will be selected for processing\n\n");
    printf("\n\n");

}
