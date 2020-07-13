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
#include "../include/fast_print.h"

unsigned int reverse(long *num) 
{ 
    unsigned int digits=0;
    long rev_num=0; 
    while((*num) > 0) 
    { 
        digits++;
        rev_num=rev_num*10+(*num)%10; 
        (*num)/=10; 
    } 
    (*num)=rev_num;
    return digits;  
} 

void fprint_f(FILE* fp, float val, int precision)
{
    unsigned int j;

    //(int)pow(10,digits);
    long unsigned power10=1;      
    for(j=0; j<(unsigned int)precision; j++) 
        power10*=10;

    //Helper Array for our conversion
    char d_arr[10]="0123456789";

    //val x power10 gives us a long representation of our float
    long long_val=val*power10; 
    //the integer part of the float (using long to prevent overflow)
    long int_part=(long)val;
    //long_val now holds just the decimals
    long_val=long_val-int_part*power10;

    //we reverse the numbers, since we are gonna fill the buffer backwards,
    //keeping the count of the digits of each part to zero fill later,
    //since the method ignores leading and trailing zeros
    unsigned int i_digits=reverse(&int_part);
    unsigned int d_digits=reverse(&long_val);
   
    //if there is no int part, print a zero
    if(int_part==0)
    { 
        fputc('0', fp);
    }
    else
    {
        //keep count of the non-zero int part to zero fill with trailing zeros later
        unsigned int i_numbers=0;
        //print the integer non-zero part
        while (int_part != 0)
        {
            fputc(d_arr[int_part%10], fp);
            int_part=int_part/10;
            i_numbers++;
        }
        //zero-fill trailing zeros
        for(j=0; j<i_digits-i_numbers; j++)
        {
            fputc('0', fp);
        } 
    }
    //Could respect locale here, but nvm for now
    fputc('.',fp);
    //if there is no decimal part, print nothing
    if(long_val == 0)
    {
        return;
    }
    else
    {
        //zero-fill leading zeros
        for(j=0; j<precision-d_digits; j++)
        {
            fputc('0', fp);
        }
        //print the decimal non-zero part
        for(j=0; j<d_digits; j++)
        {
            fputc(d_arr[long_val%10], fp);
            long_val=long_val/10;
        }
    }
}
