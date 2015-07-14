/*! \file dlplpfileconverter.cpp
 \brief
 
 */

#include <stdio.h>
#include <iostream>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include "polytope.h"
#include "main.h"

globalArgs_t globalArgs;

static const char *optString = "vh?";

//using namespace std;

static const struct option longOpts[] = {
    { "version", no_argument, NULL, 0 },
    { "test", required_argument, NULL, 0},
    { NULL, no_argument, NULL, 0 }
};

void print_help(void)
{
	printf( "no manual yet...\n");
    
}

int detectFileType(char* Buffer, FILE* file)
{
    getNextLine(Buffer, file);
    
    switch(Buffer[0])
    {
        case 'V':
            return 1;
            break;
        case 'H':
            return 2;
            break;
        case 'N':
            return 3;
    }
    
    return 0;   // file type unknown
}

int getNumberOfVertices(char* Buffer, FILE *file, int fileType)
{
    if( 1 != fileType && 3 != fileType)
    {
        std::cout << "unmatching file type for getNumberOfVertices/n";
        return 0;
    }
    getNextLine(Buffer, file);
    return atoi(Buffer);
}

int main(int argc, char* argv[])
{
    int opt = 0;
    int longIndex;
    
	globalArgs.inputFiles = NULL;
	globalArgs.numInputFiles = 0;
    
    char Buffer[INPUT_BUFFER_SIZE];
    
    opt = getopt_long(argc, argv, optString, longOpts, &longIndex );
	while( -1 != opt)
	{
		switch(opt)
		{
			/*case 'v':
				globalArgs.verbosity++;
				//printf("verbose\n");
                break;*/
			case 'h':	// intentional fall-through
			case '?':
				print_help();
				return 1;
                break;
			case 0: 	// long options without short options
				if( strcmp( "version", longOpts[longIndex].name ) == 0 ) {
                    printf("version0.0.0\n");
                }
                else if( strcmp( "test", longOpts[longIndex].name ) == 0 ) {
                    int a = atoi(optarg);
                    printf("test back: %d \n",a);
                }
                break;
            default:
                break;
		}
		opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
	}
	globalArgs.inputFiles = argv + optind;
	globalArgs.numInputFiles = argc - optind;
    
    for(int i=0; i<globalArgs.numInputFiles; ++i)
    {
        FILE *fileIn, *fileOut;
        char outputFileName[256];
        
        strcpy(outputFileName,globalArgs.inputFiles[i]);
        strcat(outputFileName, "CONCAT");
        
        std::cout << std::endl << "now converting " << globalArgs.inputFiles[i]  <<
        "to " << outputFileName << std::endl;
        
        // open file
        fileIn = fopen(globalArgs.inputFiles[i], "r");
        
        int fileType = detectFileType(Buffer, fileIn);
        if(!fileType)
        {
            std::cout << "unknown file type, can't convert this file\n";
            fclose(fileIn);
            continue;
        }
        if(2 == fileType)
        {
            std::cout << "h-type conversion not yet implemented, can't convert this file\n";
            fclose(fileIn);
            continue;
        }
        
        //
        int numberOfVertices = getNumberOfVertices(Buffer, fileIn, fileType);
        std::cout << "numV" << numberOfVertices << std::endl;
        
        
        // open output file
        //fileOut = fopen(
        
        fclose(fileIn);
	}
}

