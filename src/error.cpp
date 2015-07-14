#include "error.h"

const char* ERROR_LINE_INTRO = "***** ERROR *****";

int reportError(int num, char* error, char* descr, char* file, int line=-1)
{
	printf("%s ", ERROR_LINE_INTRO);
	if(-1 < num)	printf("nr: %d; ",num);
	if(error)		printf("%s; ", error);
	if(descr)		printf("description: %s; ", descr);
	if(file)		printf("in file: %s ", file);
	if(-1 < line)   printf("line: %d", line);
	
	return 1;
}
