#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define XY_TO_ARRAY(x, y, width) ((x * width) + y)

float * convolutionCPU(float * originalImage, int imageWidth, int imageHeight, float * originalFilter, int filterWidth, int filterHeight)
{
	float * hostImage = 0;
	float * hostFilter = 0;
	float * hostImageProcessed = 0;
	float * hostImageToProcess = 0;

	float * finalImage = 0;

		
	printf("Convolution on CPU started\n");

	//Get memory for hostFilter
	hostFilter = (float *)malloc(filterHeight * filterWidth * sizeof(float));

	if (hostFilter == 0)
	{
		printf("Error getting memory hostFilter CPU\n");
		exit(1);
	}

	for (int i = 0; i < filterHeight * filterWidth; ++i)
	{
		hostFilter[i] = originalFilter[i];
	}


	// obtener memoria para las imagenes

	// Imagen original
	hostImage = (float *)malloc(imageWidth * imageHeight * sizeof(float));


	if (hostImage == 0)
	{
		printf("Error getting memory hostImage CPU\n");
		exit(1);
	}

	for (int i = 0; i < imageWidth * imageHeight; ++i)
	{
		hostImage[i] = originalImage[i];
	}

	// Reservar memoria para la imagen resultante
	finalImage = (float *)malloc(imageWidth * imageHeight * sizeof(float));

	if (finalImage == 0)
	{
		printf("Error getting memory finalImage CPU\n");
		exit(1);
	}

	// Imagen procesada. Tiene que ser mas grande para que durante el procesado no haya problemas de indices
	// Calcular el numero de filas necesarias
	int filasNecesarias = imageHeight + (2 * (filterHeight / 2));
	if (filasNecesarias % filterHeight)
		filasNecesarias += filterHeight - (filasNecesarias % filterHeight);

	// Calcular el numero de columnas necesarias
	int columnasNecesarias = imageWidth + (2 * (filterWidth / 2));
	if (columnasNecesarias % filterWidth)
		columnasNecesarias += filterWidth - (columnasNecesarias % filterWidth);

	//reservar memoria
	hostImageToProcess = (float *)malloc(filasNecesarias * columnasNecesarias * sizeof(float));

	if (hostImageToProcess == 0)
	{
		printf("Error getting memory hostImageToProcess CPU\n");
		exit(1);
	}

	hostImageProcessed = (float *)malloc(filasNecesarias * columnasNecesarias * sizeof(float));

	if (hostImageProcessed == 0)
	{
		printf("Error getting memory hostImageProcessed CPU\n");
		exit(1);
	}

	// Set all the pixels in the image to 0
	memset(hostImageToProcess, 0, filasNecesarias * columnasNecesarias * sizeof(float));

	//Copy the image to process to a bigger array. This avoids problems when processing the image edges with the convolution matrix
	for (int i = 0; i < imageHeight; ++i)
		for (int j = 0; j < imageWidth; ++j)
			hostImageToProcess[XY_TO_ARRAY((i + (filterHeight / 2)), (j + (filterWidth / 2)), columnasNecesarias)] = hostImage[XY_TO_ARRAY(i, j, imageWidth)];

	// Convolution process
	// Loop through all the cells in the image but skipping the added rows and columns
	for(int i = filterHeight / 2; i < imageHeight + (filterHeight / 2); ++i)
	{
		for(int j = filterWidth / 2; j < imageWidth + (filterWidth / 2); ++j)
		{
			float result = 0;
			
			// Loop through all the cells in the filter
			for (int y = 0; y < filterHeight; y++)
			{
				for(int x = 0; x < filterWidth; x++)
				{
					result += hostFilter[XY_TO_ARRAY(y, x, filterWidth)] * 
						hostImageToProcess[XY_TO_ARRAY((i - ((filterHeight / 2) - y)), (j - ((filterWidth / 2) - x)), columnasNecesarias)];
				}
			}
			
			hostImageProcessed[XY_TO_ARRAY(i, j, columnasNecesarias)] = result;
			finalImage[XY_TO_ARRAY((i - (filterHeight / 2)), (j - (filterWidth / 2)), imageWidth)] = result;
		}
	}
	

	//release the memory allocated
	free(hostFilter);
	free(hostImage);
	free(hostImageToProcess);
	free(hostImageProcessed);

	printf("Convolution on CPU finished\n");

	// return the final image
	return finalImage;

}
// float * convolutionCPU(float * originalImage, float * originalFilter)


int main(int argc, char **argv)
{
	char filename[1024];
	int rows;
	int cols;

	char filterFilename[1024];
	int filterRows;
	int filterCols;


	char resultFilename[1024];

	float* image;
	float* filter;

	FILE * fileOpened;

	char rowString[30];
	
	// Check all needed parameters are provided
	if (argc != 8)
	{
		printf("Number of parameters is not correct\n");
		printf("Usage: convolutionCPU <filename> <rows> <columns> <filterFilename> <filterRows> <filterColumns> <resultFile>\n\n");

		return -1;
	}

	// Get the data for the image to process
	strcpy(filename, argv[1]);
	rows = atoi(argv[2]);
	cols = atoi(argv[3]);

	printf("Filename to process: %s, number of Rows: %d, number of Columns: %d\n", filename, rows, cols);

	// Get the data for the filter to use
	strcpy(filterFilename, argv[4]);
	filterRows = atoi(argv[5]);
	filterCols = atoi(argv[6]);

	printf("Filter to use: %s, number of Rows: %d, number of Columns: %d\n", filterFilename, filterRows, filterCols);


	// Result file
	strcpy(resultFilename, argv[7]);

	printf("Filter to use: %s, number of Rows: %d, number of Columns: %d\n", resultFilename, rows, cols);


	// reserve needed memory for the image and the filter
	image = (float *) malloc(rows * cols * sizeof(float));
	filter = (float *) malloc(rows * cols * sizeof(float));


	// open the file storing the image
	if(!(fileOpened = fopen(filename, "r")))
	{
		printf("Could not open Image file\n");
		free(image);
		free(filter);
		return -1;
	}

	// read Image from file
	int counter = 0;	
	while(!feof(fileOpened))
	{
		fgets(rowString, 30, fileOpened);
		
		image[counter++] = atof(rowString);
	}

	// Close the file
	fclose(fileOpened);

	// Display image for testing purposes
	/*for(int i = 0; i < rows * cols; i++)
		printf("%.2f\n", image[i]);*/


	// read Filter from file
	if(!(fileOpened = fopen(filterFilename, "r")))
	{
		printf("Could not open Filter file\n");
		free(image);
		free(filter);
		return -1;
	}

	counter = 0;	
	while(!feof(fileOpened))
	{
		fgets(rowString, 30, fileOpened);
		
		filter[counter++] = atof(rowString);
	}

	// Close the file
	fclose(fileOpened);

	// Display filter for testing purposes
	/*for(int i = 0; i < filterRows * filterCols; i++)
	{
		printf("%.2f\n", filter[i]);
	}*/


	// reserve memory for the result of the processing
	float * imageProcessed = convolutionCPU(image, cols, rows, filter, filterCols, filterRows);

	// Display image for testing purposes
	/*for(int i = 0; i < rows * cols; i++)
		printf("%.2f\n", imageProcessed[i]);*/
	
	//Open a file to store the result of the processing
	if((fileOpened = fopen(resultFilename, "w")))
	{
		// Save the result to file
		for(int i = 0; i < rows * cols; ++i)
		{
			fprintf(fileOpened, "%f\n", imageProcessed[i]);
		}

		// Close the file containing the result of the processing
		fclose(fileOpened);
	}
	else
	{
		// The file to store the data could not be opened. Show an error
		printf("Processed image could not be saved. Error opening file\n\n");
		
		// Release the memory
		free(image);
		free(filter);
		free(imageProcessed);
		
		return -1;
	}

	
	// Release the memory
	free(image);
	free(filter);
	free(imageProcessed);

	return 0;
}
