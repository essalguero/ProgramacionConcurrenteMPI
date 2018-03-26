#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <omp.h>
#include <mpi.h>

#define XY_TO_ARRAY(x, y, width) ((x * width) + y)

float * convolutionOpenMP(float * originalImage, int imageWidth, int imageHeight, float * originalFilter, int filterWidth, int filterHeight)
{
	float * hostImage = 0;
	float * hostFilter = 0;
	float * hostImageProcessed = 0;
	float * hostImageToProcess = 0;

	float * finalImage = 0;

		
	printf("Convolution with OpenMP started\n");

	//Get memory for hostFilter
	hostFilter = (float *)malloc(filterHeight * filterWidth * sizeof(float));

	if (hostFilter == 0)
	{
		printf("Error getting memory hostFilter CPU\n");
		exit(1);
	}

	#pragma omp parallel for
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
	
	#pragma omp parallel for
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
	#pragma omp parallel for
	for (int i = 0; i < imageHeight; ++i)
		#pragma omp parallel for
		for (int j = 0; j < imageWidth; ++j)
			hostImageToProcess[XY_TO_ARRAY((i + (filterHeight / 2)), (j + (filterWidth / 2)), columnasNecesarias)] = hostImage[XY_TO_ARRAY(i, j, imageWidth)];

	// Convolution process
	// Loop through all the cells in the image but skipping the added rows and columns
	#pragma omp parallel for
	for(int i = filterHeight / 2; i < imageHeight + (filterHeight / 2); ++i)
	{
		#pragma omp parallel for
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

	printf("Convolution with OpenMP finished\n");

	// return the final image
	return finalImage;

}
// float * convolutionCPU(float * originalImage, int imageWidth, int imageHeight, float * originalFilter, int filterWidth, int filterHeight)


void master(int argc, char** argv, int rank, int nproc)
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
		printf("Usage: readFile <filename> <rows> <columns> <filterFilename> <filterRows> <filterColumns> <resultFile>\n\n");

		//return -1;
		return;
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


	printf("Number of processes: %d\n", nproc);

	// reserve needed memory for the image and the filter
	image = (float *) malloc(rows * cols * sizeof(float));
	filter = (float *) malloc(rows * cols * sizeof(float));


	// open the file storing the image
	if(!(fileOpened = fopen(filename, "r")))
	{
		printf("Could not open Image file\n");
		free(image);
		free(filter);

		//return -1;
		return;
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

		//return -1;
		return;
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


	printf("Image has been loaded. Sending jobs to slave processes\n");


	

	// reserve memory for the result of the processing
	float * imageProcessed = (float *)malloc(rows * cols * sizeof(float));

	// Display image for testing purposes
	/*for(int i = 0; i < rows * cols; i++)
		printf("%.2f\n", imageProcessed[i]);*/
	


	// Determine the number of processes to use
	int numberSlaves = (nproc - 1) % 2;
	numberSlaves = (nproc - 1) - numberSlaves;
	if (numberSlaves > 8)
		numberSlaves = 8;

	// Image has been read from file. Now the data should be distributed to the slaves in order to process it
	
	// Divide the number of rows between the number of slave processes.
	int rowsToSlave = (float)rows / (float)numberSlaves;

	if ((numberSlaves * rowsToSlave) < rows)
		 rowsToSlave++;
	printf("Number of rows to be processed per slave: %d\n", rowsToSlave);
	
	int initialPosition;
	int endPosition;

	// split the image among the number of slaves available
	for (int slave = 1; slave <= numberSlaves; slave++)
	{
		if (1 == slave)
		{
			initialPosition = 0;
		}
		else
		{
			initialPosition = (rowsToSlave * (slave - 1) - 1) * cols;
		}
		
		if (numberSlaves == slave)
		{
			endPosition = (rows * cols) - 1;
		}
		else
		{
			endPosition = (initialPosition + ((rowsToSlave + 2) * cols)) - 1;
		}

		if (1 == slave)
			endPosition -= cols;


		int totalPositions = endPosition - initialPosition + 1;

		printf("Slave %d to process begin Position %d and endPosition %d, total of %d floats\n", slave, initialPosition, endPosition, totalPositions);

		float * imageToSlave = (float *)malloc(totalPositions * sizeof(float));

		// Copy the part of the image to be processed by the current slave
		memcpy(imageToSlave, (void*)(image + initialPosition), totalPositions * sizeof(float));
		
		/*for(int i = 0; i < totalPositions; ++i)
		{
			printf("%f\n", *(image + initialPosition + i));
		}*/
		

		/*for(int i = 0; i < totalPositions; ++i)
		{
			printf("%f\n", *(imageToSlave + i));
		}*/
		
		// Send the data to the slave
		// All the parameters needed to perform the convolution
		// number of floats of the image
		// number of columns
		// the image
		// number of rows of the filter
		// number of columns of the filter
		// the filter
		MPI_Send(&totalPositions, 1, MPI_INT, slave, 0, MPI_COMM_WORLD);
		MPI_Send(&cols, 1, MPI_INT, slave, 0, MPI_COMM_WORLD);
		MPI_Send(imageToSlave, totalPositions, MPI_FLOAT, slave, 0, MPI_COMM_WORLD);

		MPI_Send(&filterRows, 1, MPI_INT, slave, 0, MPI_COMM_WORLD);
		MPI_Send(&filterCols, 1, MPI_INT, slave, 0, MPI_COMM_WORLD);
		MPI_Send(filter, filterRows * filterCols, MPI_FLOAT, slave, 0, MPI_COMM_WORLD);

		free(imageToSlave);
	}

	// receive back the processed data;
	for(int slave=1; slave <= numberSlaves; slave++)
	{
		MPI_Status status;				
		
		int result;
		
		float* subImage;
		int totalPositions;
		int initialPosition = 0;
		int endPosition = 0;

		// Receive the result of the convolution from the slaves
		// Data received should be:
		// number of floats processed by the slave
		// subImage generated by the slave
		result = MPI_Recv(&totalPositions, 1, MPI_INT, slave, 0, MPI_COMM_WORLD, &status);
		if (result != MPI_SUCCESS)
		{
			printf("Error receiving totalPositions from slave %d. Aborting execution\n", slave);
			
			free(image);
			free(filter);
			free(imageProcessed);
			
			return;
		}
		endPosition = totalPositions - 1;
	
		subImage = (float *)malloc(totalPositions * sizeof(float));
		
		result = MPI_Recv(subImage, totalPositions, MPI_FLOAT, slave, 0, MPI_COMM_WORLD, &status);
		if (result != MPI_SUCCESS)
		{
			printf("Error receiving subImage from slave %d. Aborting execution\n", slave);
			
			free(image);
			free(filter);
			free(imageProcessed);
			free(subImage);
			
			return;
		}
		/*for (int i = 0; i < totalPositions; ++i)
		{
			printf("%f\n", *(subImage + i));
		}*/
	
		// Merge the results
		if (1 != slave)
		{
			initialPosition = cols;
		}
		
		if (slave != numberSlaves)
		{
			endPosition -= cols;
		}
		
		int positionsToCopy = (endPosition - initialPosition) + 1;
		int totalImagePosition = ((slave - 1) * rowsToSlave * cols);
		printf("Slave %d to copy begin Position %d and endPosition %d, total of %d floats, position in result image: %d\n", slave, initialPosition, endPosition, positionsToCopy, totalImagePosition);

		// Copy the part of the image to be processed by the current slave
		memcpy(imageProcessed + totalImagePosition, (void*)(subImage + initialPosition), positionsToCopy * sizeof(float));
		
		// free memory
		free(subImage);
	}



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
		
		//return -1;
		return;
	}

	
	// Release the memory
	free(image);
	free(filter);
	free(imageProcessed);

}


void esclavo(int argc, char** argv, int rank, int nproc)
{

	int totalPositions;
	int cols;
	float* image;
	
	int filterRows;
	int filterCols;
	float* filter;
	
	MPI_Status status;
	int result;

	// Receive all the parameters needed to perform the convolution
	// number of floats of the image
	// number of columns
	// the image
	// number of rows of the filter
	// number of columns of the filter
	// the filter
	result = MPI_Recv(&totalPositions, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	if (result != MPI_SUCCESS)
	{
		printf("Slave %d: Error receiving totalPositions from master. Aborting execution\n", rank);
			
		return;
	}

	result = MPI_Recv(&cols, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	if (result != MPI_SUCCESS)
	{
		printf("Slave %d: Error receiving cols from master. Aborting execution\n", rank);
			
		return;
		
	}

	image = (float *)malloc(totalPositions * sizeof(float));
	
	result = MPI_Recv(image, totalPositions, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
	if (result != MPI_SUCCESS)
	{
		printf("Slave %d: Error receiving image from master. Aborting execution\n", rank);
		free(image);
			
		return;
		
	}

	// Print the info of the image to process for debugging purposes
	/*for(int i = 0; i < totalPositions; ++i)
	{
		printf("%f\n", *(image + i));
	}*/
	
	result = MPI_Recv(&filterRows, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	if (result != MPI_SUCCESS)
	{
		printf("Slave %d: Error receiving filterRows from master. Aborting execution\n", rank);
		free(image);
			
		return;
		
	}
	
	result = MPI_Recv(&filterCols, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	if (result != MPI_SUCCESS)
	{
		printf("Slave %d: Error receiving filterCols from master. Aborting execution\n", rank);
		free(image);
			
		return;
		
	}
	
	filter = (float *)malloc(filterRows * filterCols * sizeof(float));
	
	result = MPI_Recv(filter, filterRows * filterCols, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
	if (result != MPI_SUCCESS)
	{
		printf("Slave %d: Error receiving filter from master. Aborting execution\n", rank);
		free(image);
		free(filter);
			
		return;
		
	}
	
	// Print the info of the filter for debugging purposes
	/*for(int i = 0; i < filterRows * filterCols; ++i)
	{
		printf("%f\n", *(filter + i));
	}*/
	
	// Once all the data has been received. Execute the convolution process
	float* imageProcessed = convolutionOpenMP(image, cols, totalPositions / cols, filter, filterCols, filterRows);
	
	// Send the data back to the master
	MPI_Send(&totalPositions, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	MPI_Send(imageProcessed, totalPositions, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
	
	// Print the info of the processed image for debugging purposes
	/*for(int i = 0; i < totalPositions; ++i)
	{
		printf("%f\n", *(imageProcessed + i));
	}*/
	
	// release the nenmory used
	free(image);
	free(filter);
	free(imageProcessed);


	printf("Finished slave: %d\n\n", rank);
}


int main(int argc, char **argv)
{
	int rank;
	int nproc;
	
	int returnValue = 0;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nproc);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
	
	// only start an even number of slaves less or equal to 8
	int numberSlaves = (nproc - 1) % 2;
	numberSlaves = (nproc - 1) - numberSlaves;
	if (numberSlaves > 8)
		numberSlaves = 8;

	// Check all needed parameters are provided
	if (argc == 8)
	{
		
		if (rank <= numberSlaves)
		{
			switch(rank)
			{
				case 0: printf("Number of processors: %d\n", numberSlaves);
					printf("Rank number: %d\n", rank);
					master(argc, argv, rank, numberSlaves + 1);
					break;
				default:
					esclavo(argc, argv, rank, numberSlaves + 1);
					printf("End slave: %d\n", rank);
					break;
			};
		}
	}
	else
	{
		if (0 == rank)
		{
			printf("Number of parameters is not correct\n");
			printf("Usage: mpirun convolutionMPI <filename> <rows> <columns> <filterFilename> <filterRows> <filterColumns> <resultFile>\n\n");
		}
		
		returnValue = -1;
	}
	MPI_Finalize();

	return returnValue;;

}
