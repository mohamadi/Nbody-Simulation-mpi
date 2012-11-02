/*
 	CPSC 521: Parallel Algorithms and Architectures
	Assignment 1: N-body Simulation on a Ring
	Author: Hamid Mohamadi, mohamadi@alumni.ubc.ca
*/

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>

typedef struct{
	double x,y;
	double m;
} body;

typedef struct{
	double x, y;
} fvp; // data type for force and velocity

void frcUpdt(body*, body, fvp*);
void posUpdt(body*, fvp, fvp*);
void scatter(const char *, int, int, body *);
void simulate(const int, int, int, body *);
void gather(int, int, body);

int main(int argc, char *argv[]){    
	int tmax=atoi(argv[1]); /* number of simulation rounds */
	const char *pName=argv[2];
	int rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	body bodies;
	
	/****************************************************************
	1. Initialization:
	Reading the input from file and assigning each body to a process
	****************************************************************/
	scatter(pName, rank, size, &bodies);
	
	/****************************************************************
	2. Simulation:
	Passing the data of each process to other processes and simulate
	****************************************************************/
	simulate(tmax, rank, size, &bodies);

	/****************************************************************
	3.Termination: 
	Writing the the last positions to nbody_out.txt 
	****************************************************************/
	gather(rank, size, bodies);
		
	
	/* Finalizing */
	MPI_Finalize();
	return(0);
}

void frcUpdt(body *bodies, body rcvBodies, fvp *f){
	int i,j;
	double r=0.0;
	const double  G =  0.0000000000667384;
	r = sqrt((rcvBodies.x-bodies->x)*(rcvBodies.x-bodies->x)+
			 (rcvBodies.y-bodies->y)*(rcvBodies.y-bodies->y));
	f->x += G*bodies->m*rcvBodies.m*(rcvBodies.x-bodies->x)/(r*r*r);
	f->y += G*bodies->m*rcvBodies.m*(rcvBodies.y-bodies->y)/(r*r*r);
}

void posUpdt(body *bodies, fvp f, fvp *v){
	int i;
	const double dt= 1.0;
	v->x += f.x * dt / bodies->m;
	v->y += f.y * dt / bodies->m;
	bodies->x += v->x * dt;
	bodies->y += v->y * dt;
}

void scatter(const char *pName, int rank, int size, body *bodies){
	int i;
	int sendto = (rank + 1) % size;
	int recvfrom = ((rank + size) - 1) % size;
	
	MPI_Datatype bodytype;
	MPI_Type_contiguous(3, MPI_DOUBLE, &bodytype);
	MPI_Type_commit(&bodytype);
	MPI_Status status;
	
	body outbuf;
	
	if(rank==0){
		FILE *pFile;
		pFile = fopen(pName, "rb");
		fscanf(pFile,"%lf %lf %lf", &outbuf.x, &outbuf.y, &outbuf.m);
		*bodies=outbuf;
		for(i=0; i<size-rank-1;i++){
			fscanf(pFile,"%lf %lf %lf", &outbuf.x, &outbuf.y, &outbuf.m);
			MPI_Send(&outbuf, 1, bodytype, sendto, 0, MPI_COMM_WORLD);
		}
		fclose(pFile);
	}	
	else{
		MPI_Recv(&outbuf, 1, bodytype, recvfrom, 0, MPI_COMM_WORLD, &status);
		*bodies=outbuf;
		for(i=0; i<size-rank-1;i++){
			MPI_Recv(&outbuf, 1, bodytype, recvfrom, 0, MPI_COMM_WORLD, &status);
			MPI_Send(&outbuf, 1, bodytype, sendto, 0, MPI_COMM_WORLD);
		}
	}	
}

void simulate(const int tmax, int rank, int size, body *bodies){
	int t=tmax, i, round;
	int sendto = (rank + 1) % size;
	int recvfrom = ((rank + size) - 1) % size;
	
	MPI_Datatype bodytype;
	MPI_Type_contiguous(3, MPI_DOUBLE, &bodytype);
	MPI_Type_commit(&bodytype);
	MPI_Status status;
	body inbuf, outbuf;
	fvp f,v;
	v.x=0; v.y=0;
	t=tmax;
	while(1){
		--t;
		round=size;
		outbuf=*bodies;
		f.x = 0; f.y=0;
		while (round > 1) {
			--round;			
			if (!(rank % 2)){
				MPI_Send(&outbuf, 1, bodytype, sendto, 0, MPI_COMM_WORLD);
				MPI_Recv(&inbuf, 1, bodytype, recvfrom, 0, MPI_COMM_WORLD, &status);			
			}
			else
			{
				MPI_Recv(&inbuf, 1, bodytype, recvfrom, 0, MPI_COMM_WORLD, &status);
				MPI_Send(&outbuf, 1, bodytype, sendto, 0, MPI_COMM_WORLD);
			}
			outbuf=inbuf;
			frcUpdt(bodies, inbuf, &f);
		}
		posUpdt(bodies, f, &v);
		if(t==0) break;
	}
}

void gather(int rank, int size, body bodies){
	int i;
	int sendto = (rank + 1) % size;
	int recvfrom = ((rank + size) - 1) % size;
	
	MPI_Datatype bodytype;
	MPI_Type_contiguous(3, MPI_DOUBLE, &bodytype);
	MPI_Type_commit(&bodytype);
	MPI_Status status;
	
	body outbuf;
	
	if (rank != 0){
		outbuf = bodies;
		MPI_Send(&outbuf, 1, bodytype, recvfrom, 0, MPI_COMM_WORLD);
		for(i=0;i<size-rank-1;i++){
			MPI_Recv(&outbuf, 1, bodytype, sendto, 0, MPI_COMM_WORLD, &status);
			MPI_Send(&outbuf, 1, bodytype, recvfrom, 0, MPI_COMM_WORLD);
		}
	}	
	else {
		FILE *oFile;
		oFile = fopen("nbody_out.txt", "w");
		outbuf=bodies;
		fprintf(oFile,"%15.10f %15.10f %15.10f\n", outbuf.x, outbuf.y, outbuf.m);
		for(i=0; i<size-rank-1;i++){
			MPI_Recv(&outbuf, 1, bodytype, sendto, 0, MPI_COMM_WORLD, &status);
			fprintf(oFile,"%15.10f %15.10f %15.10f\n", outbuf.x, outbuf.y, outbuf.m);	
		}
		fclose(oFile);
	}	
}
