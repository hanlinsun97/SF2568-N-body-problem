#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "list.h"



// Structure defining a node in the tree
struct node{
	// In 2D each node has up to 4 children
	struct node* children[4];
	// They are centered in : 
	double x_center, y_center;
	// They are clustered in :
	double x_min, x_max, y_min, y_max;
	// N_particle is the number of particles in the subcube
	// List Particles is a list containing the number "i" of the particles contained in the subcube
	List Particles;
	// Last Particle in the list is
	List EndParticle;
	int N_particle;
	
	
	// Pointer to the root
	struct node* root;
};

typedef struct node node;




// Function for producing a random number between two double values
double frand(double min, double max)
{
    return min + (max - min) * rand() / RAND_MAX;
}

// Function to distribute randomly the initial positions of the particles in a 2D square limited by (x_max, y_max, x_min, y_min). 
void random_position_distribution(double x_min, double x_max, double y_min, double y_max, double* x, double* y, int N){
	int i;
	for(i=0;i<N;i++){
		y[i] = frand(y_min,y_max);
		x[i] = frand(x_min,x_max);
	}
}


void tree_initialization(node *father, double* x, double* y){
	int i, j;
	double x_half, y_half, x_i, y_i,x_MAX,y_MAX,y_MIN,x_MIN;
	int N_particle = father->N_particle;
		
	List Particles = father->Particles;
	// Memory allocation
	for(i=0 ; i<4 ; i++){
		father->children[i]=calloc(1, sizeof(*(father->children[i])));	
		(father->children[i])->Particles=NULL;
	}
	node *child;
	
	// Coordinates of where we are going to cut
	x_half = (father->x_max + father->x_min)*0.5;
	y_half = (father->y_max + father->y_min)*0.5;

	// For every particles in the father cube
	for(i = 0; i < N_particle; i++){
		
		x_i = x[Particles->index];
		y_i = y[Particles->index];
		// Where is the particle ? In which child subcube ? 
		// subcube 1
		if ( (x_i < x_half) && (y_i < y_half)){
			child = father->children[0];
			child->N_particle = child->N_particle + 1; // 1 more particle
			// Now we add the particle to the list of particles included in the subcube 1
			child->EndParticle = append(Particles->index,child->EndParticle);
			child->EndParticle = (child->EndParticle)->next;
			// We need a beginning for our Particle list
			if(child->N_particle == 1) child->Particles = child->EndParticle;
		}
		// subcube 2
		if ( (x_i >= x_half) && (y_i < y_half)){
			child = father->children[1];
			child->N_particle = child->N_particle + 1; // 1 more particle
			// Now we add the particle to the list of particles included in the subcube 1
			child->EndParticle = append(Particles->index,child->EndParticle);
			child->EndParticle = (child->EndParticle)->next;
			// We need a beginning for our Particle list
			if(child->N_particle == 1) child->Particles = child->EndParticle;
		}
		// subcube 3
		if ( (x_i < x_half) && (y_i >= y_half)){
			child = father->children[2];
			child->N_particle = child->N_particle + 1; // 1 more particle
			// Now we add the particle to the list of particles included in the subcube 1
			child->EndParticle = append(Particles->index,child->EndParticle);
			child->EndParticle = (child->EndParticle)->next;
			// We need a beginning for our Particle list
			if(child->N_particle == 1) child->Particles = child->EndParticle;
		}
		// subcube 4
		if ( (x_i >= x_half) && (y_i >= y_half)){
			child = father->children[3];
			child->N_particle = child->N_particle + 1; // 1 more particle
			// Now we add the particle to the list of particles included in the subcube 1
			child->EndParticle = append(Particles->index,child->EndParticle);
			child->EndParticle = (child->EndParticle)->next;
			// We need a beginning for our Particle list
			if(child->N_particle == 1) child->Particles = child->EndParticle;
		}
		
	Particles = Particles->next;	
	}
	
	
	// Now, we need to keep dividing until all particles are alone in one cube
	// For every subcubes (children)
	for(i=0 ; i<4 ; i++){
		// We need this otherwise, in the part where we change the boundaries for the children, one child will modify the boundaries of the other child
		x_MAX = father->x_max;
		y_MAX = father->y_max;
		x_MIN = father->x_min;
		y_MIN = father->y_min;
		child = father->children[i];
		
		// If the cube is empty
		if (child->N_particle == 0){
	
			father->children[i] = NULL;
		}
		// If the cube has only particle inside, we put the right information
		if (child->N_particle == 1){
			
			// It will have no children
			for(j=0;j<4;j++){
				child->children[i] = NULL;
			}
			// Recovery of which particle it is
			Particles = child->Particles;
			// Setting the information
			child->x_center = x[Particles->index];
			child->y_center = y[Particles->index];
			// Changing the boundaries ? or make a condition when there is only 1 particles when computing the forces
			
		}
		if((child->N_particle > 1) && (child !=NULL)){
			// Rewriting the new subcubes boundaries
			if(i==0){
				x_MAX = x_half;
				y_MAX = y_half;
			}
			if(i==1){
				x_MIN = x_half;
				y_MAX = y_half;
			}
			if(i==2){
				x_MAX = x_half;
				y_MIN = y_half;
			}
			if(i==3){
				x_MIN = x_half;
				y_MIN = y_half;
			}
			
			// Filling information + changing the boundaries
			child->x_center = 0;
			child->y_center = 0;
			child->x_max = x_MAX;
			child->y_max = y_MAX;
			child->x_min = x_MIN;
			child->y_min = y_MIN;
			
			// Computing centers
			Particles = child->Particles;
			for(j=0;j<child->N_particle;j++){
				child->x_center = child->x_center + x[Particles->index];
				child->y_center = child->y_center + y[Particles->index];
				Particles = Particles->next;
			}
			// Averaging
			child->x_center = (child->x_center)/(child->N_particle);
			child->y_center = (child->y_center)/(child->N_particle);

			
			// Recursive call
			tree_initialization(father->children[i], x, y);
		}
		
	}		
}

// Visualise the constructed tree
void visualize_tree(node* father, int depth){
	printf("Depth %d : We have particles : \n",depth);
	print(father->Particles);
	printf("Center_x : %lf - Center_y : %lf \n",father->x_center,father->y_center);
	printf("\n");
	int i, child_depth;
	for(i=0;i<4;i++){
		child_depth = depth;
		if (father->children[i] != NULL){
			child_depth++;
			visualize_tree(father->children[i], child_depth);
		}
	}
}

// Compute the forces applied on the particle "index" using cluster approximation of parameter "parameter"
void compute_force(node* father, int index, double* force_x, double* force_y, double* x, double* y, double parameter){
	// If the node is not empty
	if(father != NULL){
		// Subcube length;
		double d = (father->x_max) - (father->x_min);
	
		double r_x = father->x_center - x[index];
		double r_y = father->y_center - y[index];
	
		// Distance between center and particle we are interesting into
		double r = sqrt( r_x*r_x + r_y*r_y);
		
		// If the cluster approximation is valid or it is a single particle
		if((r*parameter >= d) || (father->N_particle == 1)){
			// If there are more than one particule, or if the single particle is NOT the considered particle
			if( (father->N_particle > 1) || ((father->Particles)->index != index)){
				// Actualize forces
				force_x[index] = force_x[index] + 1./(r_x*r_x);
				force_y[index] = force_y[index] + 1./(r_y*r_y);
			}
		}
		else{
			int i;
			for(i=0;i<4;i++){
				compute_force(father->children[i], index, force_x, force_y, x, y, parameter);
			}
		}
	}
}





///////////////////////////////
//////////// MAIN /////////////
///////////////////////////////

int main(){
	// index of particle
	int t, i, N=100;
    // Memory initialization
    double *x = (double *)calloc(N, sizeof(double));
    double *y = (double *)calloc(N, sizeof(double));
    double *x_new = (double *)calloc(N, sizeof(double));
    double *y_new = (double *)calloc(N, sizeof(double));
    double *v_x = (double *)calloc(N, sizeof(double));
    double *v_y = (double *)calloc(N, sizeof(double));
    double *force_x = (double *)calloc(N, sizeof(double));
    double *force_y = (double *)calloc(N, sizeof(double));
    List Root_Particles = NULL;
    node* Root = calloc(1, sizeof(*Root));
    
    // Electron properties
 	double charge_e = 1;
 	double m_e = 1;
 	
	// Time step 
	double dt = 0.01;
	int N_t = 100;
	// Time measurement
	clock_t begin = clock();
	
	// Cluster approximation parameter
	double parameter = 0.5;
 	
	// Initial limits for the position distributon
	double x_max=10, y_max = 10;
	double x_min=-10, y_min = -10;
	
	// Position initialization 
	random_position_distribution(x_min, x_max, y_min, y_max, x, y, N);	
	
//	 Writing initial positions
//	FILE *fp = fopen("/Users/pyl/Desktop/ParallelProject/positions_initial.txt", "w");
//    for (i= 0; i < N; i++){
//        fprintf(fp, "%f, %f\n", x[i], y[i]);
//    }
//    fclose(fp);


	// Initialization of the root
	for(i=0 ; i < N ; i++){
		Root->EndParticle = append(i,Root->EndParticle);
		Root->EndParticle = (Root->EndParticle)->next;
	}
	Root->Particles = Root_Particles;
	Root->N_particle = N;
	Root->x_max = x_max;
	Root->y_max = y_max;
	Root->x_min = x_min;
	Root->y_min = y_min;
	



	
	for(t=0;t<N_t;t++){
		// Initialization of the tree
		tree_initialization(Root, x, y);
		for(i=0;i<N;i++){
			// Computing the force for parameter i with the parameter "parameter"
			compute_force(Root, i, force_x, force_y, x, y, parameter);
			
			v_x[i] = v_x[i] + force_x[i]*dt;
			v_y[i] = v_y[i] + force_y[i]*dt;
			x_new[i] = x[i] + v_x[i]*dt;
			y_new[i] = y[i] + v_y[i]*dt;
		}

		// We are going to need the max and min for boundaries!
		x_max = x[0];
		x_min = x[0];
		y_max = y[0];
		y_min = y[0];
		for(i=0;i<N;i++){
			x[i] = x_new[i];
			y[i] = y_new[i];
			force_x[i] = 0;
			force_y[i] = 0;
			
			if (x[i]> x_max) x_max = x[i];
			if (x[i]< x_min) x_min = x[i];
			if (y[i]> y_max) y_max = y[i];
			if (y[i]< y_min) y_min = y[i];
		}
		// Changing the boundaries
		Root->x_max = x_max;
		Root->y_max = y_max;
		Root->x_min = x_min;
		Root->y_min = y_min;
		
	}
	
//		 Writing initial positions
//	FILE *fp1 = fopen("/Users/pyl/Desktop/ParallelProject/positions_final.txt", "w");
//    for (i= 0; i < N; i++){
//        fprintf(fp1, "%f, %f\n", x[i], y[i]);
//    }
//    fclose(fp1);

	// Time measurement
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Time spent for N = %d particles : %lf \n",N,time_spent);
	
    //Free memory
    free(x);
    free(y);
    free(x_new);
    free(y_new);
    free(v_x);
    free(v_y);
    free(force_x);
    free(force_y);
	return 0;
}