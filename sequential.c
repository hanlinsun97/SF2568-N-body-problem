#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "list.h"



// Structure defining a node in the tree
struct node{
	// In 2D each node has up to 4 children
	struct node* children[4];
	double charge_node, x_center, y_center;
	
	// N_particle is the number of particles in the subcube
	// List Particles is a list containing the number "i" of the particles contained in the subcube
	List Particles;
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
		x[i] = frand(x_min,x_max);
		y[i] = frand(y_min,y_max);
	}
}


void tree_initialization(node *father, double x_min, double x_max, double y_min, double y_max, double* x, double* y, double charge_e){
	int i, index;
	double x_half, y_half, x_i, y_i;
	int N_particle = father->N_particle;
		
	List Particles = father->Particles;
	// Memory allocation
	for(i=0 ; i<4 ; i++){
		father->children[i]=calloc(1, sizeof(*(father->children[i])));	
		(father->children[i])->Particles=NULL;
	}
	node *child;
	
	// Coordinates of where we are going to cut
	x_half = (x_max + x_min)*0.5;
	y_half = (y_max + y_min)*0.5;
	
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
			child->Particles = append(Particles->index,child->Particles);
		}
		// subcube 2
		if ( (x_i >= x_half) && (y_i < y_half)){

			child = father->children[1];
			child->N_particle = child->N_particle + 1; // 1 more particle
			// Now we add the particle to the list of particles included in the subcube 1
			child->Particles = append(Particles->index,child->Particles);
		}
		// subcube 3
		if ( (x_i < x_half) && (y_i >= y_half)){
			
			child = father->children[2];
			child->N_particle = child->N_particle + 1; // 1 more particle
			// Now we add the particle to the list of particles included in the subcube 1
			child->Particles = append(Particles->index,child->Particles);
		}
		// subcube 4
		if ( (x_i >= x_half) && (y_i >= y_half)){
			
			child = father->children[3];
			child->N_particle = child->N_particle + 1; // 1 more particle
			// Now we add the particle to the list of particles included in the subcube 1
			child->Particles = append(Particles->index,child->Particles);
		}
		
	Particles = Particles->next;	
	}
	
	
	// Now, we need to keep dividing until all particles are alone in one cube
	// For every subcubes (children)
	for(i=0 ; i<4 ; i++){
		getchar();
		child = father->children[i];
		// If the cube is empty
		if (child->N_particle == 0){
			printf("empty \n");
			printf("Number of particle : %d \n", (father->children[i])->N_particle);
			father->children[i] = NULL;
		}
		// If the cube has only particle inside, we put the right information
		if (child->N_particle == 1){
			printf("Alone \n");
			printf("Number of particle : %d \n", (father->children[i])->N_particle);
			print((father->children[i])->Particles);
			
			int j;
			// It will have no children
			for(j=0;j<4;j++){
				child->children[i] = NULL;
			}
			// Recovery of which particle it is
			Particles = child->Particles;
			index = Particles->index;
			// Setting the information
			child->charge_node = charge_e;
			child->x_center = x[i];
			child->y_center = y[i];
			
		}
		if((child->N_particle > 1) && (child !=NULL)){
			// Rewriting the new subcubes boundaries
			if(i==0){
				x_max = x_half;
				y_max = y_half;
			}
			if(i==1){
				x_min = x_half;
				y_max = y_half;
			}
			if(i==2){
				x_max = x_half;
				y_min = y_half;
			}
			if(i==3){
				x_min = x_half;
				y_min = y_half;
			}
			printf("Dividing \n");
			printf("Number of particle : %d \n", (father->children[i])->N_particle);
			print((father->children[i])->Particles);
			printf("\n");
			// Recursive call
			tree_initialization(father->children[i], x_min, x_max, y_min, y_max, x, y, charge_e);
		}
		
	}		
}

void visualize_tree(node* Root){
	printf("Root \n");
	print(Root->Particles);	

	node* child1 = Root->children[0];
	node* child2 = Root->children[1];
	node* child3 = Root->children[2];
	node* child4 = Root->children[3];
	

}

int main(){
	// index of particle
	int N = 10;
    // Memory initialization
    double *x = (double *)malloc(N * sizeof(double));
    double *y = (double *)malloc(N * sizeof(double));
    double *v_x = (double *)malloc(N * sizeof(double));
    double *v_y = (double *)malloc(N * sizeof(double));
    double *force_x = (double *)calloc(N, sizeof(double));
    double *force_y = (double *)calloc(N, sizeof(double));
    List Root_Particles = NULL;
    node* Root = calloc(1, sizeof(*Root));
    
        // charge of 1 electron
 	double charge_e = 1;
 	
	// Initial limits for the position distributon
	double x_max, y_max = 1;
	double x_min, y_min = -1;
	
	// Position initialization 
	random_position_distribution(x_min, x_max, y_min, y_max, x, y, N);	
	// Initialization of the tree
	int i;
	for(i=0 ; i < N ; i++) Root_Particles = append(i,Root_Particles);
	
	Root->Particles = Root_Particles;
	Root->N_particle = N;
	

	tree_initialization(Root, x_min, x_max, y_min, y_max, x, y, charge_e);
	visualize_tree(Root);
	
	return 0;
}