#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <GLUT/glut.h>

struct node_t
{
    int particle;
    int has_particle;
    int has_children;
    double min_x, max_x, min_y, max_y, total_mass, c_x, c_y;
    struct node_t *children;
};

double frand(double xmin, double xmax);
void put_particle_in_tree(int new_particle, struct node_t *node);
void place_particle(int particle, struct node_t *node);
void set_node(struct node_t *node);
void free_node(struct node_t *node);

double calculate_mass(struct node_t *node);
double calculate_center_of_mass_x(struct node_t *node);
double calculate_center_of_mass_y(struct node_t *node);
void update_forces();
void update_forces_help(int particle, struct node_t *node);
void calculate_force(int particle, struct node_t *node, double r);

//Some constants and global variables
int N = 250;
const double L = 1, W = 1, dt = 1e-3, alpha = 0.25, V = 50, epsilon = 1e-1, grav = 0.04; //grav should be 100/N
double *x, *y, *u, *v, *force_x, *force_y, *mass;
struct node_t *root;

/*
 * Function for producing a random number between two double values
 */
double frand(double xmin, double xmax)
{
    return xmin + (xmax - xmin) * rand() / RAND_MAX;
}

/*
 * Prints the time between two clocks in seconds
 */
void print_time(clock_t s, clock_t e)
{
    printf("Time: %f seconds\n", (double)(e - s) / CLOCKS_PER_SEC);
}

/*
 * Updates the positions of the particles of a time step.
 */
void time_step()
{
    //Allocate memory for root
    root = malloc(sizeof(struct node_t));
    set_node(root);
    float MAX_X = x[0];
    float MAX_Y = y[0];
    float MIN_X = x[0];
    float MIN_Y = y[0];
    for (int i = 0; i < N; i++)
    {
        if (x[i] > MAX_X)
        {
            MAX_X = x[i];
        }
        if (y[i] > MAX_Y)
        {
            MAX_Y = y[i];
        }
        if (x[i] < MIN_X)
        {
            MIN_X = x[i];
        }
        if (y[i] < MIN_Y)
        {
            MIN_Y = y[i];
        }
    }
    root->min_x = MIN_X;
    root->max_x = MAX_X;
    root->min_y = MIN_Y;
    root->max_y = MAX_Y;

    //Put particles in tree
    for (int i = 0; i < N; i++)
    {
        put_particle_in_tree(i, root);
    }

    //Calculate mass and center of mass
    calculate_mass(root);
    calculate_center_of_mass_x(root);
    calculate_center_of_mass_y(root);

    //Calculate forces
    update_forces();

    //Update velocities and positions
    for (int i = 0; i < N; i++)
    {
        double ax = force_x[i] / mass[i];
        double ay = force_y[i] / mass[i];
        u[i] += ax * dt;
        v[i] += ay * dt;
        x[i] += u[i] * dt;
        y[i] += v[i] * dt;

        /* This of course doesn't make any sense physically,
         * but makes sure that the particles stay within the
         * bounds. Normally the particles won't leave the
         * area anyway.
         */
    }

    //Free memory
    free_node(root);
    free(root);
}

/*
 * Puts a particle recursively in the Barnes Hut quad-tree.
 */
void put_particle_in_tree(int new_particle, struct node_t *node)
{
    //If no particle is assigned to the node
    if (!node->has_particle)
    {
        node->particle = new_particle;
        node->has_particle = 1;
    }
    //If the node has no children
    else if (!node->has_children)
    {
        //Allocate and initiate children
        node->children = malloc(4 * sizeof(struct node_t));
        for (int i = 0; i < 4; i++)
        {
            set_node(&node->children[i]);
        }

        //Set boundaries for the children
        node->children[0].min_x = node->min_x;
        node->children[0].max_x = (node->min_x + node->max_x) / 2;
        node->children[0].min_y = node->min_y;
        node->children[0].max_y = (node->min_y + node->max_y) / 2;

        node->children[1].min_x = (node->min_x + node->max_x) / 2;
        node->children[1].max_x = node->max_x;
        node->children[1].min_y = node->min_y;
        node->children[1].max_y = (node->min_y + node->max_y) / 2;

        node->children[2].min_x = node->min_x;
        node->children[2].max_x = (node->min_x + node->max_x) / 2;
        node->children[2].min_y = (node->min_y + node->max_y) / 2;
        node->children[2].max_y = node->max_y;

        node->children[3].min_x = (node->min_x + node->max_x) / 2;
        node->children[3].max_x = node->max_x;
        node->children[3].min_y = (node->min_y + node->max_y) / 2;
        node->children[3].max_y = node->max_y;

        //Put old particle into the appropriate child
        place_particle(node->particle, node);

        //Put new particle into the appropriate child
        place_particle(new_particle, node);

        //It now has children
        node->has_children = 1;
    }
    //Add the new particle to the appropriate children
    else
    {
        //Put new particle into the appropriate child
        place_particle(new_particle, node);
    }
}

/*
 * Puts a particle in the right child of a node with children.
 */
void place_particle(int particle, struct node_t *node)
{
    if (x[particle] <= (node->min_x + node->max_x) / 2 && y[particle] <= (node->min_y + node->max_y) / 2)
    {
        put_particle_in_tree(particle, &node->children[0]);
    }
    else if (x[particle] > (node->min_x + node->max_x) / 2 && y[particle] < (node->min_y + node->max_y) / 2)
    {
        put_particle_in_tree(particle, &node->children[1]);
    }
    else if (x[particle] < (node->min_x + node->max_x) / 2 && y[particle] > (node->min_y + node->max_y) / 2)
    {
        put_particle_in_tree(particle, &node->children[2]);
    }
    else
    {
        put_particle_in_tree(particle, &node->children[3]);
    }
}

/*
 * Sets initial values for a new node
 */
void set_node(struct node_t *node)
{
    node->has_particle = 0;
    node->has_children = 0;
}

/*
 * Frees memory for a node and its children recursively.
 */
void free_node(struct node_t *node)
{
    if (node->has_children)
    {
        free_node(&node->children[0]);
        free_node(&node->children[1]);
        free_node(&node->children[2]);
        free_node(&node->children[3]);
        free(node->children);
    }
}

/*
 * Calculates the total mass for the node. It recursively updates the mass
 * of itself and all of its children.
 */
double calculate_mass(struct node_t *node)
{
    if (!node->has_particle)
    {
        node->total_mass = 0;
    }
    else if (!node->has_children)
    {
        node->total_mass = mass[node->particle];
    }
    else
    {
        node->total_mass = 0;
        for (int i = 0; i < 4; i++)
        {
            node->total_mass += calculate_mass(&node->children[i]);
        }
    }
    return node->total_mass;
}

/*
 * Calculates the x-position of the centre of mass for the
 * node. It recursively updates the position of itself and
 * all of its children.
 */
double calculate_center_of_mass_x(struct node_t *node)
{
    if (!node->has_children)
    {
        node->c_x = x[node->particle];
    }
    else
    {
        node->c_x = 0;
        double m_tot = 0;
        for (int i = 0; i < 4; i++)
        {
            if (node->children[i].has_particle)
            {
                node->c_x += node->children[i].total_mass * calculate_center_of_mass_x(&node->children[i]);
                m_tot += node->children[i].total_mass;
            }
        }
        node->c_x /= m_tot;
    }
    return node->c_x;
}

/*
 * Calculates the y-position of the centre of mass for the
 * node. It recursively updates the position of itself and
 * all of its children.
 */
double calculate_center_of_mass_y(struct node_t *node)
{
    if (!node->has_children)
    {
        node->c_y = y[node->particle];
    }
    else
    {
        node->c_y = 0;
        double m_tot = 0;
        for (int i = 0; i < 4; i++)
        {
            if (node->children[i].has_particle)
            {
                node->c_y += node->children[i].total_mass * calculate_center_of_mass_y(&node->children[i]);
                m_tot += node->children[i].total_mass;
            }
        }
        node->c_y /= m_tot;
    }
    return node->c_y;
}

/*
 * Calculates the forces in a time step of all particles in
 * the simulation using the Barnes Hut quad tree.
 */

void update_forces()
{
    for (int i = 0; i < N; i++)
    {
        force_x[i] = 0;
        force_y[i] = 0;
        update_forces_help(i, root);
    }
}

/*
 * Help function for calculating the forces recursively
 * using the Barnes Hut quad tree.
 */
void update_forces_help(int particle, struct node_t *node)
{
    //The node is a leaf node with a particle and not the particle itself
    if (!node->has_children && node->has_particle && node->particle != particle)
    {
        double r = sqrt((x[particle] - node->c_x) * (x[particle] - node->c_x) + (y[particle] - node->c_y) * (y[particle] - node->c_y));
        calculate_force(particle, node, r);
    }
    //The node has children
    else if (node->has_children)
    {
        //Calculate r and theta
        double r = sqrt((x[particle] - node->c_x) * (x[particle] - node->c_x) + (y[particle] - node->c_y) * (y[particle] - node->c_y));
        double theta = (node->max_x - node->min_x) / r;

        /* If the distance to the node's centre of mass is far enough, calculate the force,
         * otherwise traverse further down the tree
         */
        if (theta < 0.5)
        {
            calculate_force(particle, node, r);
        }
        else
        {
            update_forces_help(particle, &node->children[0]);
            update_forces_help(particle, &node->children[1]);
            update_forces_help(particle, &node->children[2]);
            update_forces_help(particle, &node->children[3]);
        }
    }
}

/*
 * Calculates and updates the force of a particle from a node.
 */
void calculate_force(int particle, struct node_t *node, double r)
{
    double temp = grav * mass[particle] * node->total_mass / ((r + epsilon) * (r + epsilon) * (r + epsilon));
    force_x[particle] += (x[particle] - node->c_x) * temp;
    force_y[particle] += (y[particle] - node->c_y) * temp;
}

/*
 * Main function.
 */

int main(int argc, char *argv[])
{

    //The second argument sets the number of time steps
    int time_steps = 1000;
    //Initiate memory for the vectors
    x = (double *)malloc(N * sizeof(double));
    y = (double *)malloc(N * sizeof(double));
    u = (double *)malloc(N * sizeof(double));
    v = (double *)malloc(N * sizeof(double));
    force_x = (double *)calloc(N, sizeof(double));
    force_y = (double *)calloc(N, sizeof(double));
    mass = (double *)malloc(N * sizeof(double));

    //Set the initial values
    for (int i = 0; i < N; i++)
    {
        mass[i] = 1;
        double R = frand(0, L / 4);
        double theta = frand(0, 2 * M_PI);
        x[i] = L / 2 + R * cos(theta);
        y[i] = W / 2 + alpha * R * sin(theta);
        double R_prim = sqrt(pow(x[i] - L / 2, 2) + pow(y[i] - W / 2, 2));
        //u[i] = -V * R_prim * sin(theta);
        //v[i] = V * R_prim * cos(theta);
        u[i] = 1000;
        v[i] = 0;
    }

    /* Run the GLUT display function if the graphics mode is on.
     * Otherwise just run the simulations without graphics
     */

    //Begin taking time
    long start = clock();

    //The main loop
    for (int i = 0; i < time_steps; i++)
    {
        time_step();
    }
    FILE *fp = fopen("/Users/iconoclast/Desktop/position.txt", "w");
    for (int j = 0; j < N; j++)
    {
        fprintf(fp, "%f, %f\n", x[j], y[j]);
    }
    fclose(fp);

    //Stop taking time and print elapsed time
    long stop = clock();
    print_time(start, stop);
    printf("Success\n");

    //Free memory
    free(x);
    free(y);
    free(u);
    free(v);
    free(force_x);
    free(force_y);
    free(mass);

    return 0;
}
