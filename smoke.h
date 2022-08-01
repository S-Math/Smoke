#include <iostream>
#include <vector>
#include <cstdlib>
#include <Eigen/Dense>
#include <glm/glm.hpp>
#include <ctime>
#include <random>
#include <chrono>

class ClassSolverNS
{
    public:
        ClassSolverNS(const int N_grid, const double dt, const double diff);
        ~ClassSolverNS();

        struct Vertex
        {
            glm::vec2 pos;
            glm::vec3 color;
        };

        void initialize_density();
        void add_source(const std::vector<Vertex> &s);
        void diffuse(const int &k_max);
        void set_boundaries_velocity();
        void set_boundaries_density();
        void advect();
        void density_step(const Eigen::VectorXd &s, const int &k_max);
        int indices(const int i, const int j);
        void reverseIndices(const int &indices, float &i, float &j);
        void comparison();
        void print();

    protected:
        // Number of point per side of the grid
        const int N;
        int get_N(){return N;};
 
        // time step
        const double dt;

        // diff constant
        const double diff;

    public:
        // Create velocity vector
        Eigen::MatrixXd velocity;

        std::vector<Vertex> vertexDensity;

};

// create destructor function (empty so far)
ClassSolverNS::~ClassSolverNS()
{

}




ClassSolverNS::ClassSolverNS(const int N_grid, const double dt, const double diff)
:N(N_grid),
dt(dt),
diff(diff),
velocity((N+2)*(N+2),2),
vertexDensity((N+2)*(N+2))
{
    /**
     * @brief Create constructor function (initialize number of point
     *  per side of the grid, time step, velocity matrix and density vector
     * 
     */
}

void ClassSolverNS::comparison()
{
    int N_test = 1000000;
    std::vector<double> stdvec(N_test);
    Eigen::VectorXd eigenvec(N_test);



    std::random_device rdevice; // used to obtain a seed
    std::mt19937 gen(rdevice()); // standard mersenne_twister engine seeded
    std::uniform_real_distribution<> dis(0.0, 1.0);

    auto start_stdvec = std::chrono::high_resolution_clock::now();

    for (size_t i =0; i<stdvec.size(); i++)
    {
        stdvec[i] = dis(gen);
    }

    auto stop_stdvec = std::chrono::high_resolution_clock::now();

    auto start_eigenvec = std::chrono::high_resolution_clock::now();

    eigenvec = Eigen::VectorXd::Random(N_test).array();

    auto stop_eigenvec = std::chrono::high_resolution_clock::now();

    auto duration_eigenvec = 
        std::chrono::duration_cast<std::chrono::milliseconds>(stop_eigenvec- start_eigenvec);

    auto duration_stdvec = 
        std::chrono::duration_cast<std::chrono::milliseconds>(stop_stdvec - start_stdvec);

    std::cout << "duration for stdvec : " << duration_stdvec.count() << std::endl;
    std::cout << "duration for eigenvec : " << duration_eigenvec.count() << std::endl;

}


void ClassSolverNS::print()
{
    std::cout << std::endl;
    std::cout << "Density :" << std::endl;
    //for (int i=0; i<vertexDensity.size(); i++)
    //{
    //    std::cout << "   " << vertexDensity[i].color[0] << std::endl;
    //}

    for (int i=0; i<get_N()+2; i++)
    {
        for (int j=0; j<get_N()+2; j++)
        {
            std::cout << vertexDensity[indices(i,j)].color[0] << " ";
        }

        std::cout << std::endl;
    }

    std::cout << "Density position :" << std::endl;
    for (int i=0; i<get_N()+2; i++)
    {
        for (int j=0; j<get_N()+2; j++)
        {
            std::cout << "(" << vertexDensity[indices(i,j)].pos[0] << "; " << vertexDensity[indices(i,j)].pos[1] <<") ";
        }

        std::cout << std::endl;
    }


}

int ClassSolverNS::indices(const int i, const int j)
{
    return i + (this->N+2)*j;
}

void ClassSolverNS::reverseIndices(const int &indices, float &i, float &j)
{
    i = indices % (this->N+2);
    j = indices - i * (this->N+2);
    std::cout << "i = " << i << std::endl;
    std::cout << "j = " << j << std::endl;
}




// initialise random density
void ClassSolverNS::initialize_density()
{
    /**
     * @brief the vector density is randomely initialised between 0 and 1
     * 
     */

    std::random_device rdevice; // used to obtain a seed
    std::mt19937 gen(rdevice()); // standard mersenne_twister engine seeded
    std::uniform_real_distribution<> dis(0.0, 1.0);


    for (size_t i = 0; i<vertexDensity.size(); i++)
    {
        double randDensity = dis(gen);
        vertexDensity[i].color = {randDensity,
                                  randDensity,
                                  randDensity};

        reverseIndices(i, 
                       vertexDensity[i].pos[0],
                       vertexDensity[i].pos[1]);

        //std::cout << "randdensity : " << randDensity << std::endl;
    }

}


// Add source to density
void ClassSolverNS::add_source(const std::vector<Vertex> &s)
{
    /**
     * @brief this function simply adds the source density to
     * the density field defined in this clas
     * 
     */

    for (size_t i = 0; i <vertexDensity.size(); i++)
    {
        vertexDensity[i].color += s[i].color;
    }

}

void ClassSolverNS::diffuse(const int &k_max)
{
    /**
     * @brief This function diffuses the density to the neighbouring cells
     * A linear system results in a sparse matrix. To solve it we use a
     * couple of Gaussian Siedel iteration to approximate the true solution
     * 
     */
    float a = this->dt * this->diff * this->N *this->N;

    std::vector<Vertex> vertexDensity_0(vertexDensity);

    for (unsigned int k=0; k<k_max; k++)
    {
        for (unsigned int i=1; i<=N; i++)
        {
            for (unsigned int j=1; j<=N; j++)
            {                
                vertexDensity[indices(i,j)].color = (vertexDensity_0[indices(i,j)].color 
                                            + a*(vertexDensity[indices(i-1,j)].color
                                                +vertexDensity[indices(i+1,j)].color
                                                +vertexDensity[indices(i,j-1)].color
                                                +vertexDensity[indices(i,j+1)].color))
                                                / (1+4*a);
            }
        }
        set_boundaries_density();
    }
}

void ClassSolverNS::set_boundaries_velocity()
{
    /**
     * @brief This function sets the boundary condition for the 2
     * components of the velocity as well as the density field
     * 
     */

    using Eigen::all;

    for (unsigned int i=1; i<=N; i++)
    {
        velocity(indices(0,i), 0) = -velocity(indices(1,i), 0); 
        velocity(indices(0,i), 1) =  velocity(indices(1,i), 1);

        velocity(indices(N+1,i), 0) = -velocity(indices(N,i), 0); 
        velocity(indices(N+1,i), 1) =  velocity(indices(N,i), 1);

        velocity(indices(i,0), 1) = -velocity(indices(i,1), 1);
        velocity(indices(i,0), 0) =  velocity(indices(i,1), 0);

        velocity(indices(i,N+1), 1) = -velocity(indices(i,N), 1);
        velocity(indices(i,N+1), 0) =  velocity(indices(i,N), 0);
    }

    velocity(indices(0,0), all) = 
            0.5 * (velocity(indices(1,0), all)+velocity(indices(0,1), all));
    velocity(indices(0,N+1), all) = 
            0.5 * (velocity(indices(1,N+1), all)+velocity(indices(0,N), all));
    velocity(indices(N+1,0), all) = 
            0.5 * (velocity(indices(N,0), all)+velocity(indices(N+1,1), all));
    velocity(indices(N+1,N+1), all) = 
            0.5 * (velocity(indices(N,N+1), all)+velocity(indices(N+1,N), all));
}

void ClassSolverNS::set_boundaries_density()
{
    for (size_t k = 0; k<glm::length(vertexDensity[0].color) ; k++ )
    {
        for (int i=1; i<=N; i++)
        {
            vertexDensity[indices(0,i)].color = vertexDensity[indices(1,i)].color;
            vertexDensity[indices(N+1,i)].color  = vertexDensity[indices(N,i)].color;
            vertexDensity[indices(i,0)].color  = vertexDensity[indices(i,1)].color;
            vertexDensity[indices(i,N+1)].color  = vertexDensity[indices(i,N)].color;
        }

        vertexDensity[indices(0,0)].color[k] = 
                    0.5 * (vertexDensity[indices(1,0)].color[k] + vertexDensity[indices(0,1)].color[k]);
        vertexDensity[indices(0,N+1)].color[k] = 
                    0.5 * (vertexDensity[indices(1,N+1)].color[k] + vertexDensity[indices(0,N)].color[k]);
        vertexDensity[indices(N+1,0)].color[k] = 
                    0.5 * (vertexDensity[indices(N,0)].color[k] + vertexDensity[indices(N+1,1)].color[k]);
        vertexDensity[indices(N+1,N+1)].color[k] = 
                    0.5 * (vertexDensity[indices(N,N+1)].color[k] + vertexDensity[indices(N+1,N)].color[k]);
    }
}

void ClassSolverNS::advect()
{
    /**
     * @brief This function advect the density along the velocity field
     * The idea is to model density as a set of particle.
     * At time t+dt we look which particle end up in the middle of the
     * cell and we trace back its position at time t
     * 
     * To calculate the density advected by this particle at one cell at time t+dt
     * we interpolate its value at time t with the four corners of the cell it belonged
     * at time t 
     */

    //Coordinate of the particle at time t
    double x,y;

    //Coordinate of the cell who host the particle at time t
    int i0, j0;

    for (unsigned int i=1; i<=N; i++)
    {
        for (unsigned int j=1; j<=N; j++)
        {
            //Calculate the position it had at time t
            x = i - dt*velocity(indices(i,j), 0); 
            y = j - dt*velocity(indices(i,j), 1);

            if (x < 0.5) x = 0.5;
            if (y < 0.5) y = 0.5;
            if (x > N+0.5) x = N+0.5;
            if (y > N+0.5) y = N+0.5;

            i0 = static_cast<int>(x);
            j0 = static_cast<int>(y);

            // at time t, interpolate the position
            for (size_t k = 0; k<glm::length(vertexDensity[0].color) ; k++ )
            {
                vertexDensity[indices(i,j)].color[k] = 
                                (1-x+i0)*
                                        ((1-y+j0) * vertexDensity[indices(i0,j0)].color[k]+
                                        (y-j0) * vertexDensity[indices(i0,j0+1)].color[k]
                                        )+
                                    (x-i0)*
                                        ((1-y+i0) * vertexDensity[indices(i0+1,j0)].color[k]+
                                        (y-j0) * vertexDensity[indices(i0+1, j0+1)].color[k]
                                        );

            }

        }
    }

    set_boundaries_velocity();
}

void ClassSolverNS::density_step(const Eigen::VectorXd &s, const int &k_max)
{
    /**
     * @brief This function realises a full step for evolving
     * the density from time t to time t+dt
     * 
     * 
     */
    //this->add_source(s);
    this->diffuse(k_max);
    this->advect();

}