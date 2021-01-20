#pragma once
/**
 *  Methods and Functions file for 3DZeroPlane
 * 
 * 
 * 
*/
#include <vector>
#include <math.h>
#include <complex>
#include <fstream>
#include <iostream>
#include <boost/function.hpp>




// Type definitions
// typedef boost::multiprecision::cpp_dec_float_50 mp_type;
namespace mp = boost::multiprecision;
typedef mp::mpfr_float_100 mp_type; // faster compile time


// template <typename T>
// using R2_function = boost::function<T(T,T)>;


// template <typename T>
// using R3_function = boost::function<T(T,T,T)>;




template <class T>
T steffensen(T (*f)(T,T), T x0, T x1, T eta, size_t iter) {
    // Implement steffensens method for root finding.
}

template <class T>
T secant(T (*f)(T,T), T x0, T x1, T eta, size_t iter) {
    // N iterations of the secant method.
    /*
        F: R2 -> R is curried along the second dimension and the root is refined
        given a proper bracket.        
    */

    T xr;
    std::vector<T> roots;
    roots.push_back(x0);
    roots.push_back(x1);

    for (size_t i = 0; i < iter; i++) {
        xr = roots[i+1] - f(roots[i+1],eta)*(roots[i+1]-roots[i])/(f(roots[i+1],eta) - f(roots[i],eta));
        roots.push_back(xr);
    }

    return roots[iter];
}



template <class T>
T secant2(T (*f)(T,T), T x, T e0, T e1, size_t iter) {
    // Secant method as above but function is curried in the first dimension.
    T er;
    std::vector<T> roots;
    roots.push_back(e1);
    roots.push_back(e0);

    for (size_t i = 0; i < iter; i++) {
        er = roots[i+1] - f(x,roots[i+1])*(roots[i+1]-roots[i])/(f(x,roots[i+1]) - f(x,roots[i]));
        roots.push_back(er);
    }

    return roots[iter];
}




template <class T>
void matrix_to_file2(std::vector<std::vector<std::vector<T> > > matrix, std::ofstream& ofile) {
    for (size_t i = 0; i < matrix[0].size(); i++) {
        for (size_t j = 0; j < matrix.size(); j++) {
            ofile << matrix[j][i][0] << " "; 
        }
        ofile << "\n ";
    }
}

template <class T>
void matrix_to_file(std::vector<std::vector<T> > matrix, std::ofstream& ofile) {
    for (size_t i = 0; i < matrix[0].size(); i++) {
        for (size_t j = 0; j < matrix.size(); j++) {
            ofile << std::setprecision(17) << matrix[j][i] << " "; 
        }
        ofile << "\n";
    }
}

template <class T>
void matrix_to_file3(std::vector<std::vector<T> > matrix, std::string file) {
    std::ofstream ofile;
    ofile.open(file);
    for (size_t i = 0; i < matrix.size(); i++) {
        for (size_t j = 0; j < matrix[i].size(); j++) {
            ofile << std::setprecision(17) << matrix[i][j] << " "; 
        }
        ofile << "\n";
    }
}



template <class T>
void vector_to_file(std::vector<T> vector, std::ofstream& ofile) {
    for (size_t i = 0; i < vector.size(); i++) {
        ofile << vector[i] << " \n";
    }
}

template <class T>
std::vector<std::vector<T> > compute_zero_set(T (*f)(T,T),std::vector<std::vector<T> > bracket, size_t secant_N){

    std::vector<T> z_zero, eta_zero, time_zero;
    std::vector<std::vector<T> > result;
    
    for (size_t i = 0; i < bracket[1].size(); i++) {
        z_zero.push_back(secant(f,bracket[0][i],bracket[2][i],bracket[1][i],secant_N));
        eta_zero.push_back(bracket[1][i]);
    }

    // Add in the zeros from the other dimension!
    // for (int i = 0; i < bracket[1].size(); ++i)
    // {
    //     /* code */
    // }



    result.push_back(z_zero);
    result.push_back(eta_zero);
    // return the bracketing set too:
    for (int i = 0; i < bracket.size(); ++i)
    {
        result.push_back(bracket[i]);
    }



    return result;
}

template <class T>
std::vector<std::vector<T> > bracket_roots(std::vector<std::vector<T> > surface, std::vector<T> z_grid, std::vector<T> eta_grid) {
    // Scanning bracket method, not very robust at all!!
    // Needs something better
    std::vector<T> bl_z, bl_eta, br_z, br_eta;
    std::vector<std::vector<T> > result;

    for (size_t i = 0; i < eta_grid.size(); ++i) {
        for (size_t j = 0; j < z_grid.size(); ++j) {
            if (j > 2 && j < z_grid.size()) {
                if (surface[j-2][i]*surface[j][i] < 0) {
                    bl_z.push_back(z_grid[j-2]);
                    bl_eta.push_back(eta_grid[i]);
                    br_z.push_back(z_grid[j]);
                    br_eta.push_back(eta_grid[i]);
                }
            }
        }
    }
    result.push_back(bl_z);
    result.push_back(bl_eta);
    result.push_back(br_z);
    result.push_back(br_eta);

    return result;
}


template <class T>
std::vector<std::vector<T> > bracket_roots_2(std::vector<std::vector<T> > surface, std::vector<T> z_grid, std::vector<T> eta_grid) {

    std::vector<T> bl_z, bl_eta, br_z, br_eta;
    std::vector<std::vector<T> > result;
    T bracket_point_left, bracket_point_right;
    // Initial bracketing set
    bracket_point_left = z_grid[0];
    bracket_point_right = z_grid[z_grid.size()-1];


    // For each slice along the eta direction
    for (int i = 0; i < eta_grid.size(); ++i)
    {
        // Scan through each eta=const slice and find the sign changes
        // relative to the first element, if it isnt nan
        for (int j = 0; j < z_grid.size(); ++j)
        {
            // Check if bracket_point_left is nan:
            if (!isnan(bracket_point_left))
            {
                if (bracket_point_left*surface[j][i] <= 0 )
                {
                    // If above is true, we found a sign change.
                    // set right bracket pt to this value:
                    bracket_point_right = surface[j][i];
                    bl_z.push_back(bracket_point_left);
                    bl_eta.push_back(eta_grid[j]);
                    br_z.push_back(bracket_point_right);
                    br_eta.push_back();
                }
            }
        }

    }

}



template <class T> 
std::vector<std::vector<T> > bracket_secant_method(T (*f)(T,T), std::vector<T> zgrid, std::vector<T> etagrid, size_t secant_N) {
    // Use a naive method to bracket N roots in the domain
    // Note: This method may miss roots if the root is a constant value
    //       along one of the grid lines.
    // BIG NOTE:
    //  Should do this from both directions and combine the zero sets together.
    size_t Nroots;
    size_t Nz = zgrid.size();
    size_t Ne = etagrid.size();
    T dz = (zgrid[Nz-1] - zgrid[0])/Nz;
    T deta = (etagrid[Ne-1] - etagrid[0])/Ne;
    T bl, br;
    std::vector<T> zleft_set,zright_set;
    std::vector<T> eleft_set,eright_set;
    std::vector<std::vector<T> > result;
    std::vector<T> roots;
    std::vector<T> z_zero,eta_zero;

    // Bracketing along the z direction first. i.e. eta is fixed, r is determined.
    std::cout << "First bracketing direction" << std::endl;
    for (int i = 0; i < etagrid.size(); ++i)
    {
        bl= zgrid[0];
        for (int j = 1; j < zgrid.size()-1; ++j)
        {
            br = zgrid[j];
            if (!isnan(f(bl,etagrid[i])) && !isnan(f(br,etagrid[i])))
            {
                if (f(bl,etagrid[i])*f(br,etagrid[i]) < 0)
                {
                    // if f_(j) brackets with f_0, then 
                    // f_(j-1) is also a bracket with f_(j) 
                    // move the next bracket to start at f_(j) on the left
                    // then continue f_(j+1) on the right and determine crossing.
                    zleft_set.push_back(br-1.1*dz);
                    zright_set.push_back(br);
                    eleft_set.push_back(etagrid[i]);
                    eright_set.push_back(etagrid[i]);
                    bl = zgrid[j];    
                }       
            }
        }
    }
    // Nroots = 0;
    Nroots = zleft_set.size();
    // Refine each bracket along this dimension with the secant method:
    for (size_t i = 0; i < zleft_set.size(); i++) {
        z_zero.push_back(secant(f,zleft_set[i],zright_set[i], eleft_set[i],secant_N));
        eta_zero.push_back(eleft_set[i]);
    }



    std::cout << "Switching bracketing direction" << std::endl;
    // Root find from the other grid direction, use a staggered grid
    // so that we arent overlapping the same roots.
    // z is fixed, eta is determined.

    
    for (int i = 0; i < zgrid.size(); ++i)
    {
        bl = etagrid[0];
        for (int j = 1; j < etagrid.size()-1; ++j)
        {
            br = etagrid[j];
            if (!isnan(f(zgrid[i],bl)) && !isnan(f(zgrid[i],br)))
            {
                if (f(zgrid[i],bl)*f(zgrid[i],br) < 0)
                {
                    zleft_set.push_back(zgrid[i]);
                    zright_set.push_back(zgrid[i]);
                    eleft_set.push_back(br - 1.1*deta);
                    eright_set.push_back(br);
                    bl = etagrid[j];

                }
            }
        }
    }
    
    for (size_t i = Nroots; i < zleft_set.size(); i++) {
        eta_zero.push_back(secant2(f,zleft_set[i],eleft_set[i], eright_set[i],secant_N));
        z_zero.push_back(zleft_set[i]);
    }

    std::cout << "Found " << Nroots << " roots along z. And " << zleft_set.size() - Nroots << " roots along eta.\n";



    // Save the roots
    result.push_back(z_zero);
    result.push_back(eta_zero);

    // Save the bracketing set 
    result.push_back(zleft_set);
    result.push_back(eleft_set);
    result.push_back(zright_set);
    result.push_back(eright_set);

    return result;
}

template <class T>
std::vector<std::vector<T> > radial_bracket(T (*f)(T,T), std::vector<T> zgrid, std::vector<T> etagrid) {
    /* Idea here:
        Take a function f(x,y) on a square grid in cartesion coordinates.
        take radial lines g(r,theta) and fix theta, then bracket the intersection
        between f and g. 
        Go through a full domain theta = [0,2pi].
        store the roots as x,y.
    

    */
}












template <class T>
std::vector<T> create_grid(T min, T max, size_t N) {
    std::vector<T> grid;
    T step = (max - min)/N;
    for (size_t i = 0; i < N; i++) {
        grid.push_back(min + i*step);
    }
    return grid;
}


template <class T>
T parameterize_interval(T a, T b, T t) {
    return a + t*t*t*exp(t-1)*(b-a);
}


/**
 *  Function that thats a function R^2 -> R, along with grid parameters and computes the surface necessary to draw
 *  in OpenGL.
 * 
 * */
template <class T>
std::vector<std::vector<T> > compute_surface( T (*f)(T,T), std::vector<T> z_grid, std::vector<T> eta_grid) {
    std::vector<T> Frame_line;
    std::vector<std::vector<T> > Figure;// Figure(z_grid.size(),std::vector<T>(eta_grid.size()));
    // Loop for surface frames. Compute f(x,y) at each x,y in grid.
    for (size_t i = 0; i < z_grid.size(); i++) {
        // Naive parallelization. Maybe update this block?
        // #pragma omp parallel for schedule(static)
        for (size_t j = 0; j < eta_grid.size(); j++) {
            // Figure[i][j] = f(z_grid[i],eta_grid[j]);  
            Frame_line.push_back(f(z_grid[i],eta_grid[j]));
        }
        Figure.push_back(Frame_line);
        Frame_line.clear();
    }
    return Figure;
}

template <class T>
std::vector<std::vector<T> > removeNAN(std::vector<std::vector<std::vector<T> > >& input, int k = 0) {
    // SPECIAL NOTE: Need to fix up the way the full solution is handled.
    //               Currently it has 1 too many dimensions in the vector.
    
    // Iterate along all solution curves for each z, stop before NaN values appear.
    // Return only the solution curve without NaN values.

    // k in the arguments picks the solution curve of the N-Dim ode to return. i.e. the nth variable solution
    std::vector<std::vector<T> > noNAN; 
    std::vector<T> section;

    for (size_t i = 0; i < input.size(); i++)
    {
        for (size_t j = 0; j < input[0].size(); j++)
        {
            if (!isnan(input[i][j][k]))
            {
                section.push_back(input[i][j][k]);
            }
        }
        noNAN.push_back(section);
        section.clear();
    }
    return noNAN;
}



class ZeroSetFinder
{
    /*
        Create an N dimensional grid,

        operator() should take a function of N variables
        -- issue doing this..

        iterate through all grids to construct a zero set

        return the zero set as a column of points in N-D space

        The domain input:
            NxM vector: N dimensions(rows) by M grid points (columns)
        
    */
    public:
        // Constructor 
        ZeroSetFinder(mp_type (*f_in)(mp_type, mp_type),std::vector<std::vector<mp_type> > Ndomain)
        :Ndomain(Ndomain) {
            N_dim = Ndomain.size();
            std::cout << "Created finder.\n";
            std::cout << "Dimension: " << N_dim << std::endl;
            f = f_in;

            // Grid size:
            std::cout << "Grid sizing: " << std::endl;
            for (int i = 0; i < N_dim; ++i)
            {
                std::cout << "N_X[" << i <<  "] -> " << Ndomain[i].size() << std::endl; 
            }
            // print domain properties to screen
            std::cout << "Domain on R^" << N_dim << std::endl;
            for (int i = 0; i < N_dim; ++i)
            { 
                int end_point = Ndomain[i].size() - 1;
                std::cout << "X[" << i << "] -> [" << Ndomain[i][0] << "," << Ndomain[i][end_point] << "]\n";
            }
            secant_N = 10;

        };

        // Destructor, clear vectors and free memory!
        ~ZeroSetFinder() {
            // Quick vector cleanup just to be sure.
            Zero_set.clear();
            Bracket_left_set.clear();
            Bracket_right_set.clear();
            Ndomain.clear();
            for (int i = 0; i < surface.size(); ++i)
            {
                surface[i].clear();
            }
        };
    
        // Functions 
        void compute_zero_set();
        void Save_zeros(std::string file);
        void print_zero_set();


        // operator overloads
        void operator()(mp_type (*f)(mp_type,mp_type)){
            // Here we want to use the domain and number of dimensions to
            // apply to a function and root find if we can.
            // essentially: bracket, then refine with secant or newton .



        }

        // Variables
        std::string output_file;
        // NB: Store zero sets as vectors of points (vector of vectors)
        //     so that this is generalizable to ND problems.
        std::vector<std::vector<mp_type> > Zero_set,Bracketing_set,Bracket_left_set,Bracket_right_set;
        std::vector<std::vector<mp_type> > Ndomain, surface;
        std::vector<double> Zero_set_trunc;
        size_t N_dim, secant_N;

        // The type of function we want to evaluate the roots on.
        // i.e. f: R^2 -> R
        boost::function<mp_type (mp_type x, mp_type y)> f;

        // Data transformation functions, i.e. transforming one variable to another
        boost::function<mp_type (mp_type x, mp_type y)> g;

        // Functions derivative, possibly. Note this depends on direction.
        // May have to modify this is include a jacobian matrix!
        boost::function<mp_type (mp_type x, mp_type y)> f_der;        



    private:
        // Private functions
        void Construct_file_stream();
        void compute_surface();
        void bracket_roots();

        // Private variables
        std::ofstream output_stream;



};


void ZeroSetFinder::compute_surface() {
    // Compute the function at every point of the domain and save
    // the result in a member variable: surface.
    std::vector<mp_type> Frame_line;
    // Loop for surface frames. Compute f(x,y) at each x,y in grid.
    for (size_t i = 0; i < Ndomain[0].size(); i++) {
        // Naive parallelization. Maybe update this block?
        // #pragma omp parallel for schedule(static)
        for (size_t j = 0; j < Ndomain[1].size(); j++) {
            // Figure[i][j] = f(z_grid[i],eta_grid[j]);  
            Frame_line.push_back(f(Ndomain[0][i],Ndomain[1][j]));
        }
        surface.push_back(Frame_line);
        Frame_line.clear();
    }
}


void ZeroSetFinder::bracket_roots() {

    std::vector<mp_type> bl_z, bl_eta, br_z, br_eta;

    for (size_t i = 0; i < Ndomain[1].size(); ++i) {
        for (size_t j = 0; j < Ndomain[0].size(); ++j) {
            if (j !=0 && j < Ndomain[0].size()-1) {
                if (surface[j-1][i]*surface[j][i] < 0) {
                    bl_z.push_back(Ndomain[0][j-1]);
                    bl_eta.push_back(Ndomain[1][i]);
                    br_z.push_back(Ndomain[0][j]);
                    br_eta.push_back(Ndomain[1][i]);
                }
            }
        }
    }
    Bracketing_set.push_back(bl_z);
    Bracketing_set.push_back(bl_eta);
    Bracketing_set.push_back(br_z);
    Bracketing_set.push_back(br_eta);
}

// void ZeroSetFinder::compute_zero_set() {
//     // 2D implementation of root finding. using secant method
//     std::vector<mp_type> z_zero, eta_zero, time_zero;
    
//     for (size_t i = 0; i < Bracketing_set[1].size(); i++) {
//         z_zero.push_back(secant(f,Bracketing_set[0][i],Bracketing_set[2][i],Bracketing_set[1][i],secant_N));
//         eta_zero.push_back(Bracketing_set[1][i]);
//     }

//     Zero_set.push_back(z_zero);
//     Zero_set.push_back(eta_zero);
// }