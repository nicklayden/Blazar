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

template <class T>
T steffensen(T (*f)(T,T), T x0, T x1, T eta, size_t iter) {
    // Implement steffensens method for root finding.
}

template <class T>
T secant(T (*f)(T,T), T x0, T x1, T eta, size_t iter) {
    // 4 iterations of the secant method.
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
    for (size_t i = 0; i < matrix.size(); i++) {
        for (size_t j = 0; j < matrix[i].size(); j++) {
            ofile << matrix[i][j] << " "; 
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

    result.push_back(z_zero);
    result.push_back(eta_zero);
    return result;
}

template <class T>
std::vector<std::vector<T> > bracket_roots(std::vector<std::vector<T> > surface, std::vector<T> z_grid, std::vector<T> eta_grid) {
    std::vector<T> bl_z, bl_eta, br_z, br_eta;
    std::vector<std::vector<T> > result;

    for (size_t i = 0; i < eta_grid.size(); ++i) {
        for (size_t j = 0; j < z_grid.size(); ++j) {
            if (j !=0 && j < z_grid.size()-1) {
                if (surface[j-1][i]*surface[j][i] < 0) {
                    bl_z.push_back(z_grid[j-1]);
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
std::vector<std::vector<T> > removeNAN(std::vector<std::vector<std::vector<T> > >& input) {
    // SPECIAL NOTE: Need to fix up the way the full solution is handled.
    //               Currently it has 1 too many dimensions in the vector.
    
    // Iterate along all solution curves for each z, stop before NaN values appear.
    // Return only the solution curve without NaN values.
    std::vector<std::vector<T> > noNAN; 
    std::vector<T> section;

    for (size_t i = 0; i < input.size(); i++)
    {
        for (size_t j = 0; j < input[0].size(); j++)
        {
            if (!isnan(input[i][j][0]))
            {
                section.push_back(input[i][j][0]);
            }
        }
        noNAN.push_back(section);
        section.clear();
    }
    return noNAN;
}