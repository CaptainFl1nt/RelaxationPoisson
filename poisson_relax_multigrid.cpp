#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

// populates the boundaries matrix with the boundary geometry
void createBoundaries(bool** bounds, int size, double x0, double x1, double y0, double y1) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            bounds[i][j] = false;
        }
    }
    for (int i = 0; i < size; i++) {
        bounds[0][i] = true;
        bounds[size-1][i] = true;
        //bounds[i][0] = true;
        //bounds[i][size-1] = true;
    }
    /*
    double x,yp,yn;
    int j1, j2;
    double h = (x1 - x0)/(size-1);
    for (int i = 0; i < size; i++) {
        x = x0 + i*h;
        if (abs(x) <= 0.5) {
            yp = 0.25*sqrt(1-x*x/0.25);
            yn = -yp;
            j1 = (int) ((yp-y0)/h);
            j2 = (int) ((yn-y0)/h);
            bounds[i][j1] = 1;
            bounds[i][j2] = 1;
        }
    }
    */
}

// enforces the boundary conditions
void imposeBCs(double** phi, int size, double x0, double x1, double y0, double y1) {
    double h = (x1 - x0)/(size-1);
    double x,y;
    
    for (int i = 0; i < size; i++) {
        phi[0][i] = 1;
        phi[size-1][i] = -1;
    }
    /*
    for (int i = 0; i < size; i++) {
        x = x0 + i*h;
        phi[i][0] = 0;
        phi[i][size-1] = 0;
    }*/
    /*
    double yp,yn;
    int j1, j2;
    for (int i = 0; i < size; i++) {
        x = x0 + i*h;
        if (abs(x) <= 0.5) {
            yp = 0.25*sqrt(1-x*x/0.25);
            yn = -yp;
            j1 = (int) ((yp-y0)/h);
            j2 = (int) ((yn-y0)/h);
            phi[i][j1] = 1;
            phi[i][j2] = 1;
        }
    }
    */
}


// right hand side of Poisson's equation
// equal to -rho/e0
double g(double x, double y) {
    return 0;
}


int main() {
    // get input information
    // Need a file of the form
    // ******************
    // start_size levels
    // tolerance
    // x0 x1
    // y0 y1
    //*******************
    ifstream input;
    input.open("input.txt");
    
    int size;
    int levels;
    double tolerance,x0,x1,y0,y1;
    input >> size;
    input >> levels;
    input >> tolerance;
    input >> x0 >> x1;
    input >> y0 >> y1;
    
    cout << levels << " " << size << endl;
    cout << tolerance << endl;
    cout << x0 << " " << x1 << endl;
    cout << y0 << " " << y1 << endl;
    
    input.close();
    
    int MAX_ITERATIONS = 1000000;
    
    int max_size = (pow(2,levels-1))*size - (pow(2,levels-1)-1);
    cout << max_size << endl;
    
    // Initialize data arrays
    bool** bounds = new bool*[max_size];
    double** phi = new double*[max_size];
    double** newPhi = new double*[max_size];
    double** save;
    for (int i = 0; i < max_size; i++) {
        bounds[i] = new bool[max_size];
        phi[i] = new double[max_size];
        newPhi[i] = new double[max_size];
        for (int j = 0; j < max_size; j++) {
            bounds[i][j] = false;
            phi[i][j] = 0;
            newPhi[i][j] = 0;
        }
    }
    
    
    double error, count, oldV, newV, diff, h, x, y;
    // Loop of each sucessive level
    for (int lvl = 0; lvl < levels; lvl++) {
        cout << "lvl: " << lvl << ", size: " << size << ", count: ";
        // set up boundary conditions.
        createBoundaries(bounds, size, x0, x1, y0, y1);
        imposeBCs(phi, size, x0, x1, y0, y1);
        
        error = 1;
        count = 0;
        // loop until convergence
        while (error > tolerance && count < MAX_ITERATIONS) {
            h = (x1 - x0)/(size - 1);
            error = 0;
            // loop over vector
            for (int i = 0; i < size; i++ ) {
                for (int j = 0; j < size; j++) {
                    // is not a boundary point
                    if (!bounds[i][j]) {
                        x = x0 + i*h;
                        y = y0 + j*h;
                        oldV = phi[i][j];
                        newV = (phi[i+1][j]+phi[i-1][j]+phi[i][j+1]+phi[i][j-1] - h*h*g(x,y))/4.0;
                        phi[i][j] = newV;
                        diff = oldV - newV;
                        error += diff*diff;
                    }
                }
            }
            count++;
            error = sqrt(error);
        }
        
        cout << count << endl;
        
        // stop at last level without resize
        if (lvl == levels-1) {
            break;
        }
        
        // interpolate data to larger grid
        for (int i = 0; i < size-1; i++) {
            for (int j = 0; j < size-1; j++) {
                newPhi[2*i][2*j] = phi[i][j];
                newPhi[2*i+1][2*j] = (phi[i][j]+phi[i+1][j])/2.0;
                newPhi[2*i][2*j+1] = (phi[i][j]+phi[i][j+1])/2.0;
                newPhi[2*i+1][2*j+1] = (phi[i][j]+phi[i+1][j]+phi[i][j+1]+phi[i+1][j+1])/4.0;
            }
        }
        for (int i = 0; i < size - 1; i++) {
            newPhi[2*i][2*size-2] = phi[i][size-1];
            newPhi[2*i+1][2*size-2] = (phi[i][size-1]+phi[i+1][size-1])/2.0;
            newPhi[2*size-2][2*i] = phi[size-1][i];
            newPhi[2*size-2][2*i+1] = (phi[size-1][i]+phi[size-1][i+1])/2.0;
        }
        newPhi[2*size-2][2*size-2] = phi[size-1][size-1];
        
        size = 2*size-1;
        
        save = phi;
        phi = newPhi;
        newPhi = save;
        
    }
    
    // output data
    ofstream output;
    output.open("mg_out.csv");
    cout << "size: " << size << endl;
    for (int i = 0; i < size; i++ ) {
        for (int j = 0; j < size-1; j++) {
            output << phi[i][j] << ", ";
        }
        output << phi[i][size-1] << endl;
    }
    output.close();
    
    return 0;
}
