# Optimizers 👀

Optimizers is a header-only library for C++ to optimize low-dimensional functions.

The following optimization algorithms are included/will be included:

 - [x] Vanilla Particle Swarm Optimization (VPSO)
 - [ ] Stochastic Gradient Descent (SGD)
 - [ ] RMSprop
 - [ ] Adam
 - [ ] Nadam
 - [ ] L-BFGS
 - [ ] Newton's Method

Note: this library is not meant for production use, but rather for testing and research purposes only.


## Requirements:
 - Eigen (only tested with Eigen 3.3.4)


Here is an example of using the `VPSO` class:

```c++
#include <iostream>
#include <Eigen/Core>
#include "optimizer.hpp"

double f(Eigen::ArrayXd x)
{
    return (x * x).sum(); // Equivalent to x^2 + y^2 (in two dimensions)
}

int main()
{
    int num_particles = 10000;
    int dims = 2;
    int iterations = 4000;
    
    VPSO optimizer = VPSO(f, num_particles, dims);
    optimizer.initialize(); // Just initializes the particles for VPSO
    for (int i = 0; i < iterations; ++i)
    {
        optimizer.score(); // Calculates the "loss" of each particle
        optimizer.move(1.0); // Updates the particle (timestep = 1.0)
    }

    std::cout << optimizer.record_val << "\t" \
            << optimizer.record_pos.transpose() << std::endl; 
    // Output lowest "loss" and corresponding input.
    
    return 0;
}
```