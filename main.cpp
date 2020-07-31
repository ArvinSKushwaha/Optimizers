#include <iostream>
#include "tqdm.h"
#include "optimizer.hpp"
#include <eigen3/Eigen/Core>

double f(Eigen::VectorXd n)
{
	// double x = n(0);
	// double y = n(1);
	// return -abs(sin(x) * cos(y) * exp(abs(1 - n.norm() / M_PI)));

	double sum = 0;
	for (size_t i = 0; i < n.size(); ++i)
	{
		sum += pow(n(i), 4) - 16 * pow(n(i), 2) + 5 * n(i);
	}
	return sum / 2.0;
}

int main()
{
	int timesteps = 400000;

	tqdm bar;
	bar.set_theme_circle();
	std::cout.precision(10);
	SwarmOptimizer optim = SwarmOptimizer(f, 10000, 10, 0.99, 0.01, 0.01, -3, 3);
	optim.initialize();
	optim.score();
	for (size_t i = 0; i < timesteps; ++i)
	{
		optim.move(0.1);
		bar.progress(i, timesteps);
		if (i % 100 == 0)
		{
			int best = std::min_element(optim.y_best, optim.y_best + optim.particle_count) - optim.y_best;
			std::cout << optim.y_best[best] << "\t" << optim.pos_best[best].mean() << std::endl;
		}
	}
	bar.finish();
	return 0;
}
