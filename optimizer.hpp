#include <eigen3/Eigen/Core>
#include <functional>
#include <float.h>

class PVSO
{
public:
    std::function<double(Eigen::ArrayXd)> f;
    Eigen::ArrayXd *pos;
    Eigen::ArrayXd *pos_best;
    Eigen::ArrayXd *vel;
    double *y;
    double *y_best;
    size_t particle_count;
    size_t dims;
    double omega = 0.95;
    double phi_p = 0.025;
    double phi_g = 0.025;
    double l_constraint = -1;
    double u_constraint = 1;

    PVSO(std::function<double(Eigen::ArrayXd)> x_, size_t n, size_t dims_)
    {
        f = x_;
        particle_count = n;
        dims = dims_;

        pos = (Eigen::ArrayXd *)malloc(particle_count * sizeof(Eigen::ArrayXd));
        pos_best = (Eigen::ArrayXd *)malloc(particle_count * sizeof(Eigen::ArrayXd));

        vel = (Eigen::ArrayXd *)malloc(particle_count * sizeof(Eigen::ArrayXd));

        y = (double *)malloc(particle_count * sizeof(double));
        y_best = (double *)malloc(particle_count * sizeof(double));
    }

    PVSO(std::function<double(Eigen::ArrayXd)> x_, size_t n, size_t dims_, double w, double pp, double pg)
    {
        f = x_;
        particle_count = n;
        dims = dims_;

        pos = (Eigen::ArrayXd *)malloc(particle_count * sizeof(Eigen::ArrayXd));
        pos_best = (Eigen::ArrayXd *)malloc(particle_count * sizeof(Eigen::ArrayXd));

        vel = (Eigen::ArrayXd *)malloc(particle_count * sizeof(Eigen::ArrayXd));

        y = (double *)malloc(particle_count * sizeof(double));
        y_best = (double *)malloc(particle_count * sizeof(double));

        omega = w;
        phi_g = pg;
        phi_p = pp;
    }

    PVSO(std::function<double(Eigen::ArrayXd)> x_, size_t n, size_t dims_, double lc, double uc)
    {
        f = x_;
        particle_count = n;
        dims = dims_;

        pos = (Eigen::ArrayXd *)malloc(particle_count * sizeof(Eigen::ArrayXd));
        pos_best = (Eigen::ArrayXd *)malloc(particle_count * sizeof(Eigen::ArrayXd));

        vel = (Eigen::ArrayXd *)malloc(particle_count * sizeof(Eigen::ArrayXd));

        y = (double *)malloc(particle_count * sizeof(double));
        y_best = (double *)malloc(particle_count * sizeof(double));

        l_constraint = lc;
        u_constraint = uc;
    }

    PVSO(std::function<double(Eigen::ArrayXd)> x_, size_t n, size_t dims_, double w, double pp, double pg, double lc, double uc)
    {
        f = x_;
        particle_count = n;
        dims = dims_;

        pos = (Eigen::ArrayXd *)malloc(particle_count * sizeof(Eigen::ArrayXd));
        pos_best = (Eigen::ArrayXd *)malloc(particle_count * sizeof(Eigen::ArrayXd));

        vel = (Eigen::ArrayXd *)malloc(particle_count * sizeof(Eigen::ArrayXd));

        y = (double *)malloc(particle_count * sizeof(double));
        y_best = (double *)malloc(particle_count * sizeof(double));

        omega = w;
        phi_g = pg;
        phi_p = pp;

        l_constraint = lc;
        u_constraint = uc;
    }

    double evaluate(Eigen::ArrayXd n)
    {
        return f(n);
    }

    void initialize()
    {
        for (size_t i = 0; i < particle_count; ++i)
        {
            pos[i] = (u_constraint - l_constraint) * (Eigen::ArrayXd::Random(dims)/2.0 + 0.5) + l_constraint;
            pos_best[i] = pos[i];
            vel[i] = Eigen::ArrayXd::Random(dims);

            y_best[i] = DBL_MAX;
        }
    }

    void score()
    {
        for (size_t i = 0; i < particle_count; ++i)
        {
            y[i] = evaluate(pos[i]);
        }
    }

    void move(double dt)
    {
        Eigen::ArrayXd best_pos;
        double best_val;
        int mindex;

        score();

        mindex = std::min_element(y_best, y_best + particle_count) - y_best;
        best_pos = pos_best[mindex];
        best_val = y_best[mindex];

        for (size_t i = 0; i < particle_count; ++i)
        {
            vel[i] = omega * vel[i] +
                     phi_p * Eigen::ArrayXd::Random(dims).abs() * (pos_best[i] - pos[i]) +
                     phi_g * Eigen::ArrayXd::Random(dims).abs() * (best_pos - pos[i]);

            pos[i] = pos[i] + dt * vel[i];

            pos_best[i] = y[i] < y_best[i] ? pos[i] : pos_best[i];
            y_best[i] = y[i] < y_best[i] ? y[i] : y_best[i];

            pos[i] = pos[i].max(-10).min(10);
        }
    }
};