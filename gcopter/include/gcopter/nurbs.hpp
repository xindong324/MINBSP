#ifndef NURBS_HPP
#define NURBS_HPP

#include "gcopter/trajectory.hpp"

#include <Eigen/Eigen>

#include <cmath>
#include <vector>
/**
 * @brief nonuniform b-spline construct
 *         given the control points and time intervals
 *         p_ is the degree of the b spline, for minimum jerk degree = 3
 *         p_ + 1 is the order of the bspline
 *         
 *         n_ + 1 is the control points number
 *         m_ is the number of the time interval number
 *         m_ + 1 is the number of the knot 
 *         For a bspline withe n_ + 1 control points and p_ + 1 order
 *         , then m_ = n_ + p_ + 1 
 * 
 */
namespace nurbs{
    // NURBS for s=3 and non-uniform time
    class NURBS_S3U{
    public:
        NURBS_S3U() = default;
        ~NURBS_S3U() {}

    private:
        int N;
        Eigen::Matrix<double, 3, 4> headPVAJ;
        Eigen::Matrix<double, 3, 4> tailPVAJ;
        // control points for B-spline with different dimensions.
        // Each row represents one single control point
        // The dimension is determined by column number
        // e.g. B-spline with N points in 3D space -> Nx3 matrix
        Eigen::Matrix3Xd control_points_;

        int p_, n_, m_;     // p degree, n+1 control points, m = n+p+1
        Eigen::VectorXd u_; // knots vector
        Eigen::VectorXd interval_;   // knot span \delta t

        Eigen::MatrixXd getDerivativeControlPoints();

        double limit_vel_, limit_acc_, limit_ratio_, feasibility_tolerance_; // physical limits and time adjustment ratio

    public:
        Eigen::MatrixXd get_control_points(void){ return control_points_; }

        inline void setParameters(const Eigen::Matrix3Xd &cpts,
                                  const Eigen::VectorXd &ts)
        {
            interval_ = ts;
            control_points_ = cpts;
        }
        
        void getKnot();
        bool getTimeSpan(double &um, double &um_p);
    
        // TODO: change into nonuniform
        inline void parameterizeToBspline(const Eigen::VectorXd &ts, const std::vector<Eigen::Vector3d> &point_set,
                                             const std::vector<Eigen::Vector3d> &start_end_derivative,
                                             Eigen::MatrixXd &ctrl_pts)
        {
            // if (ts <= 0)
            // {
            //     std::cout << "[B-spline]:time step error." << std::endl;
            //     return;
            // }

            if (point_set.size() <= 3)
            {
                std::cout << "[B-spline]:point set have only " << point_set.size() << " points." << std::endl;
                return;
            }

            if (start_end_derivative.size() != 4)
            {
                std::cout << "[B-spline]:derivatives error." << std::endl;
            }

            int K = point_set.size();

            // write A
            Eigen::Vector3d prow(3), vrow(3), arow(3);
            prow << 1, 4, 1;
            vrow << -1, 0, 1;
            arow << 1, -2, 1;

            Eigen::MatrixXd A = Eigen::MatrixXd::Zero(K + 4, K + 2);

            for (int i = 0; i < K; ++i)
            A.block(i, i, 1, 3) = (1 / 6.0) * prow.transpose();

            A.block(K, 0, 1, 3) = (1 / 2.0 / ts) * vrow.transpose();
            A.block(K + 1, K - 1, 1, 3) = (1 / 2.0 / ts) * vrow.transpose();

            A.block(K + 2, 0, 1, 3) = (1 / ts / ts) * arow.transpose();
            A.block(K + 3, K - 1, 1, 3) = (1 / ts / ts) * arow.transpose();

            //cout << "A" << endl << A << endl << endl;

            // write b
            Eigen::VectorXd bx(K + 4), by(K + 4), bz(K + 4);
            for (int i = 0; i < K; ++i)
            {
                bx(i) = point_set[i](0);
                by(i) = point_set[i](1);
                bz(i) = point_set[i](2);
            }

            for (int i = 0; i < 4; ++i)
            {
                bx(K + i) = start_end_derivative[i](0);
                by(K + i) = start_end_derivative[i](1);
                bz(K + i) = start_end_derivative[i](2);
            }

            // solve Ax = b
            Eigen::VectorXd px = A.colPivHouseholderQr().solve(bx);
            Eigen::VectorXd py = A.colPivHouseholderQr().solve(by);
            Eigen::VectorXd pz = A.colPivHouseholderQr().solve(bz);

            // convert to control pts
            ctrl_pts.resize(3, K + 2);
            ctrl_pts.row(0) = px.transpose();
            ctrl_pts.row(1) = py.transpose();
            ctrl_pts.row(2) = pz.transpose();

            // cout << "[B-spline]: parameterization ok." << endl;
        }
    };
}

#endif