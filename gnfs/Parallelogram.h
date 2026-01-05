#ifndef PARALLELOGRAM_H
#define PARALLELOGRAM_H
#include <utility>
#include <iostream>
#include <cmath>

// A parallelogram can be defined by:
// i) 3 points (the 4th point can be derived)
// or equivalently
// ii) 2 vectors and a starting point
// i.e. 6 values are enough to define it

template <class T>
struct Point
{
    T x;
    T y;
    Point()
    {}
    Point(T u, T v) : x(u), y(v)
    {}
    bool operator<(const Point<T>& pt)
    {
        return (x < pt.x);
    }
};

template <class T>
struct Vector
{
    T x;
    T y;
    Vector()
    {}
    Vector(T u, T v) : x(u), y(v)
    {}
    Vector(const Point<T>& p1, const Point<T>& p2) : x(p2.x - p1.x), y(p2.y - p1.y)
    {}
};

template <class T>
inline Point<T> operator+(const Point<T>& p, const Vector<T>& v)
{
    return Point<T>(p.x + v.x, p.y + v.y);
}

class Parallelogram
{
public:
    friend class LatticeSieverTest;
    Parallelogram()
    {}

    Parallelogram(const Point<double>& p1, const Point<double>& p2, const Point<double>& p3)
        : p1_(p1), p2_(p2), p3_(p3)
    {
        p4_ = p2_ + Vector<double>(p1_, p3_);
        sort_points();
        init();
    }

    Parallelogram(const Point<double>& p1, const Vector<double>& v1, const Vector<double>& v2)
        : p1_(p1)
    {
        p2_ = p1_ + v1;
        p3_ = p1_ + v2;
        p4_ = p2_ + v2;
        sort_points();
        init();
    }

    // create another parallelogram by transforming using basis given by b1, b2
    Parallelogram(const Parallelogram& p, const std::pair<int32_t, int32_t>& b1, const std::pair<int32_t, int32_t>& b2)
    {
        int32_t det = b1.first * b2.second - b2.first * b1.second;
        p1_ = transform(p.p1_, b1, b2, det);
        p2_ = transform(p.p2_, b1, b2, det);
        p3_ = transform(p.p3_, b1, b2, det);
        p4_ = p2_ + Vector<double>(p1_, p3_);
        sort_points();
        init();
    }

    ~Parallelogram()
    {}

    double min_x() const
    {
        return min_x_;
    }
    double max_x() const
    {
        return max_x_;
    }

    // for a given value of x, find the minimum and maximum values of y
    inline bool y_limits(double x, long int& y_min, long int& y_max)
    {
        // assumption: points are sorted so p1_.x <= p2_.x <= p3_.x <= p4_.x
        double min_y;
        double max_y;
        if (x < p1_.x)
            return false;

        if (x <= p2_.x)
        {
            if (p1_.x == p2_.x)
            {
                if (p1_.y < p2_.y)
                {
                    min_y = p1_.y;
                    max_y = p2_.y;
                }
                else
                {
                    max_y = p1_.y;
                    min_y = p2_.y;
                }
            }
            else
            {
                // we intersect p1 -> p2 and p1 -> p3
                min_y = x * y12_ + d12_;
                max_y = min_y;
                double y = x * y13_ + d13_;
                if (y < min_y)
                    min_y = y;
                else if (y > max_y)
                    max_y = y;
            }
        }
        else if (x <= p3_.x)
        {
            // we intersect p1 -> p3 and p2 -> p4
            min_y = x * y13_ + d13_;
            max_y = min_y;
            double y = x * y24_ + d24_;
            if (y < min_y)
                min_y = y;
            else if (y > max_y)
                max_y = y;
        }
        else if (x <= p4_.x)
        {
            // we intersect p2 -> p4 and p3 -> p4
            min_y = x * y24_ + d24_;
            max_y = min_y;
            double y = x * y34_ + d34_;
            if (y < min_y)
                min_y = y;
            else if (y > max_y)
                max_y = y;
        }
        else
            return false;// if (x > p4_.x) return false;

        constexpr double epsilon = 1e-10;
        double min_y_ceil = std::ceil(min_y - epsilon);
        if (min_y_ceil > max_y + epsilon)
        {
            //std::cerr << "y_limits() failed : x = " << x << ", min_y = " << min_y << ", max_y = " << max_y << std::endl;
            return false;
        }
        //std::cerr << "y_limits() : x = " << x << ", min_y = " << min_y << ", max_y = " << max_y << ", max_y - floor(max_y) = ";
        y_min = static_cast<long int>(min_y_ceil);
        y_max = static_cast<long int>(std::floor(max_y + epsilon));
        return true;
    }

    inline bool y_limits1(double x, int32_t& y_min, int32_t& y_max)
    {
        // assumption: points are sorted so p1_.x <= p2_.x <= p3_.x <= p4_.x
        double min_y;
        double max_y;
        
        // Early exit for out-of-bounds
        if (x < p1_.x || x > p4_.x)
            return false;
        
        // Calculate min_y and max_y based on x position
        if (x <= p2_.x)
        {
            // we intersect p1 -> p2 and p1 -> p3
            min_y = x * y12_ + d12_;
            double y = x * y13_ + d13_;
            max_y = (y > min_y) ? y : min_y;
            min_y = (y < min_y) ? y : min_y;
        }
        else if (x <= p3_.x)
        {
            // we intersect p1 -> p3 and p2 -> p4
            min_y = x * y13_ + d13_;
            double y = x * y24_ + d24_;
            max_y = (y > min_y) ? y : min_y;
            min_y = (y < min_y) ? y : min_y;
        }
        else // x <= p4_.x (already checked above)
        {
            // we intersect p2 -> p4 and p3 -> p4
            min_y = x * y24_ + d24_;
            double y = x * y34_ + d34_;
            max_y = (y > min_y) ? y : min_y;
            min_y = (y < min_y) ? y : min_y;
        }

        constexpr double epsilon = 1e-10;
        double min_y_ceil = std::ceil(min_y - epsilon);
        if (min_y_ceil > max_y + epsilon)
        {
            //std::cerr << "y_limits() failed : x = " << x << ", min_y = " << min_y << ", max_y = " << max_y << std::endl;
            return false;
        }
        //std::cerr << "y_limits() : x = " << x << ", min_y = " << min_y << ", max_y = " << max_y << ", max_y - floor(max_y) = ";
        //y_min = static_cast<long int>(min_y_ceil) - 1;
        //y_max = static_cast<long int>(std::floor(max_y + epsilon)) - 1;
        y_min = static_cast<int32_t>(min_y_ceil);
        y_max = static_cast<int32_t>(std::floor(max_y + epsilon));
        return true;
    }

    void display(std::ostream& os)
    {
        os << "p1_ = (" << p1_.x << "," << p1_.y << ")" << std::endl;
        os << "p2_ = (" << p2_.x << "," << p2_.y << ")" << std::endl;
        os << "p3_ = (" << p3_.x << "," << p3_.y << ")" << std::endl;
        os << "p4_ = (" << p4_.x << "," << p4_.y << ")" << std::endl;
    }

private:
    Point<double> p1_;
    Point<double> p2_;
    Point<double> p3_;
    Point<double> p4_;
    double min_x_;
    double max_x_;
    double y12_;
    double y13_;
    double y24_;
    double y34_;
    double d12_;
    double d13_;
    double d24_;
    double d34_;
    void init()
    {
        min_x_ = p1_.x;
        max_x_ = p4_.x;
        // pre-calculate stuff for calculating intersection
        double x12 = p1_.x - p2_.x;
        double x13 = p1_.x - p3_.x;
        double x24 = p2_.x - p4_.x;
        double x34 = p3_.x - p4_.x;
        y12_ = (p1_.y - p2_.y) / x12;
        y13_ = (p1_.y - p3_.y) / x13;
        y24_ = (p2_.y - p4_.y) / x24;
        y34_ = (p3_.y - p4_.y) / x34;
        d12_ = (p1_.x * p2_.y - p1_.y * p2_.x) / x12;
        d13_ = (p1_.x * p3_.y - p1_.y * p3_.x) / x13;
        d24_ = (p2_.x * p4_.y - p2_.y * p4_.x) / x24;
        d34_ = (p3_.x * p4_.y - p3_.y * p4_.x) / x34;
    }
    void sort_points()
    {
        if (p2_ < p1_) std::swap(p1_, p2_);
        if (p3_ < p2_)
        {
            std::swap(p2_, p3_);
            if (p2_ < p1_) std::swap(p1_, p2_);
        }
        if (p4_ < p3_)
        {
            std::swap(p3_, p4_);
            if (p3_ < p2_)
            {
                std::swap(p2_, p3_);
                if (p2_ < p1_) std::swap(p1_, p2_);
            }
        }
    }

public:
    static Point<double> transform(const Point<double>& pt,
                                   const std::pair<int32_t, int32_t>& b1,
                                   const std::pair<int32_t, int32_t>& b2,
                                   int32_t det)
    {
        Point<double> new_pt;
        new_pt.x = static_cast<double>((b2.second * static_cast<double>(pt.x) - b2.first * static_cast<double>(pt.y)) / static_cast<double>(det));
        new_pt.y = static_cast<double>((b1.first * static_cast<double>(pt.y) - b1.second * static_cast<double>(pt.x)) / static_cast<double>(det));
        return new_pt;
    }

};
#endif
