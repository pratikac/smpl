#ifndef __utils_h___
#define __utils_h___

#include <sys/time.h>

#define debug(x) \
    std::cout<<"DBG("<<__FILE__<<":"<<__LINE__<<") "<<x<<std::endl

#define SQ(x)   ((x)*(x))
#define RANDF (rand()/(RAND_MAX+1.0))

typedef struct tt{
    struct timeval _time;
    void tic()
    {
        gettimeofday(&_time, NULL);
    }
    double toc() const
    {
        struct timeval t1;
        gettimeofday(&t1, NULL);
        double sec = t1.tv_sec - _time.tv_sec;
        double usec = t1.tv_usec - _time.tv_usec;

        return (sec*1000.0 + usec/1000.0);
    }
}tt;

class range {
    public:
        class iterator {
            friend class range;
            public:
            long int operator *() const { return i_; }
            const iterator &operator ++() { ++i_; return *this; }
            iterator operator ++(int) { iterator copy(*this); ++i_; return copy; }

            bool operator ==(const iterator &other) const { return i_ == other.i_; }
            bool operator !=(const iterator &other) const { return i_ != other.i_; }

            protected:
            iterator(long int start) : i_ (start) { }

            private:
            unsigned long i_;
        };

        iterator begin() const { return begin_; }
        iterator end() const { return end_; }
        range(long int  begin, long int end) : begin_(begin), end_(end) {}
    private:
        iterator begin_;
        iterator end_;
};

#endif
