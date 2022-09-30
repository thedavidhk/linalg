#ifndef LINALG_VIEW_H
#define LINALG_VIEW_H

namespace LinAlg {

template<typename T>
struct Matrix_view {
    T* Mat;
    size_t size;
    
};

}

#endif
