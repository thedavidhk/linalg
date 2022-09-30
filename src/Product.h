#ifndef LINALG_PRODUCT_H
#define LINALG_PRODUCT_H

#include <cassert>

namespace LinAlg {

template<typename _Lhs, typename _Rhs>
class Product {
  public:
    Product(const _Lhs& lhs, const _Rhs& rhs) : m_lhs(lhs), m_rhs(rhs) {
        assert(lhs.cols() == rhs.rows() && "Invalid matrix product.");
    }

  private:
    _Lhs m_lhs;
    _Rhs m_rhs;
};

}

#endif
