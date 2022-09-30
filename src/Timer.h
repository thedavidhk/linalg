#ifndef LINALG_TIMER_H
#define LINALG_TIMER_H

#include <chrono>
#include <iostream>
#include <iomanip>

namespace LinAlg {

class Timer {
  public:
    Timer(const char* name)
        : m_name(name), m_t0(std::chrono::system_clock::now()) {}

    ~Timer(void) {
        m_t1 = std::chrono::system_clock::now();
        double duration = std::chrono::duration_cast<std::chrono::nanoseconds>(m_t1 - m_t0).count() / 1000000.0;
        std::cout << std::setprecision(5) << std::fixed;
        std::cout << "Timer[" << m_name << "]: " << duration << " ms" << std::endl;
    }

  private:
    std::string m_name;
    std::chrono::time_point<std::chrono::system_clock> m_t0;
    std::chrono::time_point<std::chrono::system_clock> m_t1;
};

#define LA_FUNC_TIMER Timer _functimer(__func__)

} // namespace Lin_alg
#endif
