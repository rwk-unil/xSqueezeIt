#include <chrono>

void printElapsedTime(const std::chrono::steady_clock::time_point& begin,
                      const std::chrono::steady_clock::time_point& end) {
    std::cerr << "Time elapsed = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s] "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms] "
                << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[us] " << std::endl;
}