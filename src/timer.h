
#ifndef TIMER_H_
#define TIMER_H_

#include <chrono>
#include <stack>
#include <map>
#include <sstream>
#include <iomanip>

class Timer {
private:
  struct adjusted_time_t {
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
    std::chrono::duration<double>                               adjustment;
  };
  std::stack<adjusted_time_t> _times;

public:
  void start() {
    std::chrono::duration<double> _zero;
    _times.push(adjusted_time_t {std::chrono::high_resolution_clock::now(), _zero.zero()});
  }
  std::chrono::duration<double> stop() {
    std::chrono::duration<double> _elapsed;

    // If the time stack is empty, there's nothing for us to do.
    if (_times.empty())
      return _elapsed.zero();

    // Calculate the elapsed time and remove the most recent time from the stack.
    auto _stop  = std::chrono::high_resolution_clock::now();
    auto _time  = _times.top(); _times.pop();
    _elapsed = _stop - _time.start_time;

    // Subtract the time for the current op (*before* we make any adjustment)
    // from any enclosing ops that we may be tracking.
    {
      std::stack<adjusted_time_t> _new_times;
      while (!_times.empty()) {
        auto prev_time = _times.top(); _times.pop();
        prev_time.adjustment += _elapsed;
        _new_times.push(prev_time);
      }
      while (!_new_times.empty()) {
        _times.push(_new_times.top());
        _new_times.pop();
      }
    }

    // Return the elapsed time for the current op, minus any adjustment.
    return _elapsed - _time.adjustment;
  }
};

// Add a 'to_string()' method to a (string, duration) map.
class TimeMap : public std::map<std::string, std::chrono::duration<double>> {
public:
  std::string to_string() const {
    std::stringstream ss;
    std::streamsize sz = ss.width();
    for (auto i = begin(); i != end(); ++i) {
      ss
        << ">> "
        << std::setw(20) << std::left << i->first << std::setw(sz)
        << " : "
        << std::setw(11) << std::setprecision(8) << std::right << i->second.count() << std::setw(sz)
        << std::endl;
    }
    return ss.str();
  }
};

#endif /* TIMER_H_ */
