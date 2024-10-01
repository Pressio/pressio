// Copyright(c) 2016 Alexander Dalshov.
// Distributed under the MIT License (http://opensource.org/licenses/MIT)

#ifndef QR_LOGGER_SPDLOG_SINKS_UNUSED_MSVC_SINK_H_
#define QR_LOGGER_SPDLOG_SINKS_UNUSED_MSVC_SINK_H_

#if defined(_WIN32)

#include "../details/null_mutex.h"
#include "./base_sink.h"

#include <mutex>
#include <string>


// Avoid including windows.h (https://stackoverflow.com/a/30741042)
extern "C" __declspec(dllimport) void __stdcall OutputDebugStringA(const char *lpOutputString);

namespace spdlog {
namespace sinks {
/*
 * MSVC sink (logging using OutputDebugStringA)
 */
template<typename Mutex>
class msvc_sink : public base_sink<Mutex>
{
public:
    msvc_sink() = default;

protected:
    void sink_it_(const details::log_msg &msg) override
    {
        memory_buf_t formatted;
        base_sink<Mutex>::formatter_->format(msg, formatted);
        OutputDebugStringA(fmt::to_string(formatted).c_str());
    }

    void flush_() override {}
};

using msvc_sink_mt = msvc_sink<std::mutex>;
using msvc_sink_st = msvc_sink<details::null_mutex>;

using windebug_sink_mt = msvc_sink_mt;
using windebug_sink_st = msvc_sink_st;

} // namespace sinks
} // namespace spdlog

#endif
#endif  // QR_LOGGER_SPDLOG_SINKS_UNUSED_MSVC_SINK_H_
