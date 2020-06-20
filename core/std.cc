#include <core/std.h>
#include <core/callstack.h>
#include <mutex>

std::mutex g_cout_mutex;

void Check(bool value, string_view message, const char* file, unsigned line, const char* function) {
    if (value) return;

    {
        std::unique_lock lock(g_cout_mutex);
        cout << "Check failed in " << function << " at " << file << ':' << line << " with message: " << message << endl;
        string s;
        Callstack().write(s, {});
        cout << s;
        cout.flush();
    }
    exit(0);
}

void Fail(string_view message, const char* file, unsigned line, const char* function) {
    {
        std::unique_lock lock(g_cout_mutex);
        cout << "Failed in " << function << " at " << file << ':' << line << " with message: " << message << endl;
        string s;
        Callstack().write(s, {});
        cout << s;
        cout.flush();
    }
    exit(0);
}
