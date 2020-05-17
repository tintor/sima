#include <core/string_util.h>
#include <core/util.h>
#include <cxxabi.h>

string Demangle(const char* name) {
    int status;
    char* demangled = abi::__cxa_demangle(name, 0, 0, &status);
    string out = demangled;
    free(demangled);
    return out;
}

void PrintTable(cspan<string> rows, char separator, string_view splitter) {
    vector<size_t> columns;
    for (const string& row : rows) {
        auto s = split(row, separator, false);
        if (s.size() > columns.size()) columns.resize(s.size());
        for (size_t i = 0; i < s.size(); i++) columns[i] = max(columns[i], s[i].size());
    }

    for (const string& row : rows) {
        auto s = split(row, separator, false);
        for (size_t i = 0; i < s.size(); i++) {
            if (i > 0) cout << splitter;
            cout << s[i];
            for (size_t j = s[i].size(); j < columns[i]; j++) cout << ' ';
        }
        cout << '\n';
    }
    cout.flush();
}
