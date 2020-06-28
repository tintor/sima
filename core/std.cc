#include <core/std.h>
#include <core/callstack.h>
#include <core/range.h>
#include <mutex>
#include <sstream>

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

std::ostream* fos = &cout;

std::stringstream column;
std::vector<std::string> columns;
int column_section_width = 0;
int column_count = 0;

bool extract_line(std::string& line, std::string& column) {
    int chars = 0;
    // TODO optimize
    size_t read = 0;
    while (chars < column_section_width && read < column.size()) {
        if (column.substr(read, 2) == "\033[") {
            // scan until m
            auto i = column.find('m', read + 2);
            if (i == std::string::npos) break;
            line += column.substr(read, i + 1 - read);
            read = i + 1;
        } else if (column[read] == '\n') {
            read += 1;
            break;
        } else {
            line += column[read];
            read += 1;
            chars += 1;
        }
    }
    column.erase(0, read);
    
    bool result = chars > 0;
    while (chars < column_section_width) {
        line += ' ';
        chars += 1;
    }
    line += ' ';
    return result;
}

bool is_space_only(const std::string& a) {
    for (char c : a) if (c != ' ') return false;
    return true;
}

void flush_columns() {
    std::string line;
    while(true) {
        bool empty = true, empty2 = true;
        FOR(i, column_count) {
            if (i >= columns.size()) break;
            if (columns[i].size() > 0) empty = false;
            if (extract_line(line, columns[i])) empty2 = false;
        }
        if (empty) break;
        if (!empty2) cout << line << endl;
        line.clear();
    }
    if (columns.size() <= column_count) columns.clear(); else columns.erase(columns.begin(), columns.begin() + column_count);
}

void column_section(int width, int count) {
    if (column.str().size() > 0) {
        columns.push_back(std::move(column.str()));
        column = std::stringstream();
    }
    if (columns.size() >= column_count) flush_columns();
    column_section_width = width;
    column_count = count;
    fos = column_section_width ? &column : &cout;
}

void end_column_section() {
    column_section(column_section_width, column_count);
    if (columns.size() > 0) flush_columns();
    columns.clear();
    column_section_width = 0;
    column_count = 0;
    fos = &cout;
}
