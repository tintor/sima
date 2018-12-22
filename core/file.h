#pragma once
#include <string_view>

class MappedFile {
public:
    MappedFile(std::string_view filename);
    MappedFile(const MappedFile&) = delete;
    MappedFile& operator=(const MappedFile&) = delete;
    ~MappedFile();

    std::string_view view() const { return std::string_view(m_addr, m_length); }

private:
    int m_fd;
    size_t m_length;
    const char* m_addr;
};

// TODO this can just be a string iterator (detached from file)
class FileReader {
public:
    FileReader(std::string_view filename) : m_file(filename), m_pos(m_file.view().begin()) { }

    std::string_view readline();

private:
    MappedFile m_file;
    const char* m_pos;
};
