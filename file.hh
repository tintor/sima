#ifndef __FILE_H__
#define __FILE_H__

#include <string_view>

class MappedFile {
public:
    MappedFile(const char* fname);

    MappedFile(const MappedFile&) = delete;
    MappedFile& operator=(const MappedFile&) = delete;

    ~MappedFile();

    std::string_view view() const { return std::string_view(m_addr, m_length); }

private:
    int m_fd;
    size_t m_length;
    const char* m_addr;
};

class FileReader {
public:
    FileReader(const char* fname) : m_file(fname), m_pos(m_file.view().begin()) { }

    bool getline(std::string_view& line);

private:
    MappedFile m_file;
    const char* m_pos;
};

#endif
