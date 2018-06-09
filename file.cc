#include "file.hh"
#include <sys/stat.h>
#include <iostream>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>

MappedFile::MappedFile(const char* fname) {
	m_fd = open(fname, O_RDONLY);
	if (m_fd == -1)
		throw std::runtime_error("open");

	struct stat sb;
	if (fstat(m_fd, &sb) == -1)
		throw std::runtime_error("fstat");

	m_length = sb.st_size;
	m_addr = static_cast<const char*>(mmap(NULL, m_length, PROT_READ, MAP_PRIVATE, m_fd, 0u));
	if (m_addr == MAP_FAILED)
		throw std::runtime_error("mmap");
}

MappedFile::~MappedFile() {
	munmap(const_cast<char*>(m_addr), m_length);
	close(m_fd);
}

bool FileReader::getline(std::string_view& line) {
	if (m_pos == m_file.view().end())
		return false;

	const char* b = m_pos;
	const char* e = m_file.view().end();

	while (m_pos < e) {
		if (*m_pos == '\n') {
			line = std::string_view(b, m_pos - b);
			m_pos += 1;
			return true;
		}
		m_pos += 1;
	}
	line = std::string_view(b, m_pos - b);
	return true;
}
