#ifndef __FS_HPP__
#define __FS_HPP__

#if __cplusplus >= 201703L
    #include <filesystem>
    namespace fs = std::filesystem;
    using fs::remove;
#else
    #include <stdio.h> // Has remove()
    #include <sys/stat.h>
    namespace fs {
        inline size_t file_size(const std::string& filename) {
            struct stat st;
            if (stat(filename.c_str(), &st) < 0) {
                std::cerr << "Size of file : " << filename << " could not be determined !" << std::endl;
                return 0;
            } else {
                return st.st_size;
            }
        }

        inline bool exists(const std::string& filename) {
            struct stat st;
            return stat(filename.c_str(), &st) == 0;
        }
    };
#endif

#endif /* __FS_HPP__ */