#ifndef __FS_HPP__
#define __FS_HPP__

#include <libgen.h> // Has dirname() / basename()
/*
The functions dirname() and basename() break a null-terminated pathname string
into directory and filename components. In the usual case, dirname() returns the
string up to, but not including, the final '/', and basename() returns the
component following the final '/'. Trailing '/' characters are not counted as
part of the pathname.
*/

#if __cplusplus >= 201703L
    #include <filesystem>
    namespace fs = std::filesystem;
    using fs::remove;
    // This can be used instead of basename() (more portable, if C++17 is available)
    // fs::path( "/foo/bar.txt" ).filename() => "bar.txt"
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