#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>

std::string human_readable_size(size_t size)
{
    std::string suffix[] = {"B", "KB", "MB", "GB", "TB"};
    int i = 0;
    double dblBytes = static_cast<double>(size);
    if (size > 1024)
    {
        for (i = 0; (size / 1024) > 0 && i < 5; i++, size /= 1024)
        {
            dblBytes = size / 1024.0;
        }
    }
    char output[200];
    sprintf(output, "%.02lf %s", dblBytes, suffix[i].c_str());
    return std::string(output);
}

#endif // UTILS_HPP