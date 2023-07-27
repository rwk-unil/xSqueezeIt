/*#include "bcf_traversal.hpp"
#include "CLI11.hpp"

class PPExtractTraversal : public BcfTraversal
{
public:
    float *gp_arr;
    int gp_arr_size;
    bool extract;
    std::string ofname;

    PPExtractTraversal(bool extract = false, std::string ofname = "") : gp_arr(NULL),
                                                                        gp_arr_size(0), extract(extract), ofname(ofname)
    {
    }

    virtual ~PPExtractTraversal()
    {
        // Free memory
        if (gp_arr)
        {
            free(gp_arr);
        }
    }

    virtual void handle_bcf_file_reader() override
    {
        std::cout << "Starting extraction" << std::endl;
        gp_map = std::map<Triplet, int>();
    }

    virtual void handle_bcf_line() override
    {
        auto line = bcf_fri.line;
        auto header = bcf_fri.sr->readers[0].header;

        int res = bcf_get_format_float(header, line, "GP", &gp_arr, &gp_arr_size);

        if (res <= 0)
        {
            std::cerr << "Could not extract GPs !" << std::endl;
            exit(-1);
        }

        for (size_t i = 0; i < bcf_fri.n_samples; ++i)
        {
            // std::cout << "[" << gp_arr[i * 3] << "," << gp_arr[i * 3 + 1] << "," << gp_arr[i * 3 + 2] << "]\t";
            if (!extract)
            {
                Triplet t = Triplet(gp_arr[i * 3], gp_arr[i * 3 + 1], gp_arr[i * 3 + 2]);
                if (gp_map.find(t) == gp_map.end())
                {
                    gp_map[t] = 1;
                }
                else
                {
                    gp_map[t]++;
                }
            }
        }
    }

    void show_results()
    {
        for (auto &kv : gp_map)
        {
            std::cout << kv.first << "\t" << kv.second << std::endl;
        }
        std::cout << std::endl;
    }

private:
    class Triplet
    {
    public:
        Triplet(float a, float b, float c) : a(a), b(b), c(c) {}
        float a;
        float b;
        float c;

        friend std::ostream &operator<<(std::ostream &os, const Triplet &t)
        {
            os << "[" << t.a << "," << t.b << "," << t.c << "]";
            return os;
        }

        bool operator<(const Triplet &other) const
        {
            return (a < other.a) || (a == other.a && b < other.b) || (a == other.a && b == other.b && c < other.c);
        }
    };

    std::map<Triplet, int> gp_map;
};

int main(int argc, char **argv)
{

    CLI::App app;
    std::string ifname = "-";
    std::string ofname = "-";
    bool extract = false;
    bool analyze = true;
    app.add_option("-f,--filename", ifname, "VCF/BCF Filename");
    app.add_option("-o,--output", ofname, "Output Filename");
    app.add_flag("-a, --analyze", analyze, "Analyze the file");
    app.add_flag("-e, --extract", extract, "Extract the file");
    CLI11_PARSE(app, argc, argv);

    if (analyze && extract)
    {
        std::cerr << "Cannot analyze and extract at the same time !" << std::endl;
        exit(app.exit(CLI::CallForHelp()));
    }

    if (analyze)
    {
        std::cout << "Analyzing file " << ifname << std::endl;
        PPExtractTraversal ppet;
        ppet.traverse(ifname);
        ppet.show_results();
    }
    else
    {
        std::cout << "Extracting file " << ifname << std::endl;
        PPExtractTraversal ppet(extract, ofname);
        ppet.traverse(ifname);
    }

    return 0;
}*/