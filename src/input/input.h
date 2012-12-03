/******************************************************************************
 * src/input/input.h
 *
 * Common functions to generate input strings for main memory and external
 * memory algorithms.
 *
 ******************************************************************************
 * Copyright (C) 2012 Timo Bingmann <tb@panthema.net>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include "debruijn.h"
#include "fibonacci.h"
#include "yuta-tests.h"
#include "skyline.h"

struct datasource
{
    std::string     name;
    size_t          default_size;
};

typedef std::vector<datasource> datasource_list_type;

class InputRandomInt
{
public:

    typedef int		alphabet_type;

    static const int K	= 1024;

    template <typename Container>
    bool generate(const std::string& /* datasource */, Container& string)
    {
        // Fill input vector with random data

        unsigned int xn = 2342358;
        unsigned int string_size = 400000;

        srand(xn);

        std::cout << "testing xn=" << xn << " - N=" << string_size << std::endl;

        string.resize(string_size);

        for(size_t i = 0; i < string_size; ++i)
        {
            xn = rand() % ((1 << 10) - 1);
            //xn = ((xn * 0x06E6A587) + 0xF57347EF) % (1<<10);
            string[i] = xn + 1; // don't add nulls
        }

        string.push_back(0);

        return true;
    }
};

class InputASCII
{
public:

    typedef unsigned char	alphabet_type;

    static const int K	= 256;

    template <typename Container>
    bool generate_artificial(const datasource& datasource, Container& string, size_t reducesize)
    {
        bool fixzero = false;

        string.clear();

        // check for size suffix after datasource name

        std::string dataname = datasource.name;
        size_t size = 0;

        std::string::size_type dashpos = dataname.rfind('-');
        if (dashpos != std::string::npos)
        {
            std::string sizesuffix = dataname.substr(dashpos+1);
            dataname = dataname.substr(0, dashpos);

            // read number and kb/mb/gb suffix

            char* endptr;
            size = strtoul(sizesuffix.c_str(),&endptr,10);
            if (!endptr) return false;

            sizesuffix = endptr;

            if (sizesuffix == "b")
                size *= 1;
            else if (sizesuffix == "kb")
                size *= 1024;
            else if (sizesuffix == "mb")
                size *= 1024*1024;
            else if (sizesuffix == "gb")
                size *= 1024*1024*1024;
            else
                return false;

            if (reducesize) size = reducesize;
        }
        else if (reducesize)
        {
            size = reducesize;
        }
        else
        {
            if ( datasource.default_size == 0 ) {
                std::cout << "Error: unknown size for artificial data source.\n";
                return false;
            }
            size = datasource.default_size;

            if (reducesize) size = reducesize;
        }

        // generate data

        if (dataname == "artificial/aaa")
        {
            string.resize(size);
            for (size_t i = 0; i < string.size(); ++i)
                string[i] = 'a';
        }
        else if (dataname == "artificial/abc")
        {
            string.resize(size);
            for (size_t i = 0; i < string.size(); ++i)
                string[i] = 'a' + (i % 3);
        }
        else if (dataname == "artificial/ascseq20")
        {
            string.resize(size);
            for (size_t i = 0; i < string.size(); ++i)
                string[i] = 'a' + (i % 20);
        }
        else if (dataname == "artificial/ascseq250")
        {
            string.resize(size);
            for (size_t i = 0; i < string.size(); ++i)
                string[i] = 1 + (i % 250);
        }
        else if (dataname == "artificial/descseq250")
        {
            string.resize(size);
            for (size_t i = 0; i < string.size(); ++i)
                string[i] = 250 - (i % 250);
        }
        else if (dataname == "artificial/ascseq2500")
        {
            string.resize(size);
            for (size_t i = 0; i < string.size(); ++i)
                string[i] = 1 + ((i / 10) % 250);
        }
        else if (dataname == "artificial/descseq2500")
        {
            string.resize(size);
            for (size_t i = 0; i < string.size(); ++i)
                string[i] = 250 - ((i / 10) % 250);
        }
        else if (dataname == "artificial/ascseq25000")
        {
            string.resize(size);
            for (size_t i = 0; i < string.size(); ++i)
                string[i] = 1 + ((i / 100) % 250);
        }
        else if (dataname == "artificial/descseq250000")
        {
            string.resize(size);
            for (size_t i = 0; i < string.size(); ++i)
                string[i] = 250 - ((i / 1000) % 250);
        }
        else if (dataname == "artificial/squares250")
        {
            string.resize(size);
            for (size_t i = 0; i < string.size(); ++i)
                string[i] = 1 + ((i*i) % 250);
        }
        else if (dataname == "artificial/rand250")
        {
            string.resize(size);
            srand(913427589);
            for (size_t i = 0; i < string.size(); ++i)
                string[i] = 1 + (rand() % 250);
        }
        else if (dataname == "artificial/rand256")
        {
            string.resize(size);
            srand(913427589);
            for (size_t i = 0; i < string.size(); ++i)
                string[i] = (rand() % 256);
        }
        else if (dataname == "artificial/randrepeat")
        {
            string.resize(size);
            srand(913427589);
            for (size_t i = 0; i < string.size()/2; ++i)
            {
                string[i]
                    = string[i + string.size()/2]
                    = 1 + (rand() % 250);
            }
        }
        else if (dataname == "artificial/rand10")
        {
            string.resize(size);
            srand(913427589);
            for (size_t i = 0; i < string.size(); ++i)
                string[i] = '0' + (rand() % 10);
        }
        else if (dataname == "artificial/debruijn-(10,8)")
        {
            string.resize(1e8);
            debruijn::generate(10,8,string);
            fixzero = true;
        }
        else if (dataname == "artificial/debruijn-(100,4)")
        {
            string.resize(1e8);
            debruijn::generate(100,4,string);
            fixzero = true;
        }
        else if (dataname == "artificial/fibonacci")
        {
            string.resize(size);
            fibonacci::generate(string);
        }
        else if (dataname == "artificial/yutatest1")
        {
            string.resize(size);
            yuta_tests::generate1(string);
            fixzero = true;
        }
        else if (dataname == "artificial/yutatest2")
        {
            string.resize(size);
            yuta_tests::generate2(string);
            fixzero = true;
        }
        else if (dataname == "artificial/yutatest3")
        {
            string.resize(size);
            yuta_tests::generate3(string);
            fixzero = true;
        }
        else if (dataname == "artificial/yutatest4")
        {
            string.resize(size);
            yuta_tests::generate4(string);
            fixzero = true;
        }
        else if (dataname == "artificial/yutatest5")
        {
            string.resize(size);
            yuta_tests::generate5(string);
            fixzero = true;
        }
        else if (dataname == "artificial/skyline")
        {
            string.resize(size);
            skyline::generate(string);
        }
        else
        {
            return false;
        }

        // many SACA do not support 0 in the string.
        if (fixzero)
        {
            for (size_t i = 0; i < string.size(); ++i)
            {
                if (string[i] < 255) ++string[i];
            }
        }

        //string[ string.size()-1 ] = 0;

        return true;
    }

    template <typename Container>
    bool load_corpusfile(const datasource& datasource, Container& string, size_t reducesize)
    {
        std::string dataname = datasource.name;

        std::ostringstream oss;
        if (getenv("HOME")) oss << getenv("HOME");
        else oss << ".";
        oss << "/sac-corpus/" << dataname;

        std::cout << "loading corpus " << dataname;

        if (dataname.size() >= 3 && dataname.substr(dataname.size()-3,3) == ".gz")
        {
            (std::cout << " - unzipping").flush();

            size_t filesize = 0;

            { // extract filesize from filename
                size_t i = dataname.size()-4, j = 1;
                while ( isdigit(dataname[i]) ) {
                    filesize += (dataname[i] - '0') * j;
                    j *= 10; --i;
                }
                if (filesize == 0 || dataname[i] != '.') {
                    std::cout << "\nCould not find extract size in filename " << dataname << "\n";
                    return false;
                }
            }

            if (reducesize && filesize > reducesize) filesize = reducesize;

            gzFile gz = gzopen(oss.str().c_str(), "rb");
            if (!gz) {
                std::cout << "\nCould not open file " << dataname << "\n";
                return false;
            }

            string.resize(filesize);

            const size_t buffersize = 1*1024*1024;
            unsigned char buffer[buffersize];
            size_t wp = 0, rs;

            int rb;
            while ( (rs = std::min(buffersize, filesize - wp),
                     rb = gzread(gz, buffer, rs)) > 0 )
            {
                assert(wp + rb <= string.size());
                std::copy( buffer, buffer + rb, string.begin() + wp);
                wp += rb;
            }

            if (wp == 0 || rb < 0) {
                std::cout << "\nCould not read file " << dataname << "\n";
                return false;
            }
            if (wp != filesize) {
                std::cout << "File size information is incorrect on " << dataname << "\n";
                return false;
            }

            gzclose(gz);

            std::cout << " - size = " << wp << std::endl;
        }
        else
        {
            std::ifstream infile (oss.str().c_str());

            infile.seekg(0, std::ios::end);
            if (!infile.good()) {
                std::cout << "\nCould not open file " << dataname << "\n";
                return false;
            }

            size_t filesize = infile.tellg();

            std::cout << " - size = " << filesize << std::endl;

            if (reducesize && filesize > reducesize) filesize = reducesize;

            string.resize(filesize);

            infile.seekg(0, std::ios::beg);
            if (!infile.good()) {
                std::cout << "Could not seek in file " << dataname << "\n";
                return false;
            }

            const size_t buffersize = 1*1024*1024;
            unsigned char buffer[buffersize];
            size_t wp = 0;

            while (infile.good() && wp < filesize)
            {
                size_t rs = std::min(buffersize, filesize - wp);
                infile.read((char*)buffer, rs);

                assert(wp + infile.gcount() <= string.size());
                std::copy( buffer, buffer + infile.gcount(), string.begin() + wp);

                wp += infile.gcount();
            }

            if (wp != filesize) {
                std::cout << "Could not read file " << dataname << "\n";
                return false;
            }
        }

        if (0)
        {
            for (size_t i = 0; i < string.size(); ++i)
            {
                if (string[i] == 0) {
                    string[i] = 255;
                }
            }
        }

        //string[ string.size()-1 ] = 0;

        return true;
    }

    template <typename Container>
    bool get(const datasource& datasource, Container& string)
    {
        char* endptr = NULL;
        size_t reducesize = getenv("SIZE") ? strtoul(getenv("SIZE"),&endptr,10) : 0;
        if (endptr && *endptr != '\0') reducesize = 0;

        if (datasource.name.substr(0,7) == "simple/")
        {
            std::string charstring = datasource.name.substr(7);

            string.resize( charstring.size() );
            std::copy( charstring.begin(), charstring.end(), string.begin() );
            //string.push_back(0);
        }
        else if (generate_artificial(datasource,string,reducesize))
        {
        }
        else if (load_corpusfile(datasource,string,reducesize))
        {
        }
        else
        {
            std::cout << "cannot find string " << datasource.name << "\n";
            return false;
        }

        if (reducesize && string.size() > reducesize)
        {
            string.resize(reducesize);
            string.push_back(0);	// readd zero character at end
        }

        //string[ string.size() / 2 - 1 ] = 1;

        std::cout << "testing " << datasource.name << " size " << string.size() << "\n";

        if (getenv("WRITEINPUT") && *getenv("WRITEINPUT") == '1')  // write out data into file
        {
            std::ofstream out(datasource.name.c_str());

            const size_t buffersize = 1*1024*1024;
            alphabet_type buffer[buffersize];

            std::cout << "Writing input to " << datasource.name << "\n";

            for (size_t p = 0; p < string.size(); )
            {
                size_t left = std::min<size_t>(buffersize, string.size() - p);

                std::copy( string.begin() + p, string.begin() + p + left, buffer );

                out.write((char*)buffer, left);
                p += left;
            }
        }

        return true;
    }
};

void read_datasource_directory(datasource_list_type& list, const std::string& subdir)
{
    std::string fullpath = getenv("HOME") ? getenv("HOME") : ".";
    fullpath += "/sac-corpus" + subdir;

    struct stat st;

    if (stat(fullpath.c_str(), &st) != 0) return;

    if (S_ISDIR(st.st_mode))
    {
        DIR* d = opendir(fullpath.c_str());
        std::vector<std::string> filelist;
        while( struct dirent* de = readdir(d) )
        {
            if (de->d_name[0] == '.') continue;
            filelist.push_back(de->d_name);
        }
        closedir(d);
        std::sort(filelist.begin(), filelist.end());
        for (size_t i = 0; i < filelist.size(); ++i)
            read_datasource_directory(list, subdir + "/" + filelist[i]);
    }
    else if (S_ISREG(st.st_mode))
    {
        std::string dn = subdir;

        size_t size = st.st_size;

        if (dn[0] == '/') {
            dn.erase(0,1);
        }
        if (dn.substr(dn.size()-3) == ".gz") {
            size = 0;
            size_t i = dn.size()-4, j = 1;
            while ( isdigit(dn[i]) ) {
                size += (dn[i] - '0') * j;
                j *= 10; --i;
            }
        }

        datasource ds = { dn, size };
        list.push_back(ds);
    }
}

datasource_list_type get_datasource_list()
{
    static const datasource artificial[] = {
        { "simple/mmiissiissiippimmiissiissiippi",      0 },
        { "simple/cdcbabddbbbcccabdccacdcbabddbbbcccabdcca", 0},
        // generated strings
        { "artificial/aaa",                     100000000 },
        { "artificial/abc",                     100000000 },
        { "artificial/ascseq20",		100000000 },
        { "artificial/ascseq250",		100000000 },
        { "artificial/ascseq25000",             100000000 },
        { "artificial/squares250",              100000000 },
        { "artificial/rand250",                 100000000 },
        { "artificial/randrepeat",              100000000 },
        { "artificial/rand10",                  100000000 },
        { "artificial/debruijn-(10,8)",         100000000 },
        { "artificial/debruijn-(100,4)",	100000000 },
        { "artificial/fibonacci",		100000000 },
        { "artificial/yutatest1",		100000000 },
        { "artificial/yutatest2",		100000000 },
        { "artificial/yutatest3",		100000000 },
        { "artificial/yutatest4",		100000000 },
        { "artificial/yutatest5",               100000000 },
        { "artificial/skyline",                 100000000 }
    };

    datasource_list_type list;

    for (unsigned int i = 0; i < sizeof(artificial) / sizeof(artificial[0]); ++i)
    {
        list.push_back(artificial[i]);
    }

    read_datasource_directory(list, "");

    return list;
};

bool select_datasource_list(int argc, char* argv[], datasource_list_type& datasource_list)
{
    datasource_list = get_datasource_list();        // default list

    // parse command line parameters

    if (argc == 1)
    {
        std::cout << "Available data sources:\n";

        unsigned int i = 0;
        for (datasource_list_type::const_iterator ds = datasource_list.begin();
             ds != datasource_list.end(); ++ds, ++i)
        {
            std::cout << std::setw(4) << i
                      << "  " << std::left << std::setw(60) << ds->name
                      << " " << std::right << std::setw(20) << ds->default_size << "\n";
        }
        return false;
    }
    else
    {
        char *endptr;
        int seqstart = strtol(argv[1],&endptr,10);
        if (*endptr) {
            // argument is not a number
            datasource_list.clear();
            for (int ai = 1; ai < argc; ++ai)
            {
                datasource ds = { argv[ai], 0 };
                datasource_list.push_back(ds);
            }
        }
        else if (argc == 3)
        {
            int seqend = strtol(argv[2],&endptr,10);
            if (*endptr) {
                std::cout << "Bad number for sequence end\n";
            }
            else
            {
                seqend = std::max<size_t>(seqend, datasource_list.size()-1);
                datasource_list_type newlist ( datasource_list.begin() + seqstart,
                                               datasource_list.begin() + seqend+1 );
                datasource_list.swap(newlist);
            }
        }
        return true;
    }
}
