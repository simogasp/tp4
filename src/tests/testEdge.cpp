/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "../core.hpp"
#include "../ObjModel.hpp"
#include <cstdio>
#include <map>
#include <optional>
#include <regex>
#include <sstream>
#include <unordered_map>

std::unordered_map< int, int > map;

bool parseFaceString( const std::string &toParse, face &out )
{
//  PRINTVAR( toParse );
    if ( toParse.c_str( )[0] == 'f' )
    {
        idxtype a;
        // now check the different formats: %d, %d//%d, %d/%d, %d/%d/%d
        if ( strstr( toParse.c_str( ), "//" ) )
        {
            // v//n
            return ( sscanf( toParse.c_str( ), "f %u//%u %u//%u %u//%u", &(out.v1), &a, &(out.v2), &a, &(out.v3), &a ) == 6);
        }
        else if ( sscanf( toParse.c_str( ), "f %u/%u/%u", &a, &a, &a ) == 3 )
        {
            // v/t/n
            return ( sscanf( toParse.c_str( ), "f %u/%u/%u %u/%u/%u %u/%u/%u", &(out.v1), &a, &a, &(out.v2), &a, &a, &(out.v3), &a, &a ) == 9);
        }
        else if ( sscanf( toParse.c_str( ), "f %u/%u", &a, &a ) == 2 )
        {
            // v/t .
            return ( sscanf( toParse.c_str( ), "f %u/%u %u/%u %u/%u", &(out.v1), &a, &(out.v2), &a, &(out.v3), &a ) == 6);
        }
        else
        {
            // v
//            sscanf( toParse.c_str( ), "f %u %u %u", &(out.v1), &(out.v2), &(out.v3) );

//            PRINTVAR( out );
            return ( sscanf( toParse.c_str( ), "f %u %u %u", &(out.v1), &(out.v2), &(out.v3) ) == 3);
        }
    }
    else
    {
        return false;
    }
}



int main(int, char **)
{
    edge2vertex mylist;


    mylist.insert(std::pair<edge,GLushort>(edge(6,5),6));
    mylist.insert(std::pair<edge,GLushort>(edge(1,2),1));
    mylist.insert(std::pair<edge,GLushort>(edge(2,1),4));
    mylist.insert(std::pair<edge,GLushort>(edge(1,4),5));
    mylist.insert(std::pair<edge,GLushort>(edge(3,5),6));
//     mylist[edge(3,2)] = 10;

    PRINTVAR(mylist);
    PRINTVAR(mylist.size());

    bool res = (mylist.find(edge(5,3))) != mylist.end();
    PRINTVAR(res);

    PRINTVAR(((mylist.find(edge(4,1))) != mylist.end()));
    PRINTVAR(((mylist.find(edge(3,2))) != mylist.end()));
    PRINTVAR(((mylist.find(edge(1,2))) != mylist.end()));
    PRINTVAR(((mylist.find(edge(2,1))) != mylist.end()));

    PRINTVAR(mylist[edge(4,1)]);
    PRINTVAR(mylist[edge(2,1)]);
    PRINTVAR(mylist[edge(3,5)]);
    PRINTVAR(((mylist.find(edge(5,3))) != mylist.end()));
    PRINTVAR(mylist[edge(6,5)]);
    PRINTVAR(((mylist.find(edge(5,6))) != mylist.end()));
    mylist[edge(0,1)]=100;

    PRINTVAR(mylist);

    std::string line("f 34//23 11//65 123//98");
    idxtype a;
    face out;
    std::string b;
//    std::stringstream parser (line);
//    parser >> b
//           >> out.v1 >> b >> a
//           >> out.v2 >> b >> a
//           >> out.v3 >> b >> a;
//    PRINTVAR(parser.fail());
//    std::cout << out << a << b << std::endl;
    PRINTVAR(line);
    PRINTVAR( (sscanf(line.c_str(), "f %u//%u %u//%u %u//%u", &(out.v1), &a, &(out.v2), &a, &(out.v3), &a) == 6) );
    std::cout << out << a << b << std::endl;

    PRINTVAR(v3f(1,2,4)*10.f);
    PRINTVAR(10.f*v3f(1,2,4));

    PRINTVAR((edge(3,5)==edge(3,5)));
    PRINTVAR((edge(3,5)==edge(5,3)));
    PRINTVAR((edge(1,2)==edge(3,5)));

    idxtype idx;

    PRINTVAR((face(1,2,3).containsEdge( edge(3,5), idx )));
    PRINTVAR( idx );
    PRINTVAR((face(1,2,3).containsEdge( edge(2,3), idx )));
    PRINTVAR( idx );
    PRINTVAR((face(1,2,3).containsEdge( edge(3,2), idx )));
    PRINTVAR( idx );
    PRINTVAR((face(1,2,3).containsEdge( edge(0,2), idx )));
    PRINTVAR( idx );
    PRINTVAR((face(1,2,3).containsEdge( edge(3,1), idx )));
    PRINTVAR( idx );

     std::map<std::string, face> test {{"f 12 13 1", {12, 13, 1}},
                                         {"f 12/13 13/1 1/5", {12, 13, 1}},
                                         {"f 12/13/1 13/1/5 1/5/9", {12, 13, 1}},
                                         {"f 12//1 13//5 1//9", {12, 13, 1}},
                                         {"f 12/13 1 5", {}},
                                         {"f 12/12/12 13/32/32 1/2332/332", {12,13,1},},
                                         {"f 12//15 13//302 1//3200", {12,13,1},},
                                         {"not a face", {}}
    };

//    const std::regex color_regex(R"(f(\s+(\d+)(\/\d*){0,2})((\s+\d+)(\/\d*){0,2})(\s+(\d+)(\/\d*){0,2}))");
//
//    // simple match
//    for (const auto &line : lines) {
//        std::cout << line << ": " << std::boolalpha
//                  << std::regex_search(line, color_regex) << '\n';
//    }
//    std::cout << '\n';
//
//    // show contents of marked subexpressions within each match
//    std::smatch color_match;
//    for (const auto& line : lines) {
//        if(std::regex_search(line, color_match, color_regex)) {
//            std::cout << "matches for '" << line << "'\n";
//            std::cout << "Prefix: '" << color_match.prefix() << "'\n";
//            for (size_t i = 0; i < color_match.size(); ++i)
//                std::cout << i << ": " << color_match[i] << '\n';
//            std::cout << "Suffix: '" << color_match.suffix() << "\'\n\n";
//            std::cout << "res: '" << color_match[2] << "\t" << color_match[5] << "\t" << color_match[8] << "\'\n\n";
//        }
//    }

    for(const auto& line : test)
    {
        face f1{};
        const auto ok1 = parseFaceString(line.first, f1);

        const auto f2 = parseFaceStringRegex(line.first);

        if(ok1 != f2.has_value())
        {
            std::cerr << "ERROR: " << line.first << std::endl;
            continue;
        }

        if(ok1 && f1 != f2.value())
        {
            std::cerr << "ERROR: " << line.first << std::endl;
            continue;
        }

        if(ok1 && f1 != line.second)
        {
            std::cerr << "ERROR wrt gt: " << line.first << " " << f1 << " " << line.second <<  std::endl;
            continue;
        }
    }

    std::map<std::string, std::optional<point3d>> testVertices = {
         {"v 505.000 264.000 41.356", point3d(505.000f, 264.000f, 41.356f)},
         {"v -0.078125 0.242188 0.656250", point3d{-0.078125f, 0.242188f, 0.656250}},
         {"f -0.078125 0.242188 0.656250", std::nullopt},
         {"v -10.1603 5.71902 -0.957758", point3d(-10.1603f, 5.71902f, -0.957758f)},
         {"v 2.422296 -1.510915 -0.494169", point3d(2.422296f, -1.510915f, -0.494169f)},
         {"", std::nullopt,},
         {"not a vertex", std::nullopt}
    };

    for(const auto& vline : testVertices)
    {
        const auto p1 = parseVertexStringRegex(vline.first);
        if(p1.has_value() == vline.second.has_value())
        {
            std::cout << vline.first << std::endl;
            if(p1.has_value())
            {
                std::cout << p1.value() << std::endl;
                std::cout << vline.second.value() << std::endl;
            }
            else
            {
                std::cout << "no value" << std::endl;
            }
        }
        else
        {
            std::cerr << "ERROR: " << vline.first << std::endl;
            continue;
        }
    }
}
