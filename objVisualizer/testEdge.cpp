#include "core.hpp"
#include <sstream>
#include <stdio.h>

#if HAVE_STD_UNORDERED_MAP
#include <unordered_map>
#elif HAVE_TR1_UNORDERED_MAP
#include <tr1/unordered_map>
#elif HAVE_NO_UNORDERED_MAP
#include <map>
#endif

#if HAVE_TR1_UNORDERED_MAP || HAVE_STD_UNORDERED_MAP_IN_TR1_NAMESPACE
std::tr1::unordered_map< int, int > map;
#elif HAVE_STD_UNORDERED_MAP
std::unordered_map< int, int > map;
#else
std::map< int, int > map;
#endif


int main(int argc, char **argv)
{
    _edgeList mylist;


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
	tindex out;
	std::string b;
//	std::stringstream parser (line);
//	parser >> b
//		   >> out.v1 >> b >> a
//		   >> out.v2 >> b >> a
//		   >> out.v3 >> b >> a;
//	PRINTVAR(parser.fail());
//	std::cout << out << a << b << std::endl;
	PRINTVAR(line);
	PRINTVAR( (sscanf(line.c_str(), "f %u//%u %u//%u %u//%u", &(out.v1), &a, &(out.v2), &a, &(out.v3), &a) == 6) );
	std::cout << out << a << b << std::endl;
	
	PRINTVAR(v3f(1,2,4)*10.f);
	PRINTVAR(10.f*v3f(1,2,4));
}
