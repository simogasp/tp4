#include "core.hpp"


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
}
