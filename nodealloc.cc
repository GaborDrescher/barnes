#include "nodealloc.h"

//static std::vector<node*> buffer;

void destroyNode(node *n)
{
    //buffer.push_back(n);
    delete n;
}

node* getNewNode()
{
    return new node;
    //node *no = 0;
    //
    //if(buffer.size() == 0) {
    //    no = new node;
    //}
    //else {
    //    no = buffer.back();
    //    buffer.pop_back();
    //}

    //no->reset();
    //return no;
}
