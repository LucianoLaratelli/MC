#include <vector>
#include <iostream>

typedef struct test_struct_t
{
    double a;
    int b;
    char c;
} test_struct;

using namespace std;

int main()
{
    vector<test_struct> test_Vector; //declares a variable of type variable whose name is test_Vector
    //test_Vector is an OBJECT, it has object methods associated with it like push_back

    test_struct apple;// apple is a struct with three pieces of data associated with it 
    apple[1].a=3;
    apple.a = 1.0;// the value a inside the struct apple is = 1 
    apple.b = 3;
    apple.c = 'v';

    test_Vector.push_back(apple);

    return 0;
}


