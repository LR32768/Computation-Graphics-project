#include "maths.h"
#include <iostream>
using namespace std;

void vprint(Vec v)
{
	cout << '(' << v.x << ',' << v.y << ',' << v.z << ")" << endl;
}

int main()
{
	Vec a(1, 2, 3), b(-2, 1, 4);
	vprint(a + b);
	vprint(a - b);
	vprint(a * b);
	cout << "dot: " << a.dot(b) << endl;
	vprint(a % b);
	for (int i = 0; i < 100; i++)
		cout << rand01() << endl;
	return 0;
}
