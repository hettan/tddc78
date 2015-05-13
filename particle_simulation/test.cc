#include <iostream>
#include <cstdio>
#include <sstream>

using namespace std;


int to_int(const string str){
  int i;
  istringstream (str) >> i;
  return i;
}

string to_string(const int a){
  ostringstream oss;
  oss << a;
  return oss.str();
}

int main(){
  
  int a = 5;
  int b = 10;
  cout << to_string(a) << to_string(b) << endl;

  const string c = to_string(a) + to_string(b);
  cout << to_int(c) << endl;
  
  return 0;
}
