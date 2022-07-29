#ifndef READ_FROM_FILE
#define READ_FROM_FILE

#include <map>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iterator>     // back_inserter
#include <vector>


using namespace std;

template <typename T>
void NumLinetoIntVec(std::vector<std::vector<T>>& arg_vec, 
		     std::string& arg_str) 
{
	std::istringstream iss(arg_str);
	std::vector<T> v{std::istream_iterator<T>(iss),
	std::istream_iterator<T>()};
	arg_vec.push_back(v);
}

std::vector<std::string> str_seq_to_int_vec(const std::string& arg_str) 
{
	std::vector<std::string> strings_vec;
	std::istringstream f(arg_str);
	std::string s;    
	while (getline(f, s, '\n')) {
		strings_vec.push_back(s);
	}
	return strings_vec;
}

template <typename T>
std::vector<std::vector<T>> MatrixRead(const std::string&& filename) 
{
	std::ifstream file_to_read(filename);
	const std::string str_file((
			std::istreambuf_iterator<char>(file_to_read)),
			std::istreambuf_iterator<char>());
	file_to_read.close();
	std::vector<std::string> num_string_vec = str_seq_to_int_vec(str_file);
	std::vector<std::vector<T>> vec_2d_int;
	std::for_each (num_string_vec.begin(), num_string_vec.end(), 
		       [&vec_2d_int](std::string str_to_pass) {
		NumLinetoIntVec(vec_2d_int, str_to_pass);}
	);
	for(auto it = vec_2d_int.begin(); it != vec_2d_int.end() - 1; ++it) {
		if(it->size() != (it + 1)->size()) {
			std::cout << "the input data from file is not matrix type !\n";
			exit(0);
		}
	}
	return vec_2d_int;
}

template <typename T>
void Generate_Output_File(const std::map<int, T>& X_N) 
{
	std::ofstream file_to_write("output.txt");
	std::string str_to_write = "";
	for(auto elem: X_N) {
		str_to_write = str_to_write + "  " + std::to_string(elem.second);
	}
	str_to_write += "\n";
	file_to_write << str_to_write;
}

#endif // READ_FROM_FILE
