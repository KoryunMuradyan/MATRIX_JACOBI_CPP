#ifndef	__JACOBI__HPP__
#define	__JACOBI__HPP__


#include <bits/stdc++.h>
#include <iterator>
#include "matrix.hpp"
#include "Read_from_write_into_file.hpp"
using namespace math;

template <typename T> bool 
is_Jacobi_mathod_appliable(const matrix_type<T>& raw_matrix)
{
	if(raw_matrix[0].size() - raw_matrix.size() != 1) {
		return false;
	}
	std::vector<T> row;
	for(int i = 0; i < raw_matrix.size(); ++i) {
		T x_i = raw_matrix[i][i];
		row = raw_matrix[i];
		row[i] = 0;
		row.back() = 0;
		T sum = std::transform_reduce(row.begin(), row.end(),0,
				std::plus<>{}, static_cast<T (*)(T)>(std::fabs)
			);
		if(sum > abs(x_i)) {
			return false;
		}
	}
	return true;
}

template <typename T> 
matrix_type<T> get_Jacobi_matrix(matrix_type<T> raw_matrix)
{
	for(int i = 0; i < raw_matrix.size(); ++i) {
		int tmp = raw_matrix[i][i];
		raw_matrix[i][i] = -raw_matrix[i].back();
		raw_matrix[i].back() = -tmp;
		T myconstant = -tmp;
		tmp = raw_matrix[i][i];
		raw_matrix[i].erase(raw_matrix[i].begin() + i);
		raw_matrix[i].emplace(raw_matrix[i].begin(), tmp);
		reverse(raw_matrix[i].begin() + 1, raw_matrix[i].end() - 1);
		raw_matrix[i].pop_back();
		reverse(raw_matrix[i].begin(), raw_matrix[i].end());
		std::transform(raw_matrix[i].begin(), raw_matrix[i].end(), 
				raw_matrix[i].begin(), [&](auto& c){
						return c/myconstant;
					}
		);
	}
	return raw_matrix;
}

template <typename T> std::map<int, T> 
get__X_N_T_next(matrix_type<T> raw_matrix, typename std::map<int, T>& X_N_t_1)
{
	std::map<int, T> tmp_map = X_N_t_1;
	std::map<int, T> X_N_t_2 = X_N_t_1;
	typename std::map<int, T>::iterator it;
	typename std::map<int, T>::iterator it_origin;
	for(int i = 0; i < raw_matrix.size(); ++i) {
		std::map<int, T> tmp_map_i = tmp_map;
		tmp_map_i.erase(i);
		it = tmp_map_i.begin();
		it_origin = X_N_t_2.find(i);
		it_origin->second = raw_matrix[i].back();
		for(int j = 0; j < raw_matrix[i].size() - 1; ++j) {
			it_origin->second += (it->second)*raw_matrix[i][j];
			++it;
		}
	}
	return X_N_t_2;
}

template <typename T> bool 
should_iteration_go_on(std::map<int, T> X_t_1, std::map<int, T> X_t_2, T aprox)
{ 
	int sz = X_t_1.size();
	typename std::map<int, T>::iterator it_1 = X_t_1.begin();
	typename std::map<int, T>::iterator it_2 = X_t_2.begin();
	T biggest = abs(abs(it_2->second) - abs(it_1->second));
	for(int i = 0; i < sz; ++i) {
		if(biggest < abs(abs(it_2->second) - abs(it_1->second))) {
			biggest = abs(abs(it_2->second) - abs(it_1->second));
		}
		it_1++;
		it_2++;
	}
	if(biggest < aprox){
		return false;
	}
	return true;
}

template <typename T> std::map<int, T> 
iterate_Jacobi(matrix_type<T>& B, std::map<int, T>& X_N_t, T& aprox)
{
	std::map<int, T> X_N_t_NEXT = get__X_N_T_next(B, X_N_t);
	while(should_iteration_go_on(X_N_t, X_N_t_NEXT, aprox)) {
		std::map<int, T> tmp_map = X_N_t_NEXT;
		X_N_t_NEXT = get__X_N_T_next(B, X_N_t_NEXT);
		X_N_t = tmp_map;
	}
	return X_N_t_NEXT;
}

template <typename T> std::map<int, T> 
solve_Jacobi(Matrix<T> matrix, T aprox = 0.0001)
{
	matrix_type<T> raw_matrix = matrix.get_matrix();
	if(!is_Jacobi_mathod_appliable(raw_matrix)) {
		std::cout << "Not possible to apply Jacobi method to given matrix !\n";
		exit(0);
	}
	std::vector<std::vector<T>> B = get_Jacobi_matrix(raw_matrix);
	std::map<int, T> X_N_t;
	int n = 0;
	for(auto row: B){    
		X_N_t.insert({n++, 0});
	}
	std::map<int, T> result_X_n = iterate_Jacobi(B, X_N_t, aprox);
	return result_X_n;
}

#endif // __JACOBI__HPP__
