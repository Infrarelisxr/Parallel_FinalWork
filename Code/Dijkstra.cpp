#include <iostream>
#include <cstdio>
#include <random>
#include <vector>
#include <fstream>
using namespace std;
const int MAX = 9999999;
const int MAXN = 20001;
int dist[MAXN];//存储最短路径 
bool s[MAXN];//存储已经找到最短路径的结点 
int u[MAXN][MAXN];//存储邻接矩阵 
int front[MAXN];//记录每个点最短路径上的前一个点 
int n;

void dijkstra() {
	for (int i = 1; i <= n; i++) {//预置最短路径为最大 
		dist[i] = MAX;
	}
	for (int i = 1; i <= n; i++) {
		s[i] = false;
	}
	s[1] = true;  //起始点放入 s 数组
	front[1] = 1;
	dist[1] = 0;
	int temp = 1;//中转点 
	int num = 1;
	int tempmin = MAX;
	while (num < n) {
		for (int i = 1; i <= n; i++) {
			if (u[temp][i] != MAX && !s[i]) {//中转点到某个点 i 有路径，并且点i还没有进入 s数组 
				//dist[i]=min(dist[temp]+u[temp][i],dist[i]);
				if (dist[temp] + u[temp][i] < dist[i]) {
					dist[i] = dist[temp] + u[temp][i];
					front[i] = temp;//记录一下前一个元素 
				}
			}
		}
		for (int i = 1; i <= n; i++) {
			if (!s[i] && dist[i] < tempmin) {
				temp = i;//设为中转点 
				tempmin = dist[i];
			} //找到当前最短距离的点 

		}
		s[temp] = true;//放入 s数组 
		num++;
		tempmin = MAX;
	}
}

void traceback() {//路径查找 
	int stack[MAXN];
	int top = 0;
	for (int i = 2; i <= n; i++) {
		int j = i;
		stack[top++] = i;
		while (front[j] != j) {
			stack[top++] = front[j];
			j = front[j];
		}
		while (top > 0) {
			printf("%d ", stack[--top]);
		}
		printf("\n");
	}
}
//void data(int n) {
//	random_device rd;
//	mt19937 gen(rd());
//	uniform_int_distribution<> distrib(1, 100);
//	vector<std::vector<int>> matrix(n, std::vector<int>(n));
//	// 填充矩阵
//	for (int i = 0; i < n; ++i) {
//		for (int j = 0; j < n; ++j) {
//			if (i == j)
//			{
//				matrix[i][j] = 0;
//			}
//			else
//			{
//				matrix[i][j] = distrib(gen);
//			}
//
//		}
//	}
//	ofstream file("matrix.txt");
//	// 检查文件是否成功打开
//	if (!file.is_open()) {
//		std::cerr << "Failed to open the file for writing." << std::endl;
//		return;
//	}
//
//	// 写入矩阵到文件
//	for (const auto& row : matrix) {
//		for (int val : row) {
//			file << val << " ";
//		}
//		file << std::endl; // 每行结束后换行
//	}
//
//	// 关闭文件
//	file.close();
//
//	std::cout << "Matrix has been written to matrix.txt" << std::endl;
//}
int main() {
	cin >> n;
	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= n; j++) {
			if (i == j)
			{
				u[i][j] = 0;
			}
			else
			{
				u[i][j] = 1;
			}

			
		}
	}
	struct timespec sts, ets;
	timespec_get(&sts, TIME_UTC);
	for (int i = 0; i < 100; i++)
	{
		dijkstra();
	}
	timespec_get(&ets, TIME_UTC);
	time_t dsec = ets.tv_sec - sts.tv_sec;
	long dnsec = ets.tv_nsec - sts.tv_nsec;
	if (dnsec < 0)
	{
		dsec--;
		dnsec += 1000000000ll;
	}
	printf("%lld.%09lld\n", dsec, dnsec);

	cout << "SUCEESS!" << endl;

	return 0;
}