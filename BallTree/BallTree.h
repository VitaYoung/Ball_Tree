#ifndef __BALL_TREE_H
#define __BALL_TREE_H

#include <iostream>
#include <time.h>
#include <cmath>
#include <fstream> 
#include <stdlib.h>
using namespace std;

#define N0 20

/*
    关于Ball Tree数据结构的说明
        首先，主体的数据为一个二维数组，我们知道，二维数组其实是由一维数组组成的，那么这个二维数组可以想象成一个表
            即每行（一个一维数组）代表一组数据，也就是一个点，行数即为点的数量 
            每列代表该组数据分量的具体值，列数即为点的维度，应除去第一列，因为第一列为主键id，不能作为计算距离的依据 
          
        举个例子，现在有一张表如下：
        姓名   年龄   电话   工号 
        张三   23     123     1
        李四   20     456     2
        王五   18     789     3
        
        那么用3，4, 5分别代表他们的姓名，则这张表可以抽象为如下的二维数组
        3   23   123   1
        4   20   456   2
        5   18   789   3
        
        其中每行表示一个点，每列的值相当于该点在不同维度的分量
        
        其实可以看成3个一维数组：   （维度为3） 
        num[1] = {3, 23, 123, 1}  对应点(23, 123，1) 
        num[2] = {4, 20, 456，2}  对应点(20, 456，2) 
        num[3] = {5, 18, 789，3}  对应点(18, 789，3) 
        
        这样就实现了把一组数据抽象为一个点，然后就能进行距离的判断，进而进行后面最大内积的计算与求解
*/ 


//  节点定义 ，这里我没有分叶子节点和非叶子节点了，虽然那样更合理，但实现起来就有点麻烦，偷个懒orz 
struct Node{
    int id;     //  节点的ID,任务二可能要使用 
    Node* left;     //  左节点 
    Node* right;    //  右节点 
    float** data;   //  如果点的数量小于20才用来存数据，即相当于这个"圈"里面包含的所有的点的信息 
    int count;      //  即这个"圈"里面总共有多少个点，如果小于20，则这个点为叶子节点，不再分裂并进行数据写入，否则要分裂
    float* center;    //  圆心
    float radius;     //  半径 
    int dim;        //   维度 
    int DataID;     //   数据的ID，写入数据页时用 
    int LeftID;     //  左节点的ID
    int RightID;    //  右节点的ID 
    
    Node(int i, int c, float* cen, float r, int d) {
        id = i;
        count = c;
        center = cen;
        radius = r;
        dim = d;
        left = NULL;
        right = NULL;
        data = NULL;
        DataID = 0;
        LeftID = 0;
        RightID = 0;
    }

    Node(int _id, int _LeftID, int _RightID, int _count, int _dim, int _DataID, int _radius, float* _center) {
        id = _id;
        LeftID = _LeftID;
        RightID = _RightID;
        count = _count;
        dim = _dim;
        DataID = _DataID;
        radius = _radius;
        center = _center;
        data = NULL;
        left = NULL;
        right = NULL;
    }
};

struct Query{
	int d;
	float* test;
	int bm;
	float inner;
	Query(int t, float* a) {
		d = t;
		test = new float[t+1];
		test = a;
		bm = 0;
		inner = -9999;
	}
};




class BallTree {
public:
    Node* root;
    int NodeNum;
    int DataNum;
    int CurPageID;
    int NodeslotNum;
    int DataslotNum;
    int NodeslotSize;
    int DataslotSize;
    ofstream out2;
    int CurDataID;
	BallTree() {
	    root = NULL;
	    NodeNum = 0;
	    DataNum = 0;
    }
    
    BallTree(Node* r) {
        root = r;
        NodeNum = 0;
	    DataNum = 0;
    }
    
	~BallTree() {
	    
	};

	bool buildTree(
		int n,
		int d,
		float** data);

    void MakeBallTree(Node* node, int n, int d, float** data);     // 对应文献相应函数的伪代码 
     
    void MakeBallTreeSplit(float* A, float* B, Node* node, int n, int d, float** data);   //  对应文献相应函数的伪代码 
    
    void storeInit(const char* index_path);
    void storeNode(Node* node, ofstream &out, const char* index_path);
    void storeData(Node* node, const char* index_path);
    void writeNode(Node* node, ofstream &out);
    void writeData(Node* node, ofstream &out);
	bool storeTree(
		const char* index_path);

	bool restoreTree(
		const char* index_path);
    void readNode(Node* node, const char* index_path);
    void readData(Node* node, const char* index_path);
    void restoreInit();
	
	bool isLeaf(Node* T);//判断当前节点是否是叶节点 
	
	float MIP(Query* q, Node* T);//当前节点可能最大內积 
	
	float innerPro(Query* q, float* p, int dim);//计算內积 
	
	float innerPro2(Query* q, float* p, int dim);
	
	void TreeSearch(Query* q, Node* T);
	
	
	int mipSearch(int d,float* query);  //d为维度 


	void LinearSearch(Query* q, Node* T);//当节点为叶节点的时候进行遍历 



	// optional
	bool insertData(
		int d,
		float* data);

	// optional
	bool deleteData(
		int d,
		float* data);

	// optional
	bool buildQuadTree(
		int n,
		int d,
		float** data);
};

#endif
