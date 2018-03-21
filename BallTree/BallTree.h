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
    ����Ball Tree���ݽṹ��˵��
        ���ȣ����������Ϊһ����ά���飬����֪������ά������ʵ����һά������ɵģ���ô�����ά������������һ����
            ��ÿ�У�һ��һά���飩����һ�����ݣ�Ҳ����һ���㣬������Ϊ������� 
            ÿ�д���������ݷ����ľ���ֵ��������Ϊ���ά�ȣ�Ӧ��ȥ��һ�У���Ϊ��һ��Ϊ����id��������Ϊ������������ 
          
        �ٸ����ӣ�������һ�ű����£�
        ����   ����   �绰   ���� 
        ����   23     123     1
        ����   20     456     2
        ����   18     789     3
        
        ��ô��3��4, 5�ֱ�������ǵ������������ű���Գ���Ϊ���µĶ�ά����
        3   23   123   1
        4   20   456   2
        5   18   789   3
        
        ����ÿ�б�ʾһ���㣬ÿ�е�ֵ�൱�ڸõ��ڲ�ͬά�ȵķ���
        
        ��ʵ���Կ���3��һά���飺   ��ά��Ϊ3�� 
        num[1] = {3, 23, 123, 1}  ��Ӧ��(23, 123��1) 
        num[2] = {4, 20, 456��2}  ��Ӧ��(20, 456��2) 
        num[3] = {5, 18, 789��3}  ��Ӧ��(18, 789��3) 
        
        ������ʵ���˰�һ�����ݳ���Ϊһ���㣬Ȼ����ܽ��о�����жϣ��������к�������ڻ��ļ��������
*/ 


//  �ڵ㶨�� ��������û�з�Ҷ�ӽڵ�ͷ�Ҷ�ӽڵ��ˣ���Ȼ������������ʵ���������е��鷳��͵����orz 
struct Node{
    int id;     //  �ڵ��ID,���������Ҫʹ�� 
    Node* left;     //  ��ڵ� 
    Node* right;    //  �ҽڵ� 
    float** data;   //  ����������С��20�����������ݣ����൱�����"Ȧ"������������еĵ����Ϣ 
    int count;      //  �����"Ȧ"�����ܹ��ж��ٸ��㣬���С��20���������ΪҶ�ӽڵ㣬���ٷ��Ѳ���������д�룬����Ҫ����
    float* center;    //  Բ��
    float radius;     //  �뾶 
    int dim;        //   ά�� 
    int DataID;     //   ���ݵ�ID��д������ҳʱ�� 
    int LeftID;     //  ��ڵ��ID
    int RightID;    //  �ҽڵ��ID 
    
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

    void MakeBallTree(Node* node, int n, int d, float** data);     // ��Ӧ������Ӧ������α���� 
     
    void MakeBallTreeSplit(float* A, float* B, Node* node, int n, int d, float** data);   //  ��Ӧ������Ӧ������α���� 
    
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
	
	bool isLeaf(Node* T);//�жϵ�ǰ�ڵ��Ƿ���Ҷ�ڵ� 
	
	float MIP(Query* q, Node* T);//��ǰ�ڵ�������Ȼ� 
	
	float innerPro(Query* q, float* p, int dim);//����Ȼ� 
	
	float innerPro2(Query* q, float* p, int dim);
	
	void TreeSearch(Query* q, Node* T);
	
	
	int mipSearch(int d,float* query);  //dΪά�� 


	void LinearSearch(Query* q, Node* T);//���ڵ�ΪҶ�ڵ��ʱ����б��� 



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
