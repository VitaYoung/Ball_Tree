#include "BallTree.h"

#include <iostream>
#include <time.h>
#include <cmath>
#include <fstream> 
#include <stdlib.h>

#define PAGE_SIZE 65536
using namespace std;

bool BallTree::buildTree(int n, int d, float** data) {
    //  若根节点为空则无法建树 
    if (root == NULL) {
        
      root = new Node(0, n, NULL, 0, d);    
                 
    }
    MakeBallTree(root, n, d, data);
    return true;
}

/*
    建树的基本思路：
        （1）随机选一个点；
        （2）找到离这个点最远的点A；
        （3）找到离点A最远的点B；
        （4）对所有的点，如果它离A比较近就划入A的那个圈，离B比较近就划入B的那个圈；
        （5）如果圈内的点比设定的阀值大，那么继续重复上述过程；
*/
         
void BallTree::MakeBallTree(Node* node, int n, int d, float** data) {
    NodeNum++;
    node -> id = NodeNum - 1;
    //  判断该节点内数据的个数是否小于阀值（即某个圈内点的数量是否小于阀值），小于的话就直接存数据，半径和中心点也没必要算了 ，直接return 
    if (n < N0) {
        node -> data = new float*[n];
        for (int i = 0; i < n; i++) {
            node -> data[i] = new float[d + 1];
            for (int j = 0; j < d + 1; j++) {
                node -> data[i][j] = data[i][j];
            }
        }
        
        
        float* center = new float[d];
    
        for (int j = 1; j <= d; j++) {
            float sum = 0;
            for (int i = 0; i < n; i++) {
                sum += data[i][j];
            }
            center[j-1] = sum/(n + 0.0);
        }
        
        
        //  算半径，以离圆心最远的点到圆心的距离为半径 
        float radius = 0;
        
        for (int i = 0; i < n; i++) {
            float temp = 0;
            for (int j = 1; j <= d; j++) {
                temp += (center[j-1] - data[i][j]) * (center[j-1] - data[i][j]);
            }
            if (temp > radius) {
                radius = temp;
            }
        } 
        radius = sqrt(radius);
        
        node -> radius = radius;
        node -> center = center;
        node -> DataID = DataNum;
        DataNum++;
        return;
    } else {
        node->DataID = -1;
    }
    
    
    //  随机选一个点X
     
    srand((unsigned)time(NULL));   //  这个是用来改变种子的，不然的话种子一样随机出来的数会一直都是一样的 
    int random = rand() % n;
    
    
    float* x = data[random];      //  随机选一个点x， 注意二维数组的每行（一维数组）其实就是一组数据（一个点） 
    
    
    //   找出离点X最远的点A
     
    float distance = 0;     //  记录两点之间的最大距离的平方 
    int posA = 0;           //  A点数据在总数据中的位置 
    float temp = 0;         //  暂存某点与点X的距离的平方 
    
    for (int i = 0; i < n; i++) {
        temp = 0;
        for (int j = 1; j <= d; j++) {                            //  这个循环是算到点X的距离，可类比立体几何算距离的方法 
            temp += (x[j] - data[i][j]) * (x[j] - data[i][j]);               
        }
        if (temp > distance) {                            //  如果算出来的距离比存的最大距离的平方大，则更新最大距离的平方的值和A点数据在总数据中的位置 
            distance = temp;
            posA = i;
        }
    }
    
    
    float* A = data[posA];
    
    //  找出离A点最远的点B，以下算法同上 
    
    distance = 0;
    int posB = 0;
    
    for (int i = 0; i < n; i++) {
        temp = 0;
        for (int j = 1; j <= d; j++) {
            temp += (A[j] - data[i][j]) * (A[j] - data[i][j]);
        }
        if (temp > distance) {
            distance = temp;
            posB = i;
        }
    }
    
    
    float* B = data[posB];
    
    //  找出中心点 ，这里是用每个维度上分量的平均值作为圆心点各个分量的值 
    
    float* center = new float[d];
    
    for (int j = 1; j <= d; j++) {
        float sum = 0;
        for (int i = 0; i < n; i++) {
            sum += data[i][j];
        }
        center[j-1] = sum/(n + 0.0);
    }
    
    
    //  算半径，以离圆心最远的点到圆心的距离为半径 
    float radius = 0;
    
    for (int i = 0; i < n; i++) {
        temp = 0;
        for (int j = 1; j <= d; j++) {
            temp += (center[j-1] - data[i][j]) * (center[j-1] - data[i][j]);
        }
        if (temp > radius) {
            radius = temp;
        }
    } 
    radius = sqrt(radius);
    
    
    //  赋值，因为传进来时圆心和半径的值是空的 
    node -> center = center;
    node -> radius = radius;
    
    //  因为数据个数比阀值大，所以要分裂 
    MakeBallTreeSplit(A, B, node, n, d, data);
}

//  节点分裂的函数 
void BallTree::MakeBallTreeSplit(float* A, float* B, Node* node, int n, int d, float** data) {
    int LeftCount = 0;     // 左节点的数据个数 
    int RightCount = 0;      //  右节点的数据个数 
    float** DataLeft = new float*[n];        //  存左节点数据的二维数组 
    float** DataRight = new float*[n];       //  存右节点数据的二维数组 
    
    float distanceA = 0;        //  到点A的距离 
    float distanceB = 0;        //  到点B的距离 
    
    for (int i = 0; i < n; i++) {                  //  对每个点，分别计算它到点A和点B的距离 
        distanceA = 0;
        distanceB = 0;
        for (int j = 1; j <= d; j++) {
            distanceA += (A[j] - data[i][j]) * (A[j] - data[i][j]);
            distanceB += (B[j] - data[i][j]) * (B[j] - data[i][j]);
        }
        if (distanceA <= distanceB) {                //  离哪个圈近，就划分到哪个圈内 
            DataRight[RightCount] = data[i];
            RightCount++;
        } else {
            DataLeft[LeftCount] = data[i];
            LeftCount++;
        }
    }
    
    
    //  对当前圈进行划分 
    Node* left = new Node(0, LeftCount, NULL, 0, d);
    Node* right = new Node(0, RightCount, NULL, 0, d);
    
    
    //  赋值 
    node -> left = left;
    node -> right = right;
    
    //  继续向下建树 
    MakeBallTree(left, LeftCount, d, DataLeft);
    MakeBallTree(right, RightCount, d, DataRight);
    
    node -> LeftID = (node -> left) -> id;
    node -> RightID = (node -> right) -> id;
}

bool BallTree::storeTree(const char* index_path) {
    storeInit(index_path);
    int NodePageID = root->id / NodeslotNum;
    int NodeSlotID = root->id % NodeslotNum;
    char *node_path = new char[256];
    sprintf(node_path, "%s/NodePage%d.txt", index_path, NodePageID);
    ofstream out(node_path, ios::out | ios::binary);
    CurPageID =  NodePageID; 
    CurDataID =  -1;
    storeNode(root, out, index_path);
    out.close();
    out2.close();
}
     
void BallTree::storeInit(const char* index_path) {
     /*树节点页中每个槽的大小为 ： 树节点编号(4个字节)+左孩子节点编号(4个字节)+右孩子节点编号(4个字节)+该节点数据数量(4个字节)
                        +数据维度(4个字节)+数据组的编号(4个字节)+该节点半径+该节点圆心 */
     NodeslotSize = 4 + 4 + 4 + 4 + 4 + 4 + sizeof(float) + sizeof(float) * root->dim;
     NodeslotNum = PAGE_SIZE /  NodeslotSize;
     /*数据节点页中，每个槽的大小为 ：数据组的编号(4个字节) + 数据组的数据数量 + 数据组大小  ps:因为为定长记录所以数据组大小固定为20个数据*各每个数据大小，而读时只需读实际数据数量*/ 
     DataslotSize = 4 + 4 + sizeof(float) * (root->dim + 1) * 20;
     DataslotNum = PAGE_SIZE /  DataslotSize;

} 

void BallTree::storeNode(Node* node, ofstream &out, const char* index_path) {
    int NodePageID = node->id / NodeslotNum;
    int NodeSlotID = node->id % NodeslotNum;
    if (NodePageID != CurPageID) {
        CurPageID = NodePageID;
        out.close();
        char *node_path = new char[256];
        sprintf(node_path, "%s/NodePage%d.txt", index_path, NodePageID);
        out.open(node_path, ios::binary | ios::out);
    }
    //对该节点进行写入
    writeNode(node, out); 
    //如果节点的数据量小于20，则需要将其数据写入数据页 
    if (node->count < 20 ){
        storeData(node, index_path);
    }
    if (node->left != NULL) {
        storeNode(node->left, out, index_path);
    }
    if (node->right != NULL)
    storeNode(node->right, out, index_path);
} 

void BallTree::storeData(Node* node, const char* index_path) {
    int DataPageID = node->DataID / DataslotNum;
    if (CurDataID != DataPageID) {
        CurDataID = DataPageID;
        out2.close();
        char *data_path = new char[256];
        sprintf(data_path, "%s/DataPage%d.txt", index_path, DataPageID);
        out2.open(data_path, ios::binary | ios::out);
    }
    writeData(node, out2);
}

void BallTree::writeNode(Node* node, ofstream &out) {
    out.write((char*)&node->id, sizeof(int));
    out.write((char*)&node->LeftID, sizeof(int));
    out.write((char*)&node->RightID, sizeof(int));
    out.write((char*)&node->count, sizeof(int));
    out.write((char*)&node->dim, sizeof(int));
    out.write((char*)&node->DataID, sizeof(int));
    out.write((char*)&node->radius, sizeof(float));
    
    for (int i = 0; i < node->dim; i++) {
        out.write((char*)&node->center[i], sizeof(float));
    }
    
}

void BallTree::writeData(Node* node, ofstream &out) {
    out.write((char*)&node->DataID, sizeof(int));
    out.write((char*)&node->count, sizeof(int));
    for (int i = 0; i < node->count; i++) {
        out.write((char*)node->data[i], sizeof(float)*(node->dim+1));
    }
    float* temp = 0;

    for (int i = node->count; i < 20; i++) {
        for (int j = 0; j < node->dim+1; j++) {
            out.write((char*)&temp, sizeof(float));
        }
    }
    
   
    delete []temp;
}

bool BallTree::restoreTree(const char* index_path) {
    char *node_path = new char[256];
    sprintf(node_path, "%s/NodePage%d.txt", index_path, 0);
    ifstream in(node_path, ios::binary | ios::in);
    CurPageID = 0;

    int id, LeftID, RightID, count, dim, DataID;
    float radius;
    
    in.read((char*)&id, sizeof(int));
    in.read((char*)&LeftID, sizeof(int));
    in.read((char*)&RightID, sizeof(int));
    in.read((char*)&count, sizeof(int));
    in.read((char*)&dim, sizeof(int));
    in.read((char*)&DataID, sizeof(int));
    in.read((char*)&radius, sizeof(float));

    float* center = new float[dim];
    in.read((char*)center, sizeof(float)*dim);

    root = new Node(id, LeftID, RightID, count, dim, DataID, radius, center);
    restoreInit();
    in.close();
}

void BallTree::restoreInit() {
     /*树节点页中每个槽的大小为 ： 树节点编号(4个字节)+左孩子节点编号(4个字节)+右孩子节点编号(4个字节)+该节点数据数量(4个字节)
                        +数据维度(4个字节)+数据组的编号(4个字节)+该节点半径+该节点圆心 */
     NodeslotSize = 4 + 4 + 4 + 4 + 4 + 4 + sizeof(float) + sizeof(float) * root->dim;
     NodeslotNum = PAGE_SIZE /  NodeslotSize;
     /*数据节点页中，每个槽的大小为 ：数据组的编号(4个字节) + 数据组的数据数量 + 数据组大小  ps:因为为定长记录所以数据组大小固定为20个数据*各每个数据大小，而读时只需读实际数据数量*/ 
     DataslotSize = 4 + 4 + sizeof(float) * (root->dim + 1) * 20;
     DataslotNum = PAGE_SIZE /  DataslotSize;
} 

void BallTree::readNode(Node* node, const char* index_path) {
    // 先构造node的左节点
    int NodePageID = node->LeftID / NodeslotNum;
    int NodeSlotID = node->LeftID % NodeslotNum;
    CurPageID = NodePageID;
    char *node_path = new char[256];
    sprintf(node_path, "%s/NodePage%d.txt", index_path, NodePageID);
    ifstream in(node_path, ios::binary | ios::in);

    // 寻找所在的槽做偏移
    
    char* buffer = new char[NodeslotSize];
    for (int i = 0; i < NodeSlotID; i++) {
        in.read(buffer, NodeslotSize);
    }

    // 读取数据构造节点
    int id, LeftID, RightID, count, dim, DataID;
    float radius;
    in.read((char*)&id, sizeof(int));
    in.read((char*)&LeftID, sizeof(int));
    in.read((char*)&RightID, sizeof(int));
    in.read((char*)&count, sizeof(int));
    in.read((char*)&dim, sizeof(int));
    in.read((char*)&DataID, sizeof(int));
    in.read((char*)&radius, sizeof(float));

    float* center = new float[dim];
    for (int i = 0; i < dim; i++) {
        in.read((char*)&center[i], sizeof(float));
    }
    node->left = new Node(id, LeftID, RightID, count, dim, DataID, radius, center);

    //再构造node的右节点
    NodePageID = node->RightID / NodeslotNum;
    NodeSlotID = node->RightID % NodeslotNum;
    if (1) {
        CurPageID = NodePageID;
        in.close();
        char *node_path = new char[256];
        sprintf(node_path, "%s/NodePage%d.txt", index_path, NodePageID);
        in.open(node_path, ios::binary | ios::in);
    }
    
    // 寻找所在的槽做偏移
    buffer = new char[NodeslotSize];
    for (int i = 0; i < NodeSlotID; i++) {
        in.read(buffer, NodeslotSize);
    }

    // 读取数据构造节点
    in.read((char*)&id, sizeof(int));
    in.read((char*)&LeftID, sizeof(int));
    in.read((char*)&RightID, sizeof(int));
    in.read((char*)&count, sizeof(int));
    in.read((char*)&dim, sizeof(int));
    in.read((char*)&DataID, sizeof(int));
    in.read((char*)&radius, sizeof(float));

    center = new float[dim];
    for (int i = 0; i < dim; i++) {
        in.read((char*)&center[i], sizeof(float));
    }
    
    
    
    node->right = new Node(id, LeftID, RightID, count, dim, DataID, radius, center);
    in.close();
}

void BallTree::readData(Node* node, const char* index_path) {
    int DataPageID = node->DataID / DataslotNum;
    int DataSlotID = node->DataID % DataslotNum;
    char *data_path = new char[256];
    sprintf(data_path, "%s/DataPage%d.txt", index_path, DataPageID);
    ifstream in(data_path, ios::binary | ios::in);

    char* buffer = new char[2*DataslotSize];
    for (int i = 0; i < DataSlotID; i++) {
        in.read(buffer, DataslotSize);
    }

    int DataID, count;
    in.read((char*)&DataID, sizeof(int));
    in.read((char*)&count, sizeof(int));

    node->data = new float*[node->count];
    float temp;
    for (int i = 0; i < node->count; i++) {
        node->data[i] = new float[node->dim+1];
        
        
        for (int j = 0; j < node->dim+1; j++) {
            in.read((char*)&node->data[i][j], sizeof(float));
            temp = node->data[i][j];
        }
    }
    

    in.close();
}

float BallTree::innerPro(Query* q, float* p, int dim) {
	int n = dim;
	float sum = 0.0;
	for (int i = 0; i < n; i++) {
		sum += ((q -> test[i+1])*p[i]);
	}
	return sum;
}

float BallTree::innerPro2(Query* q, float* p, int dim) {
	int n = dim;
	float sum = 0.0;
	for (int i = 0; i < n; i++) {
		sum += ((q -> test[i+1])*p[i+1]);
	
	}
	return sum;
}

float BallTree::MIP(Query* q, Node* T) {
	float* p = T -> center;
	float mode = 0.0;
	float sum = 0.0;
	float inner = 0.0;
	for (int i = 0; i < T -> dim; i++) {
		sum += ((q -> test[i+1])*(q ->test[i+1]));
	}
	mode = sqrt(sum);
	inner = innerPro(q, p, T -> dim);
	float temp = (inner+((T ->radius)*mode));
	return temp;
}

bool BallTree::isLeaf(Node* T) {
	if (T->count < 20) {           //当一个节点中的数据数小于20，则是叶子节点。 
		return true;
	}
	return false;
}

 

void BallTree::LinearSearch(Query* q, Node* T) {
	for (int i = 0 ; i < T -> count ; i++) {
		if (innerPro2(q, T -> data[i], T -> dim) > q -> inner) {
			q -> bm = T -> data[i][0];
			q -> inner = innerPro2(q, T -> data[i], T->dim);
		}
	}
}


void BallTree::TreeSearch(Query* q, Node* T) {
	if (q -> inner < MIP(q, T)) {
		if (isLeaf(T)) {
            readData(T, "Netflix/index");
			LinearSearch(q, T);
		} else {
            readNode(T, "Netflix/index");
			float I1 = MIP(q, T -> left);
			float I2 = MIP(q, T -> right);
			if (I1 <= I2) {
				TreeSearch(q, T -> right);
				TreeSearch(q, T -> left);
			} else {
				TreeSearch(q, T -> left);
				TreeSearch(q, T -> right);
			}
		}
	}
}


int BallTree::mipSearch(int d,float* query) {
	Query* q = new Query(d, query);
	TreeSearch(q, root);
	return q -> bm;
}

