#include "BallTree.h"

#include <iostream>
#include <time.h>
#include <cmath>
#include <fstream> 
#include <stdlib.h>

#define PAGE_SIZE 65536
using namespace std;

bool BallTree::buildTree(int n, int d, float** data) {
    //  �����ڵ�Ϊ�����޷����� 
    if (root == NULL) {
        
      root = new Node(0, n, NULL, 0, d);    
                 
    }
    MakeBallTree(root, n, d, data);
    return true;
}

/*
    �����Ļ���˼·��
        ��1�����ѡһ���㣻
        ��2���ҵ����������Զ�ĵ�A��
        ��3���ҵ����A��Զ�ĵ�B��
        ��4�������еĵ㣬�������A�ȽϽ��ͻ���A���Ǹ�Ȧ����B�ȽϽ��ͻ���B���Ǹ�Ȧ��
        ��5�����Ȧ�ڵĵ���趨�ķ�ֵ����ô�����ظ��������̣�
*/
         
void BallTree::MakeBallTree(Node* node, int n, int d, float** data) {
    NodeNum++;
    node -> id = NodeNum - 1;
    //  �жϸýڵ������ݵĸ����Ƿ�С�ڷ�ֵ����ĳ��Ȧ�ڵ�������Ƿ�С�ڷ�ֵ����С�ڵĻ���ֱ�Ӵ����ݣ��뾶�����ĵ�Ҳû��Ҫ���� ��ֱ��return 
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
        
        
        //  ��뾶������Բ����Զ�ĵ㵽Բ�ĵľ���Ϊ�뾶 
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
    
    
    //  ���ѡһ����X
     
    srand((unsigned)time(NULL));   //  ����������ı����ӵģ���Ȼ�Ļ�����һ���������������һֱ����һ���� 
    int random = rand() % n;
    
    
    float* x = data[random];      //  ���ѡһ����x�� ע���ά�����ÿ�У�һά���飩��ʵ����һ�����ݣ�һ���㣩 
    
    
    //   �ҳ����X��Զ�ĵ�A
     
    float distance = 0;     //  ��¼����֮����������ƽ�� 
    int posA = 0;           //  A���������������е�λ�� 
    float temp = 0;         //  �ݴ�ĳ�����X�ľ����ƽ�� 
    
    for (int i = 0; i < n; i++) {
        temp = 0;
        for (int j = 1; j <= d; j++) {                            //  ���ѭ�����㵽��X�ľ��룬��������弸�������ķ��� 
            temp += (x[j] - data[i][j]) * (x[j] - data[i][j]);               
        }
        if (temp > distance) {                            //  ���������ľ���ȴ���������ƽ����������������ƽ����ֵ��A���������������е�λ�� 
            distance = temp;
            posA = i;
        }
    }
    
    
    float* A = data[posA];
    
    //  �ҳ���A����Զ�ĵ�B�������㷨ͬ�� 
    
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
    
    //  �ҳ����ĵ� ����������ÿ��ά���Ϸ�����ƽ��ֵ��ΪԲ�ĵ����������ֵ 
    
    float* center = new float[d];
    
    for (int j = 1; j <= d; j++) {
        float sum = 0;
        for (int i = 0; i < n; i++) {
            sum += data[i][j];
        }
        center[j-1] = sum/(n + 0.0);
    }
    
    
    //  ��뾶������Բ����Զ�ĵ㵽Բ�ĵľ���Ϊ�뾶 
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
    
    
    //  ��ֵ����Ϊ������ʱԲ�ĺͰ뾶��ֵ�ǿյ� 
    node -> center = center;
    node -> radius = radius;
    
    //  ��Ϊ���ݸ����ȷ�ֵ������Ҫ���� 
    MakeBallTreeSplit(A, B, node, n, d, data);
}

//  �ڵ���ѵĺ��� 
void BallTree::MakeBallTreeSplit(float* A, float* B, Node* node, int n, int d, float** data) {
    int LeftCount = 0;     // ��ڵ�����ݸ��� 
    int RightCount = 0;      //  �ҽڵ�����ݸ��� 
    float** DataLeft = new float*[n];        //  ����ڵ����ݵĶ�ά���� 
    float** DataRight = new float*[n];       //  ���ҽڵ����ݵĶ�ά���� 
    
    float distanceA = 0;        //  ����A�ľ��� 
    float distanceB = 0;        //  ����B�ľ��� 
    
    for (int i = 0; i < n; i++) {                  //  ��ÿ���㣬�ֱ����������A�͵�B�ľ��� 
        distanceA = 0;
        distanceB = 0;
        for (int j = 1; j <= d; j++) {
            distanceA += (A[j] - data[i][j]) * (A[j] - data[i][j]);
            distanceB += (B[j] - data[i][j]) * (B[j] - data[i][j]);
        }
        if (distanceA <= distanceB) {                //  ���ĸ�Ȧ�����ͻ��ֵ��ĸ�Ȧ�� 
            DataRight[RightCount] = data[i];
            RightCount++;
        } else {
            DataLeft[LeftCount] = data[i];
            LeftCount++;
        }
    }
    
    
    //  �Ե�ǰȦ���л��� 
    Node* left = new Node(0, LeftCount, NULL, 0, d);
    Node* right = new Node(0, RightCount, NULL, 0, d);
    
    
    //  ��ֵ 
    node -> left = left;
    node -> right = right;
    
    //  �������½��� 
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
     /*���ڵ�ҳ��ÿ���۵Ĵ�СΪ �� ���ڵ���(4���ֽ�)+���ӽڵ���(4���ֽ�)+�Һ��ӽڵ���(4���ֽ�)+�ýڵ���������(4���ֽ�)
                        +����ά��(4���ֽ�)+������ı��(4���ֽ�)+�ýڵ�뾶+�ýڵ�Բ�� */
     NodeslotSize = 4 + 4 + 4 + 4 + 4 + 4 + sizeof(float) + sizeof(float) * root->dim;
     NodeslotNum = PAGE_SIZE /  NodeslotSize;
     /*���ݽڵ�ҳ�У�ÿ���۵Ĵ�СΪ ��������ı��(4���ֽ�) + ��������������� + �������С  ps:��ΪΪ������¼�����������С�̶�Ϊ20������*��ÿ�����ݴ�С������ʱֻ���ʵ����������*/ 
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
    //�Ըýڵ����д��
    writeNode(node, out); 
    //����ڵ��������С��20������Ҫ��������д������ҳ 
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
     /*���ڵ�ҳ��ÿ���۵Ĵ�СΪ �� ���ڵ���(4���ֽ�)+���ӽڵ���(4���ֽ�)+�Һ��ӽڵ���(4���ֽ�)+�ýڵ���������(4���ֽ�)
                        +����ά��(4���ֽ�)+������ı��(4���ֽ�)+�ýڵ�뾶+�ýڵ�Բ�� */
     NodeslotSize = 4 + 4 + 4 + 4 + 4 + 4 + sizeof(float) + sizeof(float) * root->dim;
     NodeslotNum = PAGE_SIZE /  NodeslotSize;
     /*���ݽڵ�ҳ�У�ÿ���۵Ĵ�СΪ ��������ı��(4���ֽ�) + ��������������� + �������С  ps:��ΪΪ������¼�����������С�̶�Ϊ20������*��ÿ�����ݴ�С������ʱֻ���ʵ����������*/ 
     DataslotSize = 4 + 4 + sizeof(float) * (root->dim + 1) * 20;
     DataslotNum = PAGE_SIZE /  DataslotSize;
} 

void BallTree::readNode(Node* node, const char* index_path) {
    // �ȹ���node����ڵ�
    int NodePageID = node->LeftID / NodeslotNum;
    int NodeSlotID = node->LeftID % NodeslotNum;
    CurPageID = NodePageID;
    char *node_path = new char[256];
    sprintf(node_path, "%s/NodePage%d.txt", index_path, NodePageID);
    ifstream in(node_path, ios::binary | ios::in);

    // Ѱ�����ڵĲ���ƫ��
    
    char* buffer = new char[NodeslotSize];
    for (int i = 0; i < NodeSlotID; i++) {
        in.read(buffer, NodeslotSize);
    }

    // ��ȡ���ݹ���ڵ�
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

    //�ٹ���node���ҽڵ�
    NodePageID = node->RightID / NodeslotNum;
    NodeSlotID = node->RightID % NodeslotNum;
    if (1) {
        CurPageID = NodePageID;
        in.close();
        char *node_path = new char[256];
        sprintf(node_path, "%s/NodePage%d.txt", index_path, NodePageID);
        in.open(node_path, ios::binary | ios::in);
    }
    
    // Ѱ�����ڵĲ���ƫ��
    buffer = new char[NodeslotSize];
    for (int i = 0; i < NodeSlotID; i++) {
        in.read(buffer, NodeslotSize);
    }

    // ��ȡ���ݹ���ڵ�
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
	if (T->count < 20) {           //��һ���ڵ��е�������С��20������Ҷ�ӽڵ㡣 
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

