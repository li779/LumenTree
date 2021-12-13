#include<nori/octree.h>

NORI_NAMESPACE_BEGIN
struct treeNode
    {
        void Init(int meshIndex_, int flag_, int nTri_, treeNode* child_){
            meshIndex = meshIndex_;
            flag = flag_;
            nTri = nTri_;
            child = child_;
        };
        int meshIndex;
        int flag;
        int nTri;
        treeNode* child;
    };
Octree::Octree(Mesh* mesh){
    bounds = mesh->getBoundingBox();
    int nums = mesh->getTriangleCount();
    std::vector<int> triNums(nums);
    for (size_t i = 0; i < nums; ++i){
        triNums[i] = i;
        tribounds.push_back(mesh->getBoundingBox(i));
    } 
    root = new treeNode();
    root->Init(0,0,nums,nullptr);
    buildTree(triNums, bounds, root, nums);

}
void Octree::buildTree(std::vector<int> meshIndex, BoundingBox3f bound, treeNode *cur, int nTri){
    // leaf condition
    if(nTri<5){
        cur->nTri = nTri;
        cur->meshIndex = meshIndices.size();
        for(size_t i = 0; i < nTri; i++){
            meshIndices.push_back(meshIndex[i]);
        }
        cur->flag = 1;
    }else{
        std::vector<BoundingBox3f> subbounds;
        Point3f halfsize =0.5*(bound.max-bound.min);
        cur->child = new treeNode[8];
        std::vector<std::vector<int>> indexArray(8);
        int index = 0;
        // get all bounding boxes for all 8 subnodes
        for (int i = 0; i < 2; i++){
            for (int j = 0; j < 2; j++){
                for (int k = 0; k < 2; k++){
                    Point3f add(i*halfsize.x(),j*halfsize.y(),k*halfsize.z());
                    subbounds.push_back(BoundingBox3f(bound.min+add,bound.min+add+halfsize));
                }
            }
        }
        // search through all mesh triangles
        for(int index = 0; index < nTri; index++){
            BoundingBox3f bbox = mesh->getBoundingBox(meshIndex[index]);
            for(int i = 0; i< 8; i++){
                if(bbox.overlaps(subbounds[i])){
                    indexArray[i].push_back(meshIndex[index]);
                }
            }
        }
        for(int i = 0; i < 8; i++){
            cur->child->Init(0,0,indexArray[i].size(),nullptr);
            buildTree(indexArray[i],subbounds[i],&(cur->child)[i],indexArray[i].size());
        }
    }
}

NORI_NAMESPACE_END