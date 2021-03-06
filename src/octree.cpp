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
Octree::Octree(std::vector<Mesh *> &mesh, BoundingBox3f &bbox){
    this->mesh = mesh;
    bounds = bbox;
    std::vector<std::pair<int,int>> triNums;
    int total_num = 0;
    for (int i = 0; i < mesh.size(); i++){
        int nums = mesh[i]->getTriangleCount();
        total_num += nums;
        std::vector<BoundingBox3f> subtribounds(nums);
        for (int j = 0; j < nums; ++j){
            triNums.push_back({i,j});
            subtribounds[j] = mesh[i]->getBoundingBox(j);
        }
        tribounds.push_back(subtribounds);
    }
     
    root = new treeNode();
    root->Init(0,0,total_num,nullptr);
    int treeDepth = buildTree(triNums, bounds, root, total_num);
    std::cout<<"Tree built: depth="<<treeDepth<<std::endl;
    tribounds.clear();
    tribounds.shrink_to_fit();
}
int Octree::buildTree(std::vector<std::pair<int,int>> meshIndex, BoundingBox3f bound, treeNode *cur, int nTri){
    int depth = 1;
    // leaf condition
    if(nTri < 8 || bound.getVolume() < 0.00001f){
        cur->nTri = nTri;
        cur->meshIndex = meshIndices.size();
        for(int i = 0; i < nTri; i++){
            meshIndices.push_back(meshIndex[i]);
        }
        cur->flag = 1;
    }
    else{
        std::vector<BoundingBox3f> subbounds;
        Point3f halfsize =0.5*(bound.max-bound.min);
        cur->child = new treeNode[8];
        std::vector<std::vector<std::pair<int,int>>> indexArray(8);
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
            for(int i = 0; i< 8; i++){
                if(tribounds[meshIndex[index].first][meshIndex[index].second].overlaps(subbounds[i])){
                    indexArray[i].push_back(meshIndex[index]);
                }
            }
        }
        for(int i = 0; i < 8; i++){
            (cur->child)[i].Init(0,0,indexArray[i].size(),nullptr);
            depth = std::max(depth,1+buildTree(indexArray[i],subbounds[i],&(cur->child)[i],indexArray[i].size()));
        }
    }
    return depth;
}
bool Octree::IntersectOctree(Ray3f &ray, Intersection &its, bool shadowRay, int &triIndex) const {
    if(!bounds.rayIntersect(ray))
        return false;
    if(shadowRay)
        return Octree::recursiveTest(*root,ray,bounds,triIndex);
    else
        return Octree::recursive(*root,ray,its,bounds,triIndex);
}

bool Octree::recursive(const treeNode &node, Ray3f &ray_, Intersection &its, BoundingBox3f bbox, int &triIndex) const {
    bool foundIntersection = false;

    // only check triangles of node and its children if ray intersects with node bbox
    if (!bbox.rayIntersect(ray_)) {
        return false;
    }
    
    // search through all triangles in node
    if (node.flag==1)
    for (uint32_t i = 0; i < node.nTri; ++i) {
        float u, v, t;
        uint32_t mesh_idx = meshIndices[node.meshIndex+i].first;
        uint32_t triangle_idx = meshIndices[node.meshIndex+i].second;
        if (mesh[mesh_idx]->rayIntersect(triangle_idx, ray_, u, v, t) && t < ray_.maxt) {
            /* An intersection was found! Can terminate
            immediately if this is a shadow ray query */
            ray_.maxt = its.t = t;
            its.uv = Point2f(u, v);
            its.mesh = mesh[mesh_idx];
            triIndex = triangle_idx;
            foundIntersection = true;
        }
    }
    else {
        std::vector<BoundingBox3f> subbounds;
        Point3f halfsize =0.5*(bbox.max-bbox.min);
        for (int i = 0; i < 2; i++){
            for (int j = 0; j < 2; j++){
                for (int k = 0; k < 2; k++){
                    Point3f add(i*halfsize.x(),j*halfsize.y(),k*halfsize.z());
                    subbounds.push_back(BoundingBox3f(bbox.min+add,bbox.min+add+halfsize));
                }
            }
        }
        float T=0;
        std::vector<std::pair<int, float>> children;
        for (int i = 0; i < 8; i++){
            float nearT, farT;
            if(subbounds[i].rayIntersect(ray_, nearT, farT) &&(ray_.mint <= farT && nearT <= ray_.maxt)){
                T = std::max(T,farT);
                children.push_back({i,nearT});
            }
        }

        std::sort(children.begin(), children.end(), [](const std::pair<int, float>& l, const std::pair<int, float>& r) {
            return l.second < r.second;
        });

        for (auto childpair: children) {
            if(T < childpair.second)
                break;
            int index = childpair.first;
            foundIntersection = Octree::recursive(node.child[index], ray_, its, subbounds[index], triIndex) || foundIntersection;
            if(foundIntersection)
                T = its.t;
        }
    }
    return foundIntersection;
}

bool Octree::recursiveTest(const treeNode &node, Ray3f &ray_, BoundingBox3f bbox, int &triIndex) const {
    bool foundIntersection = false;
    
    // search through all triangles in node
    if (node.flag==1)
    for (uint32_t i = 0; i < node.nTri; ++i) {
        float u, v, t;
        uint32_t mesh_idx = meshIndices[node.meshIndex+i].first;
        uint32_t triangle_idx = meshIndices[node.meshIndex+i].second;
        if (mesh[mesh_idx]->rayIntersect(triangle_idx, ray_, u, v, t) && t < ray_.maxt) {
            return true;
        }
    }
    else {
        std::vector<BoundingBox3f> subbounds;
        Point3f halfsize =0.5*(bbox.max-bbox.min);
        for (int i = 0; i < 2; i++){
            for (int j = 0; j < 2; j++){
                for (int k = 0; k < 2; k++){
                    Point3f add(i*halfsize.x(),j*halfsize.y(),k*halfsize.z());
                    subbounds.push_back(BoundingBox3f(bbox.min+add,bbox.min+add+halfsize));
                }
            }
        }
        std::vector<std::pair<int, float>> children;
        for (int i = 0; i < 8; i++){
            if(subbounds[i].rayIntersect(ray_))
                children.push_back({i,subbounds[i].distanceTo(ray_.o)});
        }

        std::sort(children.begin(), children.end(), [](const std::pair<int, float>& l, const std::pair<int, float>& r) {
            return l.second < r.second;
        });

        for (auto childpair: children) {
            int index = childpair.first;
            foundIntersection = Octree::recursiveTest(node.child[index], ray_, subbounds[index], triIndex);
            if (foundIntersection)
                return true;
        }
    }
    return foundIntersection;
}

NORI_NAMESPACE_END