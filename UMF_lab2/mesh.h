#ifndef MESH_H
#define MESH_H

#include <vector>
#include <cassert>

class Mesh1D {
public:
    // nodes – массив координат узлов (от 0 до L)
    Mesh1D(const std::vector<double>& nodes) : nodes_(nodes) {
        nNodes_ = nodes_.size();
        nElem_ = nNodes_ - 1;
        h_.resize(nElem_);
        for (int i = 0; i < nElem_; ++i) {
            h_[i] = nodes_[i + 1] - nodes_[i];
            assert(h_[i] > 0);
        }
    }

    int getNumNodes() const { return nNodes_; }
    int getNumElements() const { return nElem_; }
    double getNode(int i) const { return nodes_[i]; }
    double getH(int i) const { return h_[i]; }
    const std::vector<double>& getNodes() const { return nodes_; }

    // Получить глобальные номера узлов элемента
    void getElementNodes(int elem, int& i1, int& i2) const {
        i1 = elem;
        i2 = elem + 1;
    }

private:
    std::vector<double> nodes_;
    std::vector<double> h_;
    int nNodes_, nElem_;
};

#endif